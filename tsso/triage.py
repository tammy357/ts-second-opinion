"""Triage one reaction: pick the best TS guess, repair it if possible, and decide whether to trust it.

Inputs per reaction: K model predictions (here 8 GoFlow seeds) plus the gated rule checks.
  1. Run every check on every prediction.
  2. Pick the consensus prediction among those that pass all checks (fallback: all predictions).
  3. If it still breaks a stereo rule, try its mirror image (identical distances, opposite handedness).
  4. Verdict:
       red    - still breaks a rule        -> re-generate, or send straight to DFT with care
       amber  - models disagree a lot      -> run a DFT TS optimisation before trusting it
       green  - passes rules, models agree -> use as is
The amber cut-off is a budget: the most-disputed fraction of reactions (set from the disagreement
distribution alone, never from reference answers).
"""
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

from .geometry import disagreement_matrix, consensus_index, disputed_pair
from .sandbox import compile_check, run_check
from .gate import load_check_sources

STEREO_CHECKS = ('chirality', 'product_stereo')
MIRROR = np.array([1.0, 1.0, -1.0])


@dataclass
class Triage:
    reaction_id: int
    verdict: str                 # 'green' | 'amber' | 'red'
    xyz: np.ndarray              # the recommended TS guess
    chosen_model: int
    mirrored: bool
    disagreement: float          # mean D-MAE of the chosen guess to the other models (Angstrom)
    failed_rules: list           # [(name, detail)] for the recommended guess
    rejected_models: dict        # model index -> [rule names it failed]
    disputed: tuple              # ((i, j), spread in Angstrom, is the pair a forming/breaking bond)
    reasons: list = field(default_factory=list)


def load_checks():
    return {n: compile_check(s) for n, s in load_check_sources().items()}


def triage(rxn, checks, amber_cutoff, rule_text=None):
    preds = rxn.preds
    K = len(preds)
    results = {k: {n: run_check(fn, rxn, preds[k]) for n, fn in checks.items()} for k in range(K)}
    passing = np.array([all(ok is not False for ok, _ in results[k].values()) for k in range(K)])
    rejected = {k: [n for n, (ok, _) in results[k].items() if ok is False] for k in range(K) if not passing[k]}
    c = consensus_index(preds, allowed=passing)
    x, mirrored, mirror_note = preds[c].copy(), False, None
    failed = [(n, d) for n, (ok, d) in results[c].items() if ok is False]
    if failed and any(n in STEREO_CHECKS for n, _ in failed):
        mx = x * MIRROR
        mres = {n: run_check(fn, rxn, mx) for n, fn in checks.items()}
        if all(ok is not False for ok, _ in mres.values()):
            x, mirrored, failed = mx, True, []
        elif any(n in STEREO_CHECKS and ok is False for n, (ok, _) in mres.items()):
            mirror_note = ('Its mirror image was tried but also breaks a stereo rule: only some stereocentres are inverted '
                           '(the TS of a diastereomer, which has a different energy), so mirroring would invert the correct '
                           'ones too and cannot repair it.')
        else:
            mirror_note = ('Its mirror image would fix the handedness, but the guess also breaks a rule that mirroring '
                           'cannot fix (mirror images keep every distance).')
    sym = rxn.symbols
    reasons = []
    n_pass = int(passing.sum())
    reasons.append(f'{n_pass} of {K} model guesses pass every rule.')
    if K == 1:  # a single guess: rules and repair still work, but there is no disagreement signal
        if mirrored:
            reasons.append('The guess had a stereocentre the wrong way round; its mirror image fixes this '
                           'without changing any interatomic distance. The unrepaired guess has the same energy and '
                           'frequencies, so a DFT frequency check would not have caught it.')
        for n, d in failed:
            reasons.append(f'Still breaks "{n}": {d}.')
        if mirror_note:
            reasons.append(mirror_note)
        reasons.append('Only one guess was given, so model disagreement could not be measured; '
                       'supply several guesses (models or samples) for a full triage.')
        verdict = 'red' if failed else 'green'
        return Triage(rxn.id, verdict, x, c, mirrored, float('nan'), failed, rejected, ((0, 0), 0.0, False), reasons)
    M = disagreement_matrix(preds)
    u = float(M[c].sum() / (K - 1))
    pair, spread, _ = disputed_pair(preds)
    reacting = pair in set(rxn.changed_bonds)
    if rejected:
        counts = {}
        for names in rejected.values():
            for n in names:
                counts[n] = counts.get(n, 0) + 1
        reasons.append('Rejected guesses broke: ' + ', '.join(f'{n} ({v})' for n, v in sorted(counts.items())) + '.')
    if mirrored:
        reasons.append('The best guess had a stereocentre the wrong way round; its mirror image fixes this '
                       'without changing any interatomic distance. The unrepaired guess has the same energy and '
                       'frequencies, so a DFT frequency check would not have caught it.')
    for n, d in failed:
        reasons.append(f'Still breaks "{n}": {d}.')
    if mirror_note:
        reasons.append(mirror_note)
    a, b = pair
    kind = 'a forming/breaking bond' if reacting else 'not a forming/breaking bond'
    reasons.append(f'Models disagree most about {sym[a]}{a + 1}-{sym[b]}{b + 1} ({kind}): spread {spread:.2f} A.')
    reasons.append(f'Mean disagreement of the chosen guess with the other models: {u:.3f} A '
                   f'(review cut-off {amber_cutoff:.3f} A).')
    verdict = 'red' if failed else ('amber' if u >= amber_cutoff else 'green')
    return Triage(rxn.id, verdict, x, c, mirrored, u, failed, rejected, (pair, spread, reacting), reasons)


def amber_cutoff_for_budget(reactions, checks, budget=0.2):
    """Disagreement value above which the top `budget` fraction of reactions is sent for review."""
    us = []
    for r in reactions:
        M = disagreement_matrix(r.preds)
        passing = np.array([all(run_check(fn, r, p)[0] is not False for fn in checks.values()) for p in r.preds])
        c = consensus_index(r.preds, allowed=passing)
        us.append(M[c].sum() / (len(r.preds) - 1))
    return float(np.quantile(us, 1 - budget))


# ---------------------------------------------------------------- outputs
ORCA_TEMPLATE = """! wB97X D4 def2-TZVP OptTS FREQ
# TS optimisation for {rid}; method follows GoFlow's tsopt.inp (github.com/heid-lab/goflow)
# Guess: {source}

%maxcore 1500
%pal
  nprocs 2
end

%scf
  Convergence Tight
  MaxIter 500
end

%geom
  MaxIter 500
end

* xyzfile {charge} {mult} {xyzname}
"""


def write_xyz(path, symbols, xyz, comment=''):
    lines = [str(len(symbols)), comment] + [f'{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}' for s, (x, y, z) in zip(symbols, xyz)]
    Path(path).write_text('\n'.join(lines) + '\n')


def write_orca(outdir, rxn, t: Triage):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    stem = f'rxn{rxn.id}_ts_guess'
    src = f'GoFlow seed model {t.chosen_model + 1}' + (' (mirror image)' if t.mirrored else '') + f'; verdict {t.verdict}'
    write_xyz(outdir / f'{stem}.xyz', rxn.symbols, t.xyz, f'RDB7 {rxn.id} {rxn.smiles}')
    (outdir / f'{stem}.inp').write_text(ORCA_TEMPLATE.format(rid=f'RDB7 reaction {rxn.id}', source=src, xyzname=f'{stem}.xyz',
                                                              charge=0, mult=1))
    return outdir / f'{stem}.inp'
