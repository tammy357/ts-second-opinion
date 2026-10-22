"""Stereo-aware leaderboard for every TS-generator output we have (all released by the GoFlow authors).

"Success" in TS-generation papers usually means D-MAE <= 0.2 A, which cannot see handedness. This table adds the
handedness check (tsso.metrics), so a guess counts only if it is close AND correctly handed.

Usage: python scripts/leaderboard.py  -> results/leaderboard.json, results/leaderboard.md
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.metrics import evaluate
from tsso.geometry import consensus_index
from tsso.triage import load_checks, triage
from tsso.sandbox import run_check

checks = load_checks(); MIRROR = np.array([1.0, 1.0, -1.0]); STEREO = ('chirality', 'product_stereo')
cut = json.load(open('results/triage_summary.json'))['amber_cutoff_A']
KEYS = ('mean_dmae_A', 'mean_rmsd_mirror_allowed_A', 'success_dmae_only', 'handedness_correct', 'success_dmae_and_handedness')


def score(R, X):
    s = evaluate(X, [r.ts_ref_xyz for r in R], R)
    return {k: s[k] for k in KEYS}


def kabsch_align(mob, ref):
    a, b = mob - mob.mean(0), ref - ref.mean(0)
    U, S, Vt = np.linalg.svd(a.T @ b)
    if np.linalg.det(U @ Vt) < 0:
        U[:, -1] *= -1
    return a @ (U @ Vt) + ref.mean(0)


def repaired_single(r):
    """Rules + mirror repair on a single guess (no ensemble)."""
    x = r.preds[0]
    f = [n for n, fn in checks.items() if run_check(fn, r, x)[0] is False]
    if f and set(f) & set(STEREO):
        m = x * MIRROR
        if all(run_check(fn, r, m)[0] is not False for fn in checks.values()):
            return m
    return x


rows = []
R8 = load_reactions()
seeds = [score(R8, [r.preds[k] for r in R8]) for k in range(8)]
rows.append({'method': 'GoFlow, 8 seeds (mean ± s.d.)', 'size': 256, 'test_set': 'random split (1,194)',
             **{k: float(np.mean([s[k] for s in seeds])) for k in KEYS},
             'sd': {k: float(np.std([s[k] for s in seeds])) for k in KEYS}})
RS = load_reactions('data/goflow_sizes_rdb7_test.json.gz')
for k, size in enumerate((64, 128, 256, 256, 512)):
    rows.append({'method': f'GoFlow, size run {k + 1}', 'size': size, 'test_set': 'random split (1,194)',
                 **score(RS, [r.preds[k] for r in RS])})
for f, name in (('goflow_split0_test', 'reaction-core split'), ('goflow_split1_test', 'barrier-height split')):
    Rx = load_reactions(f'data/{f}.json.gz')
    rows.append({'method': 'GoFlow, 1 model', 'size': 256, 'test_set': f'{name} ({len(Rx):,})',
                 **score(Rx, [r.preds[0] for r in Rx])})
    rows.append({'method': 'GoFlow + TS Second Opinion (rules + mirror repair, no ensemble)', 'size': 256,
                 'test_set': f'{name} ({len(Rx):,})', **score(Rx, [repaired_single(r) for r in Rx])})
rows.append({'method': 'No ML: midpoint of aligned reactant and product', 'size': None, 'test_set': 'random split (1,194)',
             **score(R8, [0.5 * (r.reactant_xyz + kabsch_align(r.product_xyz, r.reactant_xyz)) for r in R8])})
rows.append({'method': 'GoFlow 8 seeds + TS Second Opinion (final answer)', 'size': 256, 'test_set': 'random split (1,194)',
             **score(R8, [triage(r, checks, cut).xyz for r in R8])})
rows.append({'method': 'GoFlow 5 size runs + TS Second Opinion (final answer)', 'size': None, 'test_set': 'random split (1,194)',
             **score(RS, [triage(r, checks, cut).xyz for r in RS])})

json.dump({'note': 'success = D-MAE <= 0.2 A vs the RDB7 reference TS; handedness = every stereocentre matches the reference',
           'rows': rows}, open('results/leaderboard.json', 'w'), indent=2)

pc = lambda v: f'{100 * v:.1f}%'
lines = ['| Method | Model size | Test set | D-MAE (Å) | RMSD, mirror allowed (Å) | Success (D-MAE ≤ 0.2 Å) | Correct handedness | **Success with correct handedness** |',
         '|---|---|---|---|---|---|---|---|']
for r in rows:
    sd = r.get('sd')
    cell = lambda k, f: f(r[k]) + (f' ± {f(sd[k])}' if sd else '')
    lines.append(f"| {r['method']} | {r['size'] or '–'} | {r['test_set']} | {cell('mean_dmae_A', lambda v: f'{v:.3f}')} | "
                 f"{cell('mean_rmsd_mirror_allowed_A', lambda v: f'{v:.3f}')} | {cell('success_dmae_only', pc)} | "
                 f"{cell('handedness_correct', pc)} | **{cell('success_dmae_and_handedness', pc)}** |")
Path('results/leaderboard.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
