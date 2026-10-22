"""Stereo-aware evaluation for any transition-state generator.

The usual scores cannot see handedness: D-MAE compares interatomic distances (identical for mirror images), and the
RMSD code shared by OA-ReactDiff, React-OT and GoFlow is called with ignore_chirality=True. These functions add the
missing check, so a TS model's "success rate" counts only guesses that are close AND correctly handed.

    from tsso.metrics import evaluate
    report = evaluate(predictions, references, reaction_smiles)   # lists of (N, 3) arrays, mapped SMILES
"""
import numpy as np

from .geometry import dmae

BAD_DMAE_A = 0.2


def rmsd(a, b, allow_mirror=False):
    """Kabsch RMSD. allow_mirror=True reproduces ignore_chirality=True (best of the guess and its mirror image)."""
    a0, b0 = a - a.mean(0), b - b.mean(0)
    U, S, Vt = np.linalg.svd(b0.T @ a0)
    d = np.sign(np.linalg.det(U @ Vt))
    proper = U.copy()
    if d < 0:
        proper[:, -1] *= -1
    best = float(np.sqrt(((b0 @ (proper @ Vt) - a0) ** 2).sum(1).mean()))
    if allow_mirror:
        best = min(best, float(np.sqrt(((b0 @ (U @ Vt) - a0) ** 2).sum(1).mean())))
    return best


def stereocentres(rxn):
    """(centre, neighbour list) for reactant stereocentres outside the reaction centre and product stereocentres."""
    nbr, nbp = rxn.neighbours('reactant'), rxn.neighbours('product')
    c = [(i, nbr[i]) for i in range(rxn.n) if rxn.reactant_chiral_tags[i] in (1, 2) and i not in rxn.reaction_centre
         and len(nbr[i]) >= 3]
    c += [(i, nbp[i]) for i in range(rxn.n) if rxn.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    return c


def _signed_volume(x, i, nb):
    v = x[nb[:3]] - x[i]
    return float(np.dot(v[0], np.cross(v[1], v[2])))


def handedness_correct(rxn, pred, ref):
    """True if every stereocentre in pred has the same handedness as in the reference TS."""
    return all(np.sign(_signed_volume(pred, i, nb)) == np.sign(_signed_volume(ref, i, nb)) for i, nb in stereocentres(rxn))


def _reaction(smiles, n):
    from .chem import parse_reaction
    from .data import Reaction
    atoms, chiral, rb, pb, pchiral = parse_reaction(smiles, with_product_chirality=True)
    if len(atoms) != n:
        raise ValueError(f'SMILES has {len(atoms)} atoms but the structure has {n}')
    nan = np.full((n, 3), np.nan)
    return Reaction(-1, smiles, float('nan'), np.array(atoms), nan, nan, nan, {(a, b): o for a, b, o in rb},
                    {(a, b): o for a, b, o in pb}, np.array(chiral), np.array(pchiral), np.zeros((0, n, 3)))


def evaluate(predictions, references, reactions, threshold=BAD_DMAE_A):
    """Score predicted TS geometries against references, with and without handedness.

    predictions, references: lists of (N, 3) arrays, atoms listed by atom-map number.
    reactions: list of atom-mapped reaction SMILES, or Reaction objects (one per prediction).
    """
    rows = []
    for p, r, rx in zip(predictions, references, reactions):
        p, r = np.asarray(p, float), np.asarray(r, float)
        rx = _reaction(rx, len(p)) if isinstance(rx, str) else rx
        d = dmae(p, r)
        has = bool(stereocentres(rx))
        ok = handedness_correct(rx, p, r)
        rows.append({'dmae_A': d, 'rmsd_mirror_allowed_A': rmsd(p, r, True), 'rmsd_A': rmsd(p, r, False),
                     'has_stereocentres': has, 'handedness_correct': ok})
    d = np.array([x['dmae_A'] for x in rows]); ok = np.array([x['handedness_correct'] for x in rows])
    has = np.array([x['has_stereocentres'] for x in rows])
    close = d <= threshold
    return {
        'n': len(rows), 'threshold_A': threshold,
        'mean_dmae_A': float(d.mean()),
        'mean_rmsd_mirror_allowed_A': float(np.mean([x['rmsd_mirror_allowed_A'] for x in rows])),
        'mean_rmsd_A': float(np.mean([x['rmsd_A'] for x in rows])),
        'success_dmae_only': float(close.mean()),
        'handedness_correct': float(ok.mean()),
        'handedness_correct_among_chiral': float(ok[has].mean()) if has.any() else None,
        'success_dmae_and_handedness': float((close & ok).mean()),
        'wrong_handed_counted_as_success': float((close & ~ok).mean()),
        'per_reaction': rows,
    }
