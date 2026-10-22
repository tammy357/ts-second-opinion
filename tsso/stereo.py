"""Per-stereocentre handedness of a TS guess, for display.

Uses the same signed-volume tests and thresholds as the two gated stereo checks (tsso/checks/chirality.py and
product_stereo.py), so what is shown always agrees with the verdict. No reference TS is needed.
  ok = True   same handedness as the reactant (spectator stereocentres) or the product (product stereocentres)
  ok = False  inverted: the gated rule fails on this centre
  ok = None   too flat to judge (the rule does not decide on it)
"""
import numpy as np


def _vol_spectator(X, c, ns):          # as in checks/chirality.py
    o = X[ns[3]] if len(ns) >= 4 else X[c]
    a, b, d = X[ns[0]] - o, X[ns[1]] - o, X[ns[2]] - o
    return float(np.dot(a, np.cross(b, d)))


def _vol_product(X, c, ns):            # as in checks/product_stereo.py
    v = [(X[k] - X[c]) / (np.linalg.norm(X[k] - X[c]) + 1e-9) for k in ns[:3]]
    return float(np.dot(v[0], np.cross(v[1], v[2])))


def stereocentre_status(rxn, xyz):
    """[{atom, label, compared_with, ok}] for every stereocentre the gated stereo rules look at."""
    out, sym = [], rxn.symbols
    nbr, nbp = rxn.neighbours('reactant'), rxn.neighbours('product')
    for c in range(rxn.n):
        if rxn.reactant_chiral_tags[c] not in (1, 2) or c in rxn.reaction_centre or len(nbr[c]) < 3:
            continue
        ns = nbr[c][:4]
        ref = _vol_spectator(rxn.reactant_xyz, c, ns)
        if abs(ref) < 0.3:
            ok = None
        else:
            now = _vol_spectator(xyz, c, ns)
            ok = True if now * ref > 0 else (False if abs(now) > 0.4 else None)
        out.append({'atom': c, 'label': f'{sym[c]}{c + 1}', 'compared_with': 'reactant', 'ok': ok})
    for c in range(rxn.n):
        if rxn.product_chiral_tags[c] not in (1, 2) or len(nbp[c]) < 3:
            continue
        ref = _vol_product(rxn.product_xyz, c, nbp[c])
        if abs(ref) < 0.1:
            ok = None
        else:
            now = _vol_product(xyz, c, nbp[c])
            thr = 0.3 if c in rxn.reaction_centre else 0.15
            ok = True if now * ref > 0 else (False if abs(now) > thr else None)
        out.append({'atom': c, 'label': f'{sym[c]}{c + 1}', 'compared_with': 'product', 'ok': ok})
    return out
