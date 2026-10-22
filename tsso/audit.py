"""Does a reaction's atom-mapped SMILES agree with its own reactant and product geometries?

A dataset (or an uploaded reaction) is flagged when a structure has two atoms at single-bond distance that the SMILES
says are not bonded, or a SMILES bond whose atoms are far apart. Thresholds are deliberately conservative, relative to
the sum of covalent radii (Cordero et al., Dalton Trans. 2008): "bonded" below 1.10x (every common single bond except
O-O fits), "not bonded" above 1.60x. Pairs with an element not in the table (most metals) are skipped.

    python -m tsso.cli audit            # the RDB7 reactions bundled with the repository
"""
import numpy as np

from .data import SYMBOL

COVALENT_RADIUS_A = {1: 0.31, 2: 0.28, 3: 1.28, 4: 0.96, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 10: 0.58,
                     11: 1.66, 12: 1.41, 13: 1.21, 14: 1.11, 15: 1.07, 16: 1.05, 17: 1.02, 18: 1.06, 19: 2.03,
                     20: 1.76, 31: 1.22, 32: 1.20, 33: 1.19, 34: 1.20, 35: 1.20, 36: 1.16, 50: 1.39, 51: 1.39,
                     52: 1.38, 53: 1.39, 54: 1.40}
BONDED_BELOW, APART_ABOVE = 1.10, 1.60


def audit_reaction(rxn):
    """[{side, kind, atoms, distance_A, ratio}] for every contradiction between the SMILES and the geometries."""
    out = []
    radii = [COVALENT_RADIUS_A.get(int(z)) for z in rxn.atoms]
    for side, xyz, bonds in (('reactant', rxn.reactant_xyz, rxn.reactant_bonds), ('product', rxn.product_xyz, rxn.product_bonds)):
        if xyz is None or not np.isfinite(xyz).all():
            continue
        d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
        for i in range(rxn.n):
            for j in range(i + 1, rxn.n):
                if radii[i] is None or radii[j] is None:
                    continue
                ratio = d[i, j] / (radii[i] + radii[j])
                said = (i, j) in bonds
                if said and ratio > APART_ABOVE:
                    kind = 'the SMILES has this bond, but the atoms are far apart'
                elif not said and ratio < BONDED_BELOW:
                    kind = 'bonded in the geometry, but not in the SMILES'
                else:
                    continue
                out.append({'side': side, 'kind': kind, 'atoms': f'{SYMBOL[int(rxn.atoms[i])]}{i + 1}-{SYMBOL[int(rxn.atoms[j])]}{j + 1}',
                            'distance_A': round(float(d[i, j]), 3), 'ratio': round(float(ratio), 2)})
    return out


def describe(issues, limit=3):
    """Plain-English warning lines for the first few issues."""
    return [f"The {x['side']} geometry contradicts the SMILES at {x['atoms']} ({x['distance_A']:.2f} A): {x['kind']}. "
            'Check the atom mapping or the structure.' for x in issues[:limit]]


def audit_dataset(reactions):
    """Rows for every flagged reaction: (reaction id, mapped SMILES, issues)."""
    return [(r.id, r.smiles, iss) for r in reactions for iss in [audit_reaction(r)] if iss]
