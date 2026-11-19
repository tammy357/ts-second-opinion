"""Read a reaction from its reactant and product geometries, so it can be checked from three .xyz files.

Without an atom-mapped SMILES, TS Second Opinion reads the bonds, bond orders and stereocentres from the reactant and
product structures and writes the atom-mapped reaction SMILES itself (atom i in the files gets map number i + 1). The
rest of the pipeline then runs exactly as it does for a SMILES you supply, and the inferred SMILES is returned so you
can check it.

  bonds          a pair is bonded when its distance is below factor x (sum of covalent radii). The factors were chosen
                 on RDB7's 10,734 reactant/product pairs (285,292 bonds); see scripts/validate_geometry_input.py.
  bond orders    RDKit's xyz2mol (rdDetermineBonds.DetermineBondOrders) for closed-shell structures; open-shell
                 structures (multiplicity > 1), or ones it cannot assign, keep single bonds, so the double-bond and
                 aromatic-ring planarity rule has nothing to check there.
  stereocentres  RDKit's assignment from the 3D structure (AssignStereochemistryFrom3D).

Needs RDKit (pip install rdkit). The in-browser app has no RDKit, so there the atom-mapped SMILES is still required.
"""
import numpy as np

from .audit import COVALENT_RADIUS_A
from .data import SYMBOL, ATOMIC_NUMBER

FACTOR_H_H = 1.30          # H2 sits at 1.20x the radii sum; non-bonded H...H pairs at 1.7x or more
FACTOR_X_H = 1.17          # heavy atom to hydrogen
FACTOR_CC_CO = 1.13        # C-C and C-O: long non-bonded contacts in strained rings start near 1.15x
FACTOR_OTHER = 1.22        # other heavy-atom pairs (N-O, N-N and C-N bonds can be long)


def _radius(z):
    r = COVALENT_RADIUS_A.get(int(z))
    if r is None:
        from rdkit import Chem
        r = Chem.GetPeriodicTable().GetRcovalent(int(z))
    return r


def _factor(za, zb):
    a, b = sorted((int(za), int(zb)))
    if a == b == 1:
        return FACTOR_H_H
    if a == 1:
        return FACTOR_X_H
    if (a, b) in ((6, 6), (6, 8)):
        return FACTOR_CC_CO
    return FACTOR_OTHER


# Most neighbours an atom can have. Where a strained cage puts two non-bonded atoms within the distance cut-off (as in
# bicyclo[1.1.1] skeletons), the atom ends up over-bonded; the relatively longest of its contacts is then dropped.
MAX_NEIGHBOURS = {1: 1, 5: 4, 6: 4, 7: 4, 8: 3, 9: 1, 14: 4, 17: 1, 35: 1, 53: 1}


def bonded_pairs(atoms, xyz):
    """[(i, j)] with i < j for every pair closer than its bonding cut-off, after the valence check above."""
    xyz = np.asarray(xyz, float)
    rad = np.array([_radius(z) for z in atoms])
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    ratio = {}
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            r = d[i, j] / (rad[i] + rad[j])
            if r < _factor(atoms[i], atoms[j]):
                ratio[(i, j)] = r
    while True:
        over = [a for a in range(len(atoms)) if sum(a in b for b in ratio) > MAX_NEIGHBOURS.get(int(atoms[a]), 99)]
        if not over:
            break
        worst = max((b for b in ratio if over[0] in b), key=ratio.get)
        del ratio[worst]
    return sorted(ratio)


def _side_smiles(atoms, xyz, charge, assign_orders):
    """Atom-mapped SMILES for one side, and whether bond orders were assigned."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdDetermineBonds
    RDLogger.DisableLog('rdApp.*')

    def build():
        m = Chem.RWMol()
        for z in atoms:
            a = Chem.Atom(int(z))
            a.SetNoImplicit(True)
            m.AddAtom(a)
        for i, j in bonded_pairs(atoms, xyz):
            m.AddBond(int(i), int(j), Chem.BondType.SINGLE)
        conf = Chem.Conformer(len(atoms))
        for k, (x, y, z) in enumerate(np.asarray(xyz, float)):
            conf.SetAtomPosition(k, (float(x), float(y), float(z)))
        m.AddConformer(conf, assignId=True)
        return m.GetMol()

    mol, orders = build(), False
    if assign_orders:
        try:
            rdDetermineBonds.DetermineBondOrders(mol, charge=int(charge))
            # xyz2mol can 'succeed' while leaving unpaired electrons on a closed-shell structure (for example
            # carbenes and some zwitterions): treat that as not assigned rather than trust a half-assigned structure
            orders = not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        except Exception:
            orders = False
        if not orders:
            mol = build()
    if not orders:
        mol.UpdatePropertyCache(strict=False)
    Chem.AssignStereochemistryFrom3D(mol)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)
        a.SetNoImplicit(True)
    return Chem.MolToSmiles(mol), orders


def reaction_smiles_from_geometry(symbols, reactant_xyz, product_xyz, charge=0, multiplicity=1):
    """{'smiles', 'bond_orders', 'notes'} for a reaction given as element symbols and two (N, 3) geometries in the
    same atom order. Raises ValueError if RDKit is missing or the inputs do not match."""
    try:
        import rdkit  # noqa: F401
    except ImportError:
        raise ValueError('reading bonds from the geometry needs RDKit (pip install rdkit); without it, give an '
                         'atom-mapped reaction SMILES')
    try:
        atoms = [ATOMIC_NUMBER[s] for s in symbols]
    except KeyError as e:
        raise ValueError(f'unknown element {e}')
    R, P = np.asarray(reactant_xyz, float), np.asarray(product_xyz, float)
    if R.shape != (len(atoms), 3) or P.shape != (len(atoms), 3):
        raise ValueError('the reactant and product must list the same atoms in the same order')
    closed_shell = int(multiplicity) == 1
    rs, ro = _side_smiles(atoms, R, charge, closed_shell)
    ps, po = _side_smiles(atoms, P, charge, closed_shell)
    notes = ['Bonds, bond orders and stereocentres were read from the reactant and product geometries (atom i in '
             'your files is map number i + 1). Check the inferred SMILES; to override it, give your own.']
    if not closed_shell:
        notes.append(f'Bond orders were not assigned because the system is open-shell (multiplicity {multiplicity}), so '
                     'the planarity rule for double bonds and aromatic rings has nothing to check.')
    elif not (ro and po):
        sides = ' and '.join(n for n, ok in (('reactant', ro), ('product', po)) if not ok)
        notes.append(f'Bond orders could not be assigned for the {sides} (no closed-shell bonding pattern fits, as for '
                     'some carbenes and zwitterions), so the planarity rule has nothing to check there. Give an '
                     'atom-mapped SMILES if you need it.')
    if not set(bonded_pairs(atoms, R)) ^ set(bonded_pairs(atoms, P)):
        notes.append('No bond forms or breaks between the reactant and product geometries as given; check that '
                     'they are different structures in the same atom order.')
    return {'smiles': f'{rs}>>{ps}', 'bond_orders': bool(ro and po), 'notes': notes}


def symbols_of(atoms):
    return [SYMBOL[int(z)] for z in atoms]
