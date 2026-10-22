"""Parse atom-mapped reaction SMILES into atoms, bonds and stereo tags indexed by atom map number - 1.

Every atom, hydrogens included, needs a map number, and the numbers must run 1..N on both sides.
Uses RDKit when it is installed, and otherwise a small pure-Python reader (tsso.smiles_lite) that gives the same
answers on these inputs; that is what lets the app run in a web browser. Set TSSO_NO_RDKIT=1 to force the fallback.
"""
import os

try:
    if os.environ.get('TSSO_NO_RDKIT'):
        raise ImportError
    from rdkit import Chem
    _ps = Chem.SmilesParserParams()
    _ps.removeHs = False
    _ps.sanitize = False
    BACKEND = 'rdkit'
except ImportError:
    Chem = None
    BACKEND = 'pure-python'
from . import smiles_lite

ORDER = {'SINGLE': 1.0, 'DOUBLE': 2.0, 'TRIPLE': 3.0, 'AROMATIC': 1.5}
_MAP_ERROR = ('every atom, including each hydrogen, needs a unique atom-map number 1..N '
              '(write hydrogens explicitly, e.g. [C:1]([H:4])([H:5])...)')


def _atoms_and_bonds(smi, radicals=False):
    """[(atomic number, map number, chiral tag, formal charge, radical electrons)], [(atom index, atom index, type)]."""
    if Chem is not None:
        m = Chem.MolFromSmiles(smi, _ps)
        if m is None:
            raise ValueError(f'could not parse SMILES: {smi}')
        if radicals:
            m.UpdatePropertyCache(strict=False)
            Chem.AssignRadicals(m)
        atoms = [(a.GetAtomicNum(), a.GetAtomMapNum(), int(a.GetChiralTag()), a.GetFormalCharge(),
                  a.GetNumRadicalElectrons()) for a in m.GetAtoms()]
        bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType())) for b in m.GetBonds()]
    else:
        m = smiles_lite.mol_from_smiles(smi)
        if radicals:
            smiles_lite.assign_radicals(m)
        atoms = [(a.z, a.map, a.chiral, a.charge, a.radicals) for a in m.atoms]
        bonds = m.bonds
    if sorted(a[1] for a in atoms) != list(range(1, len(atoms) + 1)):
        raise ValueError(_MAP_ERROR)
    return atoms, bonds


def smiles_atom_order(smi):
    """Map numbers in the order the atoms appear in the SMILES (for xyz files written in that order)."""
    return [a[1] for a in _atoms_and_bonds(smi)[0]]


def parse_mapped(smi):
    found, found_bonds = _atoms_and_bonds(smi)
    n = len(found)
    atoms, chiral = [0] * n, [0] * n
    for z, mp, tag, _, _ in found:
        atoms[mp - 1], chiral[mp - 1] = z, tag
    bonds = []
    for i, j, kind in found_bonds:
        x, y = found[i][1] - 1, found[j][1] - 1
        bonds.append([min(x, y), max(x, y), ORDER.get(kind, 1.0)])
    return atoms, chiral, sorted(bonds)


def parse_reaction(rxn_smiles, with_product_chirality=False):
    if rxn_smiles.count('>>') != 1:
        raise ValueError('expected a reaction SMILES of the form reactants>>products')
    r, p = rxn_smiles.split('>>')
    atoms, chiral, rb = parse_mapped(r)
    patoms, pchiral, pb = parse_mapped(p)
    if patoms != atoms:
        raise ValueError('reactant and product sides must contain the same mapped atoms')
    if with_product_chirality:
        return atoms, chiral, rb, pb, pchiral
    return atoms, chiral, rb, pb


def charge_and_multiplicity(rxn_smiles):
    """Total charge and spin multiplicity from the reactant side (formal charges, radical electrons; singlet if none)."""
    atoms, _ = _atoms_and_bonds(rxn_smiles.split('>>')[0], radicals=True)
    return int(sum(a[3] for a in atoms)), int(sum(a[4] for a in atoms) + 1)
