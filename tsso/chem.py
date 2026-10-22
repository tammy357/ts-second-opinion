"""Parse atom-mapped reaction SMILES into atoms, bonds and stereo tags indexed by atom map number - 1.

Every atom, hydrogens included, needs a map number, and the numbers must run 1..N on both sides.
"""
from rdkit import Chem

ORDER = {'SINGLE': 1.0, 'DOUBLE': 2.0, 'TRIPLE': 3.0, 'AROMATIC': 1.5}
_ps = Chem.SmilesParserParams()
_ps.removeHs = False
_ps.sanitize = False


def _mol(smi):
    m = Chem.MolFromSmiles(smi, _ps)
    if m is None:
        raise ValueError(f'could not parse SMILES: {smi}')
    maps = sorted(a.GetAtomMapNum() for a in m.GetAtoms())
    if maps != list(range(1, m.GetNumAtoms() + 1)):
        raise ValueError('every atom, including each hydrogen, needs a unique atom-map number 1..N '
                         '(write hydrogens explicitly, e.g. [C:1]([H:4])([H:5])...)')
    return m


def smiles_atom_order(smi):
    """Map numbers in the order the atoms appear in the SMILES (for xyz files written in that order)."""
    return [a.GetAtomMapNum() for a in _mol(smi).GetAtoms()]


def parse_mapped(smi):
    m = _mol(smi)
    n = m.GetNumAtoms()
    atoms, chiral = [0] * n, [0] * n
    for a in m.GetAtoms():
        k = a.GetAtomMapNum() - 1
        atoms[k] = a.GetAtomicNum()
        chiral[k] = int(a.GetChiralTag())
    bonds = []
    for b in m.GetBonds():
        x, y = b.GetBeginAtom().GetAtomMapNum() - 1, b.GetEndAtom().GetAtomMapNum() - 1
        bonds.append([min(x, y), max(x, y), ORDER.get(str(b.GetBondType()), 1.0)])
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
    m = _mol(rxn_smiles.split('>>')[0])
    m.UpdatePropertyCache(strict=False)
    Chem.AssignRadicals(m)
    charge = sum(a.GetFormalCharge() for a in m.GetAtoms())
    unpaired = sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())
    return int(charge), int(unpaired + 1)
