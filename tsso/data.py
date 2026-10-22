"""Load the exported GoFlow/RDB7 reactions as plain Python objects with NumPy arrays."""
import gzip, json
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np

DEFAULT = Path(__file__).resolve().parent.parent / 'data' / 'goflow8_rdb7_test.json.gz'
_ELEMENTS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr '
             'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb '
             'Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn').split()
SYMBOL = {z: s for z, s in enumerate(_ELEMENTS, start=1)}   # atomic number -> element symbol
ATOMIC_NUMBER = {s: z for z, s in SYMBOL.items()}


@dataclass
class Reaction:
    id: int
    smiles: str
    barrier_kcal_mol: float
    atoms: np.ndarray            # atomic numbers, index i = atom map number i+1
    reactant_xyz: np.ndarray     # (N, 3) Angstrom
    ts_ref_xyz: np.ndarray       # (N, 3) reference DFT transition state
    product_xyz: np.ndarray      # (N, 3)
    reactant_bonds: dict         # {(i, j): bond order}
    product_bonds: dict
    reactant_chiral_tags: np.ndarray  # RDKit ChiralType ints (1 = CW, 2 = CCW)
    product_chiral_tags: np.ndarray
    preds: np.ndarray            # (models, N, 3) GoFlow predictions
    _cache: dict = field(default_factory=dict, repr=False)

    @property
    def n(self):
        return len(self.atoms)

    @property
    def symbols(self):
        return [SYMBOL.get(int(z), str(z)) for z in self.atoms]

    @property
    def changed_bonds(self):
        """Bonds that form or break (present in only one of reactant/product)."""
        return sorted(set(self.reactant_bonds) ^ set(self.product_bonds))

    @property
    def reaction_centre(self):
        return {a for b in self.changed_bonds for a in b}

    def neighbours(self, which='reactant'):
        bonds = self.reactant_bonds if which == 'reactant' else self.product_bonds
        nb = {i: [] for i in range(self.n)}
        for i, j in bonds:
            nb[i].append(j); nb[j].append(i)
        return {i: sorted(v) for i, v in nb.items()}


def load_reactions(path=DEFAULT):
    with gzip.open(path, 'rt') as fh:
        raw = json.load(fh)
    out = []
    for r in raw['reactions']:
        out.append(Reaction(
            id=r['id'], smiles=r['smiles'], barrier_kcal_mol=r['barrier_kcal_mol'],
            atoms=np.array(r['atoms']),
            reactant_xyz=np.array(r['reactant_xyz']), ts_ref_xyz=np.array(r['ts_ref_xyz']),
            product_xyz=np.array(r['product_xyz']),
            reactant_bonds={(i, j): o for i, j, o in r['reactant_bonds']},
            product_bonds={(i, j): o for i, j, o in r['product_bonds']},
            reactant_chiral_tags=np.array(r['reactant_chiral_tags']),
            product_chiral_tags=np.array(r['product_chiral_tags']),
            preds=np.array(r['preds']),
        ))
    return out
