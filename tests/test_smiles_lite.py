"""The pure-Python SMILES reader (used when RDKit is missing, e.g. in the browser) must agree with RDKit."""
import gzip, json, random
from pathlib import Path
import pytest

Chem = pytest.importorskip('rdkit.Chem')
from tsso import chem, smiles_lite

ROOT = Path(__file__).resolve().parent.parent


def lite(fn, *args):
    saved = chem.Chem
    chem.Chem = None
    try:
        return fn(*args)
    finally:
        chem.Chem = saved


def outcome(fn, *args):
    try:
        return 'ok', fn(*args)
    except ValueError as e:
        return 'error', str(e)


def same(fn, *args):
    a, b = outcome(fn, *args), outcome(lambda *x: lite(fn, *x), *args)
    assert a == b, (fn.__name__, args, a, b)


def mapped(smiles, seed):
    """A fully atom-mapped, explicit-hydrogen SMILES of a molecule, with a random atom order."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    order = list(range(m.GetNumAtoms()))
    random.Random(seed).shuffle(order)
    for k, a in zip(order, m.GetAtoms()):
        a.SetAtomMapNum(k + 1)
    return Chem.MolToSmiles(m, doRandom=True, canonical=False, kekuleSmiles=seed % 2 == 0)


TRICKY = ['[CH3]', '[CH2]', '[OH]', '[O]', '[NH2]', '[NH]', 'C[O-]', 'C[NH3+]', '[Cl-]', 'CBr', '[O-][N+](=O)C',
          'c1ccccc1O', 'c1ccncc1', 'c1cc[nH]c1', 'C[C@H](N)C(=O)O', 'C[C@@H](O)CC', 'F[C@](Cl)(Br)I', 'C1CC2CCC1C2',
          'C/C=C/C', 'C#N', 'C=[N+]=[N-]', 'CS(=O)(=O)O', 'CP(=O)(O)O', 'C[S+](C)C', 'O=[Se]=O', '[2H]C([2H])[2H]',
          'C[Si](C)(C)C', 'B(O)(O)O', '[BH4-]', 'C[C]C', 'C=[C]', 'O=C[O]', '[N]=O', 'CI(=O)=O', '[Na+].[Cl-]',
          'C1CCCCC1C1CCCCC1C1CCCCC1C1CC1C1CC1', 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O']


def test_lite_reader_matches_rdkit_on_tricky_molecules():
    for smi in TRICKY:
        for seed in range(4):
            s = mapped(smi, seed)
            for fn in (chem.parse_mapped, chem.smiles_atom_order):
                same(fn, s)
            same(chem.charge_and_multiplicity, s + '>>' + s)


def test_lite_reader_matches_rdkit_on_every_rdb7_reaction():
    seen = set()
    for f in ('goflow8_rdb7_test.json.gz', 'rdb7_train_reference.json.gz'):
        for r in json.load(gzip.open(ROOT / 'data' / f, 'rt'))['reactions']:
            seen.add(r['smiles'])
    for f in (ROOT / 'tests' / 'fixtures').rglob('reaction.smi'):
        seen.add(f.read_text().strip())
    assert len(seen) > 10_000
    for s in sorted(seen):
        same(chem.parse_reaction, s, True)
        same(chem.charge_and_multiplicity, s)
        same(chem.smiles_atom_order, s.split('>>')[0])


def test_lite_reader_rejects_what_rdkit_rejects():
    bad = ['C1CC', 'C(C', 'CC)', '[C:1', 'C==C', '[Xx:1]', '[C:1][C:2', '[C:1]%1', 'Q', '[C:1]((C))',
           '[C:1][H:2]', '[C:1][H]', '[C:2][H:3]', '[C:1][H:1]', 'CCO', '[C:1]1[C:2]1', '[C:1].[O:2]']
    for s in bad:
        same(chem.parse_mapped, s)
