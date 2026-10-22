"""Export GoFlow's released RDB7 test-set predictions into one small, readable JSON file.

Source: github.com/heid-lab/goflow (MIT licence), logs/train_test_error_bar —
8 GoFlow models trained with seeds 1-8 on the RDB7 random split, each predicting
the 1,194 test reactions (25 samples per reaction, median-selected prediction saved).
Reference geometries and barriers come from RDB7 (Spiekermann et al., Sci. Data 2022).

Usage: python scripts/export_data.py /path/to/goflow data/goflow8_rdb7_test.json.gz [run-glob] [description]
"""
import csv, glob, gzip, json, sys
import numpy as np
from load_samples import load
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from tsso.chem import parse_reaction

root, out = sys.argv[1], sys.argv[2]
pattern = sys.argv[3] if len(sys.argv) > 3 else 'train_test_error_bar/multiruns/*/[0-7]'
runs = sorted(glob.glob(f'{root}/logs/{pattern}/test_samples/samples_all.pkl'))
assert runs, pattern
if len(sys.argv) <= 3:
    assert len(runs) == 8, runs

barrier = {}
with open(f'{root}/data/RDB7/raw_data/rdb7_full.csv') as fh:
    for row in csv.DictReader(fh):
        barrier[int(row['rxn'])] = float(row['dE0'])

ORDER = {'SINGLE': 1.0, 'DOUBLE': 2.0, 'TRIPLE': 3.0, 'AROMATIC': 1.5}

def bonds(mol):
    out = []
    for b in mol.GetBonds():
        i = b.GetBeginAtom().GetAtomMapNum() - 1
        j = b.GetEndAtom().GetAtomMapNum() - 1
        out.append([min(i, j), max(i, j), ORDER.get(str(b.GetBondType()), 1.0)])
    return sorted(out)

rxns = {}
for k, f in enumerate(runs):
    for rec in load(f):
        mp = rec._state['_store']._state['_mapping']
        rid = int(mp['rxn_index'][0])
        if rid not in rxns:
            rmol, pmol = mp['rdmol']
            rxns[rid] = {
                'id': rid,
                'smiles': mp['smiles'],
                'barrier_kcal_mol': barrier.get(rid),
                'atoms': [int(z) for z in mp['atom_type']],
                'reactant_xyz': np.round(mp['pos'][:, 0, :], 5).tolist(),
                'ts_ref_xyz': np.round(mp['pos'][:, 1, :], 5).tolist(),
                'product_xyz': np.round(mp['pos'][:, 2, :], 5).tolist(),
                'reactant_bonds': parse_reaction(mp['smiles'])[2],   # bond orders from the mapped SMILES
                'product_bonds': parse_reaction(mp['smiles'])[3],
                'reactant_chiral_tags': [int(c) for c in mp['r_chiraltype']],
                'product_chiral_tags': [int(c) for c in mp['p_chiraltype']],
                'preds': [None] * len(runs),
            }
        rxns[rid]['preds'][k] = np.round(mp['pos_gen'][-1], 5).tolist()

DEFAULT_DESC = ('GoFlow (Galustian et al., Digital Discovery 2025) predictions for the RDB7 random-split test set, '
                'from the 8 seed models (n_atom_basis 256) released in github.com/heid-lab/goflow (MIT). '
                'Index i in every array = atom map number i+1.')
data = {
    'description': sys.argv[4] if len(sys.argv) > 4 else DEFAULT_DESC,   # pass a description for other runs
    'model_seeds': list(range(1, 9)) if len(runs) == 8 else None,
    'runs': [r.split('/logs/')[1] for r in runs],
    'reactions': [rxns[r] for r in sorted(rxns)],
}
assert all(None not in r['preds'] for r in data['reactions'])
with gzip.open(out, 'wt') as fh:
    json.dump(data, fh)
print(f'wrote {len(data["reactions"])} reactions to {out}')
