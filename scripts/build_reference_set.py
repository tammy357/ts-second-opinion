"""Build the gate set: reference DFT transition states for the RDB7 *training* reactions.

Checks are gated on these 9,540 training reactions, never on the 1,194 test reactions used
for evaluation, so thresholds are not tuned on the test set.

Usage: python scripts/build_reference_set.py /path/to/goflow data/rdb7_train_reference.json.gz
"""
import csv, gzip, json, pickle, sys
import numpy as np
import sys as _s, pathlib as _pl
_s.path.insert(0, str(_pl.Path(__file__).resolve().parent.parent))
from tsso.chem import parse_mapped as parse

root, out = sys.argv[1], sys.argv[2]
raw = f'{root}/data/RDB7/raw_data'


class NoGlobals(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError('unexpected global in split file')


split = NoGlobals(open(f'{root}/data/RDB7/splits/random_split.pkl', 'rb')).load()

rows = list(csv.DictReader(open(f'{raw}/rdb7_full.csv')))
blocks, lines = [], open(f'{raw}/rdb7_full.xyz').read().splitlines()
i = 0
while i < len(lines):
    if not lines[i].strip():
        i += 1; continue
    n = int(lines[i]); xyz = []
    for ln in lines[i + 2:i + 2 + n]:
        p = ln.split(); xyz.append([float(p[1]), float(p[2]), float(p[3])])
    blocks.append(xyz); i += 2 + n
assert len(blocks) == 3 * len(rows), (len(blocks), len(rows))
rowpos = {int(r['rxn']): k for k, r in enumerate(rows)}  # reaction ids are the csv 'rxn' column, not row numbers

def build(ids):
    out = []
    for rid in ids:
        k = rowpos[rid]
        r_smi, p_smi = rows[k]['smiles'].split('>>')
        atoms, chiral, rb = parse(r_smi)
        _, pchiral, pb = parse(p_smi)
        R, T, P = blocks[3 * k:3 * k + 3]
        assert len(R) == len(atoms)
        out.append({'id': rid, 'smiles': rows[k]['smiles'], 'barrier_kcal_mol': float(rows[k]['dE0']),
                    'atoms': atoms, 'reactant_xyz': R, 'ts_ref_xyz': T, 'product_xyz': P,
                    'reactant_bonds': rb, 'product_bonds': pb, 'reactant_chiral_tags': chiral,
                    'product_chiral_tags': pchiral, 'preds': []})
    return out


data = {'description': 'RDB7 (Spiekermann et al., Sci. Data 2022) training-split reactions (GoFlow random split), '
                       'used only to gate rule checks. Index i = atom map number i+1.',
        'reactions': build(split['train'])}
with gzip.open(out, 'wt') as fh:
    json.dump(data, fh)
print(f'wrote {len(data["reactions"])} training reactions to {out}')

# consistency check against the test export, if present
try:
    test = json.load(gzip.open(out.replace('rdb7_train_reference', 'goflow8_rdb7_test'), 'rt'))['reactions']
    mine = {r['id']: r for r in build([t['id'] for t in test])}
    same_b = sum(mine[t['id']]['reactant_bonds'] == t['reactant_bonds'] and mine[t['id']]['product_bonds'] == t['product_bonds'] for t in test)
    same_c = sum(mine[t['id']]['reactant_chiral_tags'] == t['reactant_chiral_tags'] and mine[t['id']]['product_chiral_tags'] == t['product_chiral_tags'] for t in test)
    same_x = sum(np.allclose(mine[t['id']]['ts_ref_xyz'], t['ts_ref_xyz'], atol=1e-3) for t in test)
    print(f'consistency on {len(test)} test reactions: bonds {same_b}, chiral tags {same_c}, TS coordinates {same_x}')
except FileNotFoundError:
    pass
