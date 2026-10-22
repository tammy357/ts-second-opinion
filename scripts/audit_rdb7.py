"""Data-quality audit: RDB7 reactions whose reactant or product geometry contradicts the atom-mapped SMILES.

Usage: python scripts/audit_rdb7.py  -> results/rdb7_audit.csv   (same as `tsso audit --json`)
"""
import csv, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from tsso.service import audit_rdb7

res = audit_rdb7(limit=100000)
with open(ROOT / 'results' / 'rdb7_audit.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['split', 'reaction_id', 'structure', 'atoms', 'distance_A', 'distance_over_covalent_radii', 'problem'])
    for row in res['examples']:
        for i in row['issues']:
            w.writerow([row['split'], row['reaction_id'], i['side'], i['atoms'], i['distance_A'], i['ratio'], i['kind']])
print(f"{res['reactions_flagged']} of {res['reactions_checked']} reactions flagged -> results/rdb7_audit.csv")
