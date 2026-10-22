"""Triage TS guesses for one reaction from Python:  python examples/quickstart.py"""
from pathlib import Path
from tsso import second_opinion

d = str(Path(__file__).resolve().parent / 'rdb7_3509') + '/'
smiles = open(d + 'reaction.smi').read().strip()
result = second_opinion(smiles, d + 'reactant.xyz', d + 'product.xyz', d + 'guesses.xyz',
                        score_reference=d + 'reference_ts.xyz')      # the reference is optional, for scoring only

print(result['verdict'], '-', result['meaning'])
for line in result['reasons']:
    print(' -', line)
print('vs reference:', result['score_against_reference'])
open('recommended_ts.xyz', 'w').write(result['recommended_ts_xyz'] + '\n')
