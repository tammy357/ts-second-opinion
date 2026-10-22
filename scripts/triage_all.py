"""Triage every test reaction, write a table, ORCA inputs for flagged reactions, and verdict-level accuracy.

The reference TS is used only afterwards, to score the verdicts.
Usage: python scripts/triage_all.py [budget, default 0.2]
"""
import csv, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.geometry import dmae
from tsso.triage import load_checks, triage, amber_cutoff_for_budget, write_orca

budget = float(sys.argv[1]) if len(sys.argv) > 1 else 0.2
R = load_reactions(); checks = load_checks()
cut = amber_cutoff_for_budget(R, checks, budget)


def stereo_wrong(r, x):
    nbr, nbp = r.neighbours('reactant'), r.neighbours('product')
    cents = [(i, nbr) for i in range(r.n) if r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(nbr[i]) >= 3]
    cents += [(i, nbp) for i in range(r.n) if r.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    for i, nb in cents:
        v1 = r.ts_ref_xyz[nb[i][:3]] - r.ts_ref_xyz[i]; v2 = x[nb[i][:3]] - x[i]
        if np.sign(np.dot(v1[0], np.cross(v1[1], v1[2]))) != np.sign(np.dot(v2[0], np.cross(v2[1], v2[2]))):
            return True
    return False


rows, stats = [], {}
out = Path('results'); orca = out / 'orca_inputs'
for r in R:
    t = triage(r, checks, cut)
    e = dmae(t.xyz, r.ts_ref_xyz); geo = e > 0.2; st = stereo_wrong(r, t.xyz)
    rows.append({'reaction_id': r.id, 'smiles': r.smiles, 'verdict': t.verdict, 'chosen_model_seed': t.chosen_model + 1,
                 'mirrored': t.mirrored, 'disagreement_A': round(t.disagreement, 4),
                 'failed_rules': '; '.join(n for n, _ in t.failed_rules),
                 'models_rejected': len(t.rejected_models),
                 'reasons': ' '.join(t.reasons),
                 'dmae_vs_reference_A': round(e, 4), 'geometry_bad': geo, 'stereo_wrong': st})
    s = stats.setdefault(t.verdict, {'n': 0, 'geometry_bad': 0, 'stereo_wrong': 0, 'any_bad': 0, 'dmae': []})
    s['n'] += 1; s['geometry_bad'] += geo; s['stereo_wrong'] += st; s['any_bad'] += (geo or st); s['dmae'].append(e)
    if t.verdict != 'green':
        write_orca(orca, r, t)
with open(out / 'triage.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
summary = {'budget': budget, 'amber_cutoff_A': cut, 'n': len(R), 'verdicts': {}}
for v in ('green', 'amber', 'red'):
    s = stats.get(v, {'n': 0, 'geometry_bad': 0, 'stereo_wrong': 0, 'any_bad': 0, 'dmae': [0]})
    summary['verdicts'][v] = {'reactions': s['n'], 'share': s['n'] / len(R),
                              'any_bad_rate': s['any_bad'] / max(s['n'], 1), 'geometry_bad_rate': s['geometry_bad'] / max(s['n'], 1),
                              'stereo_wrong_rate': s['stereo_wrong'] / max(s['n'], 1), 'mean_dmae_A': float(np.mean(s['dmae']))}
tot_bad = sum(stats[v]['any_bad'] for v in stats)
summary['bad_guesses_in_flagged'] = sum(stats[v]['any_bad'] for v in stats if v != 'green') / max(tot_bad, 1)
summary['mirror_repairs'] = sum(r['mirrored'] for r in rows)
json.dump(summary, open(out / 'triage_summary.json', 'w'), indent=2)
print(json.dumps(summary, indent=2))
