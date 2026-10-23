"""Summarise the new DFT calculations in results/dft/*.json -> results/dft_summary.json (shown in the app and README).

Usage: python scripts/dft_summary.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    rows = [json.loads(p.read_text()) for p in sorted((ROOT / 'results' / 'dft').glob('*.json'))]
    mirror = [r for r in rows if 'unrepaired' in r and 'error' not in r['unrepaired'] and 'error' not in r.get('repaired', {'error': 1})]
    guesses = [r for r in rows if 'recommended' in r and 'error' not in r['recommended']]
    out = {'level': 'wB97X/def2-SVP (PySCF), Sella saddle search, analytic frequencies', 'mirror': [], 'by_verdict': {}, 'reactions': []}
    for r in mirror:
        u, p = r['unrepaired'], r['repaired']
        out['mirror'].append({'reaction_id': r['reaction_id'],
                              'both_first_order_saddles': u['n_imaginary'] == 1 and p['n_imaginary'] == 1,
                              'energy_difference_kcal': round(abs(u['energy_kcal'] - p['energy_kcal']), 4),
                              'imaginary_cm1': [round(u['imaginary_cm1']), round(p['imaginary_cm1'])],
                              'unrepaired_handedness_correct': u['handedness_correct'],
                              'repaired_handedness_correct': p['handedness_correct'],
                              'repaired_found_the_reference_ts': p['found_the_reference_ts']})
    for r in guesses:
        g = r['recommended']
        v = out['by_verdict'].setdefault(r['verdict'], {'n': 0, 'found_the_reference_ts': 0, 'first_order_saddle': 0})
        v['n'] += 1; v['found_the_reference_ts'] += g['found_the_reference_ts']
        v['first_order_saddle'] += g['converged'] and g['n_imaginary'] == 1
        out['reactions'].append({'reaction_id': r['reaction_id'], 'verdict': r['verdict'], 'converged': g['converged'],
                                 'n_imaginary': g['n_imaginary'], 'dmae_to_reference_ts_A': g['dmae_to_reference_ts_A'],
                                 'handedness_correct': g['handedness_correct'], 'found_the_reference_ts': g['found_the_reference_ts']})
    (ROOT / 'results' / 'dft_summary.json').write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in ('mirror', 'by_verdict')}, indent=1))


if __name__ == '__main__':
    main()
