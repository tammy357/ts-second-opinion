"""How closely does geometry-only input (no SMILES) reproduce the SMILES route? -> results/geometry_input_validation.json

For each of the 1,194 RDB7 test reactions, the atom-mapped SMILES is inferred from the reactant and product
geometries alone (tsso.infer) and compared with the dataset's own SMILES: bonds, bond orders, stereocentres, and the
triage outcome on the same 8 GoFlow guesses (verdict, chosen guess, mirror repair). Reactions flagged by the data audit
(whose geometry contradicts their own SMILES) are reported separately, since there the two routes should differ.

    python scripts/validate_geometry_input.py
"""
import json, sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from tsso.audit import audit_reaction            # noqa: E402
from tsso.chem import parse_reaction              # noqa: E402
from tsso.data import load_reactions              # noqa: E402
from tsso.geometry import dmae                    # noqa: E402
from tsso.metrics import handedness_correct       # noqa: E402
from tsso.infer import reaction_smiles_from_geometry  # noqa: E402
from tsso.service import build_reaction, amber_cutoff  # noqa: E402
from tsso.triage import load_checks, triage      # noqa: E402


def centres(tags):
    return set(np.nonzero(np.isin(np.asarray(tags), (1, 2)))[0].tolist())


def main():
    reactions, checks, cut = load_reactions(), load_checks(), amber_cutoff()
    c, differs, failed = Counter(), [], []
    final = {'smiles_route': Counter(), 'geometry_route': Counter()}   # evaluation only: against the reference TS
    for r in reactions:
        flagged = bool(audit_reaction(r))
        try:
            info = reaction_smiles_from_geometry(r.symbols, r.reactant_xyz, r.product_xyz)
            atoms, chiral, rb, pb, pchiral = parse_reaction(info['smiles'], with_product_chirality=True)
            g, err = build_reaction(info['smiles'], r.reactant_xyz, r.product_xyz, r.preds)
            if err:
                raise ValueError(err)
        except Exception as e:                    # report, never hide
            failed.append({'reaction_id': r.id, 'error': str(e)[:200]}); continue
        same = {
            'reactant_bonds': {(a, b) for a, b, _ in rb} == set(r.reactant_bonds),
            'product_bonds': {(a, b) for a, b, _ in pb} == set(r.product_bonds),
            'bond_orders': ({(a, b): o for a, b, o in rb} == r.reactant_bonds and {(a, b): o for a, b, o in pb} == r.product_bonds),
            'reactant_stereocentres': centres(chiral) == centres(r.reactant_chiral_tags),
            'product_stereocentres': centres(pchiral) == centres(r.product_chiral_tags),
        }
        t_smiles, t_geom = triage(r, checks, cut), triage(g, checks, cut)
        for route, t in (('smiles_route', t_smiles), ('geometry_route', t_geom)):
            geo_bad = dmae(t.xyz, r.ts_ref_xyz) > 0.2
            hand_bad = not handedness_correct(r, t.xyz, r.ts_ref_xyz)
            final[route]['n'] += 1; final[route]['geometry_wrong'] += geo_bad; final[route]['handedness_wrong'] += hand_bad
            final[route]['green'] += t.verdict == 'green'
            final[route]['green_but_wrong'] += t.verdict == 'green' and (geo_bad or hand_bad)
        same['verdict'] = t_smiles.verdict == t_geom.verdict
        same['recommended_guess'] = (t_smiles.chosen_model == t_geom.chosen_model and t_smiles.mirrored == t_geom.mirrored)
        same['bond_orders_assigned'] = info['bond_orders']
        for k, v in same.items():
            c[k] += bool(v)
        c['audit_flagged'] += flagged
        if not (same['verdict'] and same['recommended_guess']):
            differs.append({'reaction_id': r.id, 'audit_flagged': flagged,
                            'verdict_smiles_route': t_smiles.verdict, 'verdict_geometry_route': t_geom.verdict,
                            'same_bonds': same['reactant_bonds'] and same['product_bonds'],
                            'same_stereocentres': same['reactant_stereocentres'] and same['product_stereocentres']})
    n = len(reactions) - len(failed)
    out = {'reactions': len(reactions), 'inferred': n, 'inference_failed': failed,
           'agreement': {k: round(v / n, 4) for k, v in c.items() if k != 'audit_flagged'},
           'counts': dict(c),
           'final_answers_vs_reference_evaluation_only': {
               route: {'geometry_wrong': round(f['geometry_wrong'] / f['n'], 4),
                       'handedness_wrong': round(f['handedness_wrong'] / f['n'], 4),
                       'trusted_green': round(f['green'] / f['n'], 4),
                       'wrong_among_trusted': round(f['green_but_wrong'] / max(f['green'], 1), 4)}
               for route, f in final.items()},
           'outcome_differs': {'total': len(differs),
                               'of_which_audit_flagged': sum(d['audit_flagged'] for d in differs),
                               'reactions': differs},
           'meaning': ('Share of reactions where the SMILES inferred from the geometries gives the same result as the '
                       "dataset's own SMILES. audit_flagged reactions have a geometry that contradicts their SMILES, "
                       'so there the geometry route follows the structures and the SMILES route follows the SMILES.')}
    (ROOT / 'results' / 'geometry_input_validation.json').write_text(json.dumps(out, indent=1) + '\n')
    print(f"{n} of {len(reactions)} reactions inferred; same verdict {100 * out['agreement']['verdict']:.1f}%, "
          f"same recommended guess {100 * out['agreement']['recommended_guess']:.1f}%")
    for k, v in out['agreement'].items():
        print(f'  {k:24s} {100 * v:5.1f}%')
    print(f"  outcome differs in {len(differs)} reactions ({out['outcome_differs']['of_which_audit_flagged']} audit-flagged)")
    for route, f in out['final_answers_vs_reference_evaluation_only'].items():
        print(f"  {route}: final answers wrong-handed {100 * f['handedness_wrong']:.1f}%, geometry wrong "
              f"{100 * f['geometry_wrong']:.1f}%, trusted {100 * f['trusted_green']:.1f}%, wrong among trusted "
              f"{100 * f['wrong_among_trusted']:.1f}% (vs the reference TS; evaluation only)")


if __name__ == '__main__':
    main()
