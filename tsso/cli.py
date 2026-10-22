"""tsso: command-line interface for TS Second Opinion.

  tsso triage 3509 [--score]        triage one RDB7 test reaction
  tsso list --verdict red            list reactions by verdict
  tsso rules                         list rules and their gate results
  tsso gate "rule text" check.py     test a rule check against 9,540 real TSs
  tsso author "rule text"            let Claude write the check and iterate until the gate accepts it
  tsso check --smiles "<mapped rxn>" --reactant r.xyz --product p.xyz guesses.xyz
                                     triage your own reaction from any generator (multi-frame xyz is fine)
  tsso evaluate preds.xyz refs.xyz reactions.smi
                                     stereo-aware evaluation of any TS generator's test-set predictions
  tsso certify --target 0.05         certified trust: which guesses to trust at a guaranteed average error rate
  tsso orca 3509                     write an ORCA TS-optimisation input
  tsso summary                       headline results
  tsso mcp                           run as an MCP server for agents (Claude Code, Claude Desktop)
Add --json to any command for machine-readable output.
"""
import argparse, json, sys
from pathlib import Path

ICON = {'green': '✓ TRUST', 'amber': '! CHECK WITH DFT', 'red': '✕ RE-GENERATE'}


def _print(obj, as_json):
    if as_json:
        print(json.dumps(obj, indent=2)); return
    if 'error' in obj:
        print('error:', obj['error']); sys.exit(1)
    if 'verdict' in obj and 'reasons' in obj:
        who = 'Your reaction' if obj['reaction_id'] == 'yours' else f"Reaction {obj['reaction_id']}"
        print(f"{who}: {ICON[obj['verdict']]}  ({obj['meaning']})")
        print(f"  recommended guess: model {obj['chosen_guess']}" + ('  (mirror-repaired)' if obj['mirror_repaired'] else ''))
        for line in obj['reasons']:
            print('  -', line)
        for line in obj.get('input_warnings', []):
            print('  ! data check:', line)
        if 'score_against_reference' in obj:
            s = obj['score_against_reference']
            print(f"  vs reference (evaluation only): D-MAE {s['dmae_A']} A, stereo {'correct' if s['stereo_correct'] else 'WRONG'}")
        return
    print(json.dumps(obj, indent=2))


def main(argv=None):
    ap = argparse.ArgumentParser(prog='tsso', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--json', action='store_true')
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('triage'); p.add_argument('reaction_id', type=int); p.add_argument('--score', action='store_true')
    p.add_argument('--report', help='write a single-file HTML report card here')
    p = sub.add_parser('list'); p.add_argument('--verdict'); p.add_argument('--mirror-repaired', action='store_true')
    p.add_argument('--limit', type=int, default=20)
    sub.add_parser('rules')
    p = sub.add_parser('gate'); p.add_argument('rule'); p.add_argument('source_file')
    p = sub.add_parser('author'); p.add_argument('rule'); p.add_argument('--name'); p.add_argument('--backend', default='auto')
    p = sub.add_parser('check', help='triage your own reaction')
    p.add_argument('guesses', nargs='*', help='TS guess .xyz files (one or many frames each), or one .json spec')
    p.add_argument('--smiles', help='atom-mapped reaction SMILES, every atom mapped (or a file containing it)')
    p.add_argument('--reactant'); p.add_argument('--product')
    p.add_argument('--order', choices=['map', 'smiles'], default='map',
                   help="atoms listed by map number (default) or in SMILES order")
    p.add_argument('--reference', help='reference TS .xyz, to score the result (evaluation only)')
    p.add_argument('--write-xyz', help='write the recommended TS guess here')
    p.add_argument('--orca', help='write an ORCA OptTS+Freq input here (needs --write-xyz)')
    p.add_argument('--report', help='write a single-file HTML report card here')
    p = sub.add_parser('evaluate', help='stereo-aware evaluation of any TS generator')
    p.add_argument('predictions'); p.add_argument('references'); p.add_argument('reactions')
    p.add_argument('--order', choices=['map', 'smiles'], default='map'); p.add_argument('--threshold', type=float, default=0.2)
    p = sub.add_parser('certify'); p.add_argument('--target', type=float, default=0.05)
    p.add_argument('--count', choices=['geometry', 'geometry_or_handedness'], default='geometry')
    p.add_argument('--seed', type=int, default=0)
    p = sub.add_parser('orca'); p.add_argument('reaction_id', type=int)
    sub.add_parser('summary')
    sub.add_parser('audit', help='where RDB7 SMILES and geometries disagree')
    sub.add_parser('mcp')
    a = ap.parse_args(argv)
    if a.cmd == 'mcp':
        from .mcp_server import main as serve
        return serve()
    from . import service
    if a.cmd == 'triage':
        res = service.triage_reaction(a.reaction_id, a.score, a.report)
        _print(res, a.json)
        if not a.json and res.get('report'):
            print(f"  report card: {res['report']}")
    elif a.cmd == 'list':
        res = service.list_reactions(a.verdict, True if a.mirror_repaired else None, a.limit)
        if a.json:
            _print(res, True)
        else:
            print(f"{res['total']} reactions")
            for x in res['reactions']:
                print(f"  {x['reaction_id']:>6}  {ICON[x['verdict']]:<17} {'mirror-repaired' if x['mirror_repaired'] else ''}")
    elif a.cmd == 'rules':
        res = service.list_rules()
        if a.json:
            _print(res, True)
        else:
            for r in res['rules']:
                print(f"{'●' if r['used_in_triage'] else '○'} {r['name']}: {r['rule']}\n    {r['gate']}")
    elif a.cmd == 'gate':
        _print(service.gate_check(a.rule, Path(a.source_file).read_text()), True)
    elif a.cmd == 'author':
        from .author import author
        out = author(a.rule, a.name, a.backend)
        print('ACCEPTED -> ' + out.get('path', '') if out['accepted'] else 'REJECTED after all attempts')
    elif a.cmd == 'check':
        if len(a.guesses) == 1 and a.guesses[0].endswith('.json'):
            spec = json.loads(Path(a.guesses[0]).read_text())
            res = service.triage_structures(spec['mapped_reaction_smiles'], spec['reactant_xyz'], spec['product_xyz'],
                                            spec['guesses_xyz'], spec.get('order', 'map'))
        else:
            if not (a.smiles and a.reactant and a.product and a.guesses):
                ap.error('check needs --smiles, --reactant, --product and at least one guesses .xyz file')
            smi = Path(a.smiles).read_text().split()[0] if Path(a.smiles).is_file() else a.smiles
            res = service.second_opinion(smi, Path(a.reactant), Path(a.product), [Path(g) for g in a.guesses],
                                         a.order, Path(a.reference) if a.reference else None, a.report)
        if 'error' not in res and a.write_xyz:
            Path(a.write_xyz).write_text(res['recommended_ts_xyz'] + '\n')
            if a.orca:
                from .triage import ORCA_TEMPLATE
                Path(a.orca).write_text(ORCA_TEMPLATE.format(rid='your reaction', source=f"verdict {res['verdict']}",
                                                             xyzname=Path(a.write_xyz).name,
                                                             charge=res['charge'], mult=res['multiplicity']))
        _print(res, a.json)
        if not a.json and 'error' not in res and a.write_xyz:
            print(f"  wrote {a.write_xyz}" + (f" and {a.orca}" if a.orca else ''))
        if not a.json and res.get('report'):
            print(f"  report card: {res['report']}")
    elif a.cmd == 'evaluate':
        res = service.evaluate_files(a.predictions, a.references, a.reactions, a.order, a.threshold)
        if a.json or 'error' in res:
            _print(res, True)
        else:
            t = res['threshold_A']
            print(f"{res['n']} reactions")
            print(f"  mean D-MAE                         {res['mean_dmae_A']:.3f} A")
            print(f"  mean RMSD, mirror allowed          {res['mean_rmsd_mirror_allowed_A']:.3f} A   (what ignore_chirality=True reports)")
            print(f"  mean RMSD, proper rotations only   {res['mean_rmsd_A']:.3f} A")
            print(f"  success (D-MAE <= {t} A)            {100 * res['success_dmae_only']:.1f}%")
            print(f"  correct handedness                 {100 * res['handedness_correct']:.1f}%")
            print(f"  success with correct handedness    {100 * res['success_dmae_and_handedness']:.1f}%")
            print(f"  wrong-handed but counted a success {100 * res['wrong_handed_counted_as_success']:.1f}%")
    elif a.cmd == 'certify':
        res = service.certify_trust(a.target, a.count, a.seed)
        if a.json or 'error' in res:
            _print(res, True)
        else:
            rep = res.get('over_2000_random_splits', {})
            print(f"Certified trust at a {100 * res['target_error_rate']:.0f}% target ({res['counting'].replace('_', ' ')} errors)")
            print(f"  calibrated on {res['calibration_reactions']} reactions; {res['trusted']} of {res['new_reactions']} new "
                  f"reactions trusted ({100 * res['share_trusted']:.0f}%)")
            print(f"  actual error rate among them in this one split (evaluation only): {100 * res['actual_error_rate_among_trusted']:.1f}%")
            if rep:
                print(f"  over 2,000 random splits: {100 * rep['mean_share_trusted']:.0f}% trusted, "
                      f"{100 * rep['mean_error_among_trusted']:.1f}% of them wrong on average")
            print('  ' + res['note'])
    elif a.cmd == 'orca':
        res = service.write_orca_input(a.reaction_id)
        _print(res if a.json else {k: v for k, v in res.items() if k != 'orca_input'}, True)
    elif a.cmd == 'audit':
        res = service.audit_rdb7(limit=10 if not a.json else 1000)
        if a.json:
            _print(res, True)
        else:
            print(f"{res['reactions_flagged']} of {res['reactions_checked']} RDB7 reactions have a structure whose bonding "
                  f"contradicts its atom-mapped SMILES{' (all in the product geometry)' if res['all_on_product_side'] else ''}.")
            print('  test-set reactions:', ', '.join(map(str, res['flagged_test_ids'])))
            for row in res['examples']:
                i = row['issues'][0]
                print(f"  {row['split']} {row['reaction_id']}: {i['side']} {i['atoms']} at {i['distance_A']:.2f} A: {i['kind']}")
    elif a.cmd == 'summary':
        _print(service.benchmark_summary()['headline'], True)


if __name__ == '__main__':
    main()
