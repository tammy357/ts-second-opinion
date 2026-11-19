"""Does the DFT budget plan (tsso/plan.py) put the runs in the right place? A retrospective test on RDB7.

The plan never sees the reference TS. Here the reference is used afterwards, to count how many of the reactions it
picked really had a wrong recommended guess ("wrong" = D-MAE to the reference TS above 0.2 A; "any error" adds a
stereocentre the wrong way round). A DFT run with a reaction-path check would expose such a guess.

A. All 1,194 test reactions, 8 GoFlow seed models: wrong guesses among the first B picks, against picking at random,
   picking "check with DFT" reactions in no particular order, picking the least disputed first, and a perfect order.
   Intervals: 2,000 bootstrap resamples of the reactions.
B. Batches the size a person brings (20, 40 and 100 reactions drawn at random, 5,000 times each), budget 5.
C. The same test with a different ensemble: five GoFlow runs of four model sizes on the same reactions.
D. What is left: the error rate among answers that would be used without a check.
E. What the runs return: outcomes of the TS optimisations already in this repository, by verdict.

Usage: python scripts/plan_evaluation.py   -> results/plan_evaluation.json and results/plan_demo.json
"""
import csv, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from tsso.plan import plan_dft_budget

BUDGETS = (5, 10, 25, 50, 100)
RNG = np.random.default_rng(0)


def picks(rows, budget):
    """Indices of the reactions the plan sends to DFT, in order."""
    p = plan_dft_budget(rows, budget)
    pos = {r['id']: k for k, r in enumerate(rows)}
    return [pos[i] for i in p['run_dft']], p


def summarise(rows, wrong, label):
    """Section A for one verdict table. rows: plan inputs; wrong: bool array from the reference TS."""
    n = len(rows)
    verdict = np.array([r['verdict'] for r in rows])
    eligible = verdict != 'red'                       # the plan gives no run to re-generate reactions
    base = float(wrong[eligible].mean())
    amber = np.nonzero(verdict == 'amber')[0]
    order_all, _ = picks(rows, n)
    out = {'what': label, 'reactions': n, 'eligible_for_a_run': int(eligible.sum()),
           'wrong_among_eligible': round(base, 4), 'wrong_among_check_with_dft': round(float(wrong[amber].mean()), 4),
           'by_budget': []}
    boot = [RNG.integers(0, n, n) for _ in range(2000)]
    for B in BUDGETS:
        chosen = order_all[:B]
        hits = int(wrong[chosen].sum())
        reverse = [k for k in order_all[::-1]][:B]
        lo_hi = []
        for idx in boot:                               # resample reactions, re-plan, recount
            sub = [dict(rows[k], id=str(j)) for j, k in enumerate(idx)]
            ch, _ = picks(sub, B)
            lo_hi.append(wrong[idx][ch].mean())
        lo, hi = np.percentile(lo_hi, [2.5, 97.5])
        out['by_budget'].append({
            'budget': B, 'wrong_among_picked': hits, 'share_wrong_among_picked': round(hits / B, 3),
            'share_95pct_interval': [round(float(lo), 3), round(float(hi), 3)],
            'expected_if_picked_at_random': round(B * base, 2),
            'expected_if_any_check_with_dft_reaction': round(B * float(wrong[amber].mean()), 2),
            'least_disputed_first': int(wrong[reverse].sum()),
            'perfect_order': int(min(B, wrong[eligible].sum())),
            'times_better_than_random': round(hits / (B * base), 1)})
    return out


def batches(rows, wrong, size, budget=5, draws=5000):
    """Section B: random batches of `size` reactions."""
    n = len(rows)
    got, rand, best, tie_or_better = [], [], [], 0
    for _ in range(draws):
        idx = RNG.choice(n, size, replace=False)
        sub = [rows[k] for k in idx]
        ch, p = picks(sub, budget)
        w = wrong[idx]
        el = np.array([r['verdict'] != 'red' for r in sub])
        g = int(w[ch].sum()); r_ = budget * float(w[el].mean()) if el.any() else 0.0
        got.append(g); rand.append(r_); best.append(min(budget, int(w[el].sum())))
        tie_or_better += g >= r_
    return {'batch_size': size, 'budget': budget, 'draws': draws,
            'mean_wrong_among_picked': round(float(np.mean(got)), 2),
            'mean_if_picked_at_random': round(float(np.mean(rand)), 2),
            'mean_with_perfect_order': round(float(np.mean(best)), 2),
            'share_of_batches_at_least_as_good_as_random': round(tie_or_better / draws, 3),
            'share_of_batches_where_no_pick_was_wrong': round(float(np.mean(np.array(got) == 0)), 3)}


def second_ensemble():
    """Section C: triage the five-model-size ensemble (no reference used), then score with the reference."""
    from tsso.data import load_reactions
    from tsso.geometry import dmae
    from tsso.triage import load_checks, triage, amber_cutoff_for_budget
    R = load_reactions(ROOT / 'data' / 'goflow_sizes_rdb7_test.json.gz')
    checks = load_checks()
    cut = amber_cutoff_for_budget(R, checks, 0.2)
    rows, wrong = [], []
    for r in R:
        t = triage(r, checks, cut)
        rows.append({'id': str(r.id), 'verdict': t.verdict, 'disagreement_A': t.disagreement, 'mirror_repaired': t.mirrored})
        wrong.append(dmae(t.xyz, r.ts_ref_xyz) > 0.2)
    return rows, np.array(wrong)


def run_outcomes():
    """Section E: what the TS optimisations in this repository ended in, by verdict."""
    out = {}
    tri = {r['reaction_id']: r for r in csv.DictReader(open(ROOT / 'results' / 'triage.csv'))}
    x = json.load(open(ROOT / 'results' / 'xtb_prospective_verdicts.json'))
    tab = {}
    for r in x:
        t = tab.setdefault(r['verdict'], {'runs': 0, 'reached_the_reference_ts': 0, 'one_imaginary_frequency': 0,
                                           'one_imaginary_frequency_but_not_the_reference_ts': 0})
        t['runs'] += 1; t['reached_the_reference_ts'] += bool(r['found_the_reference_ts'])
        one = r.get('n_imaginary') == 1
        t['one_imaginary_frequency'] += one
        t['one_imaginary_frequency_but_not_the_reference_ts'] += one and not r['found_the_reference_ts']
    amber = [r for r in x if r['verdict'] == 'amber']
    bad = [r for r in amber if tri[str(r['reaction_id'])]['geometry_bad'] == 'True']
    good = [r for r in amber if tri[str(r['reaction_id'])]['geometry_bad'] != 'True']
    out['gfn2_xtb'] = {'by_verdict': tab,
                       'check_with_dft_guess_was_wrong': {'runs': len(bad), 'reached_the_reference_ts': sum(bool(r['found_the_reference_ts']) for r in bad)},
                       'check_with_dft_guess_was_right': {'runs': len(good), 'reached_the_reference_ts': sum(bool(r['found_the_reference_ts']) for r in good)},
                       'note': 'GFN2-xTB saddle searches scored against DFT reference TSs, so some misses are the method, not '
                               'the guess. First 30 reactions by id per verdict, all 6 re-generate (scripts/xtb_prospective.py).'}
    d = json.load(open(ROOT / 'results' / 'dft_summary.json'))
    out['dft'] = {'level': d['level'], 'by_verdict': d['by_verdict'],
                  'note': '15 reactions; too few for rates, shown as counts (scripts/dft_prospective.py).'}
    return out


def main():
    raw = list(csv.DictReader(open(ROOT / 'results' / 'triage.csv')))
    rows = [{'id': r['reaction_id'], 'verdict': r['verdict'], 'disagreement_A': float(r['disagreement_A']),
             'mirror_repaired': r['mirrored'] == 'True'} for r in raw]
    geo = np.array([r['geometry_bad'] == 'True' for r in raw])
    anyw = geo | np.array([r['stereo_wrong'] == 'True' for r in raw])
    out = {'question': 'Of the reactions the plan sends to DFT, how many really had a wrong recommended guess?',
           'wrong_means': 'D-MAE of the recommended guess to the RDB7 reference TS above 0.2 A (geometry); '
                          '"any error" adds a stereocentre the wrong way round. The plan never uses the reference.',
           'A_eight_seed_models_geometry': summarise(rows, geo, '8 GoFlow seed models, 1,194 RDB7 test reactions'),
           'A_eight_seed_models_any_error': summarise(rows, anyw, 'the same, counting geometry or handedness errors'),
           'B_batches': [batches(rows, geo, s) for s in (20, 40, 100)]}
    rows2, geo2 = second_ensemble()
    out['C_five_model_sizes_geometry'] = summarise(rows2, geo2, 'five GoFlow runs of four model sizes, same 1,194 reactions')
    left = []
    n_amber = sum(r['verdict'] == 'amber' for r in rows)
    for B in (0, 25, 50, 100, n_amber):
        p = plan_dft_budget(rows, B)
        act = {it['id']: it['action'] for it in p['plan']}
        a = np.array([act[r['id']] for r in rows])
        unchecked = (a == 'trust') | (a == 'hold')              # what you would be using if you ignored "hold"
        left.append({'budget': B, 'trusted': int((a == 'trust').sum()),
                     'geometry_error_among_trusted': round(float(geo[a == 'trust'].mean()), 4),
                     'any_error_among_trusted': round(float(anyw[a == 'trust'].mean()), 4),
                     'on_hold': int((a == 'hold').sum()),
                     'geometry_error_among_on_hold': round(float(geo[a == 'hold'].mean()), 4) if (a == 'hold').any() else None,
                     'geometry_error_among_everything_not_run': round(float(geo[unchecked].mean()), 4)})
    out['D_what_is_left'] = {'rows': left, 'note': 'The trusted set is fixed by the verdicts, so its error rate does not '
                             'fall with the budget; the budget shrinks the on-hold set. Errors among trusted reactions are '
                             'the ones this plan cannot find.'}
    out['E_what_the_runs_return'] = run_outcomes()
    out['limits'] = [
        'Retrospective: counts guesses a check would expose; no new DFT was run for this evaluation.',
        'One generator family (GoFlow) and one dataset (RDB7: small gas-phase H/C/N/O molecules). Both ensembles use '
        'the same test reactions, so C tests another ensemble, not other chemistry.',
        'The harder reaction-core and barrier splits in this repository have one model each, so there is no '
        'disagreement to rank by and the plan is untested there.',
        'The order has no fitted parameters, but the 0.2 A error threshold and the top-20% "check with DFT" cut-off '
        'are choices.']
    json.dump(out, open(ROOT / 'results' / 'plan_evaluation.json', 'w'), indent=1)
    demo = [{'id': r['reaction_id'], 'verdict': r['verdict'], 'disagreement_A': float(r['disagreement_A']),
             'mirror_repaired': r['mirrored'] == 'True', 'rules_broken': r['failed_rules'],
             'geometry_wrong': r['geometry_bad'] == 'True', 'handedness_wrong': r['stereo_wrong'] == 'True'} for r in raw]
    json.dump({'description': 'Triage verdicts for the 1,194 RDB7 test reactions (from results/triage.csv), for the '
                              'budget planner in the app. geometry_wrong and handedness_wrong come from the reference TS '
                              'and are for evaluation only.', 'reactions': demo},
              open(ROOT / 'results' / 'plan_demo.json', 'w'), separators=(',', ':'))
    a = out['A_eight_seed_models_geometry']['by_budget']
    for x in a:
        print(f"budget {x['budget']:>3}: {x['wrong_among_picked']} wrong among picked "
              f"({100 * x['share_wrong_among_picked']:.0f}%, 95% interval {x['share_95pct_interval']}), "
              f"random {x['expected_if_picked_at_random']}, any amber {x['expected_if_any_check_with_dft_reaction']}, "
              f"perfect {x['perfect_order']}")
    print(json.dumps(out['B_batches'], indent=1))
    print(json.dumps(out['C_five_model_sizes_geometry']['by_budget'], indent=1))
    print(json.dumps(out['D_what_is_left'], indent=1))
    print(json.dumps(out['E_what_the_runs_return'], indent=1))


if __name__ == '__main__':
    main()
