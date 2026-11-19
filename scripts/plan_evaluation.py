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
F. Two goals: "find_errors" against "finish_ts" (the picks most likely to be right).
G. Random audits: is the bound on the trusted reactions' error rate as safe as it claims, and how wide is it?
H. Fewer guesses: the same test with 2, 3, 4 and 6 of the 8 models.
I. One guess only, other reactions: the harder reaction-core and barrier splits (one GoFlow model each), ranked by
   the guess's distance from the reactant and the product.
J. The reaction-path check (tsso/irc.py) on stand-in endpoints. No IRC was run: the reactant and product geometries
   stand in for the two ends of a path.
K. Which reference-free signals rank the errors hidden among trusted reactions (and a warning about having looked).

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


def two_goals(rows, wrong):
    """Section F. wrong: geometry or handedness wrong (a run from such a guess is unlikely to finish in the right TS)."""
    pos = {r['id']: k for k, r in enumerate(rows)}
    out = []
    for B in BUDGETS:
        row = {'budget': B}
        for goal in ('find_errors', 'finish_ts'):
            ch = [pos[i] for i in plan_dft_budget(rows, B, goal=goal)['run_dft']]
            row[goal] = {'picks': len(ch), 'wrong_guesses_found': int(wrong[ch].sum()), 'picks_with_a_right_guess': int((~wrong[ch]).sum())}
        out.append(row)
    return {'by_budget': out, 'note': 'A right guess is the start most likely to finish in the right TS (in the GFN2-xTB '
            'tests, 17 of 21 right "check with DFT" guesses reached the reference TS and 1 of 9 wrong ones). finish_ts '
            'only draws on "check with DFT" reactions; trusted ones need no run to be used.'}


def random_audits(rows, wrong, draws=5000):
    """Section G: how often does the audit bound hold, and how wide is it?"""
    from tsso.plan import review_results, audits_needed
    out = []
    for size, k in ((100, 10), (100, 30), (300, 30), (300, 60)):
        held, bound, true, zero = 0, [], [], 0
        for d in range(draws):
            idx = RNG.choice(len(rows), size, replace=False)
            sub = [rows[j] for j in idx]
            p = plan_dft_budget(sub, k, audit=k, seed=d)
            truth = {rows[j]['id']: bool(wrong[j]) for j in idx}
            rev = review_results(p, {i: 'wrong' if truth[i] else 'right' for i in p['run_dft']})['trusted_reactions']
            rest = [it['id'] for it in p['plan'] if it['action'] == 'trust']
            n_wrong = sum(truth[i] for i in rest)
            held += n_wrong <= rev['at_most_wrong']; zero += rev['audits_wrong'] == 0
            bound.append(rev['at_most_wrong'] / len(rest)); true.append(n_wrong / len(rest))
        out.append({'batch_size': size, 'random_audits': k, 'draws': draws,
                    'share_of_draws_where_the_bound_held': round(held / draws, 4),
                    'mean_bound_on_error_rate': round(float(np.mean(bound)), 3),
                    'mean_true_error_rate_of_the_rest': round(float(np.mean(true)), 3),
                    'share_of_draws_with_no_wrong_audit': round(zero / draws, 3)})
    return {'rows': out, 'confidence_claimed': 0.95,
            'audits_needed_for_a_10pct_bound_if_all_come_back_right': {str(n): audits_needed(n, 0.10) for n in (30, 50, 100, 300, 1000)},
            'note': 'Wrong = geometry or handedness wrong against the RDB7 reference. The bound is exact for a random draw, '
                    'so holding in at least 95% of draws is expected, not a discovery; the point is how wide it is.'}


def fewer_guesses(n_subsets=3):
    """Section H: triage with k of the 8 seed models, then count wrong guesses among the first picks."""
    import dataclasses
    from tsso.data import load_reactions
    from tsso.geometry import dmae
    from tsso.triage import load_checks, triage
    R, checks = load_reactions(), load_checks()
    out = []
    for k in (2, 3, 4, 6):
        hits = {B: [] for B in (5, 10, 25, 50)}
        base = []
        for _ in range(n_subsets):
            keep = np.sort(RNG.choice(8, k, replace=False))
            rows, wrong = [], []
            for r in R:
                rr = dataclasses.replace(r, preds=r.preds[keep], _cache={})
                t = triage(rr, checks, float('inf'))                 # the order does not depend on the amber cut-off
                rows.append({'id': str(r.id), 'verdict': 'red' if t.verdict == 'red' else 'amber', 'disagreement_A': t.disagreement})
                wrong.append(dmae(t.xyz, r.ts_ref_xyz) > 0.2)
            wrong = np.array(wrong); pos = {x['id']: j for j, x in enumerate(rows)}
            order = [pos[i] for i in plan_dft_budget(rows, 50)['run_dft']]
            base.append(float(wrong[[x['verdict'] != 'red' for x in rows]].mean()))
            for B in hits:
                hits[B].append(int(wrong[order[:B]].sum()))
        out.append({'guesses_per_reaction': k, 'subsets_tried': n_subsets, 'wrong_among_reactions_that_can_get_a_run': round(float(np.mean(base)), 3),
                    'wrong_among_first_picks': {str(B): {'mean': round(float(np.mean(v)), 1), 'min': min(v), 'max': max(v)} for B, v in hits.items()}})
    return {'rows': out, 'note': 'Random subsets of the 8 GoFlow seed models, all 1,194 reactions. Compare with section A (8 models).'}


def single_guess_splits(test_ids):
    """Section I: one model, reactions outside the random-split test set, ranked by endpoint distance."""
    from tsso.data import load_reactions
    from tsso.geometry import dmae
    from tsso.triage import load_checks, triage
    checks, out = load_checks(), {}
    for name, f in (('reaction_core_split', 'goflow_split0_test'), ('barrier_split', 'goflow_split1_test')):
        R = [r for r in load_reactions(ROOT / 'data' / f'{f}.json.gz') if r.id not in test_ids]
        rows, wrong = [], []
        for r in R:
            t = triage(r, checks, float('inf'))
            rows.append({'id': str(r.id), 'verdict': t.verdict, 'disagreement_A': None,
                         'endpoint_distance_A': min(dmae(t.xyz, r.reactant_xyz), dmae(t.xyz, r.product_xyz))})
            wrong.append(dmae(t.xyz, r.ts_ref_xyz) > 0.2)
        wrong = np.array(wrong); pos = {x['id']: j for j, x in enumerate(rows)}
        can = np.array([x['verdict'] != 'red' for x in rows]); base = float(wrong[can].mean())
        order = [pos[i] for i in plan_dft_budget(rows, 100)['run_dft']]
        boot = [RNG.integers(0, len(rows), len(rows)) for _ in range(1000)]
        by = []
        for B in (10, 25, 50, 100):
            shares = []
            for idx in boot:
                sub = [dict(rows[j], id=str(q)) for q, j in enumerate(idx)]
                ch = [int(i) for i in plan_dft_budget(sub, B)['run_dft']]
                shares.append(wrong[idx][ch].mean())
            lo, hi = np.percentile(shares, [2.5, 97.5])
            by.append({'budget': B, 'wrong_among_picked': int(wrong[order[:B]].sum()),
                       'share_95pct_interval': [round(float(lo), 3), round(float(hi), 3)],
                       'expected_if_picked_at_random': round(B * base, 1)})
        out[name] = {'reactions': len(R), 'excluded_because_in_the_random_split_test_set': True,
                     'wrong_among_reactions_that_can_get_a_run': round(base, 3), 'by_budget': by}
    out['note'] = ('The endpoint-distance signal was chosen after looking at the random-split test reactions (section K); '
                   'these splits use other reactions and another trained model, so they are the fair test of it.')
    return out


def path_check_sanity():
    """Section J: the reaction-path check on stand-in endpoints (no IRC was run)."""
    from tsso.data import load_reactions
    from tsso.irc import irc_check
    from tsso.service import _xyz_text
    M = np.array([1.0, 1.0, -1.0])
    c = {'reactions': 0, 'true_ends_called_right': 0, 'true_ends_swapped_called_right': 0, 'both_ends_the_reactant_called_wrong': 0,
         'noisy_true_ends_called_right': 0, 'reactions_with_stereocentres_compared': 0, 'mirror_image_ends_called_wrong_with_smiles': 0,
         'mirror_image_ends_called_wrong_without_smiles': 0, 'reactant_and_product_have_the_same_bonds': 0}
    for r in load_reactions():
        t = lambda X: _xyz_text(r.symbols, X)
        R_, P_ = t(r.reactant_xyz), t(r.product_xyz)
        c['reactions'] += 1
        a = irc_check(R_, P_, R_, P_, r.smiles)
        if 'same bonds' in a.get('why', ''):
            c['reactant_and_product_have_the_same_bonds'] += 1
            continue
        c['true_ends_called_right'] += a.get('outcome') == 'right'
        c['true_ends_swapped_called_right'] += irc_check(R_, P_, P_, R_, r.smiles).get('outcome') == 'right'
        c['both_ends_the_reactant_called_wrong'] += irc_check(R_, P_, R_, R_, r.smiles).get('outcome') == 'wrong'
        noise = np.random.default_rng(r.id).normal(0, 0.03, (2,) + r.reactant_xyz.shape)
        c['noisy_true_ends_called_right'] += irc_check(R_, P_, t(r.reactant_xyz + noise[0]), t(r.product_xyz + noise[1]),
                                                       r.smiles).get('outcome') == 'right'
        if a.get('handedness') == 'right':
            c['reactions_with_stereocentres_compared'] += 1
            c['mirror_image_ends_called_wrong_with_smiles'] += irc_check(R_, P_, t(r.reactant_xyz * M), t(r.product_xyz * M),
                                                                         r.smiles).get('outcome') == 'wrong'
            c['mirror_image_ends_called_wrong_without_smiles'] += irc_check(R_, P_, t(r.reactant_xyz * M),
                                                                            t(r.product_xyz * M)).get('outcome') == 'wrong'
    c['note'] = ('Stand-ins, not IRC output: the RDB7 reactant and product geometries play the two ends of a path; noise is '
                 '0.03 A Gaussian per coordinate. Real IRC ends are less tidy (they may stop before a bond has fully formed), '
                 'so this shows the check is consistent, not how it does on real paths.')
    return c


def hidden_error_signals(raw):
    """Section K: can anything rank the wrong guesses among trusted reactions? (No reference used by the signals.)"""
    from sklearn.metrics import roc_auc_score
    from scipy.stats import rankdata
    from tsso.data import load_reactions
    from tsso.geometry import dmae, disagreement_matrix
    R = {r.id: r for r in load_reactions()}
    u, e, big = [], [], []
    for x in raw:
        r = R[int(x['reaction_id'])]; c = int(x['chosen_model_seed']) - 1
        g = r.preds[c] * (np.array([1.0, 1.0, -1.0]) if x['mirrored'] == 'True' else 1.0)
        u.append(float(x['disagreement_A'])); e.append(min(dmae(g, r.reactant_xyz), dmae(g, r.product_xyz)))
        big.append(float(disagreement_matrix(r.preds).max()))
    u, e, big = np.array(u), np.array(e), np.array(big)
    y = np.array([x['geometry_bad'] == 'True' for x in raw]); v = np.array([x['verdict'] for x in raw])
    g, nr = v == 'green', v != 'red'
    sig = {'model disagreement of the chosen guess (the rule in use)': u,
           'largest disagreement between any two guesses': big,
           'distance of the guess from the nearer of reactant and product': e}
    out = {'trusted_reactions': int(g.sum()), 'of_which_wrong': int(y[g].sum()), 'signals': {}}
    for name, s_ in sig.items():
        o = np.argsort(-s_[g], kind='stable')
        out['signals'][name] = {'auroc_among_trusted': round(float(roc_auc_score(y[g], s_[g])), 3),
                                'auroc_among_all_that_can_get_a_run': round(float(roc_auc_score(y[nr], s_[nr])), 3),
                                'wrong_in_the_top_50_trusted': int(y[g][o][:50].sum())}
    both = rankdata(u[nr]) + rankdata(e[nr]); o = np.argsort(-both, kind='stable')
    out['disagreement_and_distance_combined_by_rank'] = {
        'auroc_among_all_that_can_get_a_run': round(float(roc_auc_score(y[nr], both)), 3),
        'wrong_in_first_picks': {str(B): int(y[nr][o][:B].sum()) for B in (5, 10, 25, 50, 100)},
        'verdict': 'Not adopted: it finds fewer wrong guesses in the first 5 to 25 picks than disagreement alone '
                   '(section A), the same number in the first 50 and more only by 100.'}
    out['random_top_50_trusted_would_hold'] = round(50 * float(y[g].mean()), 1)
    out['warning'] = ('Thirteen signals were tried on these test reactions before these three were kept, so the numbers '
                      'here are optimistic for the distance signal. Its fair test is section I.')
    return out, e


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
    out['F_two_goals'] = two_goals(rows, anyw)
    out['G_random_audits'] = random_audits(rows, anyw)
    out['H_fewer_guesses'] = fewer_guesses()
    out['I_single_guess_other_splits'] = single_guess_splits({int(r['reaction_id']) for r in raw})
    out['J_reaction_path_check_on_stand_ins'] = path_check_sanity()
    out['K_signals_for_hidden_errors'], endpoint = hidden_error_signals(raw)
    out['limits'] = [
        'Retrospective: counts guesses a check would expose; no new DFT was run for this evaluation.',
        'One generator family (GoFlow) and one dataset (RDB7: small gas-phase H/C/N/O molecules). Both ensembles use '
        'the same test reactions, so C tests another ensemble, not other chemistry.',
        'The harder reaction-core and barrier splits in this repository have one model each, so model disagreement '
        'is untested there; only the weaker single-guess order is (section I).',
        'No reaction path (IRC) was computed here; section J uses stand-in endpoints.',
        'The order has no fitted parameters, but the 0.2 A error threshold and the top-20% "check with DFT" cut-off '
        'are choices.']
    json.dump(out, open(ROOT / 'results' / 'plan_evaluation.json', 'w'), indent=1)
    demo = [{'id': r['reaction_id'], 'verdict': r['verdict'], 'disagreement_A': float(r['disagreement_A']),
             'endpoint_distance_A': round(float(endpoint[k]), 4),
             'mirror_repaired': r['mirrored'] == 'True', 'rules_broken': r['failed_rules'],
             'geometry_wrong': r['geometry_bad'] == 'True', 'handedness_wrong': r['stereo_wrong'] == 'True'}
            for k, r in enumerate(raw)]
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
    for key in ('F_two_goals', 'G_random_audits', 'H_fewer_guesses', 'I_single_guess_other_splits',
                'J_reaction_path_check_on_stand_ins', 'K_signals_for_hidden_errors'):
        print(key, json.dumps(out[key]))


if __name__ == '__main__':
    main()
