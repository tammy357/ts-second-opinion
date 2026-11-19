"""Tests for the DFT budget planner (tsso/plan.py) and the planning agent (tsso/agent.py). No network, no API key."""
import asyncio, csv, json, re
from pathlib import Path

import numpy as np

from tsso import agent, service
from tsso.plan import (audits_needed, explain_plan, plan_dft_budget, read_verdict_table, review_results, table_rows,
                       tiles_html, upper_bound_wrong)

ROOT = Path(__file__).resolve().parent.parent


def row(i, verdict, d, **more):
    return dict({'id': str(i), 'verdict': verdict, 'disagreement_A': d}, **more)


ROWS = [row('a', 'green', 0.03), row('b', 'amber', 0.20), row('c', 'amber', 0.15), row('d', 'red', 0.40, rules_broken='clash'),
        row('e', 'green', 0.09), row('f', 'amber', 0.12), row('g', 'green', None), row('h', '', None, error='missing product')]


def actions(plan):
    return {it['id']: it['action'] for it in plan['plan']}


def test_the_most_disputed_reactions_that_pass_the_rules_get_the_runs():
    p = plan_dft_budget(ROWS, 2)
    assert p['run_dft'] == ['b', 'c']                                   # most disputed first; never the red one
    assert actions(p) == {'b': 'run_dft', 'c': 'run_dft', 'f': 'hold', 'a': 'trust', 'e': 'trust', 'g': 'trust',
                          'd': 'regenerate', 'h': 'not_checked'}
    assert p['counts'] == {'run_dft': 2, 'trust': 3, 'hold': 1, 'regenerate': 1, 'not_checked': 1}
    assert any('on hold' in c for c in p['caveats']) and any('IRC' in c for c in p['caveats'])
    assert any('single guess' in c for c in p['caveats'])


def test_budget_edges():
    assert plan_dft_budget(ROWS, 0)['run_dft'] == [] and actions(plan_dft_budget(ROWS, 0))['b'] == 'hold'
    big = plan_dft_budget(ROWS, 50)                                     # more runs than rankable reactions
    assert big['run_dft'] == ['b', 'c', 'f', 'e', 'a'] and big['budget_unused'] == 45
    assert actions(big)['d'] == 'regenerate' and actions(big)['g'] == 'trust'      # still no run for red or unranked
    for bad in (-1, 2.5, '3', True, None):
        assert 'error' in plan_dft_budget(ROWS, bad)
    assert 'error' in plan_dft_budget(ROWS + [row('a', 'green', 0.01)], 1)          # duplicate id
    assert plan_dft_budget([], 3)['counts']['run_dft'] == 0


def test_same_input_same_plan_and_ties_break_by_id():
    tied = [row('z', 'amber', 0.2), row('m', 'amber', 0.2), row('b', 'amber', 0.2)]
    assert plan_dft_budget(tied, 2)['run_dft'] == ['b', 'm']
    assert plan_dft_budget(tied[::-1], 2)['run_dft'] == ['b', 'm']                  # input order does not matter
    rng = np.random.default_rng(0)
    shuffled = [ROWS[k] for k in rng.permutation(len(ROWS))]
    assert plan_dft_budget(shuffled, 2) == plan_dft_budget(ROWS, 2)


def test_plain_words_only_use_numbers_from_the_plan():
    p = plan_dft_budget(ROWS, 2)
    text = explain_plan(p)
    assert 'b, c' in text and 'Hold 1' in text
    assert agent.check_numbers(text, [p]) == []


def test_benchmark_plan_matches_the_recorded_evaluation_and_ignores_the_reference():
    plain = service.plan_dft_budget(10)
    scored = service.plan_dft_budget(10, evaluate=True)
    assert plain['run_dft'] == scored['run_dft'] and 'reference_check' not in plain
    ev = json.load(open(ROOT / 'results' / 'plan_evaluation.json'))['A_eight_seed_models_geometry']
    rec = next(x for x in ev['by_budget'] if x['budget'] == 10)
    assert scored['reference_check']['run_dft']['geometry_wrong'] == rec['wrong_among_picked']
    assert rec['wrong_among_picked'] > 3 * rec['expected_if_picked_at_random']       # the claim made in the docs
    assert plain['counts']['regenerate'] == 6 and not set(plain['run_dft']) & {it['id'] for it in plain['plan']
                                                                                if it['action'] == 'regenerate'}


def test_benchmark_sample_certified_trust_and_errors():
    p = service.plan_dft_budget(5, sample=40, seed=1, target_error_rate=0.05)
    assert p['reactions'] == 40 and len(p['run_dft']) == 5
    c = p['certified_trust']
    assert 1194 - 40 - 6 <= c['calibration_reactions'] <= 1194 - 40 and 0 < c['certified'] <= 40   # the others, minus rule-breakers
    assert service.plan_dft_budget(5, sample=40, seed=1)['run_dft'] == p['run_dft']
    assert 'error' in service.plan_dft_budget(5, reaction_ids=[2, 99999])
    assert 'error' in service.plan_dft_budget(5, sample=1)
    assert 'A_eight_seed_models_geometry' in service.plan_evidence()


def test_plan_from_a_verdict_table_on_disk(tmp_path):
    cols = ['id', 'verdict', 'mirror_repaired', 'disagreement_A', 'rules_broken', 'error', 'charge', 'multiplicity']
    with open(tmp_path / 'verdicts.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in ROWS:
            w.writerow({'id': r['id'], 'verdict': r['verdict'], 'disagreement_A': '' if r['disagreement_A'] is None else r['disagreement_A'],
                        'mirror_repaired': 'False', 'rules_broken': r.get('rules_broken', ''), 'error': r.get('error', ''),
                        'charge': 0 if r['verdict'] else '', 'multiplicity': 1 if r['verdict'] else ''})
    (tmp_path / 'b').mkdir(); (tmp_path / 'b' / 'recommended_ts.xyz').write_text('1\n\nH 0 0 0\n')
    assert [r['id'] for r in read_verdict_table(tmp_path / 'verdicts.csv')] == list('abcdefgh')
    for source in (tmp_path, tmp_path / 'verdicts.csv'):
        p = service.plan_dft_budget(2, str(source))
        assert p['run_dft'] == ['b', 'c'] and list(p['orca_inputs']) == ['b']
        assert any('No ORCA input was written for: c' in c for c in p['caveats'])
    inp = (tmp_path / 'b' / 'orca_optts_irc.inp').read_text()
    assert 'OptTS FREQ IRC' in inp and '%irc' in inp and '* xyzfile 0 1 recommended_ts.xyz' in inp
    rows = list(csv.DictReader(open(tmp_path / 'dft_plan.csv')))
    assert [r['id'] for r in rows][:2] == ['b', 'c'] and rows[0]['action'] == 'run_dft' and rows[0]['dft_order'] == '1'
    # your own calibration set: certified trust for your own reactions
    with open(tmp_path / 'cal.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['id', 'verdict', 'disagreement_A', 'wrong']); w.writeheader()
        for k in range(60):
            w.writerow({'id': k, 'verdict': 'green', 'disagreement_A': 0.02 + 0.004 * k, 'wrong': k > 52})
    c = service.plan_dft_budget(2, str(tmp_path), target_error_rate=0.2, calibration=str(tmp_path / 'cal.csv'))['certified_trust']
    assert c['calibration_reactions'] == 60 and 'a' in c['certified_ids']
    assert 'd' not in c['certified_ids']                                         # still breaks a rule: never certified
    assert 'error' in service.plan_dft_budget(2, str(tmp_path), target_error_rate=0.2, calibration=str(tmp_path / 'verdicts.csv'))
    assert any('needs --calibration' in x for x in service.plan_dft_budget(2, str(tmp_path), target_error_rate=0.2)['caveats'])


def test_goal_finish_ts_picks_the_least_disputed_reactions_that_need_a_check():
    p = plan_dft_budget(ROWS, 2, goal='finish_ts')
    assert p['run_dft'] == ['f', 'c'] and actions(p)['b'] == 'hold' and actions(p)['e'] == 'trust'
    assert plan_dft_budget(ROWS, 9, goal='finish_ts')['budget_unused'] == 6          # trusted reactions need no run
    assert 'error' in plan_dft_budget(ROWS, 2, goal='cheapest')
    ev = json.load(open(ROOT / 'results' / 'plan_evaluation.json'))['F_two_goals']['by_budget']
    for row in ev:                                                               # the trade-off the docs describe
        assert row['finish_ts']['picks_with_a_right_guess'] > row['find_errors']['picks_with_a_right_guess']
        assert row['find_errors']['wrong_guesses_found'] > row['finish_ts']['wrong_guesses_found']


def test_single_guesses_are_ranked_by_endpoint_distance_after_everything_else():
    rows = [row('s1', 'green', None, endpoint_distance_A=0.30), row('s2', 'green', None, endpoint_distance_A=0.45),
            row('s3', 'green', None), row('m', 'amber', 0.2), row('r', 'red', None, endpoint_distance_A=0.9)]
    p = plan_dft_budget(rows, 2)
    assert p['run_dft'] == ['m', 's2'] and 'weaker sign' in p['plan'][1]['why']
    assert plan_dft_budget(rows, 9)['run_dft'] == ['m', 's2', 's1']                # no distance, no rank; red never
    ev = json.load(open(ROOT / 'results' / 'plan_evaluation.json'))['I_single_guess_other_splits']
    for split in ('reaction_core_split', 'barrier_split'):                       # modest, but better than chance
        for b in ev[split]['by_budget']:
            assert b['wrong_among_picked'] > 1.5 * b['expected_if_picked_at_random']


def test_the_audit_bound_is_the_exact_hypergeometric_one():
    from scipy.stats import hypergeom
    for N, k, x in ((32, 5, 0), (100, 10, 1), (47, 2, 0), (200, 40, 3), (10, 10, 2), (9, 0, 0)):
        ub = upper_bound_wrong(N, k, x)
        M = ub + x                                                               # wrong in the whole population
        assert hypergeom.cdf(x, N, M, k) > 0.05                                  # still compatible with what was seen
        assert M == N - (k - x) or hypergeom.cdf(x, N, M + 1, k) <= 0.05         # and one more would not be
    assert upper_bound_wrong(10, 10, 2) == 0 and upper_bound_wrong(9, 0, 0) == 9
    assert audits_needed(100, 0.10) == 34 and upper_bound_wrong(100, 34, 0) <= 0.10 * 66 < upper_bound_wrong(100, 33, 0) + 1
    rng = np.random.default_rng(1)                                               # coverage, by simulation
    held = 0
    for _ in range(2000):
        pop = rng.random(80) < 0.15
        idx = rng.choice(80, 12, replace=False)
        rest = np.delete(pop, idx)
        held += rest.sum() <= upper_bound_wrong(80, 12, int(pop[idx].sum()))
    assert held / 2000 >= 0.95
    ev = json.load(open(ROOT / 'results' / 'plan_evaluation.json'))['G_random_audits']
    assert all(r['share_of_draws_where_the_bound_held'] >= 0.95 for r in ev['rows'])
    assert all(r['mean_bound_on_error_rate'] > 2 * r['mean_true_error_rate_of_the_rest'] for r in ev['rows'])   # wide


def test_random_audits_are_reproducible_and_reviewed_separately():
    demo = json.load(open(ROOT / 'results' / 'plan_demo.json'))['reactions'][:120]
    p = plan_dft_budget(demo, 12, audit=8, seed=3)
    assert len(p['run_dft']) == 12 and len(p['audited']) == 8 and p['of_which_random_audits'] == 8
    assert plan_dft_budget(demo[::-1], 12, audit=8, seed=3)['audited'] == p['audited']        # same draw, any input order
    assert plan_dft_budget(demo, 12, audit=8, seed=4)['audited'] != p['audited']
    trusted_before = {it['id'] for it in plan_dft_budget(demo, 4)['plan'] if it['action'] == 'trust'}
    assert set(p['audited']) <= trusted_before                                   # audits are drawn from trusted reactions only
    assert all(it['verdict'] == 'green' for it in p['plan'] if it.get('audit'))
    assert 'error' in plan_dft_budget(demo, 3, audit=4) and 'error' in plan_dft_budget(demo, 3, audit=-1)
    targeted = [i for i in p['run_dft'] if i not in p['audited']]
    part = review_results(p, {targeted[0]: 'wrong', p['audited'][0]: 'right'})
    assert part['targeted_runs']['to_regenerate'] == [targeted[0]] and part['trusted_reactions']['error_rate_bound'] is None
    full = review_results(p, {i: 'wrong' if i in (targeted[0], p['audited'][0]) else 'right' for i in p['run_dft']})
    t = full['trusted_reactions']
    n_rest = p['counts']['trust']
    assert t['audits'] == 8 and t['audits_wrong'] == 1 and t['at_most_wrong'] == upper_bound_wrong(n_rest + 8, 8, 1)
    assert full['targeted_runs'] == {'reported': 4, 'right': 3, 'wrong': 1, 'to_regenerate': [targeted[0]]}
    assert review_results(plan_dft_budget(demo, 4), {})['trusted_reactions']['error_rate_bound'] is None   # no audit, no rate
    assert 'error' in review_results(p, {'not-a-reaction': 'right'}) and 'error' in review_results(p, {targeted[0]: 'maybe'})
    via_service = service.review_dft_results({i: 'right' for i in service.plan_dft_budget(6, sample=80, seed=2, audit=3)['run_dft']},
                                             6, sample=80, seed=2, audit=3)
    assert via_service['trusted_reactions']['audits_wrong'] == 0 and via_service['trusted_reactions']['at_most_wrong'] > 0
    text = explain_plan(p)
    assert 'drawn at random' in text and agent.check_numbers(text, [p]) == []


def test_reaction_path_check_sees_bonds_and_handedness():
    from tsso.data import load_reactions
    from tsso.irc import irc_check, ORCA_IRC_TEMPLATE
    r = next(x for x in load_reactions() if x.id == 3509)
    t = lambda X: service._xyz_text(r.symbols, X)
    R_, P_, M = t(r.reactant_xyz), t(r.product_xyz), np.array([1.0, 1.0, -1.0])
    ok = irc_check(R_, P_, P_, R_, r.smiles)                                      # either branch may lead to either end
    assert ok['outcome'] == 'right' and ok['handedness'] == 'right' and ok['ends'] == {'forward': 'product', 'backward': 'reactant'}
    mirror = irc_check(R_, P_, t(r.reactant_xyz * M), t(r.product_xyz * M), r.smiles)
    assert mirror['outcome'] == 'wrong' and mirror['connects'] and mirror['wrong_centres']
    blind = irc_check(R_, P_, t(r.reactant_xyz * M), t(r.product_xyz * M))        # bonds alone cannot see a mirror image
    assert blind['outcome'] == 'right' and blind['handedness'] == 'not checked'
    back = irc_check(R_, P_, R_, R_, r.smiles)
    assert back['outcome'] == 'wrong' and not back['connects']
    assert irc_check(R_, P_, R_ + '\n' + P_, R_, r.smiles)['outcome'] == 'right'   # a trajectory: the last frame counts
    assert 'error' in irc_check(R_, P_, '2\n\nH 0 0 0\nH 0 0 0.7', P_)
    assert 'OptTS FREQ IRC' in ORCA_IRC_TEMPLATE
    ev = json.load(open(ROOT / 'results' / 'plan_evaluation.json'))['J_reaction_path_check_on_stand_ins']
    assert ev['true_ends_called_right'] == ev['reactions'] == ev['both_ends_the_reactant_called_wrong']
    assert ev['mirror_image_ends_called_wrong_with_smiles'] == ev['reactions_with_stereocentres_compared'] > 0
    assert ev['mirror_image_ends_called_wrong_without_smiles'] == 0


# ------------------------------------------------------------------ the agent
def fake_api(script):
    """A stand-in for the Messages API: replies with the scripted content blocks, and keeps what it was sent."""
    sent = []

    def send(body):
        sent.append(json.loads(json.dumps(body)))                       # must be JSON-serialisable, like a real call
        content = script[len(sent) - 1]
        return {'type': 'message', 'model': 'fake', 'content': content,
                'stop_reason': 'tool_use' if any(b['type'] == 'tool_use' for b in content) else 'end_turn'}
    return send, sent


def use(name, args, i='t1'):
    return {'type': 'tool_use', 'id': i, 'name': name, 'input': args}


def test_agent_calls_the_planner_and_never_sees_structures():
    send, sent = fake_api([[use('plan_dft_budget', {'budget': 5, 'sample': 40, 'seed': 1}), use('triage_reaction', {'reaction_id': 3509}, 't2')],
                           [{'type': 'text', 'text': 'Run DFT on 5 of your 40 reactions; trust 32. About 99% are fine.'}]])
    out = agent.run('I have 40 reactions and 5 DFT runs.', send)
    assert [s['tool'] for s in out['steps']] == ['plan_dft_budget', 'triage_reaction']
    assert out['steps'][0]['output']['run_dft'] == service.plan_dft_budget(5, sample=40, seed=1)['run_dft']
    assert out['numbers_not_from_tools'] == ['99']                      # the made-up number is flagged, the rest traced
    assert out['numbers']['typed_and_in_the_persons_messages'] == ['5', '40'] and '32' in out['numbers']['typed_and_found_in_a_tool_result']
    results = sent[1]['messages'][-1]['content']
    assert [json.loads(r['content'])['call'] for r in results] == [1, 2]            # call numbers, for references
    assert [r['tool_use_id'] for r in results] == ['t1', 't2'] and all(r['type'] == 'tool_result' for r in results)
    text = json.dumps(sent)
    assert 'recommended_ts_xyz' not in text and '[C:1]' not in text and 'smiles' not in text
    assert not re.search(r'[A-Z][a-z]? -?\d+\.\d{3,} -?\d+\.\d{3,} -?\d+\.\d{3,}', text)      # no coordinate lines
    assert sent[0]['system'] == agent.SYSTEM and {t['name'] for t in sent[0]['tools']} >= {'plan_dft_budget', 'plan_evidence'}


def test_agent_survives_bad_tool_calls_and_stops():
    send, sent = fake_api([[use('delete_everything', {})], [use('plan_dft_budget', {'budget': 3, 'nonsense': 1})],
                           [{'type': 'text', 'text': 'The tool failed.'}]])
    out = agent.run('plan it', send)
    assert 'no tool called' in out['steps'][0]['output']['error'] and 'TypeError' in out['steps'][1]['output']['error']
    assert sent[1]['messages'][-1]['content'][0]['is_error'] is True and out['answer'] == 'The tool failed.'
    loop, _ = fake_api([[use('plan_evidence', {})]] * 20)
    assert agent.run('loop', loop, max_steps=3)['stopped'] == 'max_steps'
    try:
        agent.run('x', lambda body: {'type': 'error', 'error': {'message': 'invalid x-api-key'}})
        assert False
    except RuntimeError as e:
        assert 'invalid x-api-key' in str(e)


def test_agent_async_and_custom_tools():
    tools = {'echo': ({'name': 'echo', 'description': 'd', 'input_schema': {'type': 'object', 'properties': {}}},
                      lambda: {'value': 0.117, 'xyz': 'C 0 0 0'})}
    send, sent = fake_api([[use('echo', {})], [{'type': 'text', 'text': 'It is 11.7%.'}]])

    async def asend(body):
        return send(body)
    out = asyncio.run(agent.run_async('q', asend, tools=tools, history=[{'role': 'user', 'content': 'hi'}, {'role': 'assistant', 'content': 'hello'}]))
    assert out['answer'] == 'It is 11.7%.' and out['numbers_not_from_tools'] == []
    assert out['steps'][0]['output'] == {'value': 0.117} and len(sent[0]['messages']) == 3


def test_numbers_are_filled_in_by_reference_not_typed():
    plan = service.plan_dft_budget(5, sample=40, seed=1)
    raw = ('Run DFT on {{1:run_dft}}. {{1:counts.hold}} are on hold and {{ 1:counts.trust }} trusted; the first disagrees by '
           '{{1:plan.0.disagreement_A}} A. On the benchmark {{2:D_what_is_left.rows.0.geometry_error_among_trusted|%}} of trusted '
           'guesses were wrong. Roughly 87% of plans work, {{3:nothing}} and {{1:counts}} and {{1:plan.99.id}}.')
    send, sent = fake_api([[use('plan_dft_budget', {'budget': 5, 'sample': 40, 'seed': 1}), use('plan_evidence', {}, 't2')],
                           [{'type': 'text', 'text': raw}]])
    out = agent.run('I have 40 reactions and 5 DFT runs.', send)
    assert out['raw_answer'] == raw
    assert f"Run DFT on {', '.join(plan['run_dft'])}." in out['answer']             # ids copied by the program
    assert f"{plan['counts']['hold']} are on hold and {plan['counts']['trust']} trusted" in out['answer']
    assert '4.4% of trusted guesses were wrong' in out['answer']
    n = out['numbers']
    assert len(n['filled_in_from_tool_results']) == 5
    # every number Claude typed itself is reported as typed, whether or not it happens to match something
    assert '87' in n['typed_and_found_nowhere'] + n['typed_and_found_in_a_tool_result']
    assert n['typed_and_in_the_persons_messages'] == []
    assert n['references_that_point_at_nothing'] == ['{{3:nothing}}', '{{1:counts}}', '{{1:plan.99.id}}']
    assert out['answer'].count('[no such value:') == 3
    note = agent.describe_numbers(out)
    assert '5 values above were copied from tool results' in note and '87' in note and 'do not exist' in note
    made_up = agent.finish({'messages': [{'role': 'user', 'content': 'plan it'}], 'steps': out['steps'], 'model': None,
                            'stopped': None, 'answer': 'About 61.7% succeed.'})
    assert made_up['numbers']['typed_and_found_nowhere'] == ['61.7'] and 'Do not rely' in agent.describe_numbers(made_up)
    assert '{{N:path}}' in agent.SYSTEM
    assert agent.render('{{0:x}} {{1:a.b}}', [{'output': {'a': {'b': [1, 2]}}}])[0] == '[no such value: 0:x] 1, 2'


def test_number_tracing():
    src = [{'rate': 0.0442, 'n': 1194, 'ids': ['rxn_1093', '24', '28'], 'note': 'cut-off 0.109 A; 9 of 10'}]
    ok = 'Of 1,194 reactions, 4.4% (4%) were wrong; cut-off 0.109 A; 9 of 10; reaction 1093; 3 ids.\n1. first\n2. second'
    assert agent.check_numbers(ok, src) == []
    assert agent.check_numbers('About 87% pass and 12 fail, at wB97X/def2-SVP with C3-H7.', src) == ['87', '12']
    assert agent.check_numbers('You said 5 runs for reactions 5, 6 and 7.', [], 'I can afford 5 runs') == ['6', '7']


def test_app_helpers_and_batch_tools():
    demo = json.load(open(ROOT / 'results' / 'plan_demo.json'))['reactions']
    assert len(demo) == 1194
    batch = demo[:40]
    tools = agent.batch_tools(batch, json.load(open(ROOT / 'results' / 'plan_evaluation.json')))
    p = tools['plan_dft_budget'][1](budget=5)
    assert 'geometry_wrong' not in json.dumps(p) and len(p['run_dft']) == 5      # reference labels never reach the agent
    assert len(tools['plan_dft_budget'][1](budget=5, goal='finish_ts', audit=2)['audited']) == 2
    assert all(it['action'] != 'trust' for it in p['plan']) and p['counts']['trust'] > 0
    ev = tools['plan_evidence'][1]()
    assert ev['B_batches'][1]['batch_size'] == 40 and ev['limits']
    full = plan_dft_budget(batch, 5)
    truth = {r['id']: r for r in batch}
    rows = table_rows(full, truth)
    assert len(rows) == 40 and rows[0]['What to do'] == 'Run DFT' and rows[0]['DFT order'] == 1
    assert set(r['Reference TS says (evaluation only)'] for r in rows) <= {'right', 'wrong geometry', 'wrong handedness'}
    assert 'Reference TS says (evaluation only)' not in table_rows(full)[0]
    html = tiles_html(full)
    assert html.count('border-top:4px') == 4 and 'Fix the input' not in html


def test_number_tracing_has_a_measured_blind_spot():
    """How often an arbitrary number would pass as traced. Recorded so the docs can state it; not a quality bar."""
    demo = json.load(open(ROOT / 'results' / 'plan_demo.json'))['reactions']
    tools = agent.batch_tools(demo[:40], json.load(open(ROOT / 'results' / 'plan_evaluation.json')))
    plan_only = [tools['plan_dft_budget'][1](budget=5)]
    with_evidence = plan_only + [tools['plan_evidence'][1]()]
    slip = lambda outs, texts: np.mean([agent.check_numbers(t, outs) == [] for t in texts])
    whole = [f'about {n} of them' for n in range(101)]
    decimals = [f'a rate of {n / 10:.1f}% here' for n in range(0, 1001, 7)]
    assert slip(plan_only, decimals) < 0.1 and slip(with_evidence, decimals) < 0.25
    assert slip(plan_only, whole) < 0.25 and 0.1 < slip(with_evidence, whole) < 0.5      # weak for small whole numbers


def test_mcp_and_cli_expose_the_planner(capsys):
    from tsso.cli import main
    main(['plan', '--budget', '5', '--sample', '40', '--seed', '1', '--evaluate'])
    text = capsys.readouterr().out
    assert 'Run DFT' in text and 'reference check (evaluation only)' in text and 'IRC' in text
    main(['plan', '--budget', '6', '--sample', '80', '--seed', '2', '--audit', '3', '--goal', 'finish_ts'])
    assert 'drawn at random' in capsys.readouterr().out
    from tsso.mcp_server import mcp
    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {'plan_dft_budget', 'plan_evidence', 'review_dft_results', 'irc_check'} <= names
