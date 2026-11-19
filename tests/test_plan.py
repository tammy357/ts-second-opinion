"""Tests for the DFT budget planner (tsso/plan.py) and the planning agent (tsso/agent.py). No network, no API key."""
import asyncio, csv, json, re
from pathlib import Path

import numpy as np

from tsso import agent, service
from tsso.plan import explain_plan, plan_dft_budget, read_verdict_table, table_rows, tiles_html

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
    assert c['calibration_reactions'] == 1194 - 40 and 0 < c['certified'] <= 40      # calibrated on the other reactions
    assert service.plan_dft_budget(5, sample=40, seed=1)['run_dft'] == p['run_dft']
    assert 'error' in service.plan_dft_budget(5, reaction_ids=[2, 99999])
    assert 'error' in service.plan_dft_budget(5, sample=1)
    assert 'A_eight_seed_models_geometry' in service.plan_evidence()


def test_plan_from_a_verdict_table_on_disk(tmp_path):
    cols = ['id', 'verdict', 'meaning', 'chosen_guess', 'mirror_repaired', 'disagreement_A', 'rules_broken', 'error']
    with open(tmp_path / 'verdicts.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in ROWS:
            w.writerow({'id': r['id'], 'verdict': r['verdict'], 'disagreement_A': '' if r['disagreement_A'] is None else r['disagreement_A'],
                        'mirror_repaired': 'False', 'rules_broken': r.get('rules_broken', ''), 'error': r.get('error', '')})
    (tmp_path / 'b').mkdir(); (tmp_path / 'b' / 'orca_optts.inp').write_text('! OptTS\n')
    assert [r['id'] for r in read_verdict_table(tmp_path / 'verdicts.csv')] == list('abcdefgh')
    for source in (tmp_path, tmp_path / 'verdicts.csv'):
        p = service.plan_dft_budget(2, str(source))
        assert p['run_dft'] == ['b', 'c'] and list(p['orca_inputs']) == ['b']
        assert any('No ORCA input was found for: c' in c for c in p['caveats'])
    rows = list(csv.DictReader(open(tmp_path / 'dft_plan.csv')))
    assert [r['id'] for r in rows][:2] == ['b', 'c'] and rows[0]['action'] == 'run_dft' and rows[0]['dft_order'] == '1'


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
    results = sent[1]['messages'][-1]['content']
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
    from tsso.mcp_server import mcp
    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {'plan_dft_budget', 'plan_evidence'} <= names
