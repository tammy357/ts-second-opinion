"""The planner page (docs/planner/) carries the planning rule in JavaScript. These tests keep it equal to tsso/plan.py."""
import hashlib, json, random, shutil, subprocess
from pathlib import Path

import pytest

from tsso.plan import audits_needed, explain_plan, plan_dft_budget, review_results, upper_bound_wrong

ROOT = Path(__file__).resolve().parent.parent
JS = r'''
const P = require(process.argv[2]), d = require(process.argv[3]); let bad = [];
const eq = (a, b, w) => { if (JSON.stringify(a) !== JSON.stringify(b)) bad.push(w + ': ' + JSON.stringify(a).slice(0, 200) + ' != ' + JSON.stringify(b).slice(0, 200)); };
d.sha.forEach(([s, h]) => eq(P.sha256hex(s), h, 'sha256'));
d.bounds.forEach(([N, k, x, u]) => eq(P.upperBoundWrong(N, k, x), u, 'bound'));
d.need.forEach(([n, k]) => eq(P.auditsNeeded(n, 0.1), k, 'audits needed'));
d.cases.forEach((c, i) => { const p = P.plan(c.rows, c.budget, {goal: c.goal, audit: c.audit, seed: c.seed}), e = d.expected[i];
  eq(p.run_dft, e.run_dft, 'run_dft ' + i); eq(p.audited, e.audited, 'audited ' + i); eq(p.counts, e.counts, 'counts ' + i);
  eq(p.plan.map(it => [it.id, it.action, it.dft_order, it.why]), e.actions, 'plan ' + i);
  eq([p.budget_unused, p.of_which_optional_spot_checks], e.extra, 'extra ' + i); eq(P.explain(p), e.text, 'words ' + i);
  eq(P.review(p, e.outcomes).at_most_wrong, e.at_most_wrong, 'review ' + i); });
eq(P.parseCsv('id,verdict,disagreement_A,mirror_repaired,rules_broken\n"a,1",green,0.05,True,"x ""y"""\nb,red,,False,\n')
    .map(r => [r.id, r.verdict, r.disagreement_A, r.mirror_repaired, r.rules_broken]), [['a,1', 'green', 0.05, true, 'x "y"'], ['b', 'red', null, false, '']], 'csv');
const r = P.render('{{1:run_dft}} and {{1:counts.hold}} and {{1:rate|%}} and {{2:x}}', [{run_dft: ['4', '7'], counts: {hold: 2}, rate: 0.0442}]);
eq([r.text, r.filled.length, r.unresolved], ['4, 7 and 2 and 4.4% and [no such value: 2:x]', 3, ['{{2:x}}']], 'render');
eq(P.typedNumbers('Run {{1:run_dft}}; about 87% of your 5 runs, C3-H7.', 'I have 5 runs'), ['87'], 'typed numbers');
console.log(JSON.stringify(bad.slice(0, 5))); process.exit(bad.length ? 1 : 0);
'''


@pytest.mark.skipif(shutil.which('node') is None, reason='needs node')
def test_the_page_makes_the_same_plan_as_python(tmp_path):
    demo = json.load(open(ROOT / 'results' / 'plan_demo.json'))['reactions']
    keep = ('id', 'verdict', 'disagreement_A', 'endpoint_distance_A', 'mirror_repaired', 'rules_broken')
    rng, cases, expected = random.Random(0), [], []
    for c in range(120):
        rows = [{k: r[k] for k in keep} for r in rng.sample(demo, rng.choice([5, 12, 40, 90, 200]))]
        if c % 5 == 0:
            for r in rows:
                r['disagreement_A'] = None                       # one guess per reaction
        if c % 7 == 0:
            rows[0].update(verdict='', error='missing product')
        if c % 9 == 0:
            for r in rows[:6]:
                r['disagreement_A'] = 0.0625                     # ties, and a value exactly half way between 3-decimal numbers
        budget = rng.choice([0, 1, 3, 5, 10, 25, 60])
        cases.append({'rows': rows, 'budget': budget, 'goal': rng.choice(['find_errors', 'finish_ts']),
                      'audit': min(rng.choice([0, 0, 1, 3, 10]), budget), 'seed': rng.randrange(50)})
    for c in cases:
        p = plan_dft_budget(c['rows'], c['budget'], goal=c['goal'], audit=c['audit'], seed=c['seed'])
        out = {i: 'right' if int(hashlib.md5(i.encode()).hexdigest(), 16) % 3 else 'wrong' for i in p['run_dft']}
        expected.append({'run_dft': p['run_dft'], 'audited': p['audited'], 'counts': p['counts'],
                         'actions': [[it['id'], it['action'], it['dft_order'], it['why']] for it in p['plan']],
                         'extra': [p['budget_unused'], p['of_which_optional_spot_checks']], 'text': explain_plan(p), 'outcomes': out,
                         'at_most_wrong': review_results(p, out)['trusted_reactions'].get('at_most_wrong')})
    data = {'cases': cases, 'expected': expected,
            'bounds': [[N, k, x, upper_bound_wrong(N, k, x)] for N, k, x in ((32, 5, 0), (100, 10, 1), (951, 60, 2), (10, 10, 2), (9, 0, 0))],
            'need': [[n, audits_needed(n, 0.1)] for n in (30, 100, 300)],
            'sha': [[s, hashlib.sha256(s.encode()).hexdigest()] for s in ('', '0:432', '7:rxn_3509', 'x' * 55, 'y' * 56, 'z' * 200, 'é:1')]}
    (tmp_path / 'cases.json').write_text(json.dumps(data))
    (tmp_path / 'run.js').write_text(JS)
    res = subprocess.run(['node', str(tmp_path / 'run.js'), str(ROOT / 'apps' / 'planner' / 'planner.js'), str(tmp_path / 'cases.json')],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stdout + res.stderr


def test_the_built_page_is_current_and_self_contained():
    import sys
    sys.path.insert(0, str(ROOT / 'scripts'))
    import build_planner
    built = (ROOT / 'docs' / 'planner' / 'index.html').read_text()
    assert built == build_planner.page(), 'run python scripts/build_planner.py'
    assert '__PLANNER_JS__' not in built and '__DATA__' not in built and built.count('<script') == 2
    assert 'src=' not in built and 'https://api.anthropic.com/v1/messages' in built       # no outside scripts; one outside call
    assert '[C:1]' not in built and 'recommended_ts_xyz' not in built                     # no structures in the page
    from tsso.agent import SYSTEM
    assert json.dumps(SYSTEM)[1:40] in built                                              # the same instructions as tsso agent
