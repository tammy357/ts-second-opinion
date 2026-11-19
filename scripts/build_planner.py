"""Build the DFT budget planner page for GitHub Pages: docs/planner/index.html, one self-contained file.

The page is plain HTML and JavaScript (no Python in the browser), so it opens at once on a phone. It carries the
planning rule (apps/planner/planner.js, a port of tsso/plan.py that tests/test_planner_page.py compares with the
Python), the verdicts for the 1,194 benchmark reactions and the evaluation results.

Usage: python scripts/build_planner.py
"""
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
SRC, OUT = ROOT / 'apps' / 'planner', ROOT / 'docs' / 'planner'
EVIDENCE_KEYS = ('question', 'wrong_means', 'A_eight_seed_models_geometry', 'B_batches', 'D_what_is_left', 'E_what_the_runs_return',
                 'F_two_goals', 'G_random_audits', 'H_fewer_guesses', 'I_single_guess_other_splits',
                 'J_reaction_path_check_on_stand_ins', 'limits')


def page():
    from tsso.agent import SYSTEM
    from tsso.llm import API_MODEL
    demo = json.load(open(ROOT / 'results' / 'plan_demo.json'))['reactions']
    code = {'green': 0, 'amber': 1, 'red': 2}
    data = [[r['id'], code[r['verdict']], r['disagreement_A'], r['endpoint_distance_A'], int(r['mirror_repaired']),
             int(r['geometry_wrong']), int(r['handedness_wrong']), r['rules_broken']] for r in demo]
    ev = json.load(open(ROOT / 'results' / 'plan_evaluation.json'))
    ev = {k: ev[k] for k in EVIDENCE_KEYS}
    inline = lambda x: json.dumps(x, separators=(',', ':')).replace('</', '<\\/')
    html = (SRC / 'template.html').read_text()
    for mark, value in (('/*__PLANNER_JS__*/', (SRC / 'planner.js').read_text().replace('</', '<\\/')),
                        ('/*__DATA__*/null', inline(data)), ('/*__EVIDENCE__*/null', inline(ev)),
                        ("/*__SYSTEM__*/''", inline(SYSTEM)), ("/*__MODEL__*/''", inline(API_MODEL))):
        if html.count(mark) != 1:
            raise SystemExit(f'the template should contain {mark} exactly once')
        html = html.replace(mark, value)
    return html


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    html = page()
    (OUT / 'index.html').write_text(html)
    print(f'wrote {(OUT / "index.html").relative_to(ROOT)}: {len(html) / 1e3:.0f} kB')


if __name__ == '__main__':
    main()
