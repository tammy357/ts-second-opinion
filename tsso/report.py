"""A single-file HTML report card for one reaction: verdict, reasons, rule results, stereocentres, linked 3D views of
the raw and recommended TS guesses, the recommended geometry and an ORCA input. Opens offline in any browser."""
import datetime, html, json
from pathlib import Path

import numpy as np

from . import __version__
from .sandbox import run_check
from .stereo import stereocentre_status
from .triage import ORCA_TEMPLATE

ASSETS = Path(__file__).resolve().parent / 'assets'
LABEL = {'green': ('✓', 'Trust', 'Passes every rule and the models agree.'),
         'amber': ('!', 'Check with DFT', 'Passes the rules, but the models disagree a lot.'),
         'red': ('✕', 'Re-generate', 'Still breaks a rule after repair.')}


def _xyz(symbols, x, comment=''):
    return '\n'.join([str(len(symbols)), comment] + [f'{s} {a:.5f} {b:.5f} {c:.5f}' for s, (a, b, c) in zip(symbols, x)])


def _mark(status):
    return {True: '<span class="ok">✓ right</span>', False: '<span class="bad">✕ inverted</span>',
            None: '<span class="muted">flat, not judged</span>'}[status]


def report_html(rxn, t, checks, charge=0, mult=1, title='Your reaction', score=None):
    """rxn: tsso.data.Reaction; t: tsso.triage.Triage; checks: {name: compiled check}; score: optional dict from
    scoring against a reference TS (shown as evaluation only)."""
    raw = rxn.preds[t.chosen_model]
    rec = t.xyz
    st_raw, st_rec = stereocentre_status(rxn, raw), stereocentre_status(rxn, rec)
    rule_rows = []
    for name, fn in checks.items():
        ok, detail = run_check(fn, rxn, rec)
        res = {True: '<span class="ok">✓ passes</span>', False: '<span class="bad">✕ fails</span>',
               None: '<span class="muted">could not run</span>'}[ok]
        if ok and str(detail).startswith('rule not applicable'):
            res, detail = '<span class="muted">– n/a</span>', 'nothing in this reaction for this rule to check'
        rule_rows.append(f'<tr><td><code>{html.escape(name)}</code></td><td>{res}</td><td>{html.escape(str(detail))}</td></tr>')
    stereo_rows = ''.join(
        f'<tr><td>{a["label"]}</td><td>{a["compared_with"]}</td><td>{_mark(a["ok"])}</td><td>{_mark(b["ok"])}</td></tr>'
        for a, b in zip(st_raw, st_rec)) or '<tr><td colspan="4" class="muted">No stereocentres in this reaction.</td></tr>'

    def centres(status, x):
        return [{'pos': {'x': float(x[s['atom']][0]), 'y': float(x[s['atom']][1]), 'z': float(x[s['atom']][2])},
                 'ok': s['ok'], 'text': s['label']} for s in status]

    icon, word, meaning = LABEL[t.verdict]
    src = f'guess {t.chosen_model + 1}' + (' (mirror image)' if t.mirrored else '') + f'; verdict {t.verdict}'
    rec_xyz = _xyz(rxn.symbols, rec, f'TS Second Opinion recommended guess; {src}')
    orca = ORCA_TEMPLATE.format(rid=title, source=src, xyzname='recommended_ts.xyz', charge=charge, mult=mult)
    facts = [f'{sum(1 for k in range(len(rxn.preds)) if k not in t.rejected_models)} of {len(rxn.preds)} raw guesses pass every rule',
             'mirror-image repair applied' if t.mirrored else 'no repair needed' if t.verdict != 'red' else 'repair could not fix it',
             'model disagreement not measurable (one guess)' if np.isnan(t.disagreement) else f'model disagreement {t.disagreement:.3f} Å',
             f'charge {charge}, multiplicity {mult}']
    score_html = ''
    if score:
        score_html = (f'<section><h2>Against the reference TS <span class="tag">evaluation only</span></h2><p>D-MAE '
                      f'{score["dmae_A"]:.3f} Å; stereochemistry {"correct" if score["stereo_correct"] else "wrong"}. '
                      'The reference was not used to make the decision.</p></section>')
    js = (ASSETS / '3Dmol-min.js').read_text()
    data = json.dumps({'raw': _xyz(rxn.symbols, raw), 'rec': _xyz(rxn.symbols, rec),
                       'craw': centres(st_raw, raw), 'crec': centres(st_rec, rec)})
    page = TEMPLATE
    for k, v in {
        '__TITLE__': html.escape(title), '__SMILES__': html.escape(rxn.smiles), '__VERDICT__': t.verdict,
        '__ICON__': icon, '__WORD__': word, '__MEANING__': meaning,
        '__FACTS__': ''.join(f'<li>{html.escape(f)}</li>' for f in facts),
        '__REASONS__': ''.join(f'<li>{html.escape(r)}</li>' for r in t.reasons),
        '__RULES__': ''.join(rule_rows), '__STEREO__': stereo_rows, '__SCORE__': score_html,
        '__RAWLABEL__': f'Best raw guess (guess {t.chosen_model + 1})',
        '__RECLABEL__': 'Recommended (mirror-repaired)' if t.mirrored else 'Recommended',
        '__RECXYZ__': html.escape(rec_xyz), '__ORCA__': html.escape(orca),
        '__DATE__': datetime.date.today().isoformat(), '__VERSION__': __version__,
        '__DATA__': data, '__JS__': js,
    }.items():
        page = page.replace(k, v)
    return page


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>TS report card</title>
<style>
:root { --bg:#f7f7f5; --surface:#ffffff; --ink:#111110; --ink2:#55534e; --line:#e3e2de; --code:#f1f0ec;
        --green:#137333; --green-bg:#e6f4ea; --amber:#9a5b00; --amber-bg:#fdf1dc; --red:#b3261e; --red-bg:#fbe9e7; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#141413; --surface:#1d1d1b; --ink:#f2f1ed; --ink2:#b5b3ac; --line:#34332f; --code:#262522;
  --green:#7fd68f; --green-bg:#173a20; --amber:#f2c063; --amber-bg:#3d2e10; --red:#f28b82; --red-bg:#43201d; } }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink); font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:1040px; margin:0 auto; padding:28px 16px 48px; }
header .kicker { color:var(--ink2); font-size:13px; letter-spacing:.04em; text-transform:uppercase; }
h1 { font-size:26px; margin:4px 0 6px; } h2 { font-size:17px; margin:0 0 10px; }
.smiles { font:12px/1.4 ui-monospace,Menlo,monospace; color:var(--ink2); word-break:break-all; }
section { background:var(--surface); border:1px solid var(--line); border-radius:12px; padding:18px; margin-top:16px; }
.verdict { display:flex; gap:16px; align-items:center; border-width:2px; }
.verdict .badge { font-size:30px; font-weight:800; width:62px; height:62px; border-radius:14px; display:grid; place-items:center; flex:none; }
.verdict.green { border-color:var(--green); } .verdict.green .badge { background:var(--green-bg); color:var(--green); }
.verdict.amber { border-color:var(--amber); } .verdict.amber .badge { background:var(--amber-bg); color:var(--amber); }
.verdict.red { border-color:var(--red); } .verdict.red .badge { background:var(--red-bg); color:var(--red); }
.verdict .word { font-size:22px; font-weight:750; } .verdict p { margin:2px 0 0; color:var(--ink2); }
ul.facts { display:flex; flex-wrap:wrap; gap:8px; list-style:none; padding:0; margin:12px 0 0; }
ul.facts li { background:var(--code); border-radius:999px; padding:3px 10px; font-size:13px; color:var(--ink2); }
.views { display:grid; grid-template-columns:1fr 1fr; gap:12px; }
.view { border:1px solid var(--line); border-radius:10px; overflow:hidden; }
.view .cap { font-size:13px; color:var(--ink2); padding:8px 10px; border-bottom:1px solid var(--line); }
.viewer { width:100%; height:300px; position:relative; }
table { width:100%; border-collapse:collapse; font-size:14px; }
th, td { text-align:left; padding:7px 8px; border-bottom:1px solid var(--line); vertical-align:top; }
th { color:var(--ink2); font-weight:600; font-size:13px; }
.ok { color:var(--green); font-weight:600; } .bad { color:var(--red); font-weight:600; } .muted { color:var(--ink2); }
.tag { font-size:12px; font-weight:500; color:var(--ink2); border:1px solid var(--line); border-radius:999px; padding:1px 8px; margin-left:6px; }
pre { background:var(--code); border-radius:8px; padding:10px; overflow:auto; font-size:12px; margin:8px 0 0; }
details summary { cursor:pointer; color:var(--ink2); }
.legend { font-size:13px; color:var(--ink2); margin-top:8px; }
.dot { display:inline-block; width:10px; height:10px; border-radius:50%; vertical-align:middle; margin:0 4px 0 10px; }
footer { color:var(--ink2); font-size:12px; margin-top:20px; }
@media (max-width:720px) { .views { grid-template-columns:1fr; } }
</style></head>
<body><main>
<header><div class="kicker">TS Second Opinion · report card</div><h1>__TITLE__</h1><div class="smiles">__SMILES__</div></header>

<section class="verdict __VERDICT__"><div class="badge">__ICON__</div>
<div><div class="word">__WORD__</div><p>__MEANING__</p><ul class="facts">__FACTS__</ul></div></section>

<section><h2>Why</h2><ul>__REASONS__</ul></section>

<section><h2>Raw guess and recommended guess <span class="tag">drag to rotate; the two views move together</span></h2>
<div class="views">
  <div class="view"><div class="cap">__RAWLABEL__</div><div id="v1" class="viewer"></div></div>
  <div class="view"><div class="cap">__RECLABEL__</div><div id="v2" class="viewer"></div></div>
</div>
<div class="legend">Stereocentres:<span class="dot" style="background:#1e8e3e"></span>right handedness
<span class="dot" style="background:#d93025"></span>inverted<span class="dot" style="background:#9aa0a6"></span>too flat to judge</div>
</section>

<section><h2>Stereocentres</h2>
<table><tr><th>Atom</th><th>Compared with</th><th>Raw guess</th><th>Recommended</th></tr>__STEREO__</table></section>

<section><h2>Rules on the recommended guess</h2>
<table><tr><th>Rule</th><th>Result</th><th>Detail</th></tr>__RULES__</table>
<p class="muted" style="font-size:13px;margin:10px 0 0">Each rule was written by Claude from a plain-English sentence and kept only after
accepting at least 99% of 9,540 real DFT transition states (RDB7: H, C, N, O). For very different chemistry, test the rules
on your own reference structures with <code>tsso gate</code>.</p></section>
__SCORE__
<section><h2>Next step</h2>
<details><summary>Recommended TS geometry (.xyz)</summary><pre>__RECXYZ__</pre></details>
<details><summary>ORCA input: ωB97X-D4/def2-TZVP OptTS + Freq</summary><pre>__ORCA__</pre></details></section>

<footer>TS Second Opinion __VERSION__ · __DATE__ · github.com/tammy357/ts-second-opinion · Results are screening
advice, not a substitute for a TS optimisation. Please cite TS Second Opinion, GoFlow and RDB7.</footer>
</main>
<script>__JS__</script>
<script>
const D = __DATA__;
const dark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
function show(id, xyz, centres) {
  const v = $3Dmol.createViewer(id, {backgroundColor: dark ? '#1d1d1b' : '#ffffff'});
  v.addModel(xyz, 'xyz');
  v.setStyle({}, {stick: {radius: 0.12}, sphere: {scale: 0.2}});
  for (const c of centres) {
    const col = c.ok === true ? '#1e8e3e' : c.ok === false ? '#d93025' : '#9aa0a6';
    v.addSphere({center: c.pos, radius: 0.6, color: col, opacity: 0.45});
    v.addLabel(c.text, {position: c.pos, backgroundColor: col, fontColor: 'white', fontSize: 12, backgroundOpacity: 0.9});
  }
  v.zoomTo(); v.zoom(1.5); v.render(); return v;
}
const a = show('v1', D.raw, D.craw), b = show('v2', D.rec, D.crec);
if (a.linkViewer) { a.linkViewer(b); b.linkViewer(a); }
</script>
</body></html>
"""


def write_report(path, *args, **kwargs):
    Path(path).write_text(report_html(*args, **kwargs))
    return str(path)
