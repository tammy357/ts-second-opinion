"""Build the static showcase site in docs/ (served by GitHub Pages): the five worked examples with verdicts, linked 3D
views and report cards, the Rule Studio replays, the trust slider and the leaderboard. Everything is precomputed from
the same code as the app, so the page needs no server. Checking your own reaction needs the package (see the README).

Usage: python scripts/build_site.py   -> docs/index.html, docs/assets/, docs/reports/
"""
import json, shutil, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.geometry import dmae
from tsso.metrics import handedness_correct
from tsso.report import report_html
from tsso.sandbox import run_check
from tsso.stereo import stereocentre_status
from tsso.triage import load_checks, triage, ORCA_TEMPLATE
from tsso.viewer import align, COLOURS

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / 'docs'
EXAMPLES = {
    3509: ('Every model got the handedness wrong',
           'All 8 GoFlow guesses put C2 the wrong way round. Their interatomic distances are almost perfect, so D-MAE and '
           'chirality-blind RMSD call them accurate, but each one is the TS of the mirror-image reaction.'),
    5877: ('Only part of the molecule is flipped',
           'One of three stereocentres is inverted: the guess is the TS of a diastereomer, whose energy really differs. '
           'Mirroring would invert the two correct centres, so it cannot be repaired.'),
    1292: ('The models disagree',
           'Two guesses pass every rule, but the models disagree a lot about the geometry. This is where a DFT check pays off.'),
    2577: ('A guess you can trust', 'Every guess passes every rule, both stereocentres are right and the models agree closely.'),
    1093: ('A problem in the dataset',
           'Every guess breaks the reacting-bond rule. Looking closer, the atom-mapped SMILES says the N5–O6 bond breaks, '
           'but the reactant, TS and product geometries all keep it near 1.4 Å, so this RDB7 entry looks inconsistent.'),
}
WORD = {'green': ('✓', 'Trust', 'Passes every rule and the models agree.'),
        'amber': ('!', 'Check with DFT', 'Passes the rules, but the models disagree a lot.'),
        'red': ('✕', 'Re-generate', 'Still breaks a rule after repair.')}


def xyz(symbols, x, comment=''):
    return '\n'.join([str(len(symbols)), comment] + [f'{s} {a:.5f} {b:.5f} {c:.5f}' for s, (a, b, c) in zip(symbols, x)])


def centres(status, x):
    return [{'pos': {'x': float(x[s['atom']][0]), 'y': float(x[s['atom']][1]), 'z': float(x[s['atom']][2])},
             'color': COLOURS[s['ok']], 'text': s['label']} for s in status]


def main():
    (DOCS / 'assets').mkdir(parents=True, exist_ok=True); (DOCS / 'reports').mkdir(exist_ok=True)
    for f in ('3Dmol-min.js', '3Dmol-LICENSE.txt'):
        shutil.copy(ROOT / 'tsso' / 'assets' / f, DOCS / 'assets' / f)
    for f in ('fig_hero_mirror.png', 'fig_final_answer.png', 'fig_certify.png', 'fig_risk_coverage.png'):
        shutil.copy(ROOT / 'results' / f, DOCS / 'assets' / f)
    (DOCS / '.nojekyll').write_text('')

    R = {r.id: r for r in load_reactions()}
    checks = load_checks()
    cut = json.load(open(ROOT / 'results' / 'triage_summary.json'))['amber_cutoff_A']
    expl = json.load(open(ROOT / 'results' / 'explanations.json'))['examples']
    ex = []
    for rid, (title, note) in EXAMPLES.items():
        r = R[rid]; t = triage(r, checks, cut)
        raw = r.preds[t.chosen_model]
        st_raw, st_rec = stereocentre_status(r, raw), stereocentre_status(r, t.xyz)
        rules = []
        for name, fn in checks.items():
            ok, detail = run_check(fn, r, t.xyz)
            if ok and str(detail).startswith('rule not applicable'):
                ok, detail = 'na', 'nothing in this reaction for this rule to check'
            rules.append({'name': name, 'ok': ok, 'detail': str(detail)})
        score = {'dmae_A': dmae(t.xyz, r.ts_ref_xyz), 'stereo_correct': handedness_correct(r, t.xyz, r.ts_ref_xyz)}
        (DOCS / 'reports' / f'rdb7_{rid}.html').write_text(report_html(r, t, checks, 0, 1, f'RDB7 reaction {rid}', score))
        n = len(r.preds)
        icon, word, meaning = WORD[t.verdict]
        src = f"GoFlow guess {t.chosen_model + 1}{' (mirror image)' if t.mirrored else ''}; verdict {t.verdict}"
        ex.append({
            'id': rid, 'title': title, 'note': note, 'verdict': t.verdict, 'icon': icon, 'word': word, 'meaning': meaning,
            'facts': [f'{n - len(t.rejected_models)} of {n} raw guesses pass every rule',
                      'mirror-image repair applied' if t.mirrored else ('no repair needed' if t.verdict != 'red' else 'cannot be repaired'),
                      f'model disagreement {t.disagreement:.3f} Å (cut-off {cut:.3f} Å)',
                      f'vs reference (evaluation only): D-MAE {score["dmae_A"]:.3f} Å, handedness {"right" if score["stereo_correct"] else "wrong"}'],
            'reasons': t.reasons, 'explanation': expl[str(rid)]['explanation'], 'smiles': r.smiles,
            'raw_label': f'Best raw guess (GoFlow model {t.chosen_model + 1})',
            'rec_label': 'Recommended (mirror-repaired)' if t.mirrored else 'Recommended',
            'raw': xyz(r.symbols, raw), 'rec': xyz(r.symbols, t.xyz, f'TS Second Opinion recommended guess, RDB7 {rid}'),
            'ref': xyz(r.symbols, align(r.ts_ref_xyz, t.xyz)),
            'craw': centres(st_raw, raw), 'crec': centres(st_rec, t.xyz),
            'stereo': [{'label': a['label'], 'vs': a['compared_with'], 'raw': a['ok'], 'rec': b['ok']} for a, b in zip(st_raw, st_rec)],
            'rules': rules,
            'orca': ORCA_TEMPLATE.format(rid=f'RDB7 reaction {rid}', source=src, xyzname=f'rdb7_{rid}_ts.xyz', charge=0, mult=1),
        })
    cert = json.load(open(ROOT / 'results' / 'certify.json'))
    certd = {k: [{'t': d['target'], 'share': d['mean_share_trusted'], 'err': d['mean_error_among_trusted'],
                  'p95': d['error_among_trusted_95th_percentile']} for d in cert[k]['conformal']]
             for k in ('geometry', 'geometry_or_handedness')}
    replays = []
    for f in ('authoring_log.jsonl', 'authoring_log_run1.jsonl'):
        for line in open(ROOT / 'results' / f):
            e = json.loads(line)
            if e['name'] in ('h_transfer', 'planarity') and len(e['attempts']) > 1 and not any(x['name'] == e['name'] for x in replays):
                replays.append({'name': e['name'], 'rule': e['rule'], 'attempts': [
                    {k: a[k] for k in ('attempt', 'pass_rate', 'crash_rate', 'rejects_mirror_images', 'rejects_noisy_geometries',
                                       'accepted', 'code')} for a in e['attempts']]})
    replays.sort(key=lambda r: r['name'] != 'h_transfer')       # the clearest story first
    lead = json.load(open(ROOT / 'results' / 'leaderboard.json'))['rows']
    m = json.load(open(ROOT / 'results' / 'metrics.json'))
    fa = m['final_answer']
    stats = {'caught': m['stereo_checks_recall_of_stereo_errors'], 'hand_before': fa['single_model_avg']['stereo_bad_rate'],
             'hand_after': fa['check_aware_plus_mirror']['stereo_bad_rate'],
             'trusted5': next(d for d in cert['geometry']['conformal'] if d['target'] == 0.05)['mean_share_trusted']}
    data = json.dumps({'examples': ex, 'certify': certd, 'replays': replays, 'leaderboard': lead, 'stats': stats})
    (DOCS / 'index.html').write_text(PAGE.replace('__DATA__', data))
    print(f"wrote docs/index.html ({(DOCS / 'index.html').stat().st_size // 1024} KB), {len(ex)} report cards")


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>TS Second Opinion</title>
<meta name="description" content="Check, repair and triage machine-learned transition-state guesses before you spend DFT.">
<style>
:root { --bg:#f7f7f5; --surface:#ffffff; --ink:#151513; --ink2:#55534e; --line:#e3e2de; --soft:#f1f0ec; --accent:#2a6fd6;
  --green:#137333; --green-bg:#e6f4ea; --amber:#9a5b00; --amber-bg:#fdf1dc; --red:#b3261e; --red-bg:#fbe9e7; --viewer:#ffffff; }
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg:#131312; --surface:#1c1c1a; --ink:#f2f1ed; --ink2:#b4b2ab; --line:#33322e; --soft:#262522; --accent:#7fb0ff;
  --green:#7fd68f; --green-bg:#173a20; --amber:#f2c063; --amber-bg:#3d2e10; --red:#f28b82; --red-bg:#43201d; --viewer:#1c1c1a; } }
:root[data-theme="dark"] { --bg:#131312; --surface:#1c1c1a; --ink:#f2f1ed; --ink2:#b4b2ab; --line:#33322e; --soft:#262522;
  --accent:#7fb0ff; --green:#7fd68f; --green-bg:#173a20; --amber:#f2c063; --amber-bg:#3d2e10; --red:#f28b82; --red-bg:#43201d; --viewer:#1c1c1a; }
* { box-sizing:border-box; }
body { margin:0; background:var(--bg); color:var(--ink); font:16px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif; }
main { max-width:1080px; margin:0 auto; padding:28px 16px 64px; }
a { color:var(--accent); }
.kicker { font-size:13px; letter-spacing:.06em; text-transform:uppercase; color:var(--ink2); }
h1 { font-size:clamp(30px,5vw,44px); margin:6px 0 8px; letter-spacing:-.01em; }
h2 { font-size:24px; margin:44px 0 6px; } h3 { font-size:17px; margin:0 0 8px; }
.lede { font-size:18px; max-width:780px; margin:0 0 14px; }
.links { display:flex; gap:10px; flex-wrap:wrap; margin:6px 0 18px; }
.btn { display:inline-block; padding:8px 14px; border-radius:999px; border:1px solid var(--line); background:var(--surface);
  color:var(--ink); text-decoration:none; font-size:14px; cursor:pointer; }
.btn.primary { background:var(--ink); color:var(--bg); border-color:var(--ink); }
.tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:12px; }
.tile, .card { background:var(--surface); border:1px solid var(--line); border-radius:14px; padding:16px 18px; }
.tile .v { font-size:30px; font-weight:750; line-height:1.1; } .tile .l { color:var(--ink2); font-size:14px; margin-top:6px; }
.muted { color:var(--ink2); }
.pills { display:flex; gap:8px; flex-wrap:wrap; margin:14px 0; }
.pill { padding:8px 14px; border-radius:999px; border:1px solid var(--line); background:var(--surface); cursor:pointer; font-size:14px; color:var(--ink); }
.pill[aria-pressed="true"] { background:var(--ink); color:var(--bg); border-color:var(--ink); }
.verdict { display:flex; gap:16px; align-items:center; border-width:2px; }
.badge { flex:none; width:64px; height:64px; border-radius:14px; display:grid; place-items:center; font-size:32px; font-weight:800; }
.green { border-color:var(--green)!important; } .green .badge { background:var(--green-bg); color:var(--green); } .green .word { color:var(--green); }
.amber { border-color:var(--amber)!important; } .amber .badge { background:var(--amber-bg); color:var(--amber); } .amber .word { color:var(--amber); }
.red { border-color:var(--red)!important; } .red .badge { background:var(--red-bg); color:var(--red); } .red .word { color:var(--red); }
.word { font-size:26px; font-weight:750; }
.chips { display:flex; flex-wrap:wrap; gap:8px; margin-top:12px; }
.chip { background:var(--soft); border-radius:999px; padding:3px 10px; font-size:13px; color:var(--ink2); }
.grid2 { display:grid; grid-template-columns:1fr 1fr; gap:12px; }
.view { border:1px solid var(--line); border-radius:12px; overflow:hidden; background:var(--viewer); }
.view .cap { padding:8px 12px; border-bottom:1px solid var(--line); font-size:14px; color:var(--ink2); background:var(--surface); }
.viewer { width:100%; height:330px; position:relative; }
.legend { font-size:13px; color:var(--ink2); margin-top:8px; display:flex; flex-wrap:wrap; gap:4px 14px; align-items:center; }
.dot { display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:5px; vertical-align:middle; }
table { width:100%; border-collapse:collapse; font-size:14px; }
th, td { text-align:left; padding:7px 8px; border-bottom:1px solid var(--line); vertical-align:top; }
th { color:var(--ink2); font-weight:600; font-size:13px; } td.num, th.num { text-align:right; }
tr.ours td { background:var(--green-bg); }
.ok { color:var(--green); font-weight:600; } .bad { color:var(--red); font-weight:600; }
details { margin-top:8px; } summary { cursor:pointer; color:var(--ink2); }
pre { background:var(--soft); border-radius:10px; padding:12px; overflow:auto; font-size:12.5px; line-height:1.45; }
code { font-family:ui-monospace,Menlo,monospace; }
.explain { border-left:4px solid var(--accent); }
.attempt { border-left:4px solid var(--line); padding:6px 12px; margin:10px 0; }
.attempt.pass { border-color:var(--green); } .attempt.fail { border-color:var(--red); }
.big { font-size:44px; font-weight:750; line-height:1.1; margin:4px 0; }
input[type=range] { width:min(360px,100%); }
.fig { width:100%; border-radius:12px; border:1px solid var(--line); background:#fff; }
.stack > * + * { margin-top:12px; }
footer { margin-top:56px; color:var(--ink2); font-size:14px; }
@media (max-width:760px) { .grid2 { grid-template-columns:1fr; } .viewer { height:280px; } }
</style></head>
<body><main>
<div class="kicker">Open Scientific Intelligence Hackathon 2026</div>
<h1>TS Second Opinion</h1>
<p class="lede"><b>Check, repair and triage machine-learned transition-state guesses before you spend DFT.</b>
About one in six GoFlow guesses has a stereocentre the wrong way round, and most of those are perfect mirror images that pass
distance scores, chirality-blind RMSD and even a DFT frequency check. This page shows real results; the tool runs on any
generator's output, with any elements.</p>
<div class="links">
  <a class="btn primary" href="app/">Open the live app</a>
  <a class="btn" href="https://github.com/tammy357/ts-second-opinion">Code on GitHub</a>
  <a class="btn" href="#examples">Try the examples</a>
  <a class="btn" href="#own">Check your own reaction</a>
  <a class="btn" href="https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md">Detailed results</a>
</div>
<div class="tiles" id="tiles"></div>

<h2 id="examples">1 · Five real cases</h2>
<p class="muted">Each case is GoFlow's real prediction (8 models) for an RDB7 test reaction, checked by TS Second Opinion.
Drag a molecule to rotate it; the two views move together.</p>
<div class="pills" id="pills" role="group" aria-label="examples"></div>
<div class="stack">
  <div class="card verdict" id="verdict"><div class="badge" id="vicon"></div>
    <div><div class="muted" id="vtitle"></div><div class="word" id="vword"></div><div class="muted" id="vmeaning"></div>
      <div class="chips" id="vfacts"></div></div></div>
  <p id="vnote"></p>
  <div class="grid2">
    <div class="view"><div class="cap" id="rawcap"></div><div id="v1" class="viewer"></div></div>
    <div class="view"><div class="cap" id="reccap"></div><div id="v2" class="viewer"></div></div>
  </div>
  <div class="legend"><span><span class="dot" style="background:#1e8e3e"></span>stereocentre right</span>
    <span><span class="dot" style="background:#d93025"></span>inverted</span><span><span class="dot" style="background:#9aa0a6"></span>too flat to judge</span>
    <label><input type="checkbox" id="ref"> overlay the reference DFT TS (orange; never used for the verdict)</label></div>
  <div class="card"><h3>Why</h3><ul id="reasons"></ul></div>
  <div class="card explain"><h3>Claude's explanation</h3><p id="explain" style="margin:0"></p>
    <p class="muted" style="font-size:13px;margin:8px 0 0">Generated in advance with the Claude Code CLI from the facts above.</p></div>
  <details class="card"><summary>Stereocentres and every rule</summary><div id="tables"></div></details>
  <div class="links"><a class="btn primary" id="report" target="_blank">Open the full report card</a>
    <button class="btn" id="dlxyz">Download recommended TS (.xyz)</button><button class="btn" id="dlorca">Download ORCA input</button></div>
</div>

<h2 id="own">2 · Check your own reaction</h2>
<div class="card"><p style="margin-top:0">Any TS generator, any elements. You need an atom-mapped reaction SMILES (every atom mapped),
the reactant and product as <code>.xyz</code>, and one or more TS guesses.</p>
<div class="links" style="margin:0 0 12px"><a class="btn primary" href="app/">Upload them in the live app</a></div>
<p class="muted" style="margin-top:0">The live app is a marimo notebook running in your browser with Pyodide: nothing to install,
and your files never leave your computer. The first load takes a little while, as Python starts in your browser. To install it instead,
three commands, no GPU:</p>
<pre><code>git clone https://github.com/tammy357/ts-second-opinion &amp;&amp; cd ts-second-opinion
pip install -e .
marimo run apps/demo.py      # the app, with an upload tab
# or: tsso check --smiles rxn.smi --reactant r.xyz --product p.xyz guesses.xyz --report report.html</code></pre>
<p class="muted" style="margin-bottom:0">Agents can use it too: it is an MCP server with Claude Code commands
(<code>/second-opinion</code>, <code>/new-rule</code>).</p></div>

<h2>3 · Rule Studio: plain English in, tested code out</h2>
<p class="muted">A chemist writes a rule as a sentence. Claude writes the Python, and the gate keeps it only if it accepts at
least 99% of 9,540 real DFT transition states, never crashes, and catches at least 1% of mirrored or distorted ones.
Otherwise Claude gets the failures back and rewrites it. Two real runs:</p>
<div class="pills" id="rpills"></div><div class="card" id="replay"></div>

<h2>4 · How much can you trust?</h2>
<p class="muted">Conformal selection picks which reactions to trust so that, on average, at most the chosen share of them are
wrong. Calibrated on half of the 1,194 test reactions and tested on the other half, 2,000 times.</p>
<div class="card"><label>Target error rate: <b id="tval"></b><br><input type="range" id="target" min="0" max="9" step="1" value="3"></label>
  <div style="margin-top:8px"><label><input type="radio" name="count" value="geometry" checked> geometry errors (what matters for the barrier)</label>
  <label style="margin-left:12px"><input type="radio" name="count" value="geometry_or_handedness"> geometry or handedness</label></div>
  <div class="big" id="share"></div><div class="muted" id="sharetext"></div></div>

<h2>5 · Stereo-aware leaderboard</h2>
<p class="muted">"Success" usually means D-MAE ≤ 0.2 Å, which cannot see handedness. The last column counts only guesses that
are close <b>and</b> correctly handed. Handedness stays at 82–85% correct for every released GoFlow model size.</p>
<div class="card" style="overflow-x:auto"><table id="lead"></table></div>

<h2>6 · Results</h2>
<div class="grid2"><img class="fig" src="assets/fig_final_answer.png" alt="final answers wrong, by kind">
<img class="fig" src="assets/fig_certify.png" alt="certified trust"></div>
<p class="muted">Retrospective: GoFlow's released predictions for the 1,194 RDB7 test reactions; rules gated only on the 9,540
training reactions; no new DFT. <a href="https://arxiv.org/abs/2509.01543">FK-Flow</a> (from the GoFlow group) first measured
GoFlow's chirality errors and fixed them by steering generation; TS Second Opinion works afterwards, on any model's output.
Full results and limitations: <a href="https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md">docs/RESULTS.md</a>.</p>

<footer>Built with GoFlow's released predictions (MIT) and RDB7 (Spiekermann, Pattanaik &amp; Green, Sci. Data 2022). Rules
written by Claude and gated on real data. 3D viewer: 3Dmol.js (BSD-3). Code: MIT licence.
<a href="https://github.com/tammy357/ts-second-opinion">github.com/tammy357/ts-second-opinion</a></footer>
</main>
<script src="assets/3Dmol-min.js"></script>
<script>
const D = __DATA__;
const $ = id => document.getElementById(id);
const pc = (v, d = 0) => (100 * v).toFixed(d) + '%';
const el = (tag, cls, text) => { const e = document.createElement(tag); if (cls) e.className = cls; if (text !== undefined) e.textContent = text; return e; };
const dark = () => getComputedStyle(document.documentElement).getPropertyValue('--viewer').trim();

// tiles
[[pc(D.stats.caught), 'of wrong-handed GoFlow guesses caught by two LLM-written rules (0.2% false alarms)'],
 [pc(D.stats.hand_before) + ' → ' + pc(D.stats.hand_after), 'final answers with the wrong handedness, before → after checking and mirror repair'],
 [pc(D.stats.trusted5), 'of reactions trusted at a guaranteed 5% average geometry-error rate']].forEach(([v, l]) => {
  const t = el('div', 'tile'); t.append(el('div', 'v', v), el('div', 'l', l)); $('tiles').append(t); });

// viewers
const V1 = $3Dmol.createViewer('v1', {backgroundColor: dark() || 'white'});
const V2 = $3Dmol.createViewer('v2', {backgroundColor: dark() || 'white'});
if (V1.linkViewer) { V1.linkViewer(V2); V2.linkViewer(V1); }
function draw(v, xyz, cs, ref) {
  v.removeAllModels(); v.removeAllShapes(); v.removeAllLabels();
  v.addModel(xyz, 'xyz'); v.setStyle({model: 0}, {stick: {radius: 0.12}, sphere: {scale: 0.2}});
  if (ref) { v.addModel(ref, 'xyz'); v.setStyle({model: 1}, {stick: {radius: 0.06, color: '#eb6834', opacity: 0.75}}); }
  for (const c of cs) { v.addSphere({center: c.pos, radius: 0.6, color: c.color, opacity: 0.45});
    v.addLabel(c.text, {position: c.pos, backgroundColor: c.color, fontColor: 'white', fontSize: 12, backgroundOpacity: 0.9}); }
  v.zoomTo(); v.zoom(1.5); v.render();
}
let cur = 0;
const mark = ok => ok === true ? ['ok', '✓ right'] : ok === false ? ['bad', '✕ inverted'] : ['muted', 'too flat to judge'];
function show(i) {
  cur = i; const e = D.examples[i];
  [...$('pills').children].forEach((b, k) => b.setAttribute('aria-pressed', k === i));
  $('verdict').className = 'card verdict ' + e.verdict;
  $('vicon').textContent = e.icon; $('vtitle').textContent = 'RDB7 reaction ' + e.id; $('vword').textContent = e.word;
  $('vmeaning').textContent = e.meaning; $('vfacts').replaceChildren(...e.facts.map(f => el('span', 'chip', f)));
  $('vnote').textContent = e.note; $('rawcap').textContent = e.raw_label; $('reccap').textContent = e.rec_label;
  $('reasons').replaceChildren(...e.reasons.map(r => el('li', '', r)));
  $('explain').textContent = e.explanation;
  const st = el('table'); st.innerHTML = '<tr><th>Stereocentre</th><th>Compared with</th><th>Raw guess</th><th>Recommended</th></tr>';
  if (!e.stereo.length) { const tr = el('tr'); const td = el('td', 'muted', 'No stereocentres in this reaction.'); td.colSpan = 4; tr.append(td); st.append(tr); }
  e.stereo.forEach(s => { const tr = el('tr'); tr.append(el('td', '', s.label), el('td', '', s.vs));
    for (const ok of [s.raw, s.rec]) { const [c, t] = mark(ok); tr.append(el('td', c, t)); } st.append(tr); });
  const rt = el('table'); rt.innerHTML = '<tr><th>Rule</th><th>On the recommended guess</th><th>Detail</th></tr>';
  e.rules.forEach(r => { const tr = el('tr'); const [c, t] = r.ok === 'na' ? ['muted', '– n/a'] : r.ok ? ['ok', '✓ passes'] : ['bad', '✕ fails'];
    tr.append(el('td', '', r.name), el('td', c, t), el('td', '', r.detail)); rt.append(tr); });
  $('tables').replaceChildren(el('h3', '', 'Stereocentres'), st, el('h3', '', 'Rules'), rt);
  $('report').href = 'reports/rdb7_' + e.id + '.html';
  render();
}
function render() { const e = D.examples[cur]; draw(V1, e.raw, e.craw, null); draw(V2, e.rec, e.crec, $('ref').checked ? e.ref : null); }
D.examples.forEach((e, i) => { const b = el('button', 'pill', e.title); b.onclick = () => show(i); $('pills').append(b); });
$('ref').onchange = render;
function save(name, text) { const a = el('a'); a.href = URL.createObjectURL(new Blob([text], {type: 'text/plain'})); a.download = name; a.click(); }
$('dlxyz').onclick = () => save('rdb7_' + D.examples[cur].id + '_ts.xyz', D.examples[cur].rec + '\n');
$('dlorca').onclick = () => save('rdb7_' + D.examples[cur].id + '_ts.inp', D.examples[cur].orca);
show(0);

// rule studio
function replay(i) {
  [...$('rpills').children].forEach((b, k) => b.setAttribute('aria-pressed', k === i));
  const r = D.replays[i]; const box = $('replay'); box.replaceChildren(el('p', '', '“' + r.rule + '”'));
  r.attempts.forEach(a => {
    let v = 'Accepted: kept';
    if (!a.accepted) v = a.crash_rate > 0 ? `Rejected: crashed on ${(100 * a.crash_rate).toFixed(2)}% of structures`
      : a.pass_rate < 0.99 ? `Rejected: too strict, it called ${(100 * (1 - a.pass_rate)).toFixed(1)}% of real TSs wrong`
      : 'Rejected: it almost never fires, so it would add nothing';
    const d = el('div', 'attempt ' + (a.accepted ? 'pass' : 'fail'));
    d.append(el('b', '', 'Attempt ' + a.attempt), document.createTextNode(` · accepts ${(100 * a.pass_rate).toFixed(1)}% of real TSs (needs ≥ 99%) · catches ${(100 * a.rejects_mirror_images).toFixed(1)}% of mirror images and ${(100 * a.rejects_noisy_geometries).toFixed(1)}% of distorted TSs`),
      el('div', a.accepted ? 'ok' : 'bad', v));
    const det = el('details'); det.append(el('summary', '', 'Code Claude wrote')); const pre = el('pre'); pre.append(el('code', '', a.code)); det.append(pre);
    d.append(det); box.append(d); });
}
D.replays.forEach((r, i) => { const b = el('button', 'pill', r.name + ': ' + r.attempts.length + ' attempts'); b.onclick = () => replay(i); $('rpills').append(b); });
replay(0);

// trust slider
function trust() {
  const k = document.querySelector('input[name=count]:checked').value; const d = D.certify[k][+$('target').value];
  $('tval').textContent = Math.round(100 * d.t) + '%';
  $('share').textContent = pc(d.share) + ' of reactions';
  $('sharetext').textContent = `can be trusted; on average ${(100 * d.err).toFixed(1)}% of those were wrong (in 95% of splits at most ${(100 * d.p95).toFixed(1)}%). Recalibrate for a new generator or new chemistry.`;
}
$('target').oninput = trust; document.querySelectorAll('input[name=count]').forEach(r => r.onchange = trust); trust();

// leaderboard
const L = $('lead'); L.innerHTML = '<tr><th>Method</th><th>Test set</th><th class="num">Success (D-MAE ≤ 0.2 Å)</th><th class="num">Correct handedness</th><th class="num">Success with correct handedness</th></tr>';
D.leaderboard.forEach(r => { const tr = el('tr', r.method.includes('TS Second Opinion') ? 'ours' : '');
  const sd = r.sd; const c = k => pc(r[k], 1) + (sd ? ' ± ' + pc(sd[k], 1) : '');
  tr.append(el('td', '', r.method), el('td', '', r.test_set), el('td', 'num', c('success_dmae_only')), el('td', 'num', c('handedness_correct')));
  const last = el('td', 'num'); last.append(el('b', '', c('success_dmae_and_handedness'))); tr.append(last); L.append(tr); });
</script>
</body></html>
"""

if __name__ == '__main__':
    main()
