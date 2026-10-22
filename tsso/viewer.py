"""Side-by-side, linked 3D views of the raw and recommended TS guesses with stereocentres highlighted (3Dmol.js,
inlined so it works offline and inside an iframe)."""
import json
from pathlib import Path

import numpy as np

ASSETS = Path(__file__).resolve().parent / 'assets'
COLOURS = {True: '#1e8e3e', False: '#d93025', None: '#9aa0a6'}


def align(mobile, target):
    """Proper-rotation Kabsch alignment of mobile onto target (no reflection)."""
    a, b = mobile - mobile.mean(0), target - target.mean(0)
    U, S, Vt = np.linalg.svd(a.T @ b)
    if np.linalg.det(U @ Vt) < 0:
        U[:, -1] *= -1
    return a @ (U @ Vt) + target.mean(0)


def _xyz(symbols, x):
    return '\n'.join([str(len(symbols)), ''] + [f'{s} {a:.5f} {b:.5f} {c:.5f}' for s, (a, b, c) in zip(symbols, x)])


def _centres(status, x):
    return [{'pos': {'x': float(x[s['atom']][0]), 'y': float(x[s['atom']][1]), 'z': float(x[s['atom']][2])},
             'color': COLOURS[s['ok']], 'text': s['label']} for s in status]


def compare_html(rxn, raw, rec, st_raw, st_rec, labels=('Best raw guess', 'Recommended'), reference=None, height=320):
    data = json.dumps({'raw': _xyz(rxn.symbols, raw), 'rec': _xyz(rxn.symbols, rec), 'craw': _centres(st_raw, raw),
                       'crec': _centres(st_rec, rec),
                       'ref': None if reference is None else _xyz(rxn.symbols, align(reference, rec))})
    js = (ASSETS / '3Dmol-min.js').read_text()
    return f"""<!doctype html><html><head><meta charset="utf-8"><style>
body {{ margin:0; font:13px system-ui,-apple-system,sans-serif; color:#55534e; background:transparent; }}
.grid {{ display:grid; grid-template-columns:1fr 1fr; gap:10px; }}
.cell {{ border:1px solid #e3e2de; border-radius:10px; overflow:hidden; background:#fff; }}
.cap {{ padding:7px 10px; border-bottom:1px solid #e3e2de; }}
.v {{ width:100%; height:{height}px; position:relative; }}
.legend {{ margin-top:8px; }} .dot {{ display:inline-block; width:10px; height:10px; border-radius:50%; margin:0 4px 0 12px; vertical-align:middle; }}
@media (max-width:640px) {{ .grid {{ grid-template-columns:1fr; }} }}
</style></head><body>
<div class="grid"><div class="cell"><div class="cap">{labels[0]}</div><div id="a" class="v"></div></div>
<div class="cell"><div class="cap">{labels[1]}</div><div id="b" class="v"></div></div></div>
<div class="legend">Stereocentres:<span class="dot" style="background:#1e8e3e"></span>right handedness
<span class="dot" style="background:#d93025"></span>inverted<span class="dot" style="background:#9aa0a6"></span>too flat to judge
{'<span class="dot" style="background:#eb6834"></span>reference DFT TS (thin, orange)' if reference is not None else ''}
· drag to rotate; the two views move together</div>
<script>{js}</script><script>
const D = {data};
function show(id, xyz, cs, ref) {{
  const v = $3Dmol.createViewer(id, {{backgroundColor: '#ffffff'}});
  v.addModel(xyz, 'xyz'); v.setStyle({{model: 0}}, {{stick: {{radius: 0.12}}, sphere: {{scale: 0.2}}}});
  if (ref) {{ v.addModel(ref, 'xyz'); v.setStyle({{model: 1}}, {{stick: {{radius: 0.06, color: '#eb6834', opacity: 0.7}}}}); }}
  for (const c of cs) {{ v.addSphere({{center: c.pos, radius: 0.6, color: c.color, opacity: 0.45}});
    v.addLabel(c.text, {{position: c.pos, backgroundColor: c.color, fontColor: 'white', fontSize: 12, backgroundOpacity: 0.9}}); }}
  v.zoomTo(); v.zoom(1.5); v.render(); return v;
}}
const A = show('a', D.raw, D.craw, null), B = show('b', D.rec, D.crec, D.ref);
if (A.linkViewer) {{ A.linkViewer(B); B.linkViewer(A); }}
</script></body></html>"""
