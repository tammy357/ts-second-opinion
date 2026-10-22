"""Build the in-browser (WebAssembly) app for GitHub Pages: docs/app/index.html + docs/app/bundle.zip (about 1 MB).

The page runs apps/demo.py with Pyodide, so nothing is installed and no server is needed. The bundle carries the tsso
package, the six gated rules, the worked examples and the results the app shows; RDKit is replaced by the pure-Python
reader in tsso/smiles_lite.py, which tests/test_smiles_lite.py checks against RDKit.

Usage: python scripts/build_wasm.py [--execute]     (--execute embeds a pre-rendered preview; needs uv)
"""
import gzip, io, json, shutil, subprocess, sys, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'docs' / 'app'
EXAMPLES = (3509, 5877, 1292, 2577, 1093)
FILES = (['results/' + f for f in ('metrics.json', 'triage_summary.json', 'certify.json', 'explanations.json',
                                   'leaderboard.json', 'authoring_log.jsonl', 'authoring_log_run1.jsonl',
                                   'fig_final_answer.png', 'fig_certify.png', 'fig_risk_coverage.png')]
         + ['examples/rdb7_3509/' + f for f in ('reaction.smi', 'reactant.xyz', 'product.xyz', 'guesses.xyz')])


def bundle():
    buf = io.BytesIO()
    stamp = (2026, 1, 1, 0, 0, 0)                       # fixed timestamps: the zip only changes when its contents do

    def add(name, data):
        zi = zipfile.ZipInfo(name, stamp)
        zi.compress_type = zipfile.ZIP_DEFLATED
        zf.writestr(zi, data)

    with zipfile.ZipFile(buf, 'w') as zf:
        for p in sorted((ROOT / 'tsso').rglob('*')):
            rel = p.relative_to(ROOT).as_posix()
            if p.is_file() and '__pycache__' not in rel and not rel.startswith(('tsso/checks_draft', 'tsso/checks_discovered')):
                add(rel, p.read_bytes())
        for f in FILES:
            add(f, (ROOT / f).read_bytes())
        full = json.load(gzip.open(ROOT / 'data' / 'goflow8_rdb7_test.json.gz', 'rt'))
        full['reactions'] = [r for r in full['reactions'] if r['id'] in EXAMPLES]
        full['description'] = 'The worked examples only (for the in-browser app). ' + full.get('description', '')
        add('data/app_examples.json.gz', gzip.compress(json.dumps(full).encode(), mtime=0))
    return buf.getvalue()


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    cmd = ['marimo', 'export', 'html-wasm', str(ROOT / 'apps' / 'demo.py'), '-o', str(OUT / 'index.html'), '--mode', 'run',
           '--single-file', '-f']                     # one HTML file; marimo and Python load from public CDNs
    if '--execute' in sys.argv:
        cmd.append('--execute')
    subprocess.run(cmd, check=True)
    data = bundle()
    (OUT / 'bundle.zip').write_bytes(data)
    total = sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
    print(f'wrote {OUT.relative_to(ROOT)}: bundle {len(data) / 1e6:.1f} MB, site {total / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
