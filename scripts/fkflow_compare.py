"""Compare FK-Flow's published predictions with and without TS Second Opinion, on the same RDB7 test reactions.

    python scripts/fkflow_compare.py fkflow_export.zip   -> results/fkflow_comparison.json, results/fkflow_comparison.md

Each FK-Flow run gives one guess per reaction, so TS Second Opinion runs in single-guess mode (rules + mirror repair, no
model-disagreement signal). Reactant, product and reference geometries come from RDB7 via data/goflow8_rdb7_test.json.gz.
"""
import copy, gzip, io, json, sys, zipfile
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.metrics import evaluate
from tsso.triage import load_checks, triage

KEYS = ('mean_dmae_A', 'success_dmae_only', 'handedness_correct', 'success_dmae_and_handedness')


def runs(export):
    z = zipfile.ZipFile(export)
    for n in sorted(z.namelist()):
        if Path(n).name.startswith('preds_') and n.endswith('.json.gz'):
            yield Path(n).name, json.load(gzip.open(io.BytesIO(z.read(n)), 'rt'))


def main(export):
    R = {r.id: r for r in load_reactions()}
    checks = load_checks()
    cut = json.load(open('results/triage_summary.json'))['amber_cutoff_A']
    rows = []
    for name, d in runs(export):
        ids = [int(k) for k in d['preds'] if int(k) in R]
        same_smiles = sum(d.get('smiles', {}).get(str(i)) == R[i].smiles for i in ids)
        Rs, raw, fixed, repaired = [], [], [], 0
        for i in ids:
            x = np.array(d['preds'][str(i)])
            if x.shape != R[i].ts_ref_xyz.shape:
                continue
            r = copy.copy(R[i]); r.preds = x[None]
            t = triage(r, checks, cut)
            Rs.append(r); raw.append(x); fixed.append(t.xyz); repaired += t.mirrored
        if not Rs:
            print(f'{name}: no reactions matched the RDB7 test set'); continue
        a = evaluate(raw, [r.ts_ref_xyz for r in Rs], Rs); b = evaluate(fixed, [r.ts_ref_xyz for r in Rs], Rs)
        label = f"FK-Flow release: {d['method']}{' + ' + d['potentials'] if d.get('potentials') else ''}"
        rows.append({'run': name, 'method': label, 'n': len(Rs), 'smiles_match': same_smiles, 'mirror_repaired': repaired,
                     'raw': {k: a[k] for k in KEYS}, 'with_ts_second_opinion': {k: b[k] for k in KEYS}})
        print(f"{label} [{name}]: n={len(Rs)} (SMILES match {same_smiles}); success with correct handedness "
              f"{100 * a['success_dmae_and_handedness']:.1f}% -> {100 * b['success_dmae_and_handedness']:.1f}% "
              f"(handedness {100 * a['handedness_correct']:.1f}% -> {100 * b['handedness_correct']:.1f}%, {repaired} repaired)")
    json.dump({'source': 'FK-Flow, Zenodo 10.5281/zenodo.23416442 (CC BY 4.0)', 'rows': rows},
              open('results/fkflow_comparison.json', 'w'), indent=2)
    pc = lambda v: f'{100 * v:.1f}%'
    md = ['| Run | Reactions | Success (D-MAE ≤ 0.2 Å) | Correct handedness | Success with correct handedness |', '|---|---|---|---|---|']
    for r in rows:
        md.append(f"| {r['method']} | {r['n']:,} | {pc(r['raw']['success_dmae_only'])} | {pc(r['raw']['handedness_correct'])} | "
                  f"{pc(r['raw']['success_dmae_and_handedness'])} |")
        md.append(f"| ↳ + TS Second Opinion (rules + mirror repair) | {r['n']:,} | {pc(r['with_ts_second_opinion']['success_dmae_only'])} | "
                  f"{pc(r['with_ts_second_opinion']['handedness_correct'])} | **{pc(r['with_ts_second_opinion']['success_dmae_and_handedness'])}** |")
    Path('results/fkflow_comparison.md').write_text('\n'.join(md) + '\n')
    print('\n'.join(md))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'fkflow_export.zip')
