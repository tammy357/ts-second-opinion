"""Does it matter who writes the check? Claude-authored checks (tsso/checks) vs drafts written directly as code
during development (tsso/checks_draft), for the same six rules.

Both versions go through the same gate on the 9,540 training references, then are compared on GoFlow's test guesses.
Usage: python scripts/ablation_authors.py  -> results/ablation_authors.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.author import evaluate_source, train_set
from tsso.data import load_reactions
from tsso.geometry import dmae, consensus_index
from tsso.metrics import handedness_correct
from tsso.sandbox import compile_check, run_check

R = load_reactions(); TRAIN = train_set(); MIRROR = np.array([1.0, 1.0, -1.0]); STEREO = ('chirality', 'product_stereo')
geo_bad = np.array([[dmae(p, r.ts_ref_xyz) > 0.2 for p in r.preds] for r in R])
hand_bad = np.array([[not handedness_correct(r, p, r.ts_ref_xyz) for p in r.preds] for r in R])
any_bad = geo_bad | hand_bad


def load(folder):
    return {p.stem: p.read_text() for p in sorted(Path(folder).glob('*.py'))}


def per_rule(src):
    g = evaluate_source(src, TRAIN)
    fn = compile_check(src)
    F = np.array([[run_check(fn, r, p)[0] is False for p in r.preds] for r in R])
    return {'gate_accepts_real_TS': g['pass_rate'], 'gate_rejects_mirror': g['rejects_mirror_images'],
            'gate_rejects_noisy': g['rejects_noisy_geometries'], 'gate_accepted': g['accepted'],
            'test_flag_rate': float(F.mean()), 'test_precision_any_bad': float(any_bad[F].mean()) if F.any() else None,
            'test_recall_handedness_errors': float(F[hand_bad].mean())}, F


def pipeline(flags):
    """Final answers when these checks pick the consensus among passing guesses and drive the mirror repair."""
    names = list(flags)
    geo, hand = [], []
    for a, r in enumerate(R):
        passing = ~np.any([flags[n][a] for n in names], axis=0)
        c = consensus_index(r.preds, allowed=passing)
        x = r.preds[c]
        failed = [n for n in names if flags[n][a][c]]
        if failed and set(failed) & set(STEREO):
            m = x * MIRROR
            if not any(run_check(fns[n], r, m)[0] is False for n in names):
                x = m
        geo.append(dmae(x, r.ts_ref_xyz) > 0.2); hand.append(not handedness_correct(r, x, r.ts_ref_xyz))
    return {'final_geometry_wrong': float(np.mean(geo)), 'final_handedness_wrong': float(np.mean(hand))}


out = {}
for label, folder in (('claude_authored', 'tsso/checks'), ('direct_drafts', 'tsso/checks_draft')):
    srcs = load(folder)
    fns = {n: compile_check(s) for n, s in srcs.items()}
    rules, flags = {}, {}
    for n, s in srcs.items():
        rules[n], flags[n] = per_rule(s)
    out[label] = {'rules': rules, 'pipeline': pipeline(flags)}
    print(label, json.dumps(out[label]['pipeline']), {n: (round(v['gate_accepts_real_TS'], 3), v['gate_accepted'],
          round(v['test_flag_rate'], 3), round(v['test_recall_handedness_errors'], 3)) for n, v in rules.items()})
json.dump(out, open('results/ablation_authors.json', 'w'), indent=2)
