"""Does the approach hold beyond one setup?  (no new training: all predictions were released by the GoFlow authors)

A. Five GoFlow runs of four model sizes (n_atom_basis 64, 128, 256, 256, 512; random split) as a heterogeneous ensemble.
B. One GoFlow model on the harder reaction-core and barrier-height splits (rules only, no ensemble).
C. A naive guess with no ML at all: the midpoint of the aligned reactant and product.

Usage: python scripts/generality.py  -> results/generality.json
"""
import json, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.geometry import dmae, disagreement_matrix, consensus_index
from tsso.triage import load_checks
from tsso.sandbox import run_check

checks = load_checks(); BAD = 0.2; MIRROR = np.array([1.0, 1.0, -1.0]); STEREO = ('chirality', 'product_stereo')


def stereo_wrong(r, x):
    nbr, nbp = r.neighbours('reactant'), r.neighbours('product')
    cents = [(i, nbr) for i in range(r.n) if r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(nbr[i]) >= 3]
    cents += [(i, nbp) for i in range(r.n) if r.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    for i, nb in cents:
        v1 = r.ts_ref_xyz[nb[i][:3]] - r.ts_ref_xyz[i]; v2 = x[nb[i][:3]] - x[i]
        if np.sign(np.dot(v1[0], np.cross(v1[1], v1[2]))) != np.sign(np.dot(v2[0], np.cross(v2[1], v2[2]))):
            return True
    return False


def fails(r, x):
    return [n for n, fn in checks.items() if run_check(fn, r, x)[0] is False]


def repair(r, x):
    f = fails(r, x)
    if f and any(n in STEREO for n in f):
        m = x * MIRROR
        if not fails(r, m):
            return m, []
    return x, f


def bad(r, x):
    return dmae(x, r.ts_ref_xyz) > BAD or stereo_wrong(r, x)


def ensemble_report(R):
    K = R[0].preds.shape[0]
    E = np.array([[dmae(p, r.ts_ref_xyz) for p in r.preds] for r in R])
    S = np.array([[stereo_wrong(r, p) for p in r.preds] for r in R])
    F = np.array([[bool(set(fails(r, p)) & set(STEREO)) for p in r.preds] for r in R])
    U = np.array([disagreement_matrix(r.preds).sum(1) / (K - 1) for r in R])
    anyF = np.array([[bool(fails(r, p)) for p in r.preds] for r in R])
    any_bad = (E > BAD) | S
    final = []
    for a, r in enumerate(R):
        c = consensus_index(r.preds, allowed=~anyF[a]); x, _ = repair(r, r.preds[c]); final.append(bad(r, x))
    au_d = np.mean([roc_auc_score(any_bad[:, k], U[:, k]) for k in range(K)])
    au_c = np.mean([roc_auc_score(any_bad[:, k], U[:, k] + 10 * anyF[:, k]) for k in range(K)])
    return {'n_reactions': len(R), 'n_models': K, 'per_model_dmae_A': [round(float(v), 3) for v in E.mean(0)],
            'per_model_stereo_wrong': [round(float(v), 3) for v in S.mean(0)],
            'stereo_rule_recall': float(F[S].mean()), 'stereo_rule_false_alarm': float(F[~S].mean()),
            'wrong_single_model_avg': float(any_bad.mean()), 'wrong_after_pipeline': float(np.mean(final)),
            'auroc_any_bad': {'disagreement': float(au_d), 'disagreement+rules': float(au_c)}}


def single_report(R):
    E = np.array([dmae(r.preds[0], r.ts_ref_xyz) for r in R]); S = np.array([stereo_wrong(r, r.preds[0]) for r in R])
    F = np.array([bool(set(fails(r, r.preds[0])) & set(STEREO)) for r in R])
    after, flagged_after = [], []
    for r in R:
        x, f = repair(r, r.preds[0]); after.append(bad(r, x)); flagged_after.append(bool(f))
    after, flagged_after = np.array(after), np.array(flagged_after)
    ref_pass = np.mean([not fails(r, r.ts_ref_xyz) for r in R])
    return {'n_reactions': len(R), 'dmae_A': float(E.mean()), 'stereo_wrong': float(S.mean()),
            'stereo_rule_recall': float(F[S].mean()) if S.any() else None, 'stereo_rule_false_alarm': float(F[~S].mean()),
            'wrong_before': float(((E > BAD) | S).mean()), 'wrong_after_rules_and_mirror': float(after.mean()),
            'still_flagged_after_repair': float(flagged_after.mean()),
            'reference_TS_pass_all_rules': float(ref_pass)}


def kabsch_align(mob, ref):
    a, b = mob - mob.mean(0), ref - ref.mean(0)
    U, S, Vt = np.linalg.svd(a.T @ b)
    if np.linalg.det(U @ Vt) < 0:
        U[:, -1] *= -1
    return a @ (U @ Vt) + ref.mean(0)


out = {}
out['A_five_model_sizes'] = ensemble_report(load_reactions('data/goflow_sizes_rdb7_test.json.gz'))
out['B_reaction_core_split'] = single_report(load_reactions('data/goflow_split0_test.json.gz'))
out['B_barrier_split'] = single_report(load_reactions('data/goflow_split1_test.json.gz'))
R = load_reactions()
mid = []
for r in R:
    x = 0.5 * (r.reactant_xyz + kabsch_align(r.product_xyz, r.reactant_xyz))
    f = fails(r, x)
    mid.append((dmae(x, r.ts_ref_xyz), stereo_wrong(r, x), bool(f), bad(r, x)))
mid = np.array(mid, dtype=float)
out['C_naive_midpoint'] = {'dmae_A': float(mid[:, 0].mean()), 'stereo_wrong': float(mid[:, 1].mean()),
                           'flagged_by_rules': float(mid[:, 2].mean()), 'wrong': float(mid[:, 3].mean()),
                           'wrong_among_flagged': float(mid[mid[:, 2] == 1, 3].mean()),
                           'wrong_among_unflagged': float(mid[mid[:, 2] == 0, 3].mean())}
json.dump(out, open('results/generality.json', 'w'), indent=2)
print(json.dumps(out, indent=2))
