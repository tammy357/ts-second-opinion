"""Stereo-aware benchmark for TS generators.

1. Mirror blindness: the usual TS metrics give a mirror-image TS a perfect score.
   D-MAE compares interatomic distances (identical for mirror images), and the RMSD code shared by
   OA-ReactDiff, React-OT and GoFlow is called with ignore_chirality=True (min over the mirror image).
2. How many wrong-handed GoFlow guesses would count as "accurate" under the usual threshold.
3. Selective prediction: risk-coverage curves and AURC ("how wrong are the answers you choose to trust?").

Usage: python scripts/benchmark.py   -> results/benchmark.json, results/fig_risk_coverage.png
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.geometry import dmae, disagreement_matrix
from tsso.triage import load_checks, triage
from tsso.sandbox import run_check

R = load_reactions(); checks = load_checks(); K = R[0].preds.shape[0]
MIRROR = np.array([1.0, 1.0, -1.0]); BAD = 0.2


def kabsch_rmsd(a, b, allow_reflection=False):
    a0, b0 = a - a.mean(0), b - b.mean(0)
    U, S, Vt = np.linalg.svd(b0.T @ a0)
    d = np.sign(np.linalg.det(U @ Vt))
    if not allow_reflection and d < 0:
        U[:, -1] *= -1
    rot = U @ Vt
    return float(np.sqrt(((b0 @ rot - a0) ** 2).sum(1).mean()))


def stereo_centres(r):
    nbr, nbp = r.neighbours('reactant'), r.neighbours('product')
    c = [(i, nbr) for i in range(r.n) if r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(nbr[i]) >= 3]
    c += [(i, nbp) for i in range(r.n) if r.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    return c


def stereo_wrong(r, x):
    for i, nb in stereo_centres(r):
        v1 = r.ts_ref_xyz[nb[i][:3]] - r.ts_ref_xyz[i]; v2 = x[nb[i][:3]] - x[i]
        if np.sign(np.dot(v1[0], np.cross(v1[1], v1[2]))) != np.sign(np.dot(v2[0], np.cross(v2[1], v2[2]))):
            return True
    return False


out = {}
# 1. mirror blindness on the reference TSs themselves
chiral = [r for r in R if stereo_centres(r)]
mir = [(dmae(r.ts_ref_xyz * MIRROR, r.ts_ref_xyz), kabsch_rmsd(r.ts_ref_xyz * MIRROR, r.ts_ref_xyz, True),
        kabsch_rmsd(r.ts_ref_xyz * MIRROR, r.ts_ref_xyz, False)) for r in chiral]
mir = np.array(mir)
out['mirror_blindness'] = {
    'n_reactions_with_stereocentres': len(chiral),
    'mirror_of_reference_dmae_max_A': float(mir[:, 0].max()),
    'mirror_of_reference_rmsd_reflection_allowed_max_A': float(mir[:, 1].max()),
    'mirror_of_reference_rmsd_proper_rotation_median_A': float(np.median(mir[:, 2])),
}

# 2. wrong-handed GoFlow guesses that look accurate
E = np.array([[dmae(p, r.ts_ref_xyz) for p in r.preds] for r in R])
S = np.array([[stereo_wrong(r, p) for p in r.preds] for r in R])
out['wrong_handed_guesses'] = {
    'share_of_all_guesses': float(S.mean()),
    'share_of_wrong_handed_with_dmae_below_0.2A': float((E[S] <= BAD).mean()),
    'median_dmae_wrong_handed_A': float(np.median(E[S])), 'median_dmae_right_handed_A': float(np.median(E[~S])),
    'success_rate_dmae_only': float((E <= BAD).mean()),
    'success_rate_dmae_and_stereo': float(((E <= BAD) & ~S).mean()),
}
# Whole mirror image (every checked stereocentre flipped, so the mirror image has the right handedness) or a
# partial flip? A whole mirror image has exactly the same energy, forces and frequencies as its mirror image: it is
# a TS for the mirror-image reaction, so no energy-based check (xTB, MLIP or DFT) can tell it apart.
# A partial flip corresponds to a diastereomer, whose energy genuinely differs.
W = np.array([[bool(S[a, k]) and not stereo_wrong(r, p * MIRROR) for k, p in enumerate(r.preds)] for a, r in enumerate(R)])
out['wrong_handed_guesses'].update({
    'n_wrong_handed': int(S.sum()),
    'n_whole_mirror_images': int(W.sum()),
    'share_whole_mirror_images': float(W.sum() / S.sum()),
    'share_partial_flips': float(1 - W.sum() / S.sum()),
    'share_whole_mirror_with_dmae_below_0.2A': float((E[W] <= BAD).mean()),
})

# 3. risk-coverage: single model + disagreement vs full pipeline (rules + repair + disagreement)
cut = 0.10873146814426746
single_risk, single_score = [], []
for k in range(K):
    rk, sk = [], []
    for a, r in enumerate(R):
        M = disagreement_matrix(r.preds)
        rk.append(E[a, k] > BAD or S[a, k]); sk.append(M[k].sum() / (K - 1))
    single_risk.append(np.array(rk)); single_score.append(np.array(sk))
pipe_risk, pipe_geo, pipe_score = [], [], []
for r in R:
    t = triage(r, checks, cut)
    g = dmae(t.xyz, r.ts_ref_xyz) > BAD
    pipe_geo.append(g); pipe_risk.append(g or stereo_wrong(r, t.xyz))
    pipe_score.append(t.disagreement + (10.0 if t.failed_rules else 0.0))
pipe_risk, pipe_geo, pipe_score = np.array(pipe_risk), np.array(pipe_geo), np.array(pipe_score)
cov = np.linspace(0.05, 1.0, 96)


def curve(risk, score):
    o = np.argsort(score, kind='stable')          # most confident first
    r = risk[o].astype(float)
    return np.array([r[:max(1, int(round(c * len(r))))].mean() for c in cov])


single_curve = np.mean([curve(single_risk[k], single_score[k]) for k in range(K)], axis=0)
pipe_curve = curve(pipe_risk, pipe_score)
oracle_curve = curve(pipe_risk, pipe_risk.astype(float))
aurc = lambda c: float(np.trapezoid(c, cov) / (cov[-1] - cov[0]))
out['risk_coverage'] = {
    'aurc_single_model_with_disagreement': aurc(single_curve),
    'aurc_pipeline': aurc(pipe_curve), 'aurc_pipeline_oracle_order': aurc(oracle_curve),
    'risk_at_50pct_coverage': {'single_model_with_disagreement': float(single_curve[np.argmin(abs(cov - 0.5))]),
                               'pipeline': float(pipe_curve[np.argmin(abs(cov - 0.5))])},
    'risk_at_80pct_coverage': {'single_model_with_disagreement': float(single_curve[np.argmin(abs(cov - 0.8))]),
                               'pipeline': float(pipe_curve[np.argmin(abs(cov - 0.8))])},
    'risk_at_full_coverage': {'single_model': float(single_curve[-1]), 'pipeline': float(pipe_curve[-1])},
}
# The same, counting geometry errors only (D-MAE > 0.2 A), i.e. ignoring handedness, which a single model's
# disagreement cannot see. This is the fairer comparison for anything that only needs the barrier.
single_geo = np.mean([curve(E[:, k] > BAD, single_score[k]) for k in range(K)], axis=0)
pipe_geo_curve = curve(pipe_geo, pipe_score)
at = lambda c, x: float(c[np.argmin(abs(cov - x))])
out['risk_coverage_geometry_only'] = {
    'aurc_single_model_with_disagreement': aurc(single_geo), 'aurc_pipeline': aurc(pipe_geo_curve),
    'risk_at_50pct_coverage': {'single_model_with_disagreement': at(single_geo, 0.5), 'pipeline': at(pipe_geo_curve, 0.5)},
    'risk_at_80pct_coverage': {'single_model_with_disagreement': at(single_geo, 0.8), 'pipeline': at(pipe_geo_curve, 0.8)},
    'risk_at_full_coverage': {'single_model': float(single_geo[-1]), 'pipeline': float(pipe_geo_curve[-1])},
}
json.dump(out, open('results/benchmark.json', 'w'), indent=2)
print(json.dumps(out, indent=2))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), facecolor=SURF, sharey=True)
panels = [(axes[0], single_geo, pipe_geo_curve, None, 'Geometry errors only',
           'D-MAE > 0.2 Å (what matters for the barrier)'),
          (axes[1], single_curve, pipe_curve, oracle_curve, 'Geometry or handedness',
           'adds a stereocentre the wrong way round')]
for ax, sc, pc, oc, title, sub in panels:
    ax.set_facecolor(SURF)
    ax.plot(cov * 100, sc * 100, color='#2a78d6', lw=2)
    ax.plot(cov * 100, pc * 100, color='#eb6834', lw=2)
    if oc is not None:
        ax.plot(cov * 100, oc * 100, color=INK2, lw=1.2, ls=(0, (3, 3)))
    for s_ in ('top', 'right', 'left'):
        ax.spines[s_].set_visible(False)
    ax.spines['bottom'].set_color(GRID); ax.yaxis.grid(True, color=GRID, lw=0.8); ax.tick_params(length=0, colors=INK2)
    ax.set_xlim(5, 100); ax.set_ylim(0, 32)
    ax.set_xlabel('reactions you choose to trust (%)', color=INK2)
    ax.set_title(title, loc='left', color=INK, fontsize=11.5, fontweight='bold', pad=22)
    ax.text(0, 1.02, sub, transform=ax.transAxes, color=INK2, fontsize=9, va='bottom')
axes[0].set_ylabel('of those, share that are wrong (%)', color=INK2)
axes[0].legend(handles=[plt.Line2D([], [], color='#2a78d6', lw=2), plt.Line2D([], [], color='#eb6834', lw=2)],
               labels=['one model, ranked by disagreement', 'TS Second Opinion'], frameon=False, loc='upper left',
               labelcolor=INK2, fontsize=9)
axes[1].text(101, single_curve[-1] * 100, 'one model,\nranked by\ndisagreement', color=INK2, fontsize=8.5, va='center')
axes[1].text(101, pipe_curve[-1] * 100, 'TS Second\nOpinion', color=INK, fontsize=8.5, va='center', fontweight='bold')
axes[1].text(22, 1.0, 'dashed: perfect ranking', color=INK2, fontsize=8)
fig.suptitle('Knowing when to trust a TS guess (lower is better)', x=0.01, ha='left', color=INK, fontsize=13,
             fontweight='bold')
fig.text(0.01, 0.885, '1,194 RDB7 test reactions, most confident first. One model ranked by disagreement cannot see handedness errors.',
         color=INK2, fontsize=9.5)
fig.subplots_adjust(left=0.07, right=0.89, top=0.74, bottom=0.13, wspace=0.12)
fig.savefig('results/fig_risk_coverage.png', dpi=200, facecolor=SURF)
print('wrote results/fig_risk_coverage.png')
