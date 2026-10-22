"""Certified trust: how many TS guesses can be trusted at a guaranteed average error rate?

Conformal selection (tsso/certify.py): calibrate on reactions with known answers, then trust new reactions
chosen by Benjamini-Hochberg on conformal p-values. Tested by splitting the 1,194 RDB7 test reactions at random
into calibration and "new" halves 2,000 times.

Also:
  - a naive plug-in rule (trust everything scoring above the cut-off where the calibration error rate is still
    under the target), to show why the correction matters;
  - a deliberate change of models: calibrate on the 8-seed ensemble, then certify guesses from the 5 size runs
    ensemble on the other half of the reactions (the guarantee assumes nothing changed; does it survive?).

Usage: python scripts/certify.py  -> results/certify.json, results/fig_certify.png
"""
import csv, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.certify import certify, trust_score
from tsso.data import load_reactions
from tsso.geometry import dmae
from tsso.triage import load_checks, triage

BAD = 0.2
TARGETS = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20]
N_SPLITS = 2000
rng = np.random.default_rng(0)

rows = list(csv.DictReader(open('results/triage.csv')))
ids = np.array([int(r['reaction_id']) for r in rows])
score = np.array([trust_score(r['verdict'], float(r['disagreement_A'])) for r in rows])
wrong = {'geometry': np.array([r['geometry_bad'] == 'True' for r in rows]),
         'geometry_or_handedness': np.array([r['geometry_bad'] == 'True' or r['stereo_wrong'] == 'True' for r in rows])}
n = len(rows)


def naive(score_cal, wrong_cal, score_new, alpha):
    o = np.argsort(-score_cal, kind='stable')
    run = np.cumsum(wrong_cal[o]) / np.arange(1, len(o) + 1)
    ok = np.nonzero(run <= alpha)[0]
    if not len(ok):
        return np.zeros(len(score_new), bool)
    return score_new >= score_cal[o][ok.max()]


def repeated(method, w, alpha):
    err, share = [], []
    for _ in range(N_SPLITS):
        perm = rng.permutation(n); c, t = perm[:n // 2], perm[n // 2:]
        sel = method(score[c], w[c], score[t], alpha)
        err.append(w[t][sel].mean() if sel.any() else 0.0); share.append(sel.mean())
    err, share = np.array(err), np.array(share)
    return {'target': alpha, 'mean_error_among_trusted': float(err.mean()),
            'error_among_trusted_95th_percentile': float(np.quantile(err, 0.95)),
            'share_of_draws_above_target': float((err > alpha).mean()), 'mean_share_trusted': float(share.mean())}


out = {'method': 'conformal selection (Jin & Candes, JMLR 2023): BH on conformal p-values; trust score = -(model '
                 'disagreement + 10 if a rule is still broken)',
       'n_reactions': n, 'calibration_share': 0.5, 'n_random_splits': N_SPLITS,
       'unconditional_error_rate': {k: float(v.mean()) for k, v in wrong.items()}}
for k, w in wrong.items():
    out[k] = {'conformal': [repeated(certify, w, a) for a in TARGETS],
              'naive_plug_in': [repeated(naive, w, a) for a in TARGETS]}
    print(k, [(d['target'], round(d['mean_error_among_trusted'], 3), round(d['mean_share_trusted'], 3)) for d in out[k]['conformal']])

# ---- change of models: calibrate on the 8-seed ensemble, certify the 5-size ensemble (different reactions)
checks = load_checks()
cut = json.load(open('results/triage_summary.json'))['amber_cutoff_A']
S5 = {r.id: r for r in load_reactions('data/goflow_sizes_rdb7_test.json.gz')}
from tsso.service import _stereo_wrong as stereo_wrong  # same stereo test as scripts/triage_all.py
s5, w5 = [], {'geometry': [], 'geometry_or_handedness': []}
for rid in ids:
    r = S5[rid]
    t = triage(r, checks, cut)
    g = dmae(t.xyz, r.ts_ref_xyz) > BAD
    s5.append(trust_score(t.verdict, t.disagreement))
    w5['geometry'].append(g); w5['geometry_or_handedness'].append(g or stereo_wrong(r, t.xyz))
s5 = np.array(s5); w5 = {k: np.array(v) for k, v in w5.items()}
shift = {}
for k in wrong:
    res = {'calibrated_on_8_seeds': [], 'recalibrated_on_5_sizes': []}
    for a in (0.05, 0.10):
        e1, e2, sh1, sh2 = [], [], [], []
        for _ in range(N_SPLITS):
            perm = rng.permutation(n); c, t = perm[:n // 2], perm[n // 2:]
            sel = certify(score[c], wrong[k][c], s5[t], a); e1.append(w5[k][t][sel].mean() if sel.any() else 0.0); sh1.append(sel.mean())
            sel = certify(s5[c], w5[k][c], s5[t], a); e2.append(w5[k][t][sel].mean() if sel.any() else 0.0); sh2.append(sel.mean())
        res['calibrated_on_8_seeds'].append({'target': a, 'mean_error_among_trusted': float(np.mean(e1)), 'mean_share_trusted': float(np.mean(sh1))})
        res['recalibrated_on_5_sizes'].append({'target': a, 'mean_error_among_trusted': float(np.mean(e2)), 'mean_share_trusted': float(np.mean(sh2))})
    shift[k] = res
out['change_of_models'] = {'description': 'calibrate on the 8-seed ensemble (half of the reactions), certify the guesses '
                                          'of the 5 GoFlow size runs (four model sizes) on the other half; then recalibrate on that ensemble',
                           **shift}
print(json.dumps(out['change_of_models'], indent=1))
json.dump(out, open('results/certify.json', 'w'), indent=2)

# ---- figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE = '#2a78d6', '#eb6834'
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6), facecolor=SURF)
T = np.array(TARGETS) * 100
for ax in axes:
    ax.set_facecolor(SURF)
    for s_ in ('top', 'right', 'left'):
        ax.spines[s_].set_visible(False)
    ax.spines['bottom'].set_color(GRID); ax.yaxis.grid(True, color=GRID, lw=0.8); ax.tick_params(length=0, colors=INK2)
    ax.set_xlabel('target error rate among trusted reactions (%)', color=INK2)
    ax.set_xlim(1, 21); ax.set_xticks(range(0, 21, 5))
ax = axes[0]
ax.plot(T, T, color=INK2, lw=1.2, ls=(0, (3, 3)))
for k, c, lab in (('geometry', BLUE, 'geometry errors'), ('geometry_or_handedness', ORANGE, 'geometry or handedness')):
    conf = out[k]['conformal']
    ax.plot(T, [100 * d['mean_error_among_trusted'] for d in conf], color=c, lw=2, marker='o', ms=5, mec=SURF, mew=1.5, label=lab)
ax.set_ylim(0, 21); ax.set_yticks(range(0, 21, 5)); ax.set_ylabel('actual error rate among trusted (%)', color=INK2)
ax.text(14.5, 16.2, 'target', color=INK2, fontsize=8.5, rotation=0)
ax.legend(frameon=False, loc='upper left', labelcolor=INK2, fontsize=9)
ax.set_title('The guarantee holds', loc='left', color=INK, fontsize=11.5, fontweight='bold', pad=22)
ax.text(0, 1.02, 'average over 2,000 random calibration splits', transform=ax.transAxes, color=INK2, fontsize=9, va='bottom')
ax = axes[1]
for k, c in (('geometry', BLUE), ('geometry_or_handedness', ORANGE)):
    conf = out[k]['conformal']
    ax.plot(T, [100 * d['mean_share_trusted'] for d in conf], color=c, lw=2, marker='o', ms=5, mec=SURF, mew=1.5)
g5 = next(d for d in out['geometry']['conformal'] if d['target'] == 0.05)
ax.annotate(f"{100 * g5['mean_share_trusted']:.0f}% trusted at a 5% target (geometry errors)", xy=(5, 100 * g5['mean_share_trusted']),
            xytext=(1.6, 101), color=INK, fontsize=9, va='bottom', arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
ax.set_ylim(0, 112); ax.set_yticks(range(0, 101, 20)); ax.set_ylabel('reactions trusted (%)', color=INK2)
ax.set_title('How much you can trust', loc='left', color=INK, fontsize=11.5, fontweight='bold', pad=22)
ax.text(0, 1.02, 'share of new reactions certified', transform=ax.transAxes, color=INK2, fontsize=9, va='bottom')
fig.suptitle('Certified trust: a guaranteed average error rate among the TS guesses you trust', x=0.01, ha='left',
             color=INK, fontsize=13, fontweight='bold')
fig.text(0.01, 0.885, 'Conformal selection on 1,194 RDB7 test reactions (8 GoFlow models), half used for calibration. No DFT at decision time.',
         color=INK2, fontsize=9.5)
fig.subplots_adjust(left=0.07, right=0.98, top=0.74, bottom=0.13, wspace=0.22)
fig.savefig('results/fig_certify.png', dpi=200, facecolor=SURF)
print('wrote results/fig_certify.png')
