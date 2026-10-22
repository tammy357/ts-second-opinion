"""Physics check with GFN2-xTB (tblite): does a cheap quantum-chemistry single point add anything?

For every GoFlow guess we compute the GFN2-xTB energy relative to the reactant ("xTB barrier") and the
RMS gradient (a true TS is a stationary point, so its forces should be small). We test them
(a) as signals for flagging bad guesses and (b) for choosing among the 8 guesses. No fitting: signals are
combined by averaging ranks.

Usage: python scripts/xtb_check.py  -> results/xtb.json
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score
from tblite.interface import Calculator
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.geometry import dmae, disagreement_matrix, consensus_index

ANG = 1 / 0.52917721067; HARTREE_KCAL = 627.509
R = load_reactions(); K = R[0].preds.shape[0]; BAD = 0.2


def xtb(atoms, xyz):
    calc = Calculator('GFN2-xTB', atoms, xyz * ANG)
    calc.set('verbosity', 0); calc.set('max-iter', 300)
    try:
        res = calc.singlepoint()
        g = res.get('gradient')
        return res.get('energy') * HARTREE_KCAL, float(np.sqrt((g ** 2).sum(1).mean())) * HARTREE_KCAL * ANG
    except Exception:
        return np.nan, np.nan


t0 = time.time()
Eb = np.zeros((len(R), K)); G = np.zeros((len(R), K)); ref_b = np.zeros(len(R)); err = np.zeros((len(R), K)); U = np.zeros((len(R), K))
for a, r in enumerate(R):
    e_r, _ = xtb(r.atoms, r.reactant_xyz)
    e_ts, _ = xtb(r.atoms, r.ts_ref_xyz); ref_b[a] = e_ts - e_r
    M = disagreement_matrix(r.preds)
    for k in range(K):
        e, g = xtb(r.atoms, r.preds[k]); Eb[a, k] = e - e_r; G[a, k] = g
        err[a, k] = dmae(r.preds[k], r.ts_ref_xyz); U[a, k] = M[k].sum() / (K - 1)
print(f'{len(R) * (K + 2)} single points in {time.time() - t0:.0f}s; failed: {int(np.isnan(Eb).sum())}')
barrier_rdb7 = np.array([r.barrier_kcal_mol for r in R])
bad = err > BAD
ok = ~np.isnan(Eb).any(1)


def auroc(score):
    return float(np.mean([roc_auc_score(bad[ok, k], score[ok, k]) for k in range(K)]))


rank = lambda s: np.column_stack([rankdata(s[:, k]) for k in range(K)])
out = {
    'sanity_xtb_vs_dft_barrier_of_reference_TS_spearman': float(spearmanr(ref_b, barrier_rdb7, nan_policy='omit').correlation),
    'flag_geometry_bad_auroc': {
        'disagreement': auroc(U), 'xtb_barrier': auroc(Eb), 'xtb_rms_gradient': auroc(G),
        'disagreement+xtb_barrier (rank average)': auroc(rank(U) + rank(Eb)),
        'disagreement+xtb_gradient (rank average)': auroc(rank(U) + rank(G)),
    },
}
# choosing among the 8 guesses
pick = {'consensus': [], 'lowest_xtb_energy': [], 'lowest_xtb_gradient': [], 'consensus_rank+energy_rank': []}
for a, r in enumerate(R):
    if not ok[a]:
        continue
    c = consensus_index(r.preds)
    D = np.array([[np.abs(np.linalg.norm(p[:, None] - p[None], axis=-1) - np.linalg.norm(q[:, None] - q[None], axis=-1)).mean()
                   for q in r.preds] for p in r.preds])
    cons_score = D.sum(1)
    pick['consensus'].append(err[a, c]); pick['lowest_xtb_energy'].append(err[a, np.argmin(Eb[a])])
    pick['lowest_xtb_gradient'].append(err[a, np.argmin(G[a])])
    pick['consensus_rank+energy_rank'].append(err[a, np.argmin(rankdata(cons_score) + rankdata(Eb[a]))])
out['choose_one_of_8'] = {k: {'mean_dmae_A': float(np.mean(v)), 'geometry_bad_rate': float(np.mean(np.array(v) > BAD))}
                          for k, v in pick.items()}
out['note'] = 'GFN2-xTB via tblite, neutral singlets assumed; energies relative to the RDB7 reactant geometry.'
json.dump(out, open('results/xtb.json', 'w'), indent=2)
print(json.dumps(out, indent=2))
