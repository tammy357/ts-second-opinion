"""Decision-point experiment: do gated rule checks add to model disagreement for triaging TS guesses?

Retrospective, on GoFlow's released predictions for the 1,194 RDB7 test reactions (8 seed models).
Checks were gated and thresholded on the 9,540 *training* references only.

Usage: python scripts/evaluate.py  [writes results/metrics.json]
"""
import json, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.gate import load_check_sources
from tsso.sandbox import compile_check, run_check
from tsso.geometry import dmae, disagreement_matrix, consensus_index

BAD = 0.2      # Angstrom D-MAE vs reference TS; about twice the average error
BUDGETS = (0.1, 0.2, 0.3)
R = load_reactions()
checks = {n: compile_check(s) for n, s in load_check_sources().items()}
STEREO = ('chirality', 'product_stereo')
K = R[0].preds.shape[0]

# Per reaction, per model: error, chirality truth, disagreement, check results
err = np.zeros((len(R), K)); chi_bad = np.zeros((len(R), K), bool); fails = {n: np.zeros((len(R), K), bool) for n in checks}
unc = np.zeros((len(R), K)); details = {}
for a, r in enumerate(R):
    M = disagreement_matrix(r.preds)
    for k in range(K):
        err[a, k] = dmae(r.preds[k], r.ts_ref_xyz)
        unc[a, k] = M[k].sum() / (K - 1)
        for n, fn in checks.items():
            ok, d = run_check(fn, r, r.preds[k])
            fails[n][a, k] = (ok is False)
    # stereo label against the reference TS: reactant spectator centres and product-defined centres
    nbr, nbp = r.neighbours('reactant'), r.neighbours('product')
    cents = [(i, nbr) for i in range(r.n) if r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(nbr[i]) >= 3]
    cents += [(i, nbp) for i in range(r.n) if r.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    for k in range(K):
        for i, nb in cents:
            v1 = r.ts_ref_xyz[nb[i][:3]] - r.ts_ref_xyz[i]; v2 = r.preds[k][nb[i][:3]] - r.preds[k][i]
            if np.sign(np.dot(v1[0], np.cross(v1[1], v1[2]))) != np.sign(np.dot(v2[0], np.cross(v2[1], v2[2]))):
                chi_bad[a, k] = True

geo_bad = err > BAD
any_bad = geo_bad | chi_bad
any_fail = np.zeros_like(geo_bad)
for n in checks:
    any_fail |= fails[n]
shape_fail = np.zeros_like(geo_bad)
for n in checks:
    if n not in STEREO:
        shape_fail |= fails[n]


def catch(label, score, budget):
    out = []
    for k in range(K):
        o = np.argsort(-score[:, k], kind='stable')[:int(budget * len(R))]
        out.append(label[o, k].sum() / max(label[:, k].sum(), 1))
    return float(np.mean(out))


def auroc(label, score):
    return float(np.mean([roc_auc_score(label[:, k], score[:, k]) for k in range(K)]))


combo = unc + 10.0 * any_fail          # failing a check puts a reaction at the front of the queue
combo_shape = unc + 10.0 * shape_fail
m = {'n_reactions': len(R), 'n_models': K, 'bad_threshold_A': BAD,
     'mean_dmae_single_model': float(err.mean()),
     'rate_geometry_bad': float(geo_bad.mean()), 'rate_stereo_bad': float(chi_bad.mean()),
     'rate_any_bad': float(any_bad.mean()),
     'reactions_with_spectator_stereocentre': int(sum(1 for r in R if any(
         r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(r.neighbours()[i]) >= 3 for i in range(r.n)))),
     'flag_rate_per_check': {n: float(f.mean()) for n, f in fails.items()},
     'precision_per_check_any_bad': {n: float(any_bad[f].mean()) if f.any() else None for n, f in fails.items()},
     'stereo_checks_recall_of_stereo_errors': float((fails['chirality'] | fails['product_stereo'])[chi_bad].mean()),
     'stereo_checks_false_alarm_rate': float((fails['chirality'] | fails['product_stereo'])[~chi_bad].mean()),
     'auroc': {
         'geometry_bad': {'disagreement': auroc(geo_bad, unc), 'disagreement+shape_checks': auroc(geo_bad, combo_shape)},
         'any_bad': {'disagreement': auroc(any_bad, unc), 'disagreement+all_checks': auroc(any_bad, combo)},
     },
     'catch_rate': {}}
for b in BUDGETS:
    m['catch_rate'][f'{int(b * 100)}%'] = {
        'geometry_bad': {'disagreement': catch(geo_bad, unc, b), 'disagreement+shape_checks': catch(geo_bad, combo_shape, b)},
        'any_bad': {'disagreement': catch(any_bad, unc, b), 'disagreement+all_checks': catch(any_bad, combo, b)},
    }
m['flagged_by_checks_fraction'] = float(any_fail.mean())

# Final answer: plain consensus vs consensus among predictions that pass every check
final = {'consensus': [], 'check_aware_consensus': [], 'check_aware_plus_mirror': [], 'oracle_best': []}
labels = {'consensus': [], 'check_aware_consensus': [], 'check_aware_plus_mirror': []}
n_mirrored = 0


def stereo_wrong(r, x):
    nbr, nbp = r.neighbours('reactant'), r.neighbours('product')
    cents = [(i, nbr) for i in range(r.n) if r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(nbr[i]) >= 3]
    cents += [(i, nbp) for i in range(r.n) if r.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    for i, nb in cents:
        v1 = r.ts_ref_xyz[nb[i][:3]] - r.ts_ref_xyz[i]; v2 = x[nb[i][:3]] - x[i]
        if np.sign(np.dot(v1[0], np.cross(v1[1], v1[2]))) != np.sign(np.dot(v2[0], np.cross(v2[1], v2[2]))):
            return True
    return False

for a, r in enumerate(R):
    passing = ~np.array([any(fails[n][a, k] for n in checks) for k in range(K)])
    c = consensus_index(r.preds); c2 = consensus_index(r.preds, allowed=passing)
    final['consensus'].append(err[a, c]); final['check_aware_consensus'].append(err[a, c2]); final['oracle_best'].append(err[a].min())
    labels['consensus'].append((geo_bad[a, c], chi_bad[a, c])); labels['check_aware_consensus'].append((geo_bad[a, c2], chi_bad[a, c2]))
    # repair: if the chosen guess still fails a stereo rule, try its mirror image (distances, and so D-MAE, unchanged)
    x = r.preds[c2]
    if not passing[c2] and any(fails[n][a, c2] for n in STEREO):
        mx = x * np.array([1.0, 1.0, -1.0])
        if all(run_check(fn, r, mx)[0] for fn in checks.values()):
            x = mx; n_mirrored += 1
    final['check_aware_plus_mirror'].append(err[a, c2])
    labels['check_aware_plus_mirror'].append((geo_bad[a, c2], stereo_wrong(r, x)))
m['final_answer'] = {}
for key in ('consensus', 'check_aware_consensus', 'check_aware_plus_mirror'):
    L = np.array(labels[key])
    m['final_answer'][key] = {'mean_dmae': float(np.mean(final[key])), 'geometry_bad_rate': float(L[:, 0].mean()),
                              'stereo_bad_rate': float(L[:, 1].mean()), 'any_bad_rate': float((L[:, 0] | L[:, 1]).mean())}
m['final_answer']['single_model_avg'] = {'mean_dmae': float(err.mean()), 'geometry_bad_rate': float(geo_bad.mean()),
                                         'stereo_bad_rate': float(chi_bad.mean()), 'any_bad_rate': float(any_bad.mean())}
m['final_answer']['oracle_best_of_8_mean_dmae'] = float(np.mean(final['oracle_best']))
m['final_answer']['mirror_repairs_applied'] = n_mirrored

Path('results').mkdir(exist_ok=True)
json.dump(m, open('results/metrics.json', 'w'), indent=2)
np.savez('results/per_reaction.npz', ids=np.array([r.id for r in R]), err=err, unc=unc, geo_bad=geo_bad,
         chi_bad=chi_bad, any_fail=any_fail, **{f'fail_{n}': f for n, f in fails.items()})
print(json.dumps(m, indent=2))
