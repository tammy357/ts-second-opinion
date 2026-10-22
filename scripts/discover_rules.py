"""Let Claude propose new rules, then keep only those that help on held-out data.

1. Claude proposes N new plain-English rules for a believable TS (it sees the existing rules, not the test data).
2. Each rule goes through the normal authoring loop and gate (9,540 training references).
3. Selection half (even reaction ids): keep a rule if adding it raises the AUROC for flagging wrong guesses.
4. Report half (odd ids): AUROC of disagreement + current rules, with and without the kept rules.

Usage: python scripts/discover_rules.py [N=8]  -> results/discovery.json, tsso/checks_discovered/*.py
"""
import json, re, subprocess, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.author import author
from tsso.data import load_reactions
from tsso.geometry import dmae, disagreement_matrix
from tsso.gate import load_check_sources
from tsso.sandbox import compile_check, run_check

N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
ROOT = Path(__file__).resolve().parent.parent
OUTDIR = ROOT / 'tsso' / 'checks_discovered'
existing = [re.search(r'^RULE = "(.*)"$', s, re.M).group(1) for s in load_check_sources().values()]

prompt = f"""You are a computational chemist. A generative model predicts transition-state (TS) geometries for small
gas-phase organic reactions (H, C, N, O; at most 7 heavy atoms) from the reactant and product. We check each predicted TS
with simple geometric rules that may only use: element types, the optimised reactant and product geometries, their bond
lists with bond orders, which bonds form or break, and tetrahedral stereo tags. The true TS is not available.

Rules we already have:
{chr(10).join('- ' + e for e in existing)}

Typical model errors are subtle: atoms 0.1-0.3 Angstrom away from the true TS, often around the forming and breaking bonds.
Propose {N} NEW rules, different from the ones above, that a real TS should satisfy and that could catch such errors.
Each rule must be checkable from the information above. Reply with exactly {N} lines, each starting with "- ", plain English, one sentence each."""
res = subprocess.run(['claude', '-p'], input=prompt, capture_output=True, text=True, timeout=300)
rules = [ln[2:].strip() for ln in res.stdout.splitlines() if ln.startswith('- ')][:N]
print('Claude proposed:'); [print(' -', r) for r in rules]

# baseline signals on the test set
R = load_reactions(); K = R[0].preds.shape[0]
base_checks = {n: compile_check(s) for n, s in load_check_sources().items()}


def stereo_wrong(r, x):
    nbr, nbp = r.neighbours('reactant'), r.neighbours('product')
    cents = [(i, nbr) for i in range(r.n) if r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(nbr[i]) >= 3]
    cents += [(i, nbp) for i in range(r.n) if r.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    for i, nb in cents:
        v1 = r.ts_ref_xyz[nb[i][:3]] - r.ts_ref_xyz[i]; v2 = x[nb[i][:3]] - x[i]
        if np.sign(np.dot(v1[0], np.cross(v1[1], v1[2]))) != np.sign(np.dot(v2[0], np.cross(v2[1], v2[2]))):
            return True
    return False


bad = np.array([[dmae(p, r.ts_ref_xyz) > 0.2 or stereo_wrong(r, p) for p in r.preds] for r in R])
U = np.array([disagreement_matrix(r.preds).sum(1) / (K - 1) for r in R])
F0 = np.array([[any(run_check(fn, r, p)[0] is False for fn in base_checks.values()) for p in r.preds] for r in R])
sel = np.array([r.id % 2 == 0 for r in R]); rep = ~sel


def auroc(mask, flags):
    s = U + 10.0 * flags
    return float(np.mean([roc_auc_score(bad[mask, k], s[mask, k]) for k in range(K)]))


results, kept = [], []
for i, rule in enumerate(rules):
    print(f'\n[{i + 1}] {rule}')
    out = author(rule, name=f'discovered_{i + 1}', backend='claude', verbose=True, save_dir=OUTDIR, max_attempts=3)
    entry = {'rule': rule, 'accepted_by_gate': out['accepted'], 'attempts': out['attempts'],
             'reference_pass_rate': out['stats']['pass_rate']}
    if out['accepted']:
        fn = compile_check(out['source'])
        Fi = np.array([[run_check(fn, r, p)[0] is False for p in r.preds] for r in R])
        entry['flag_rate_on_guesses'] = float(Fi.mean())
        entry['precision_on_selection_half'] = float(bad[sel][Fi[sel]].mean()) if Fi[sel].any() else None
        entry['auroc_gain_selection_half'] = auroc(sel, F0 | Fi) - auroc(sel, F0)
        if entry['auroc_gain_selection_half'] > 0.002:
            kept.append((rule, Fi))
    results.append(entry)

F1 = F0.copy()
for _, Fi in kept:
    F1 |= Fi
summary = {'proposed': len(rules), 'passed_gate': sum(e['accepted_by_gate'] for e in results), 'kept_on_selection_half': len(kept),
           'report_half_auroc': {'disagreement': float(np.mean([roc_auc_score(bad[rep, k], U[rep, k]) for k in range(K)])),
                                 'disagreement+original_rules': auroc(rep, F0), 'disagreement+original+discovered_rules': auroc(rep, F1)},
           'rules': results}
json.dump(summary, open(ROOT / 'results' / 'discovery.json', 'w'), indent=2)
print(json.dumps({k: v for k, v in summary.items() if k != 'rules'}, indent=2))
