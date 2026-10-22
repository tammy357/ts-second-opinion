"""Certified trust: decide which TS guesses to trust with a guaranteed average error rate.

This is conformal selection (Jin & Candes, "Selection by prediction with conformal p-values", JMLR 2023) with a
yes/no label. You need a calibration set: reactions processed the same way whose answers you know (for example
from earlier DFT). Each reaction gets a trust score (higher = more trustworthy). For a new reaction j,

    p_j = (1 + number of calibration reactions that were wrong AND score at least as high as j) / (n_cal + 1)

and the reactions picked by Benjamini-Hochberg at level alpha are trusted. If the calibration and new reactions
are exchangeable (drawn the same way), the expected share of trusted reactions that are wrong is at most alpha,
for any model, any score and any sample size. Three conditions matter in practice: the score must not be tuned
on these reactions; the guarantee is an average over calibration draws, not a promise for every draw; and it
does not survive a change in the data or the models without recalibration (see scripts/certify.py).
"""
import numpy as np


def trust_score(verdict, disagreement):
    """Higher = more trustworthy. Reactions still breaking a rule rank last; then lower model disagreement wins."""
    d = 0.0 if disagreement is None or (isinstance(disagreement, float) and np.isnan(disagreement)) else float(disagreement)
    return -(d + (10.0 if verdict == 'red' else 0.0))


def conformal_pvalues(score_cal, wrong_cal, score_new):
    score_cal, wrong_cal = np.asarray(score_cal, float), np.asarray(wrong_cal, bool)
    bad = np.sort(score_cal[wrong_cal])
    n_at_least = len(bad) - np.searchsorted(bad, np.asarray(score_new, float), side='left')  # ties count: conservative
    return (1 + n_at_least) / (len(score_cal) + 1)


def benjamini_hochberg(p, alpha):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p, kind='stable')
    ok = np.nonzero(p[order] <= alpha * np.arange(1, m + 1) / m)[0]
    chosen = np.zeros(m, bool)
    if len(ok):
        chosen[order[:ok.max() + 1]] = True
    return chosen


def certify(score_cal, wrong_cal, score_new, target_error_rate):
    """Boolean mask over the new reactions: True = trust (certified at the target average error rate)."""
    return benjamini_hochberg(conformal_pvalues(score_cal, wrong_cal, score_new), target_error_rate)
