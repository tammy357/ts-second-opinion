"""Plan a DFT budget: given triage verdicts for many reactions and the number of DFT runs you can afford, say which
reactions get a run, which to trust, which to hold back and which to re-generate.

The plan is a fixed function of its inputs: no randomness, no language model, and no reference answers. The rule is

  1. re-generate   every reaction whose recommended guess still breaks a rule (red). No DFT run goes to these.
  2. run DFT       the `budget` reactions whose guesses disagree most, among those that pass the rules.
  3. hold          reactions marked "check with DFT" (amber) that the budget did not reach: not verified, so not to
                   be used as they are.
  4. trust         the rest (green and not picked).

Why "most disputed first": it puts the runs where a wrong answer is most likely, so the answers you go on to use
without a check contain as few errors as possible. Tested afterwards against the RDB7 reference TSs (8 GoFlow guesses
per reaction): in random batches of 40 reactions with 5 runs, 2.25 of the 5 picks had a wrong geometry on average,
against 0.58 when picking at random and 3.98 with a perfect order (scripts/plan_evaluation.py, written to
results/plan_evaluation.json, with intervals, other budgets and a second ensemble). Using disagreement between
models to decide where an expensive label goes is borrowed from active learning: "query by committee" (Seung, Opper
& Sompolinsky, COLT 1992; for learned potentials, Schran, Brezina & Marsalek, J. Chem. Phys. 153, 104105 (2020)).
There it picks training data; here it picks which answers to check.

What the rule does not do, and the plan says so in `caveats`:
  - It does not maximise the number of runs that end in a usable TS. A run started from a wrong guess usually ends
    somewhere else, so a picked reaction is often found to need re-generating; that is the point of checking it.
  - It cannot see an error that every model shares. Model disagreement measures spread, not bias.
  - A TS optimisation that ends with one imaginary frequency is not yet a confirmed TS for this reaction: follow the
    reaction path (IRC) to both ends.
  - With a single guess per reaction there is no disagreement to rank by.

Optional: with a calibration set (reactions processed the same way whose answers are known), `certified` marks the
reactions that conformal selection trusts at a target average error rate (tsso.certify; Jin & Candes, JMLR 2023).
"""
import csv
import math
from pathlib import Path

ACTIONS = ('run_dft', 'trust', 'hold', 'regenerate', 'not_checked')
ACTION_LABEL = {'run_dft': 'Run DFT', 'trust': 'Trust', 'hold': 'Hold (not verified)', 'regenerate': 'Re-generate',
                'not_checked': 'Fix the input'}
VERIFY = ('After each TS optimisation, check for exactly one imaginary frequency and then follow the reaction path '
          '(IRC) to confirm it joins this reactant to this product. One imaginary frequency alone is not enough.')


def _number(x):
    """A finite float, or None for missing values ('', None, NaN)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _truthy(x):
    return x is True or str(x).strip().lower() in ('true', '1', 'yes')


def rows_from_csv_text(text):
    """Rows of a verdict table given as CSV text (see read_verdict_table)."""
    import io
    raw = list(csv.DictReader(io.StringIO(text.lstrip('\ufeff'))))
    rows = []
    for r in raw:
        r = {(k or '').strip().lower(): (v or '').strip() for k, v in r.items()}
        rows.append({'id': r.get('id') or r.get('reaction_id') or f'row_{len(rows) + 1}',
                     'verdict': r.get('verdict', ''), 'disagreement_A': _number(r.get('disagreement_a')),
                     'mirror_repaired': _truthy(r.get('mirror_repaired', r.get('mirrored', ''))),
                     'rules_broken': r.get('rules_broken', r.get('failed_rules', '')), 'error': r.get('error', '')})
    return rows


def read_verdict_table(path):
    """Rows of a verdict table: the verdicts.csv written by `tsso batch`, or results/triage.csv in the repository."""
    return rows_from_csv_text(Path(path).read_text(encoding='utf-8-sig'))


def plan_dft_budget(rows, budget, calibration=None, target_error_rate=0.05):
    """Split reactions into run DFT / trust / hold / re-generate for a budget of `budget` DFT runs.

    rows: one dict per reaction with 'id', 'verdict' ('green', 'amber' or 'red') and 'disagreement_A' (model
    disagreement in Angstrom, or None when only one guess was given); 'mirror_repaired', 'rules_broken' and 'error'
    are optional. This is what triage_batch writes to verdicts.csv (see read_verdict_table).
    budget: how many DFT TS optimisations you can afford (a whole number, 0 or more).
    calibration: optional {'verdict': [...], 'disagreement_A': [...], 'wrong': [...]} for reactions with known answers
    that were processed the same way; adds a 'certified' flag at `target_error_rate` (conformal selection).

    Returns a JSON-serialisable dict. The same inputs always give the same plan; ties in disagreement are broken by id.
    """
    if isinstance(budget, bool) or not isinstance(budget, int) or budget < 0:
        return {'error': 'budget must be a whole number of DFT runs, 0 or more'}
    ids = [str(r.get('id')) for r in rows]
    if len(set(ids)) != len(ids):
        return {'error': 'reaction ids must be unique'}
    items, ranked = [], []
    for r in rows:
        v = str(r.get('verdict') or '').strip().lower()
        it = {'id': str(r.get('id')), 'verdict': v or None, 'disagreement_A': _number(r.get('disagreement_A')),
              'mirror_repaired': _truthy(r.get('mirror_repaired')), 'action': None, 'dft_order': None, 'why': ''}
        if v not in ('green', 'amber', 'red'):
            it['action'] = 'not_checked'
            it['why'] = (r.get('error') or 'no verdict for this reaction').strip()
        elif v == 'red':
            it['action'] = 'regenerate'
            broke = str(r.get('rules_broken') or '').strip()
            it['why'] = 'still breaks a rule after repair' + (f' ({broke[:120]})' if broke else '')
        elif it['disagreement_A'] is None:
            it['action'] = 'trust'
            it['why'] = 'passes every rule, but only one guess was given, so model disagreement is unknown'
        else:
            ranked.append(it)
        items.append(it)
    ranked.sort(key=lambda it: (-it['disagreement_A'], it['id']))       # most disputed first; ties by id
    for k, it in enumerate(ranked):
        d = f"{it['disagreement_A']:.3f} A"
        if k < budget:
            it['action'], it['dft_order'] = 'run_dft', k + 1
            it['why'] = (f'number {k + 1} by model disagreement ({d})' + ('' if it['verdict'] == 'amber' else
                         '; a spot check: every "check with DFT" reaction is already covered, so this run is optional'))
        elif it['verdict'] == 'amber':
            it['action'] = 'hold'
            it['why'] = f'the models disagree ({d}) and the budget ran out before it; do not use it unverified'
        else:
            it['action'] = 'trust'
            it['why'] = f'passes every rule and the models agree ({d})'
        if it['mirror_repaired']:
            it['why'] += '; wrong-handed guess repaired by mirroring'
    order = {a: k for k, a in enumerate(ACTIONS)}
    items.sort(key=lambda it: (order[it['action']], it['dft_order'] or 0,
                               -(it['disagreement_A'] if it['disagreement_A'] is not None else -1.0), it['id']))
    counts = {a: sum(it['action'] == a for it in items) for a in ACTIONS}
    n_amber = sum(it['verdict'] == 'amber' for it in items)
    single = sum(it['verdict'] == 'green' and it['disagreement_A'] is None for it in items)
    caveats = [VERIFY,
               'Model disagreement cannot see an error that every model makes; it measures spread, not bias.',
               'The order was tested on GoFlow guesses for small gas-phase H/C/N/O molecules (RDB7). For other models or '
               'chemistry, check it on reactions whose TS you already know before relying on it.']
    if counts['hold']:
        caveats.insert(0, f"{counts['hold']} reaction(s) marked \"check with DFT\" are beyond this budget. They are on "
                          f"hold, not trusted: {n_amber} runs would cover all of them.")
    if single:
        caveats.append(f'{single} trusted reaction(s) had a single guess, so nothing but the rules vouches for them. '
                       f'Supply several guesses (models or samples) to rank them.')
    if counts['regenerate']:
        caveats.append('No DFT run is planned for the re-generate reactions: in the DFT and GFN2-xTB tests in this '
                       'repository none of 6 such guesses optimised to the reference TS (results/dft_summary.json).')
    out = {'budget': budget, 'reactions': len(items), 'counts': counts, 'dft_runs_planned': counts['run_dft'],
           'of_which_optional_spot_checks': sum(it['action'] == 'run_dft' and it['verdict'] != 'amber' for it in items),
           'budget_unused': budget - counts['run_dft'],
           'policy': ('Re-generate what still breaks a rule; give the DFT runs to the reactions whose guesses disagree '
                      'most; hold the other "check with DFT" reactions; trust the rest.'),
           'run_dft': [it['id'] for it in items if it['action'] == 'run_dft'],
           'plan': items, 'caveats': caveats}
    if calibration is not None:
        out['certified_trust'] = _certify(items, calibration, target_error_rate)
    return out


def _certify(items, calibration, target):
    """Mark which reactions conformal selection trusts at the target average error rate (see tsso.certify)."""
    import numpy as np
    from .certify import certify, trust_score
    cal_v, cal_d, cal_w = calibration['verdict'], calibration['disagreement_A'], calibration['wrong']
    if not (len(cal_v) == len(cal_d) == len(cal_w)) or len(cal_v) < 20:
        return {'error': 'calibration needs matching verdict, disagreement_A and wrong lists of at least 20 reactions'}
    if not 0 < float(target) < 1:
        return {'error': 'target_error_rate must be between 0 and 1'}
    judged = [it for it in items if it['verdict'] in ('green', 'amber', 'red')]
    s_cal = np.array([trust_score(v, _number(d)) for v, d in zip(cal_v, cal_d)])
    s_new = np.array([trust_score(it['verdict'], it['disagreement_A']) for it in judged])
    chosen = certify(s_cal, np.array([_truthy(w) for w in cal_w]), s_new, float(target)) if judged else []
    for it, ok in zip(judged, chosen):
        it['certified'] = bool(ok)
    ids = [it['id'] for it in judged if it.get('certified')]
    return {'target_error_rate': float(target), 'calibration_reactions': len(cal_v), 'certified': len(ids),
            'certified_ids': ids,
            'certified_and_not_sent_to_dft': sum(it.get('certified', False) and it['action'] == 'trust' for it in judged),
            'note': ('On average over calibration draws, at most the target share of the certified reactions is wrong, '
                     'if these reactions are drawn like the calibration ones (same generator, same kind of chemistry). '
                     'The guarantee is for the certified set as a whole, not for each reaction and not for a subset.')}


def explain_plan(plan):
    """The plan in plain words. Written by this template from the plan itself, not by a language model."""
    if 'error' in plan:
        return f"Could not make a plan: {plan['error']}."
    c, n, b = plan['counts'], plan['reactions'], plan['budget']
    runs = [it for it in plan['plan'] if it['action'] == 'run_dft']
    lines = [f"You have {n} reaction{'s' if n != 1 else ''} and {b} DFT run{'s' if b != 1 else ''}."]
    if runs:
        spread = (f"{runs[0]['disagreement_A']:.3f} A down to {runs[-1]['disagreement_A']:.3f} A" if len(runs) > 1
                  else f"{runs[0]['disagreement_A']:.3f} A")
        lines.append(f"Run DFT on {len(runs)}: {', '.join(it['id'] for it in runs)}, in order of how much the guesses "
                     f"disagree ({spread}).")
        spot = sum(it['verdict'] != 'amber' for it in runs)
        if spot:
            lines.append(('One of these is an optional spot check on a reaction' if spot == 1 else
                          f'{spot} of these are optional spot checks on reactions') + ' that would otherwise be trusted: '
                         'every reaction marked "check with DFT" is already covered.')
    elif b:
        lines.append('No reaction has a disagreement value to rank, so no run is planned.')
    if plan['budget_unused'] > 0 and runs:
        lines.append(f"{plan['budget_unused']} run(s) are left over: nothing else could be ranked.")
    if c['trust']:
        lines.append(f"Trust {c['trust']}: they pass every rule and their guesses agree.")
    if c['hold']:
        lines.append(f"Hold {c['hold']}: marked \"check with DFT\" but beyond this budget, so do not use them unverified.")
    if c['regenerate']:
        lines.append(f"Re-generate {c['regenerate']}: they still break a rule, and a DFT run from such a guess is "
                     f"unlikely to reach the right TS.")
    if c['not_checked']:
        lines.append(f"{c['not_checked']} could not be checked; fix their input files first.")
    cert = plan.get('certified_trust')
    if cert and 'error' not in cert:
        lines.append(f"{cert['certified']} reaction(s) are certified at a {100 * cert['target_error_rate']:g}% target "
                     f"average error rate, which holds only if these reactions resemble the calibration set.")
    lines.append(VERIFY)
    return ' '.join(lines)


def write_plan_csv(plan, path):
    """Write the plan as a table (one row per reaction, in the order to act on)."""
    cols = ['id', 'action', 'dft_order', 'verdict', 'disagreement_A', 'mirror_repaired', 'certified', 'why']
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for it in plan['plan']:
            w.writerow({k: ('' if it.get(k) is None else it.get(k)) for k in cols})
    return str(Path(path))


# ------------------------------------------------------------------ for the app: the plan as a table and as tiles
_TILE = {'run_dft': ('#1d4ed8', 'Run DFT', 'most disputed first'), 'trust': ('#15803d', 'Trust', 'pass the rules, models agree'),
         'hold': ('#b45309', 'Hold', 'need a check, beyond the budget'), 'regenerate': ('#b91c1c', 'Re-generate', 'still break a rule'),
         'not_checked': ('#55534e', 'Fix the input', 'could not be checked')}


def tiles_html(plan):
    """One tile per action with its count, as HTML."""
    import html
    tiles = []
    for a in ACTIONS:
        n = plan['counts'][a]
        if a == 'not_checked' and not n:
            continue
        colour, name, hint = _TILE[a]
        tiles.append(f"<div style='flex:1;min-width:150px;padding:12px 14px;border:1px solid #e3e2de;border-top:4px solid {colour};"
                     f"border-radius:12px;background:#fff'><div style='font-size:13px;color:#55534e'>{html.escape(name)}</div>"
                     f"<div style='font-size:30px;font-weight:750;line-height:1.15;color:{colour}'>{n}</div>"
                     f"<div style='font-size:12px;color:#55534e'>{html.escape(hint)}</div></div>")
    return "<div style='display:flex;gap:12px;flex-wrap:wrap'>" + ''.join(tiles) + '</div>'


def table_rows(plan, truth=None):
    """The plan as rows for a table. truth: optional {id: {'geometry_wrong': bool, 'handedness_wrong': bool}} from a
    reference TS, shown in an extra column (evaluation only; the plan was made without it)."""
    rows = []
    for it in plan['plan']:
        r = {'What to do': ACTION_LABEL[it['action']], 'Reaction': it['id'],
             'DFT order': it['dft_order'] if it['dft_order'] else '',
             'Model disagreement (A)': '' if it['disagreement_A'] is None else round(it['disagreement_A'], 3),
             'Why': it['why']}
        if truth is not None:
            t = truth.get(it['id'])
            r['Reference TS says (evaluation only)'] = ('' if t is None else 'wrong geometry' if t['geometry_wrong'] else
                                                        'wrong handedness' if t['handedness_wrong'] else 'right')
        rows.append(r)
    return rows
