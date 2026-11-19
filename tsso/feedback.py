"""A 'Was this verdict right?' link that opens a pre-filled GitHub issue form.

Nothing is sent automatically: the link only opens the form in a browser. It carries the verdict, the rules that the
recommended guess broke, whether it was mirror-repaired, how many guesses there were and the version, but never any
structures; people paste those only if they want to share them.
"""
from urllib.parse import urlencode

from . import __version__

ISSUES = 'https://github.com/tammy357/ts-second-opinion/issues/new'


def feedback_url(verdict, rules_broken=(), mirror_repaired=False, n_guesses=None, smiles_source=None):
    lines = [f'Verdict: {verdict}',
             'Rules broken by the recommended guess: ' + (', '.join(rules_broken) if rules_broken else 'none'),
             f"Mirror-image repair: {'yes' if mirror_repaired else 'no'}"]
    if n_guesses is not None:
        lines.append(f'TS guesses given: {n_guesses}')
    if smiles_source:
        lines.append(f'SMILES: {smiles_source}')
    return ISSUES + '?' + urlencode({'template': 'verdict_feedback.yml', 'title': f'Verdict feedback ({verdict})',
                                     'verdict': verdict, 'details': '\n'.join(lines), 'version': __version__})
