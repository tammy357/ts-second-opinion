"""Plain-English rule -> LLM-written check -> gate on reference TSs -> rewrite until accepted.

    python -m tsso.author "Bonds that form or break have lengths between their reactant and product values." --name reacting_bonds

The gate set is the 9,540 RDB7 *training* reference transition states. A check is accepted only if it
accepts at least 99% of them without crashing. We also report how often it rejects deliberately
corrupted references (mirror images and noisy geometries), so a check that passes everything is visible.
"""
import argparse, datetime, json, re
from pathlib import Path
import numpy as np

from .data import load_reactions
from .gate import GATE_PASS_RATE

MIN_SENSITIVITY = 0.01   # must reject at least 1% of mirror-image or noisy (0.3 A) reference TSs
from .llm import ask, build_prompt, lookup_cached, CHECK_DIR
from .sandbox import compile_check, run_check

ROOT = Path(__file__).resolve().parent.parent
TRAIN = ROOT / 'data' / 'rdb7_train_reference.json.gz'
LOG = ROOT / 'results' / 'authoring_log.jsonl'
_train = None


def train_set():
    global _train
    if _train is None:
        _train = load_reactions(TRAIN)
    return _train


def evaluate_source(source, reactions, rng_seed=0):
    """Gate statistics plus sensitivity to corrupted references."""
    fn = compile_check(source)
    rng = np.random.default_rng(rng_seed)
    ok = crash = 0; failures = []
    rej_mirror = rej_noise = 0
    for r in reactions:
        passed, detail = run_check(fn, r, r.ts_ref_xyz)
        if passed is None:
            crash += 1; failures.append((r.id, detail))
        elif passed:
            ok += 1
        else:
            failures.append((r.id, detail))
        m, _ = run_check(fn, r, r.ts_ref_xyz * np.array([1.0, 1.0, -1.0]))
        z, _ = run_check(fn, r, r.ts_ref_xyz + rng.normal(0, 0.3, r.ts_ref_xyz.shape))
        rej_mirror += (m is False); rej_noise += (z is False)
    n = len(reactions)
    return {'pass_rate': ok / n, 'crash_rate': crash / n, 'failures': failures[:8], 'n_failures': len(failures),
            'rejects_mirror_images': rej_mirror / n, 'rejects_noisy_geometries': rej_noise / n,
            'accepted': ok / n >= GATE_PASS_RATE and crash == 0 and max(rej_mirror, rej_noise) / n >= MIN_SENSITIVITY}


def feedback_text(stats):
    if stats['pass_rate'] >= GATE_PASS_RATE and stats['crash_rate'] == 0:
        return (f'it accepted {100 * stats["pass_rate"]:.1f}% of real transition states, which is good, but it almost never fires: '
                f'it rejected only {100 * stats["rejects_mirror_images"]:.1f}% of mirror-image TSs and '
                f'{100 * stats["rejects_noisy_geometries"]:.1f}% of TSs with 0.3 Angstrom random noise added. A useful check must '
                f'still catch clear violations of the rule. Tighten it (for example a smaller tolerance) without rejecting real TSs. '
                'Do not wrap the function in try/except.')
    ex = '\n'.join(f'  - reaction {rid}: {d}' for rid, d in stats['failures'][:6])
    msg = (f'it rejected {100 * (1 - stats["pass_rate"]):.1f}% of {len(train_set())} real DFT transition states '
           f'(and crashed on {100 * stats["crash_rate"]:.1f}%). It must accept at least {100 * GATE_PASS_RATE:.0f}% of real '
           f'transition states and never crash. Examples of real TSs it wrongly rejected:\n{ex}\n'
           'Revise the check so that real transition states pass, while it still catches clear violations of the rule. '
           'Do not make it accept everything.')
    return msg


def save_check(name, rule, source, stats, provenance, save_dir=CHECK_DIR):
    body = re.sub(r'^RULE = .*\n', '', source, flags=re.M)
    header = (f'RULE = {json.dumps(rule)}\n'
              f'# Written by {provenance["model"]} via {provenance["backend"]} on {provenance["date"]}, '
              f'attempt {provenance["attempt"]}.\n'
              f'# Gate (RDB7 training references): passes {100 * stats["pass_rate"]:.1f}%; '
              f'rejects {100 * stats["rejects_mirror_images"]:.0f}% of mirror images and '
              f'{100 * stats["rejects_noisy_geometries"]:.0f}% of noisy geometries.\n\n')
    save_dir = Path(save_dir); save_dir.mkdir(parents=True, exist_ok=True)
    path = save_dir / f'{name}.py'
    path.write_text(header + body)
    return path


def author(rule, name=None, backend='auto', max_attempts=4, verbose=True, save_dir=CHECK_DIR):
    """Returns dict with accepted flag, source, stats and the attempts made."""
    name = name or re.sub(r'[^a-z0-9]+', '_', rule.lower()).strip('_')[:40]
    if backend == 'cache':
        src, stem = lookup_cached(rule)
        if src is None:
            raise LookupError('rule not in cache; use backend "claude" or "api"')
        stats = evaluate_source(src, train_set())
        return {'name': stem, 'accepted': stats['accepted'], 'source': src, 'stats': stats, 'attempts': 0, 'backend': 'cache'}
    attempts, previous, fb = [], None, None
    for attempt in range(1, max_attempts + 1):
        code, used, model = ask(build_prompt(rule, fb, previous), backend)
        try:
            stats = evaluate_source(code, train_set())
        except Exception as e:
            stats = {'pass_rate': 0.0, 'crash_rate': 1.0, 'failures': [(None, f'could not compile: {e}')], 'n_failures': 0,
                     'rejects_mirror_images': 0, 'rejects_noisy_geometries': 0, 'accepted': False}
        attempts.append({'attempt': attempt, 'pass_rate': stats['pass_rate'], 'crash_rate': stats['crash_rate'],
                         'accepted': stats['accepted'], 'rejects_mirror_images': stats['rejects_mirror_images'],
                         'rejects_noisy_geometries': stats['rejects_noisy_geometries'], 'code': code})
        if verbose:
            print(f'  attempt {attempt}: passes {100 * stats["pass_rate"]:.1f}% of references, '
                  f'rejects {100 * stats["rejects_mirror_images"]:.0f}% mirrors / {100 * stats["rejects_noisy_geometries"]:.0f}% noisy'
                  f' -> {"ACCEPTED" if stats["accepted"] else "rewrite"}')
        if stats['accepted']:
            prov = {'backend': used, 'model': model, 'date': datetime.date.today().isoformat(), 'attempt': attempt}
            path = save_check(name, rule, code, stats, prov, save_dir)
            _log(rule, name, used, model, attempts, True)
            return {'name': name, 'accepted': True, 'source': path.read_text(), 'stats': stats,
                    'attempts': attempt, 'backend': used, 'path': str(path)}
        previous, fb = code, feedback_text(stats)
    _log(rule, name, used, model, attempts, False)
    return {'name': name, 'accepted': False, 'source': code, 'stats': stats, 'attempts': max_attempts, 'backend': used}


async def author_async(rule, ask_fn, name=None, max_attempts=3, save_dir=CHECK_DIR, on_attempt=None):
    """Same loop as author(), for callers whose LLM call is async (the in-browser app). ask_fn(prompt) -> (code, backend, model)."""
    name = name or re.sub(r'[^a-z0-9]+', '_', rule.lower()).strip('_')[:40]
    attempts, previous, fb = [], None, None
    for attempt in range(1, max_attempts + 1):
        code, used, model = await ask_fn(build_prompt(rule, fb, previous))
        try:
            stats = evaluate_source(code, train_set())
        except Exception as e:
            stats = {'pass_rate': 0.0, 'crash_rate': 1.0, 'failures': [(None, f'could not compile: {e}')], 'n_failures': 0,
                     'rejects_mirror_images': 0, 'rejects_noisy_geometries': 0, 'accepted': False}
        attempts.append({'attempt': attempt, 'pass_rate': stats['pass_rate'], 'crash_rate': stats['crash_rate'],
                         'accepted': stats['accepted'], 'rejects_mirror_images': stats['rejects_mirror_images'],
                         'rejects_noisy_geometries': stats['rejects_noisy_geometries'], 'code': code})
        if on_attempt:
            on_attempt(attempts[-1])
        if stats['accepted']:
            prov = {'backend': used, 'model': model, 'date': datetime.date.today().isoformat(), 'attempt': attempt}
            path = save_check(name, rule, code, stats, prov, save_dir)
            return {'name': name, 'accepted': True, 'source': path.read_text(), 'stats': stats, 'attempts': attempt,
                    'attempt_log': attempts, 'backend': used, 'model': model}
        previous, fb = code, feedback_text(stats)
    return {'name': name, 'accepted': False, 'source': code, 'stats': stats, 'attempts': max_attempts,
            'attempt_log': attempts, 'backend': used, 'model': model}


def _log(rule, name, backend, model, attempts, accepted):
    LOG.parent.mkdir(exist_ok=True)
    with open(LOG, 'a') as fh:
        fh.write(json.dumps({'time': datetime.datetime.now().isoformat(timespec='seconds'), 'rule': rule, 'name': name,
                             'backend': backend, 'model': model, 'accepted': accepted, 'attempts': attempts}) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rule'); ap.add_argument('--name'); ap.add_argument('--backend', default='auto')
    ap.add_argument('--max-attempts', type=int, default=4)
    ap.add_argument('--save-dir', default=str(CHECK_DIR.parent / 'checks_user'),
                    help='where accepted checks go (default tsso/checks_user; move a file into tsso/checks to use it in triage)')
    a = ap.parse_args()
    out = author(a.rule, a.name, a.backend, a.max_attempts, save_dir=a.save_dir)
    print(('ACCEPTED -> ' + out.get('path', '')) if out['accepted'] else 'REJECTED after all attempts')
