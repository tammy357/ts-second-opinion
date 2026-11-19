"""Plain functions behind the CLI and the MCP server. Every function returns JSON-serialisable data."""
import csv, datetime, json, re
from functools import lru_cache
from pathlib import Path
import numpy as np

from .data import load_reactions, Reaction, SYMBOL, ATOMIC_NUMBER
from .chem import parse_reaction, smiles_atom_order, charge_and_multiplicity
from .geometry import dmae
from .sandbox import compile_check
from .triage import load_checks, triage, write_xyz, ORCA_TEMPLATE, amber_cutoff_for_budget
from .gate import GATE_PASS_RATE
from .audit import audit_reaction, audit_dataset, describe
from .feedback import feedback_url
from .paths import data_file, results_file, output_dir, shown

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
CHECKS = ROOT / 'tsso' / 'checks'
USER_CHECKS = ROOT / 'tsso' / 'checks_user'
VERDICT_MEANING = {'green': 'trust: passes every rule and the models agree',
                   'amber': 'check with DFT: passes the rules but the models disagree a lot',
                   'red': 're-generate: still breaks a rule after repair'}


@lru_cache(maxsize=1)
def reactions():
    return {r.id: r for r in load_reactions()}


DEFAULT_REVIEW_CUTOFF_A = 0.10873   # top-20% model disagreement on GoFlow/RDB7 (results/triage_summary.json)


@lru_cache(maxsize=1)
def amber_cutoff():
    f = RESULTS / 'triage_summary.json'
    if f.exists():
        return float(json.load(open(f))['amber_cutoff_A'])
    return DEFAULT_REVIEW_CUTOFF_A


def _xyz_text(symbols, xyz, comment=''):
    return '\n'.join([str(len(symbols)), comment] + [f'{s} {a:.5f} {b:.5f} {c:.5f}' for s, (a, b, c) in zip(symbols, xyz)])


def _stereo_wrong(r, x):
    nbr, nbp = r.neighbours('reactant'), r.neighbours('product')
    cents = [(i, nbr) for i in range(r.n) if r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre and len(nbr[i]) >= 3]
    cents += [(i, nbp) for i in range(r.n) if r.product_chiral_tags[i] in (1, 2) and len(nbp[i]) >= 3]
    for i, nb in cents:
        v1 = r.ts_ref_xyz[nb[i][:3]] - r.ts_ref_xyz[i]; v2 = x[nb[i][:3]] - x[i]
        if np.sign(np.dot(v1[0], np.cross(v1[1], v1[2]))) != np.sign(np.dot(v2[0], np.cross(v2[1], v2[2]))):
            return True
    return False


def _triage_dict(r, t, score=False):
    out = {'reaction_id': r.id, 'smiles': r.smiles, 'verdict': t.verdict, 'meaning': VERDICT_MEANING[t.verdict],
           'chosen_guess': int(t.chosen_model) + 1, 'mirror_repaired': bool(t.mirrored),
           'disagreement_A': None if np.isnan(t.disagreement) else round(t.disagreement, 4),
           'review_cutoff_A': round(amber_cutoff(), 4),
           'rules_broken_by_recommended_guess': [{'rule': n, 'detail': d} for n, d in t.failed_rules],
           'guesses_rejected_by_rules': {str(k + 1): v for k, v in t.rejected_models.items()},
           'reasons': t.reasons,
           'recommended_ts_xyz': _xyz_text(r.symbols, t.xyz, f'recommended TS guess for reaction {r.id}'),
           'feedback_url': feedback_url(t.verdict, [n for n, _ in t.failed_rules], t.mirrored, len(r.preds))}
    if score and r.ts_ref_xyz is not None and not np.isnan(r.ts_ref_xyz).any():
        out['score_against_reference'] = {'dmae_A': round(dmae(t.xyz, r.ts_ref_xyz), 4),
                                          'stereo_correct': not _stereo_wrong(r, t.xyz),
                                          'note': 'uses the reference DFT TS; for evaluation only'}
    return out


# ------------------------------------------------------------------ reactions in the benchmark set
def list_reactions(verdict=None, mirror_repaired=None, limit=20, offset=0):
    f = results_file('triage.csv')
    rows = list(csv.DictReader(open(f))) if f else []
    if verdict:
        rows = [x for x in rows if x['verdict'] == verdict]
    if mirror_repaired is not None:
        rows = [x for x in rows if (x['mirrored'] == 'True') == bool(mirror_repaired)]
    return {'total': len(rows), 'reactions': [{'reaction_id': int(x['reaction_id']), 'verdict': x['verdict'],
                                               'mirror_repaired': x['mirrored'] == 'True', 'smiles': x['smiles']}
                                              for x in rows[offset:offset + limit]]}


def triage_reaction(reaction_id, score_against_reference=False, report=None):
    r = reactions().get(int(reaction_id))
    if r is None:
        return {'error': f'reaction {reaction_id} is not in the 1,194-reaction RDB7 test set'}
    checks = load_checks()
    t = triage(r, checks, amber_cutoff())
    out = _triage_dict(r, t, score_against_reference)
    if report:
        from .report import write_report
        out['report'] = write_report(report, r, t, checks, 0, 1, f'RDB7 reaction {r.id}',
                                     out.get('score_against_reference'))
    return out


# ------------------------------------------------------------------ user-supplied structures
def _symbol(tok):
    if tok.isdigit():
        return SYMBOL[int(tok)]
    t = tok.strip().capitalize()
    if t not in ATOMIC_NUMBER:
        raise ValueError(f'unknown element "{tok}" in xyz')
    return t


def _parse_body(body, expected=None):
    sym, xyz = [], []
    for ln in body:
        p = ln.split()
        if len(p) < 4:
            continue
        try:
            coords = [float(p[1]), float(p[2]), float(p[3])]
        except ValueError:
            raise ValueError(f'could not read coordinates in the line "{ln.strip()[:60]}"')
        sym.append(_symbol(p[0])); xyz.append(coords)
    xyz = np.array(xyz, float).reshape(-1, 3)
    if not np.isfinite(xyz).all():
        raise ValueError('coordinates must be finite numbers (found NaN or infinity)')
    if expected is not None and len(sym) != expected:
        raise ValueError(f'a frame says it has {expected} atoms but {len(sym)} coordinate lines follow')
    return sym, xyz


def read_xyz_frames(text):
    """Every frame of an xyz or multi-frame xyz text, as [(symbols, (N, 3) array in Angstrom)]."""
    lines = text.replace('\ufeff', '').strip().splitlines()      # tolerate a byte-order mark (Windows editors)
    frames, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1; continue
        try:
            n = int(lines[i].split()[0])
        except ValueError:                      # a bare coordinate block without the count/comment header
            frames.append(_parse_body(lines[i:])); break
        frames.append(_parse_body(lines[i + 2:i + 2 + n], expected=n)); i += 2 + n
    return frames


def parse_xyz(text):
    return read_xyz_frames(text)[0]


def _frames(x):
    """xyz text, a path to an .xyz file, an (N, 3) or (K, N, 3) array, or a list of these -> [(symbols|None, xyz)]."""
    if isinstance(x, (list, tuple)):
        return [f for item in x for f in _frames(item)]
    if isinstance(x, np.ndarray):
        return [(None, a) for a in (x if x.ndim == 3 else x[None])]
    if isinstance(x, Path) or (isinstance(x, str) and '\n' not in x and x.strip().endswith('.xyz')):
        return read_xyz_frames(Path(x).read_bytes().decode('utf-8-sig', errors='replace'))
    return read_xyz_frames(str(x))


def _reorder(frames, order_maps):
    """Atoms listed in SMILES order -> atoms listed by map number."""
    out = []
    for sym, xyz in frames:
        idx = np.argsort(np.array(order_maps) - 1)     # position in file of map number 1, 2, ...
        out.append((None if sym is None else [sym[k] for k in idx], xyz[idx]))
    return out


def build_reaction(mapped_reaction_smiles, reactant, product, guesses, order='map', score_reference=None):
    """(Reaction, None) for any reaction and TS guesses, or (None, error message). See second_opinion for inputs."""
    try:
        atoms, chiral, rb, pb, pchiral = parse_reaction(mapped_reaction_smiles, with_product_chirality=True)
    except ValueError as e:
        return None, str(e)
    try:
        R, P, G = _frames(reactant), _frames(product), _frames(guesses)
    except (ValueError, OSError) as e:
        return None, f'could not read the xyz input: {e}'
    if not R or not P or not G:
        return None, 'need a reactant, a product and at least one TS guess'
    if order == 'smiles':
        rside, pside = mapped_reaction_smiles.split('>>')
        R, G = _reorder(R, smiles_atom_order(rside)), _reorder(G, smiles_atom_order(rside))
        P = _reorder(P, smiles_atom_order(pside))
    elif order != 'map':
        return None, "order must be 'map' or 'smiles'"
    want = [SYMBOL[z] for z in atoms]
    for name, (sym, xyz) in [('reactant', R[0]), ('product', P[0])] + [(f'guess {k + 1}', g) for k, g in enumerate(G)]:
        if xyz.shape != (len(atoms), 3):
            return None, f'{name} has {len(xyz)} atoms; the SMILES has {len(atoms)}'
        if sym is not None and sym != want:
            hint = ''
            if order == 'map':
                rside = mapped_reaction_smiles.split('>>')[0]
                if sym == [want[m - 1] for m in smiles_atom_order(rside)]:
                    hint = (" It matches the SMILES atom order: choose 'in the order of the SMILES' in the app, "
                            "order='smiles' in Python, or --order smiles on the command line.")
            return None, f'{name}: atoms are {"".join(sym)}, expected {"".join(want)} (listed by atom-map number).{hint}'
    ref = np.full((len(atoms), 3), np.nan)
    if score_reference is not None:
        F = _frames(score_reference)
        if order == 'smiles':
            F = _reorder(F, smiles_atom_order(mapped_reaction_smiles.split('>>')[0]))
        ref = F[0][1]
    r = Reaction(-1, mapped_reaction_smiles, float('nan'), np.array(atoms), R[0][1], ref, P[0][1],
                 {(a, b): o for a, b, o in rb}, {(a, b): o for a, b, o in pb}, np.array(chiral), np.array(pchiral),
                 np.array([g[1] for g in G]))
    return r, None


def infer_smiles(reactant, product, guesses, charge=0, multiplicity=1):
    """{'smiles', 'bond_orders', 'notes'} read from the reactant and product geometries (see tsso.infer), or
    {'error': ...}. The reactant, the product and every guess must list the same atoms in the same order."""
    from .infer import reaction_smiles_from_geometry
    try:
        R, P, G = _frames(reactant), _frames(product), _frames(guesses)
    except (ValueError, OSError) as e:
        return {'error': f'could not read the xyz input: {e}'}
    if not R or not P or not G:
        return {'error': 'need a reactant, a product and at least one TS guess'}
    sym = R[0][0]
    if sym is None or P[0][0] is None:
        return {'error': 'without a SMILES, give the reactant and product as xyz text or .xyz files, so the elements are known'}
    for name, (s, x) in [('product', P[0])] + [(f'guess {k + 1}', g) for k, g in enumerate(G)]:
        if len(x) != len(sym):
            return {'error': f'{name} has {len(x)} atoms; the reactant has {len(sym)}'}
        if s is not None and s != sym:
            return {'error': f'{name} lists its atoms as {"".join(s)} but the reactant as {"".join(sym)}; without a '
                             f'SMILES, every file must list the same atoms in the same order'}
    try:
        return reaction_smiles_from_geometry(sym, R[0][1], P[0][1], charge, multiplicity)
    except ValueError as e:
        return {'error': str(e)}


def second_opinion(mapped_reaction_smiles, reactant, product, guesses, order='map', score_reference=None, report=None,
                   title='Your reaction', charge=None, multiplicity=None):
    """Triage TS guesses for any reaction, from any generator.

    mapped_reaction_smiles: reactants>>products with every atom (hydrogens too) atom-mapped 1..N. Optional: give None
    or '' to read bonds, bond orders and stereocentres from the reactant and product geometries instead (needs RDKit;
    then every file must list the same atoms in the same order, and the inferred SMILES is returned for checking).
    reactant, product: one structure each; guesses: one or more TS guesses (several enable model disagreement).
    Each can be xyz text, a path to an .xyz file (multi-frame allowed), or NumPy arrays (only with a SMILES).
    order: 'map' if atoms are listed by atom-map number (atom i = map number i+1), or 'smiles' if they are
    listed in the order atoms appear in the SMILES (reactant side for the reactant and the guesses).
    charge, multiplicity: read from the SMILES when given; without a SMILES they default to 0 and 1.
    report: optional path for a single-file HTML report card.
    """
    inferred = None
    if not (mapped_reaction_smiles or '').strip():
        inferred = infer_smiles(reactant, product, guesses, charge or 0, multiplicity or 1)
        if 'error' in inferred:
            return {'error': inferred['error']}
        mapped_reaction_smiles, order = inferred['smiles'], 'map'
    r, err = build_reaction(mapped_reaction_smiles, reactant, product, guesses, order, score_reference)
    if err:
        return {'error': err}
    checks = load_checks()
    t = triage(r, checks, amber_cutoff())
    out = _triage_dict(r, t, score=score_reference is not None)
    out['reaction_id'] = 'yours'
    if inferred:
        out['charge'], out['multiplicity'] = int(charge or 0), int(multiplicity or 1)
    else:
        c, m = charge_and_multiplicity(mapped_reaction_smiles)
        out['charge'] = c if charge is None else int(charge)
        out['multiplicity'] = m if multiplicity is None else int(multiplicity)
    out['smiles_source'] = 'read from the reactant and product geometries' if inferred else 'given'
    if inferred:
        out['inference_notes'] = inferred['notes']
        out['bond_orders_assigned'] = inferred['bond_orders']
    out['feedback_url'] = feedback_url(t.verdict, [n for n, _ in t.failed_rules], t.mirrored, len(r.preds),
                                       out['smiles_source'])
    if report:
        from .report import write_report
        out['report'] = write_report(report, r, t, checks, out['charge'], out['multiplicity'], title,
                                     out.get('score_against_reference'), notes=inferred['notes'] if inferred else None,
                                     feedback=out['feedback_url'])
    out['n_guesses'] = len(r.preds)
    out['input_warnings'] = describe(audit_reaction(r))
    out['note'] = ('The review cut-off was calibrated on GoFlow guesses for RDB7. For another generator, compare '
                   'disagreement values across your own batch, or calibrate with tsso.certify.')
    return out


def audit_rdb7(limit=20):
    """Where RDB7's atom-mapped SMILES and its own reactant/product geometries disagree (see tsso.audit)."""
    from .data import load_reactions
    test = audit_dataset(load_reactions())
    train = audit_dataset(load_reactions(data_file('rdb7_train_reference.json.gz')))
    rows = [{'split': s, 'reaction_id': rid, 'issues': iss} for s, part in (('test', test), ('train', train)) for rid, _, iss in part]
    return {'reactions_checked': 10734, 'reactions_flagged': len(rows), 'flagged_test_ids': [rid for rid, _, _ in test],
            'all_on_product_side': all(i['side'] == 'product' for r in rows for i in r['issues']),
            'examples': rows[:limit],
            'meaning': 'A flagged reaction has a structure whose bonding contradicts its atom-mapped SMILES (atoms at '
                       'single-bond distance that the SMILES does not bond, or a SMILES bond with the atoms far apart).'}


def audit_files(reactions, reactants, products, order='map'):
    """Audit your own reaction dataset: one atom-mapped SMILES per line, and multi-frame reactant and product .xyz (one
    frame per reaction, same order). Returns every reaction whose geometry contradicts its SMILES (see tsso.audit)."""
    smiles = [ln.split()[0] for ln in (reactions.splitlines() if isinstance(reactions, str) else reactions) if ln.strip()]
    R, P = _frames(reactants), _frames(products)
    if not (len(R) == len(P) == len(smiles)):
        return {'error': f'{len(smiles)} reactions, {len(R)} reactant frames and {len(P)} product frames: they must match'}
    flagged, skipped = [], []
    for k, (smi, (rs, rx), (ps, px)) in enumerate(zip(smiles, R, P)):
        rxn, err = build_reaction(smi, _xyz_text(rs, rx), _xyz_text(ps, px), [_xyz_text(rs, rx)], order)
        if err:
            skipped.append({'reaction': k + 1, 'reason': err}); continue
        issues = audit_reaction(rxn)
        if issues:
            flagged.append({'reaction': k + 1, 'smiles': smi, 'issues': issues})
    return {'reactions_checked': len(smiles) - len(skipped), 'reactions_flagged': len(flagged), 'flagged': flagged,
            'skipped': skipped}


def triage_structures(mapped_reaction_smiles, reactant_xyz, product_xyz, guesses_xyz, order='map', charge=None,
                      multiplicity=None):
    """MCP/CLI wrapper around second_opinion for xyz text. mapped_reaction_smiles may be empty (read from geometry)."""
    if isinstance(guesses_xyz, str):
        guesses_xyz = [guesses_xyz]
    return second_opinion(mapped_reaction_smiles, reactant_xyz, product_xyz, list(guesses_xyz or []), order,
                          charge=charge, multiplicity=multiplicity)


def evaluate_files(predictions, references, reactions, order='map', threshold=0.2):
    """Stereo-aware evaluation of any generator's test-set predictions.

    predictions, references: multi-frame .xyz files (or text) with one frame per reaction, in the same order.
    reactions: a file with one atom-mapped reaction SMILES per line (same order), or a list of SMILES.
    """
    from .metrics import evaluate
    smiles = list(reactions) if isinstance(reactions, (list, tuple)) else \
        [ln.split()[0] for ln in Path(reactions).read_text().splitlines() if ln.strip()]
    P, Rf = _frames(predictions), _frames(references)
    if not (len(P) == len(Rf) == len(smiles)):
        return {'error': f'{len(P)} predictions, {len(Rf)} references and {len(smiles)} reactions: they must match'}
    if order == 'smiles':
        P = [_reorder([f], smiles_atom_order(s.split('>>')[0]))[0] for f, s in zip(P, smiles)]
        Rf = [_reorder([f], smiles_atom_order(s.split('>>')[0]))[0] for f, s in zip(Rf, smiles)]
    try:
        out = evaluate([f[1] for f in P], [f[1] for f in Rf], smiles, threshold)
    except ValueError as e:
        return {'error': str(e)}
    out.pop('per_reaction')
    return out


# ------------------------------------------------------------------ rules
def _header(src):
    rule = re.search(r'^RULE = (".*")$', src, re.M)
    gate = re.search(r'# Gate.*$', src, re.M)
    return (json.loads(rule.group(1)) if rule else None), (gate.group(0)[2:] if gate else None)


def list_rules():
    out = []
    for folder, used in ((CHECKS, True), (USER_CHECKS, False)):
        for p in sorted(folder.glob('*.py')) if folder.exists() else []:
            rule, gate = _header(p.read_text())
            out.append({'name': p.stem, 'rule': rule, 'gate': gate, 'used_in_triage': used})
    return {'rules': out}


def get_rule_source(name):
    for folder in (CHECKS, USER_CHECKS):
        p = folder / f'{name}.py'
        if p.exists():
            return {'name': name, 'source': p.read_text()}
    return {'error': f'no rule named {name}'}


def gate_check(rule_text, python_source):
    """Test a candidate check on the 9,540 RDB7 training reference TSs (never the test set)."""
    from .author import evaluate_source, train_set, feedback_text, MIN_SENSITIVITY
    try:
        compile_check(python_source)
    except Exception as e:
        return {'accepted': False, 'error': f'rejected before testing: {e}'}
    s = evaluate_source(python_source, train_set())
    fires = None
    try:
        fn = compile_check(python_source)
        from .sandbox import run_check
        R = list(reactions().values())
        fires = float(np.mean([[run_check(fn, r, p)[0] is False for p in r.preds] for r in R]))
    except Exception:
        pass
    return {'rule': rule_text, 'accepted': s['accepted'],
            'accepts_real_TS': round(s['pass_rate'], 4), 'required': GATE_PASS_RATE,
            'crash_rate': round(s['crash_rate'], 4),
            'rejects_mirror_images': round(s['rejects_mirror_images'], 4),
            'rejects_noisy_TS': round(s['rejects_noisy_geometries'], 4), 'min_sensitivity': MIN_SENSITIVITY,
            'examples_of_real_TS_rejected': [{'reaction': rid, 'detail': d} for rid, d in s['failures'][:6]],
            'share_of_goflow_guesses_it_flags': fires,
            'feedback': None if s['accepted'] else feedback_text(s)}


def save_rule(name, rule_text, python_source, use_in_triage=False):
    g = gate_check(rule_text, python_source)
    if not g['accepted']:
        return {'saved': False, 'reason': 'the check did not pass the gate', 'gate': g}
    name = re.sub(r'[^a-z0-9_]+', '_', name.lower()).strip('_') or 'rule'
    folder = CHECKS if use_in_triage else USER_CHECKS
    folder.mkdir(parents=True, exist_ok=True)
    body = re.sub(r'^RULE = .*\n', '', python_source, flags=re.M)
    head = (f'RULE = {json.dumps(rule_text)}\n# Written by an agent via the ts-second-opinion MCP server on '
            f'{datetime.date.today().isoformat()}.\n# Gate (RDB7 training references): passes {100 * g["accepts_real_TS"]:.1f}%; '
            f'rejects {100 * g["rejects_mirror_images"]:.0f}% of mirror images and {100 * g["rejects_noisy_TS"]:.0f}% of noisy geometries.\n\n')
    (folder / f'{name}.py').write_text(head + body)
    if use_in_triage:
        load_checks.cache_clear() if hasattr(load_checks, 'cache_clear') else None
    return {'saved': True, 'path': str((folder / f'{name}.py').relative_to(ROOT)), 'used_in_triage': use_in_triage, 'gate': g}


# ------------------------------------------------------------------ outputs
def write_orca_input(reaction_id):
    r = reactions().get(int(reaction_id))
    if r is None:
        return {'error': f'reaction {reaction_id} is not in the test set'}
    t = triage(r, load_checks(), amber_cutoff())
    stem = f'rxn{r.id}_ts_guess'
    outdir = output_dir() / 'orca_inputs'; outdir.mkdir(parents=True, exist_ok=True)
    src = f'GoFlow seed model {t.chosen_model + 1}' + (' (mirror image)' if t.mirrored else '') + f'; verdict {t.verdict}'
    write_xyz(outdir / f'{stem}.xyz', r.symbols, t.xyz, f'RDB7 {r.id} {r.smiles}')
    inp = ORCA_TEMPLATE.format(rid=f'RDB7 reaction {r.id}', source=src, xyzname=f'{stem}.xyz', charge=0, mult=1)
    (outdir / f'{stem}.inp').write_text(inp)
    return {'verdict': t.verdict, 'files': [shown(outdir / f'{stem}.inp'), shown(outdir / f'{stem}.xyz')], 'orca_input': inp}


# ------------------------------------------------------------------ planning a DFT budget
def _benchmark_rows():
    """Verdict rows for the 1,194 RDB7 test reactions, with the reference-based labels kept apart (evaluation only)."""
    f = results_file('triage.csv')
    if f is None:
        return None, None
    raw = list(csv.DictReader(open(f)))
    rows = [{'id': x['reaction_id'], 'verdict': x['verdict'], 'disagreement_A': float(x['disagreement_A']),
             'mirror_repaired': x['mirrored'] == 'True', 'rules_broken': x['failed_rules']} for x in raw]
    truth = {x['reaction_id']: {'geometry_wrong': x['geometry_bad'] == 'True', 'handedness_wrong': x['stereo_wrong'] == 'True'}
             for x in raw}
    return rows, truth


def score_plan(plan, truth):
    """What the reference TS says about a plan, by action. For evaluation only: a plan never uses it."""
    out = {}
    for action in ('run_dft', 'trust', 'hold', 'regenerate'):
        ids = [it['id'] for it in plan['plan'] if it['action'] == action and it['id'] in truth]
        out[action] = {'reactions': len(ids),
                       'geometry_wrong': sum(truth[i]['geometry_wrong'] for i in ids),
                       'handedness_wrong': sum(truth[i]['handedness_wrong'] for i in ids),
                       'geometry_wrong_ids': [i for i in ids if truth[i]['geometry_wrong']]}
    can = [it['id'] for it in plan['plan'] if it['action'] in ('run_dft', 'trust', 'hold') and it['id'] in truth]
    base = sum(truth[i]['geometry_wrong'] for i in can) / len(can) if can else 0.0
    out['expected_geometry_wrong_if_runs_were_picked_at_random'] = round(out['run_dft']['reactions'] * base, 2)
    out['note'] = ('Uses the RDB7 reference TS (wrong geometry = D-MAE above 0.2 A). For evaluation only; '
                   'the plan was made without it.')
    return out


def plan_dft_budget(budget, path='', reaction_ids=None, sample=None, seed=0, evaluate=False, target_error_rate=None,
                    out_dir=None):
    """Plan where `budget` DFT runs go (see tsso.plan). Either your own reactions or the RDB7 benchmark.

    path: a verdicts.csv written by `tsso batch`, the folder holding it, or a folder / CSV of reactions (triaged first
    into out_dir). Empty: the benchmark, restricted to reaction_ids, or to `sample` reactions drawn with `seed`.
    evaluate (benchmark only): add what the reference TS says about the picks. target_error_rate (benchmark only): add
    certified trust, calibrated on the benchmark reactions that are not in the batch.
    """
    from .plan import plan_dft_budget as make, read_verdict_table, write_plan_csv, explain_plan
    if isinstance(budget, float) and budget.is_integer():
        budget = int(budget)
    if path:
        src = Path(path)
        table = src / 'verdicts.csv' if src.is_dir() and (src / 'verdicts.csv').is_file() else src
        if not (table.is_file() and table.suffix.lower() == '.csv' and
                'verdict' in open(table, encoding='utf-8-sig').readline().lower()):
            from .batch import triage_batch
            res = triage_batch(src, out_dir or 'tsso_batch')
            if 'error' in res:
                return res
            table = Path(res['table'])
        plan = make(read_verdict_table(table), budget)
        if 'error' in plan:
            return plan
        plan['source'] = shown(table)
        plan['plan_table'] = shown(write_plan_csv(plan, table.parent / 'dft_plan.csv'))
        inputs = {it['id']: table.parent / it['id'] / 'orca_optts.inp' for it in plan['plan'] if it['action'] == 'run_dft'}
        plan['orca_inputs'] = {i: shown(p) for i, p in inputs.items() if p.is_file()}
        missing = [i for i, p in inputs.items() if not p.is_file()]
        if missing:
            plan['caveats'].append('No ORCA input was found for: ' + ', '.join(missing) + '. `tsso batch` writes one for '
                                   'every "check with DFT" reaction; for others use `tsso check --orca`.')
        if evaluate or target_error_rate is not None:
            plan['caveats'].append('evaluate and target_error_rate apply to the benchmark only: your reactions have no '
                                   'reference TS here, and certified trust needs a calibration set (tsso.certify).')
    else:
        rows, truth = _benchmark_rows()
        if rows is None:
            return {'error': 'results/triage.csv is missing; run scripts/triage_all.py first'}
        everything = list(rows)
        if reaction_ids:
            want = {str(int(i)) for i in reaction_ids}
            rows = [r for r in rows if r['id'] in want]
            unknown = sorted(want - {r['id'] for r in rows})
            if unknown:
                return {'error': 'not in the 1,194-reaction RDB7 test set: ' + ', '.join(unknown[:10])}
        elif sample:
            if not 2 <= int(sample) <= len(rows):
                return {'error': f'sample must be between 2 and {len(rows)}'}
            pick = np.random.default_rng(int(seed)).choice(len(rows), int(sample), replace=False)
            rows = [rows[k] for k in sorted(pick)]
        cal = None
        if target_error_rate is not None:
            used = {r['id'] for r in rows}
            rest = [r for r in everything if r['id'] not in used]
            cal = {'verdict': [r['verdict'] for r in rest], 'disagreement_A': [r['disagreement_A'] for r in rest],
                   'wrong': [truth[r['id']]['geometry_wrong'] for r in rest]}
        plan = make(rows, budget, cal, target_error_rate if target_error_rate is not None else 0.05)
        if 'error' in plan:
            return plan
        plan['source'] = f'RDB7 benchmark: {len(rows)} test reactions, 8 GoFlow guesses each'
        if evaluate:
            plan['reference_check'] = score_plan(plan, truth)
    plan['in_plain_words'] = explain_plan(plan)
    return plan


def plan_evidence():
    """How the planning rule did against reference TSs on RDB7 (results/plan_evaluation.json)."""
    f = results_file('plan_evaluation.json')
    if f is None:
        return {'error': 'results/plan_evaluation.json is missing; run scripts/plan_evaluation.py'}
    return json.load(open(f))


def certify_trust(target_error_rate=0.05, count='geometry', seed=0, limit=30):
    """Certified trust on the RDB7 test set: calibrate on a random half, certify the other half (conformal selection)."""
    from .certify import certify, trust_score
    if count not in ('geometry', 'geometry_or_handedness'):
        return {'error': "count must be 'geometry' or 'geometry_or_handedness'"}
    f = results_file('triage.csv')
    if f is None:
        return {'error': 'results/triage.csv is missing; run scripts/triage_all.py first'}
    rows = list(csv.DictReader(open(f)))
    ids = np.array([int(x['reaction_id']) for x in rows])
    score = np.array([trust_score(x['verdict'], float(x['disagreement_A'])) for x in rows])
    wrong = np.array([x['geometry_bad'] == 'True' or (count != 'geometry' and x['stereo_wrong'] == 'True') for x in rows])
    perm = np.random.default_rng(seed).permutation(len(rows)); c, t = perm[:len(rows) // 2], perm[len(rows) // 2:]
    sel = certify(score[c], wrong[c], score[t], float(target_error_rate))
    trusted = sorted(int(i) for i in ids[t][sel])
    out = {'target_error_rate': float(target_error_rate), 'counting': count,
           'calibration_reactions': int(len(c)), 'new_reactions': int(len(t)),
           'trusted': len(trusted), 'share_trusted': round(len(trusted) / len(t), 4),
           'trusted_reaction_ids': trusted[:limit],
           'actual_error_rate_among_trusted': round(float(wrong[t][sel].mean()), 4) if sel.any() else None,
           'note': ('Guarantee: on average over calibration draws, the share of trusted reactions that are wrong is at most '
                    'the target, if new reactions are drawn like the calibration ones. Recalibrate when the generator, '
                    'the data or the rules change. actual_error_rate_among_trusted uses the reference TS (evaluation only).')}
    f = results_file('certify.json')
    if f is not None:
        rep = json.load(open(f))[count]['conformal']
        best = min(rep, key=lambda d: abs(d['target'] - float(target_error_rate)))
        out['over_2000_random_splits'] = best
    return out


def benchmark_summary():
    out = {}
    for name in ('metrics', 'triage_summary', 'benchmark', 'certify', 'generality', 'xtb', 'discovery'):
        f = results_file(f'{name}.json')
        if f is not None:
            out[name] = json.load(open(f))
    m, b, c = out.get('metrics', {}), out.get('benchmark', {}), out.get('certify', {})
    fa = m.get('final_answer', {})
    one, ours = fa.get('single_model_avg', {}), fa.get('check_aware_plus_mirror', {})
    wh = b.get('wrong_handed_guesses', {})
    g5 = next((d for d in c.get('geometry', {}).get('conformal', []) if d['target'] == 0.05), None)
    head = {
        'geometry_errors_in_final_answers': {'one_model': one.get('geometry_bad_rate'), 'ts_second_opinion': ours.get('geometry_bad_rate'),
                                             'note': 'the rules do not fix geometry; they and model disagreement say which guesses to trust'},
        'handedness_errors_in_final_answers': {'one_model': one.get('stereo_bad_rate'), 'ts_second_opinion': ours.get('stereo_bad_rate')},
        'any_error_in_final_answers': {'one_model': one.get('any_bad_rate'), 'ts_second_opinion': ours.get('any_bad_rate')},
        'wrong_handed_guesses_caught_by_rules': m.get('stereo_checks_recall_of_stereo_errors'),
        'wrong_handed_guesses_that_are_whole_mirror_images': wh.get('share_whole_mirror_images'),
        'why_that_matters': ('a whole mirror image has the same energy and frequencies as the right TS: it is the TS of the '
                             'mirror-image reaction, so D-MAE, RMSD (ignore_chirality=True) and a DFT frequency check all '
                             'pass it; only a stereo check against the reaction as written catches it'),
        'wrong_handed_guesses_that_pass_dmae_0.2A': wh.get('share_of_wrong_handed_with_dmae_below_0.2A'),
        'certified_trust_geometry_at_5pct_target': g5,
        'risk_coverage_geometry_only': b.get('risk_coverage_geometry_only'),
        'risk_coverage_geometry_or_handedness': b.get('risk_coverage'),
        'definition': 'geometry wrong = D-MAE > 0.2 A; handedness wrong = a stereocentre the wrong way round; both vs the '
                      'RDB7 reference TS; 1,194 test reactions, 8 GoFlow models',
    }
    return {'headline': head, 'files': sorted(out)}
