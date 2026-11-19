"""Triage many reactions at once and write one verdict table.

Input, either:
  a folder with one sub-folder per reaction, each holding reactant.xyz, product.xyz and the TS guesses (any other
  .xyz files, or .xyz files in a guesses/ sub-folder; multi-frame files are fine), plus an optional reaction.smi with
  the atom-mapped SMILES. Without reaction.smi, bonds and stereocentres are read from the geometries (needs RDKit).
or
  a CSV file with columns id, reactant, product, guesses ('; '-separated paths for several files) and, optionally,
  smiles, charge, multiplicity and order ('map' or 'smiles'). Paths are relative to the CSV file.

Output (in out_dir, default ./tsso_batch):
  verdicts.csv              one row per reaction: verdict, chosen guess, mirror repair, model disagreement and its
                            rank within this batch, rules broken, input warnings, any error
  <id>/recommended_ts.xyz   the recommended TS guess for every reaction that could be checked
  <id>/orca_optts.inp       an ORCA OptTS + Freq input for reactions to check with DFT or re-generate
  <id>/report.html          a report card per reaction, with reports=True

    tsso batch my_reactions/ --out results_folder
"""
import csv
from pathlib import Path

import numpy as np

from .triage import ORCA_TEMPLATE

ENDPOINT_NAMES = ('reactant', 'product')
SKIP_PREFIXES = ('reference', 'recommended')     # never treat a reference TS or our own output as a guess
COLUMNS = ['id', 'verdict', 'meaning', 'chosen_guess', 'mirror_repaired', 'disagreement_A', 'batch_disagreement_rank',
           'rules_broken', 'smiles_source', 'n_guesses', 'input_warnings', 'error']


def _xyz_files(folder):
    return sorted(p for p in folder.iterdir() if p.is_file() and p.suffix.lower() == '.xyz')


def _from_folder(folder):
    jobs = []
    for sub in sorted(p for p in folder.iterdir() if p.is_dir()):
        files = _xyz_files(sub)
        named = {p.stem.lower(): p for p in files}
        guesses = [p for p in files if p.stem.lower() not in ENDPOINT_NAMES and not p.stem.lower().startswith(SKIP_PREFIXES)]
        if (sub / 'guesses').is_dir():
            guesses += _xyz_files(sub / 'guesses')
        smi = sub / 'reaction.smi'
        jobs.append({'id': sub.name, 'reactant': named.get('reactant'), 'product': named.get('product'),
                     'guesses': guesses, 'smiles': smi.read_text().split()[0] if smi.is_file() and smi.read_text().split() else ''})
    return jobs


def _from_csv(path):
    base, jobs = path.parent, []
    with open(path, newline='', encoding='utf-8-sig') as fh:
        for k, row in enumerate(csv.DictReader(fh)):
            row = {(key or '').strip().lower(): (v or '').strip() for key, v in row.items()}
            g = [base / x.strip() for x in row.get('guesses', '').replace(',', ';').split(';') if x.strip()]
            job = {'id': row.get('id') or f'reaction_{k + 1}', 'smiles': row.get('smiles', ''), 'guesses': g,
                   'reactant': base / row['reactant'] if row.get('reactant') else None,
                   'product': base / row['product'] if row.get('product') else None,
                   'order': row.get('order') or 'map'}
            for key in ('charge', 'multiplicity'):
                if row.get(key):
                    job[key] = int(row[key])
            jobs.append(job)
    return jobs


def triage_batch(source, out_dir='tsso_batch', charge=0, multiplicity=1, reports=False):
    """Triage every reaction in a folder or CSV (see the module docstring). Returns a JSON-serialisable summary."""
    from .service import second_opinion
    source = Path(source)
    if source.is_dir():
        jobs = _from_folder(source)
    elif source.is_file() and source.suffix.lower() == '.csv':
        jobs = _from_csv(source)
    else:
        return {'error': f'{source} is neither a folder of reaction folders nor a .csv file'}
    if not jobs:
        return {'error': f'no reactions found in {source}: expected one sub-folder per reaction, or CSV rows'}
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows, results = [], []
    for job in jobs:
        rid = str(job['id'])
        row = {'id': rid}
        missing = [k for k in ('reactant', 'product') if not job.get(k) or not Path(job[k]).is_file()]
        if missing or not job['guesses']:
            row['error'] = ('missing ' + ', '.join(missing + ([] if job['guesses'] else ['TS guesses'])))
            rows.append(row); results.append(None); continue
        folder = out / rid
        if reports:
            folder.mkdir(parents=True, exist_ok=True)
        given = bool(job.get('smiles'))           # with a SMILES, charge and multiplicity come from it unless set per row
        res = second_opinion(job.get('smiles') or None, Path(job['reactant']), Path(job['product']),
                             [Path(g) for g in job['guesses']], job.get('order', 'map'),
                             report=str(folder / 'report.html') if reports else None, title=rid,
                             charge=job.get('charge', None if given else charge),
                             multiplicity=job.get('multiplicity', None if given else multiplicity))
        if 'error' in res:
            row['error'] = res['error']
            rows.append(row); results.append(None); continue
        folder.mkdir(parents=True, exist_ok=True)
        (folder / 'recommended_ts.xyz').write_text(res['recommended_ts_xyz'] + '\n')
        if res['verdict'] in ('amber', 'red'):
            (folder / 'orca_optts.inp').write_text(ORCA_TEMPLATE.format(
                rid=rid, source=f"verdict {res['verdict']}", xyzname='recommended_ts.xyz',
                charge=res['charge'], mult=res['multiplicity']))
        row.update({'verdict': res['verdict'], 'meaning': res['meaning'], 'chosen_guess': res['chosen_guess'],
                    'mirror_repaired': res['mirror_repaired'], 'disagreement_A': res['disagreement_A'],
                    'rules_broken': '; '.join(f"{x['rule']}: {x['detail']}" for x in res['rules_broken_by_recommended_guess']),
                    'smiles_source': res['smiles_source'], 'n_guesses': res['n_guesses'],
                    'input_warnings': ' | '.join(res['input_warnings'])})
        rows.append(row); results.append(res)
    d = np.array([r['disagreement_A'] if r and r.get('disagreement_A') is not None else np.nan for r in results], float)
    ok = np.isfinite(d)
    for row, value in zip(rows, d):
        if np.isfinite(value) and ok.sum() > 1:
            row['batch_disagreement_rank'] = round(100 * float((d[ok] < value).mean()))   # % of this batch below it
    with open(out / 'verdicts.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, '') for k in COLUMNS})
    counts = {v: sum(r.get('verdict') == v for r in rows) for v in ('green', 'amber', 'red')}
    return {'reactions': len(rows), 'checked': sum(counts.values()), 'errors': sum(1 for r in rows if r.get('error')),
            'verdicts': counts, 'out_dir': str(out), 'table': str(out / 'verdicts.csv'),
            'error_rows': [{'id': r['id'], 'error': r['error']} for r in rows if r.get('error')][:20],
            'note': ('Amber uses the review cut-off calibrated on GoFlow guesses for RDB7. For another generator, '
                     'batch_disagreement_rank (the % of this batch with lower model disagreement) is the fairer guide: '
                     'review the top-ranked reactions first.')}
