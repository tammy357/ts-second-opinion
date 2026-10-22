"""Fetch FK-Flow's published predictions (Mark, Galustian, Kovar & Heid; Zenodo 10.5281/zenodo.23416442, CC BY 4.0)
and convert them to a small file TS Second Opinion can read.

Run from the repository folder (needs internet, about 225 MB of downloads, no PyTorch needed):
    python scripts/fkflow_import.py
It writes fkflow_export.zip (a few MB) with:
    inventory.txt      what is inside the Zenodo zips (so the conversion can be fixed if the format differs)
    runs.csv           every sampling run found, with its settings
    preds_*.json.gz    final TS predictions for the RDB7 test reactions of the unsteered and FK-steered runs
    results/           small Sella / IRC / DFT result tables, if found
Then compare with:  python scripts/fkflow_compare.py fkflow_export.zip
"""
import csv, gzip, io, json, re, sys, urllib.request, zipfile
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from load_samples import OuterUnpickler  # torch-free reader for GoFlow-style sample pickles

BASE = 'https://zenodo.org/records/23416442/files/'
FILES = ['README.md', 'predictions.zip', 'results.zip', 'dft_analysis.zip']
DL = Path.home() / 'Downloads' / 'fkflow_zenodo'
OUT = Path('fkflow_export')


def fetch(name):
    DL.mkdir(parents=True, exist_ok=True)
    p = DL / name
    if p.exists() and p.stat().st_size > 0:
        print(f'  have {name} ({p.stat().st_size / 1e6:.1f} MB)'); return p
    print(f'  downloading {name} ...', flush=True)
    urllib.request.urlretrieve(BASE + name + '?download=1', p)
    print(f'  got {name} ({p.stat().st_size / 1e6:.1f} MB)')
    return p


def overrides_for(zf, member, names):
    """Hydra settings of the run that produced this samples file (overrides.yaml or config.yaml nearby)."""
    parts = member.split('/')
    for up in range(1, 5):
        base = '/'.join(parts[:-up])
        for cand in (f'{base}/.hydra/overrides.yaml', f'{base}/.hydra/config.yaml'):
            if cand in names:
                return cand, zf.read(cand).decode(errors='replace')
    return None, ''


def label(text):
    t = text.lower()
    method = re.search(r'inference_sampling_method[=:]\s*([a-z_]+)', t)
    method = method.group(1) if method else ('fk_steering' if 'fk_potentials' in t else 'median')
    pot = re.search(r'fk_potentials[=:]\s*([a-z_]+)', t)
    ns = re.search(r'num_samples[=:]\s*(\d+)', t)
    temp = re.search(r'fk_steering_temperature[=:]\s*([\d.]+)', t)
    split = re.search(r'(random_split|rxn_core_split|barrier_split)', t)
    return {'method': method, 'potentials': pot.group(1) if pot else '', 'num_samples': ns.group(1) if ns else '',
            'temperature': temp.group(1) if temp else '', 'split': split.group(1) if split else ''}


def final_positions(rec):
    mp = rec._state['_store']._state['_mapping']
    pg = np.asarray(mp['pos_gen'])
    x = pg[-1] if pg.ndim == 3 else pg
    return int(np.asarray(mp['rxn_index']).ravel()[0]), np.round(x, 5).tolist(), mp.get('smiles')


def main():
    print('1/4 download'); paths = {f: fetch(f) for f in FILES}
    OUT.mkdir(exist_ok=True); (OUT / 'results').mkdir(exist_ok=True)
    (OUT / 'README_zenodo.md').write_text(paths['README.md'].read_text(errors='replace'))

    print('2/4 inventory')
    inv = []
    for f in ('predictions.zip', 'results.zip', 'dft_analysis.zip'):
        zf = zipfile.ZipFile(paths[f]); names = zf.namelist()
        ext = {}
        for n in names:
            e = Path(n).suffix or '(none)'; ext[e] = ext.get(e, 0) + 1
        inv.append(f'== {f}: {len(names)} entries; by extension: {ext}')
        inv += ['  ' + n for n in names[:400]]
    (OUT / 'inventory.txt').write_text('\n'.join(inv))
    print('  wrote inventory.txt')

    print('3/4 predictions')
    zf = zipfile.ZipFile(paths['predictions.zip']); names = zf.namelist()
    samples = [n for n in names if n.endswith('.pkl') and 'sample' in n.lower()]
    print(f'  {len(samples)} sample files')
    rows, chosen = [], {}
    for n in samples:
        cfg_path, cfg = overrides_for(zf, n, set(names))
        lab = label(cfg + ' ' + n)
        rows.append({'file': n, 'config': cfg_path or '', **lab})
        key = (lab['method'], lab['potentials'], lab['split'])
        if lab['split'] in ('', 'random_split') and len(chosen.get(key, [])) < 3:
            chosen.setdefault(key, []).append(n)
    with open(OUT / 'runs.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['file', 'config', 'method', 'potentials', 'num_samples', 'temperature', 'split'])
        w.writeheader(); w.writerows(rows)
    print(f'  wrote runs.csv ({len(rows)} runs); converting {sum(len(v) for v in chosen.values())} of them')
    k = 0
    for key, files in chosen.items():
        for n in files:
            try:
                recs = OuterUnpickler(io.BytesIO(zf.read(n))).load()
                preds, smiles = {}, {}
                for rec in recs:
                    rid, x, smi = final_positions(rec)
                    preds[rid] = x; smiles[rid] = smi
                k += 1
                tag = re.sub(r'[^a-z0-9]+', '_', f'{key[0]}_{key[1]}_{k}'.lower()).strip('_')
                with gzip.open(OUT / f'preds_{tag}.json.gz', 'wt') as fh:
                    json.dump({'source_file': n, 'method': key[0], 'potentials': key[1], 'split': key[2],
                               'preds': preds, 'smiles': smiles}, fh)
                print(f'  converted {n} -> preds_{tag}.json.gz ({len(preds)} reactions)')
            except Exception as e:  # keep going; the inventory shows what to fix
                print(f'  could not convert {n}: {type(e).__name__}: {e}')

    print('4/4 result tables')
    for f in ('results.zip', 'dft_analysis.zip'):
        zf = zipfile.ZipFile(paths[f])
        for n in zf.namelist():
            if Path(n).suffix.lower() in ('.csv', '.json', '.txt', '.tsv') and zf.getinfo(n).file_size < 5_000_000 and \
                    any(w in n.lower() for w in ('sella', 'irc', 'dft', 'result', 'summary', 'metric')):
                dest = OUT / 'results' / n.replace('/', '__')
                dest.write_bytes(zf.read(n))
    with zipfile.ZipFile('fkflow_export.zip', 'w', zipfile.ZIP_DEFLATED) as z:
        for p in OUT.rglob('*'):
            if p.is_file():
                z.write(p, p.relative_to(OUT.parent))
    print(f'\nDone: fkflow_export.zip ({Path("fkflow_export.zip").stat().st_size / 1e6:.1f} MB). Send this file back.')


if __name__ == '__main__':
    main()
