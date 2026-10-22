"""New calculations: optimise TS guesses to saddle points with GFN2-xTB (tblite) and Sella, then check frequencies.

A. Mirror check: for reactions whose best guess was mirror-repaired, optimise both the unrepaired (wrong-handed) and
   the repaired guess. Does the wrong-handed one also reach a first-order saddle point, with the same energy and
   imaginary frequency (so a frequency check passes it), while keeping the wrong handedness?
B. Verdict check: optimise the recommended guess for reactions marked trust / check with DFT / re-generate. How often
   does it reach a first-order saddle point that is the reference TS (D-MAE <= 0.2 A to the DFT TS, right handedness)?

GFN2-xTB is a semi-empirical level, not DFT; the reference TSs are DFT. Usage:
    python scripts/xtb_prospective.py mirror 20      -> results/xtb_prospective_mirror.json
    python scripts/xtb_prospective.py verdicts 20    -> results/xtb_prospective_verdicts.json
"""
import csv, json, os, sys, tempfile, time, warnings
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.vibrations import Vibrations
from sella import Sella
from tblite.ase import TBLite

warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from tsso.data import load_reactions, SYMBOL
from tsso.geometry import dmae
from tsso.metrics import handedness_correct
from tsso.triage import load_checks, triage, MIRROR
from tsso.service import amber_cutoff

KCAL = 23.0605478   # kcal/mol per eV


def optimise(symbols, xyz, steps=400):
    """Sella saddle search (order 1) at GFN2-xTB, then a numerical Hessian. Returns a dict."""
    atoms = Atoms(symbols, positions=np.array(xyz, float))
    atoms.calc = TBLite(method='GFN2-xTB', verbosity=0)
    with tempfile.TemporaryDirectory() as tmp:
        cwd = os.getcwd(); os.chdir(tmp)
        try:
            opt = Sella(atoms, order=1, internal=True, logfile=None, trajectory=None)
            converged = bool(opt.run(fmax=0.03, steps=steps))
            energy = float(atoms.get_potential_energy()) * KCAL
            vib = Vibrations(atoms, name='vib', delta=0.005)
            vib.run()
            freqs = vib.get_frequencies()
            vib.clean()
        except Exception as e:
            os.chdir(cwd)
            return {'converged': False, 'error': f'{type(e).__name__}: {e}'}
        os.chdir(cwd)
    imag = [float(abs(f.imag)) for f in freqs if abs(f.imag) > 30.0]          # cm-1, ignore numerical noise
    return {'converged': converged, 'energy_kcal': energy, 'n_imaginary': len(imag),
            'imaginary_cm1': max(imag) if imag else 0.0, 'xyz': atoms.get_positions()}


def main():
    part, n = sys.argv[1], int(sys.argv[2])
    R = {r.id: r for r in load_reactions()}
    tri = {int(x['reaction_id']): x for x in csv.DictReader(open(ROOT / 'results' / 'triage.csv'))}
    C, cut = load_checks(), amber_cutoff()
    rows, t0 = [], time.time()
    if part == 'mirror':
        ids = sorted(i for i, x in tri.items() if x['mirrored'] == 'True' and x['verdict'] == 'green')[:n]
        for i in ids:
            r = R[i]; t = triage(r, C, cut); sym = [SYMBOL[z] for z in r.atoms]
            raw, rep = optimise(sym, t.xyz * MIRROR), optimise(sym, t.xyz)
            row = {'reaction_id': i}
            for name, res in (('unrepaired', raw), ('repaired', rep)):
                row[name] = {k: v for k, v in res.items() if k != 'xyz'}
                if 'xyz' in res:
                    row[name]['handedness_correct'] = bool(handedness_correct(r, res['xyz'], r.ts_ref_xyz))
                    row[name]['dmae_to_dft_ts_A'] = round(float(dmae(res['xyz'], r.ts_ref_xyz)), 4)
            rows.append(row)
            print(f"{i}: unrepaired {row['unrepaired'].get('n_imaginary')} imag, repaired {row['repaired'].get('n_imaginary')} imag, "
                  f"dE {abs(raw.get('energy_kcal', np.nan) - rep.get('energy_kcal', np.nan)):.4f} kcal/mol ({time.time() - t0:.0f}s)", flush=True)
    else:
        by = {}
        for i, x in sorted(tri.items()):
            if x['mirrored'] == 'True' and x['verdict'] == 'green':
                continue                                   # part A covers mirror repairs
            by.setdefault(x['verdict'], []).append(i)
        ids = [(v, i) for v in ('green', 'amber', 'red') for i in by.get(v, [])[:n]]
        for v, i in ids:
            r = R[i]; t = triage(r, C, cut); sym = [SYMBOL[z] for z in r.atoms]
            res = optimise(sym, t.xyz)
            row = {'reaction_id': i, 'verdict': v, **{k: val for k, val in res.items() if k != 'xyz'}}
            if 'xyz' in res:
                row['handedness_correct'] = bool(handedness_correct(r, res['xyz'], r.ts_ref_xyz))
                row['dmae_to_dft_ts_A'] = round(float(dmae(res['xyz'], r.ts_ref_xyz)), 4)
                row['found_the_reference_ts'] = bool(row['converged'] and row['n_imaginary'] == 1 and
                                                     row['dmae_to_dft_ts_A'] <= 0.2 and row['handedness_correct'])
            rows.append(row)
            print(f"{v} {i}: {row.get('n_imaginary')} imag, D-MAE {row.get('dmae_to_dft_ts_A')}, found {row.get('found_the_reference_ts')} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    (ROOT / 'results' / f'xtb_prospective_{part}.json').write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
