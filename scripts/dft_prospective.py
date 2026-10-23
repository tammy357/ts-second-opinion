"""New DFT calculations: optimise TS guesses to saddle points with DFT (PySCF, wB97X/def2-SVP) and Sella, then compute
analytic frequencies. Same functional family as RDB7's references (wB97X-D3/def2-TZVP), smaller basis, no dispersion.

Jobs (one JSON file each, written as soon as the job finishes, in results/dft/):
  mirror:<id>   optimise both the unrepaired (wrong-handed) and the repaired recommended guess
  guess:<id>    optimise the recommended guess

Usage: python scripts/dft_prospective.py mirror:3509 guess:1093 ...
"""
import json, os, sys, tempfile, time, warnings
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Bohr, Hartree
from pyscf import dft, gto
from pyscf.hessian import thermo
from sella import Sella

warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from tsso.chem import charge_and_multiplicity
from tsso.data import load_reactions, SYMBOL
from tsso.geometry import dmae
from tsso.metrics import handedness_correct
from tsso.service import amber_cutoff
from tsso.triage import MIRROR, load_checks, triage

XC, BASIS = 'wb97x', 'def2-svp'
KCAL = 627.509474


class PySCFCalc(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, symbols, charge, spin, **kw):
        super().__init__(**kw)
        self.symbols, self.charge, self.spin, self.dm, self.mf = symbols, charge, spin, None, None

    def build(self, positions):
        mol = gto.M(atom=[(s, tuple(p)) for s, p in zip(self.symbols, positions)], basis=BASIS, charge=self.charge,
                    spin=self.spin, verbose=0)
        mf = (dft.RKS(mol) if self.spin == 0 else dft.UKS(mol)).density_fit()
        mf.xc, mf.conv_tol, mf.max_cycle = XC, 1e-9, 200
        return mf

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        mf = self.build(self.atoms.get_positions())
        e = mf.kernel(dm0=self.dm)
        self.dm, self.mf = mf.make_rdm1(), mf
        g = mf.nuc_grad_method().kernel()
        self.results = {'energy': e * Hartree, 'forces': -g * Hartree / Bohr}


def optimise(symbols, xyz, charge, spin, steps=80, fmax=0.05):
    atoms = Atoms(symbols, positions=np.array(xyz, float))
    calc = PySCFCalc(symbols, charge, spin)
    atoms.calc = calc
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        cwd = os.getcwd(); os.chdir(tmp)
        try:
            opt = Sella(atoms, order=1, internal=True, logfile=None, trajectory=None)
            converged = bool(opt.run(fmax=fmax, steps=steps))
            n_steps = opt.nsteps
        finally:
            os.chdir(cwd)
    mf = calc.build(atoms.get_positions())
    energy = mf.kernel(dm0=calc.dm)
    hess = mf.Hessian().kernel()
    freqs = thermo.harmonic_analysis(mf.mol, hess)['freq_wavenumber']
    imag = [float(abs(f.imag)) if np.iscomplexobj(f) else float(-f) for f in freqs
            if (np.iscomplexobj(f) and abs(f.imag) > 30) or (not np.iscomplexobj(f) and f < -30)]
    return {'converged': converged, 'steps': int(n_steps), 'energy_kcal': float(energy * KCAL), 'n_imaginary': len(imag),
            'imaginary_cm1': max(imag) if imag else 0.0, 'minutes': round((time.time() - t0) / 60, 1),
            'xyz': atoms.get_positions().round(6).tolist()}


def main():
    R = {r.id: r for r in load_reactions()}
    C, cut = load_checks(), amber_cutoff()
    out_dir = ROOT / 'results' / 'dft'; out_dir.mkdir(parents=True, exist_ok=True)
    for job in sys.argv[1:]:
        kind, rid = job.split(':'); rid = int(rid)
        dest = out_dir / f'{kind}_{rid}.json'
        if dest.exists():
            continue
        r = R[rid]; t = triage(r, C, cut); sym = [SYMBOL[z] for z in r.atoms]
        charge, mult = charge_and_multiplicity(r.smiles)
        starts = {'unrepaired': t.xyz * MIRROR, 'repaired': t.xyz} if kind == 'mirror' else {'recommended': t.xyz}
        row = {'reaction_id': rid, 'verdict': t.verdict, 'mirror_repaired': bool(t.mirrored), 'level': f'{XC}/{BASIS}',
               'charge': charge, 'multiplicity': mult}
        for name, x in starts.items():
            try:
                res = optimise(sym, x, charge, mult - 1)
            except Exception as e:                      # record the failure and move on to the next job
                row[name] = {'error': f'{type(e).__name__}: {e}'[:300]}
                print(f'{job} {name}: failed {e}', flush=True)
                continue
            res['handedness_correct'] = bool(handedness_correct(r, np.array(res['xyz']), r.ts_ref_xyz))
            res['dmae_to_reference_ts_A'] = round(float(dmae(np.array(res['xyz']), r.ts_ref_xyz)), 4)
            res['found_the_reference_ts'] = bool(res['converged'] and res['n_imaginary'] == 1 and
                                                 res['dmae_to_reference_ts_A'] <= 0.2 and res['handedness_correct'])
            row[name] = res
            print(f"{job} {name}: converged {res['converged']} in {res['steps']} steps, {res['n_imaginary']} imag "
                  f"({res['imaginary_cm1']:.0f} cm-1), D-MAE {res['dmae_to_reference_ts_A']}, hand {res['handedness_correct']}, "
                  f"{res['minutes']} min", flush=True)
        dest.write_text(json.dumps(row, indent=1))


if __name__ == '__main__':
    main()
