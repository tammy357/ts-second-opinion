"""Package the stereo-aware TS benchmark as plain files (CSV + xyz + dataset card), ready for Hugging Face Datasets or Zenodo.

Usage: python scripts/build_benchmark.py  -> benchmark/
"""
import csv, json, shutil, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.geometry import dmae
from tsso.metrics import handedness_correct, rmsd, stereocentres

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'benchmark'
MIRROR = np.array([1.0, 1.0, -1.0])


def frame(r, x, comment):
    return '\n'.join([str(r.n), comment] + [f'{s} {a:.5f} {b:.5f} {c:.5f}' for s, (a, b, c) in zip(r.symbols, x)])


def main():
    OUT.mkdir(exist_ok=True)
    R = load_reactions()
    tri = {int(x['reaction_id']): x for x in csv.DictReader(open(ROOT / 'results' / 'triage.csv'))}
    with open(OUT / 'goflow_guesses.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['reaction_id', 'model_seed', 'n_stereocentres', 'dmae_A', 'rmsd_mirror_allowed_A', 'rmsd_A',
                    'success_dmae_0.2A', 'handedness_correct', 'whole_mirror_image', 'success_dmae_and_handedness'])
        for r in R:
            nst = len(stereocentres(r))
            for k, p in enumerate(r.preds):
                d = dmae(p, r.ts_ref_xyz); ok = handedness_correct(r, p, r.ts_ref_xyz)
                whole = (not ok) and handedness_correct(r, p * MIRROR, r.ts_ref_xyz)
                w.writerow([r.id, k + 1, nst, f'{d:.4f}', f'{rmsd(p, r.ts_ref_xyz, True):.4f}', f'{rmsd(p, r.ts_ref_xyz):.4f}',
                            int(d <= 0.2), int(ok), int(whole), int(d <= 0.2 and ok)])
    with open(OUT / 'reactions.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['reaction_id', 'mapped_reaction_smiles', 'barrier_kcal_mol', 'n_atoms', 'n_stereocentres',
                    'ts_second_opinion_verdict', 'mirror_repaired', 'final_dmae_A', 'final_handedness_correct'])
        for r in R:
            t = tri[r.id]
            w.writerow([r.id, r.smiles, r.barrier_kcal_mol, r.n, len(stereocentres(r)), t['verdict'], int(t['mirrored'] == 'True'),
                        t['dmae_vs_reference_A'], int(t['stereo_wrong'] != 'True')])
    chiral = [r for r in R if stereocentres(r)]
    (OUT / 'mirror_test_references.xyz').write_text(
        '\n'.join(frame(r, r.ts_ref_xyz * MIRROR, f'RDB7 {r.id} reference TS, mirrored (wrong handedness)') for r in chiral) + '\n')
    (OUT / 'mirror_test_reactions.smi').write_text('\n'.join(r.smiles for r in chiral) + '\n')
    (OUT / 'mirror_test_true_references.xyz').write_text(
        '\n'.join(frame(r, r.ts_ref_xyz, f'RDB7 {r.id} reference TS') for r in chiral) + '\n')
    shutil.copy(ROOT / 'results' / 'leaderboard.md', OUT / 'leaderboard.md')
    print(f'wrote benchmark/: {len(R)} reactions, {len(R) * 8} guesses, {len(chiral)} mirror-test reactions')
    (OUT / 'README.md').write_text(CARD.replace('__NCHIRAL__', str(len(chiral))))


CARD = """---
license: mit
pretty_name: Stereo-aware TS benchmark (RDB7 + GoFlow)
tags: [chemistry, transition-state, stereochemistry, benchmark]
size_categories: [1K<n<10K]
---

# Stereo-aware transition-state benchmark (RDB7 + GoFlow)

Labels that make TS-generation benchmarks see handedness. The usual scores (D-MAE, and RMSD computed with the mirror image
allowed) give a wrong-handed TS a perfect score; these files say, for every guess, whether it is close **and** correctly
handed. From [TS Second Opinion](https://github.com/tammy357/ts-second-opinion) (OSI Hackathon 2026).

| File | Contents |
|---|---|
| `reactions.csv` | the 1,194 RDB7 test reactions (random split): atom-mapped SMILES, barrier, stereocentre count, and TS Second Opinion's verdict and final result |
| `goflow_guesses.csv` | 9,552 GoFlow guesses (8 seed models): D-MAE, RMSD with and without the mirror allowed, handedness correct, whole mirror image, success with and without handedness |
| `mirror_test_references.xyz` + `mirror_test_reactions.smi` | a unit test for TS metrics: the __NCHIRAL__ reference TSs with stereocentres, mirrored. A metric that gives these a perfect score cannot see handedness. `mirror_test_true_references.xyz` holds the unmirrored references in the same order |
| `leaderboard.md` | stereo-aware leaderboard of every released GoFlow model |

Score your own model with `tsso evaluate predictions.xyz references.xyz reactions.smi` (`pip install` the repository).

**Definitions.** Handedness is judged by the sign of the signed volume at every stereocentre that does not take part in
the reaction (compared with the reactant) and every product stereocentre, against the reference DFT TS. Success means
D-MAE at most 0.2 Å. Atom order everywhere follows the atom-map numbers.

**Sources and credit.** Reference geometries and barriers: RDB7 (Spiekermann, Pattanaik & Green, *Sci. Data* 2022).
Predictions: GoFlow (Galustian, Mark, Karwounopoulos, Kovar & Heid, *Digital Discovery* 2025; MIT), as released by the
authors. FK-Flow (Mark, Galustian, Kovar & Heid, 2025) first measured GoFlow's chirality errors with its own definitions.
Please cite RDB7, GoFlow and TS Second Opinion. The labels here are MIT-licensed; the underlying data keep their own terms.
"""

if __name__ == '__main__':
    main()
