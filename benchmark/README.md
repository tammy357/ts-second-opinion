---
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
| `mirror_test_references.xyz` + `mirror_test_reactions.smi` | a unit test for TS metrics: the 447 reference TSs with stereocentres, mirrored. A metric that gives these a perfect score cannot see handedness. `mirror_test_true_references.xyz` holds the unmirrored references in the same order |
| `leaderboard.md` | stereo-aware leaderboard of every released GoFlow model |

Score your own model with `tsso evaluate predictions.xyz references.xyz reactions.smi` (`pip install` the repository).

**Definitions.** Handedness is judged by the sign of the signed volume at every stereocentre that does not take part in
the reaction (compared with the reactant) and every product stereocentre, against the reference DFT TS. Success means
D-MAE at most 0.2 Å. Atom order everywhere follows the atom-map numbers.

**Sources and credit.** Reference geometries and barriers: RDB7 (Spiekermann, Pattanaik & Green, *Sci. Data* 2022).
Predictions: GoFlow (Galustian, Mark, Karwounopoulos, Kovar & Heid, *Digital Discovery* 2025; MIT), as released by the
authors. FK-Flow (Mark, Galustian, Kovar & Heid, 2025) first measured GoFlow's chirality errors with its own definitions.
Please cite RDB7, GoFlow and TS Second Opinion. The labels here are MIT-licensed; the underlying data keep their own terms.
