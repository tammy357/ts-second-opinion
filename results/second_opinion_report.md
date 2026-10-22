# Second opinion: RDB7 reactions 3509, 5877 and 1093

Run by Claude with the TS Second Opinion triage (`tsso.service`, the functions behind its MCP tools) and the
**Valency MCP Bond** server for literature search. The verdicts use only the TS guesses, the reactant, the product and the gated rules; the reference
DFT TS is used only in the "evaluation only" table at the end.

## 3509: ✓ trust (mirror-repaired)

All 8 GoFlow guesses broke the `chirality` rule: stereocentre C2, which takes no part in the reaction, was the wrong way
round. The best guess (model 8) was mirrored; mirroring keeps every interatomic distance, so it fixes the handedness
without touching the geometry. The repaired guess passes every rule and the models agree closely (0.024 Å, cut-off
0.109 Å). The unrepaired guess has exactly the same energy and frequencies as the repaired one: it is the TS of the
mirror-image reaction, so a DFT frequency check would not have caught it. Only the stereo rule did.

## 5877: ✕ re-generate (cannot be repaired)

All 8 guesses broke `product_stereo`: C2 has the opposite handedness to the product. Mirroring was tried but breaks the
two correct stereocentres, because only one of three is inverted: this is the TS of a diastereomer, whose energy really
differs. The models agree with each other (0.032 Å), but they agree on the wrong diastereomer. ORCA input:
`results/orca_inputs/rxn5877_ts_guess.inp` (treat the geometry as a poor start; re-generating is better).

## 1093: ✕ re-generate, and check the data

All 8 guesses broke `reacting_bonds`: the SMILES says the N5–O6 bond breaks, and the guesses stretch it to 2.78 Å,
against 1.38 Å in the reactant and 1.45 Å in the product. The product geometry itself keeps N5–O6 at a bonding
1.45 Å, so the RDB7 entry contradicts its own SMILES (`tsso audit` flags it). ORCA input:
`results/orca_inputs/rxn1093_ts_guess.inp`.

## Related reading (Valency MCP Bond)

Two `semantic_search_papers` queries. These are related reading, not evidence about these specific reactions.

- *Stereochemistry errors in machine-learned transition states* (for 3509 and 5877):
  [Right into the Saddle: Stereochemistry-Aware Generation of Molecular Transition States](https://doi.org/10.26434/chemrxiv.15001681/v1)
  (Nikitin, Anstine & Isayev, ChemRxiv 2026) makes generation itself stereo-aware by encoding product chirality;
  TS Second Opinion instead checks and repairs any generator's output afterwards.
  [Are we measuring the right thing? Evaluation practice for machine-learned transition states and reaction barriers](https://doi.org/10.26434/chemrxiv.15009215/v1)
  (Sain, ChemRxiv 2026) asks what reported TS accuracies certify, the same question the mirror test asks of metrics.
- *Inconsistencies in reaction datasets* (for 1093):
  [Reactants, Products, and Transition States of Elementary Chemical Reactions Based on Quantum Chemistry](https://doi.org/10.26434/chemrxiv.11400240.v2)
  (Grambow, Pattanaik & Green, ChemRxiv 2020) is the dataset RDB7 grew from, with its atom-mapped reaction SMILES.

| reaction id | verdict | action |
|---|---|---|
| 3509 | green | Trust model 8 (mirror-repaired) |
| 5877 | red | Re-generate; a partly flipped guess cannot be mirror-repaired |
| 1093 | red | Re-generate, and check the RDB7 entry (its product geometry contradicts its SMILES) |

**Evaluation only** (uses the reference DFT TS; played no part in the verdicts):

| reaction id | D-MAE (Å) | handedness correct | matches the reference? |
|---|---|---|---|
| 3509 | 0.019 | yes | yes |
| 5877 | 0.181 | no | no (wrong diastereomer, despite D-MAE ≤ 0.2 Å) |
| 1093 | 0.389 | yes | no (above 0.2 Å) |

The triage trusted the one guess that matched the reference and flagged both that did not. Note 5877: by the usual
D-MAE ≤ 0.2 Å score it would count as a success. This is three reactions, so it says nothing about overall reliability;
see docs/RESULTS.md for that.
