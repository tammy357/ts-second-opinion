# Second opinion: reactions 3509, 1093, 239

Triage used rules and model disagreement only. The reference TS was not used for any decision.

## 3509 (green)
All 8 GoFlow guesses failed the `chirality` rule. The best guess (model 8) had a stereocentre the wrong way round, and a mirror-image repair was applied. A mirror image leaves every interatomic distance unchanged, so the geometry is the same and only the handedness changes. The repaired guess breaks no rule. Its mean disagreement with the other models is 0.024 A, well under the 0.109 A review cut-off. The largest spread (0.09 A, H8-H14) is not at a forming or breaking bond. Trust it.

## 1093 (red)
All 8 guesses failed `reacting_bonds`. In the recommended guess (model 6) the N5-O6 distance is 2.78 A, but the reactant has 1.38 A and the product 1.45 A (tolerance 0.35 A). The N-O bond that breaks in this isoxazole ring-opening is therefore too stretched to be a plausible TS. No mirror repair was applied. Model disagreement is low (0.074 A, under the cut-off), so the models agree with each other on a geometry that fails a rule. Agreement does not make it correct. Re-generate the guess or send it to DFT. A TS-optimisation input was written, but treat it as a poor starting point.
- Files: `results/orca_inputs/rxn1093_ts_guess.inp`, `results/orca_inputs/rxn1093_ts_guess.xyz`
- Method: wB97X D4 def2-TZVP, OptTS FREQ.

## 239 (green)
All 8 guesses failed `product_stereo`. The best guess (model 6) was mirror-repaired, which leaves all distances unchanged, and it then passes every rule. Disagreement is 0.066 A, under the 0.109 A cut-off. The largest spread (0.33 A, H7-H9) is between hydrogens that are not at a reacting bond. Trust it.

## Related reading
The Valency / `semantic_search_papers` tools are not available in this session, so no literature search was run.

## Summary

| reaction id | verdict | action |
|---|---|---|
| 3509 | green | Trust model 8 (mirror-repaired) |
| 1093 | red | Re-generate, or DFT from the written ORCA input |
| 239 | green | Trust model 6 (mirror-repaired) |

## Evaluation only
This section uses the reference DFT TS. It did not influence any decision above. Wrong is defined as D-MAE > 0.2 A or wrong stereo.

| reaction id | recommended guess | D-MAE (A) | stereo correct | matched reference? |
|---|---|---|---|---|
| 3509 | 8 | 0.0189 | yes | yes |
| 1093 | 6 | 0.3886 | yes | no (D-MAE above 0.2 A) |
| 239 | 6 | 0.0761 | yes | yes |

The triage flagged the one wrong guess (1093) as red and trusted the two correct ones. This is a retrospective check on three reactions. It is not a new DFT result and not a statistic.
