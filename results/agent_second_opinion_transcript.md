Recorded on 22 October 2026: the `/second-opinion 3509 5877 1093` workflow, run by Claude with the Valency MCP Bond
server (https://mcp.valency.io/) connected. The triage steps called the functions behind the `ts-second-opinion` MCP
tools directly (`tsso.service.triage_reaction`, `write_orca_input`); the literature steps were Valency MCP tool calls.
Report: `results/second_opinion_report.md`.

Steps, in order:
1. `triage_reaction(3509)`, `triage_reaction(5877)`, `triage_reaction(1093)`: verdicts green, red, red.
2. `write_orca_input(5877)`, `write_orca_input(1093)`: `results/orca_inputs/rxn5877_ts_guess.inp`, `rxn1093_ts_guess.inp`.
3. Valency `semantic_search_papers`, query "stereochemistry and chirality errors in machine-learned transition state
   structures; mirror images, enantiomers, RMSD that ignores chirality" (from 2023, all sources): 8 results; cited
   RitS (ChemRxiv 2026) and Sain (ChemRxiv 2026).
4. Valency `semantic_search_papers`, query "errors and inconsistencies in reaction datasets: atom mapping disagrees with
   optimized geometries, quality control of transition state databases" (source chemrxiv): 8 results; cited Grambow,
   Pattanaik & Green (ChemRxiv 2020), the dataset RDB7 grew from.
5. `triage_reaction(..., score_against_reference=True)` for each reaction, reported separately as evaluation only.

Summary given to the user:

- **3509 (green):** all 8 guesses broke `chirality`; model 8 was mirror-repaired and passes every rule; the models agree
  (0.024 Å). The unrepaired guess has the same energy and frequencies, so only the stereo rule could catch it. Trust it.
- **5877 (red):** all 8 guesses have C2 inverted relative to the product; mirroring would break the two correct centres
  (a diastereomer), so it cannot be repaired. Re-generate.
- **1093 (red):** all 8 guesses stretch N5–O6 to 2.78 Å because the SMILES says it breaks, but the RDB7 product geometry
  keeps it at 1.45 Å: the dataset entry contradicts itself. Re-generate and check the data.
- **Related reading via Valency:** stereo-aware generation already exists (RitS), which is complementary: TS Second
  Opinion checks any generator's output afterwards. A 2026 review questions what TS accuracy numbers certify.

Evaluation only: 3509 matches the reference (D-MAE 0.019 Å, correct handedness); 5877 does not (D-MAE 0.181 Å but the
wrong diastereomer, which the usual score would call a success); 1093 does not (0.389 Å).
