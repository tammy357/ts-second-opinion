# TS Second Opinion (agent notes)

- Purpose: decide which machine-learned transition-state (TS) guesses to trust, repair or send to DFT.
- MCP server: `tsso mcp` (configured in .mcp.json; needs `pip install -e .`). Tools: triage_structures (any reaction, any
  elements, from xyz text), evaluate_predictions (stereo-aware scoring), certify_trust, list_reactions, triage_reaction,
  list_rules, get_rule_source, gate_check, save_rule, write_orca_input, benchmark_summary, audit_rdb7_dataset.
- CLI equivalents: tsso check / evaluate / certify / triage / gate / author / audit. Python: `from tsso import second_opinion`.
- Example inputs live in examples/. `tsso check ... --report r.html` writes a one-file HTML report card.
- The app (apps/demo.py) is verdict-first: examples, upload, linked 3D views, Rule Studio replay, trust slider.
  It also runs in the browser (docs/app, built by scripts/build_wasm.py); there tsso/smiles_lite.py stands in for RDKit.
- Extras: tsso/mirror_test.py (does a metric see handedness?), tsso/stereo_edit.py (experimental local inversion of
  partly flipped guesses; not in the pipeline), tsso/audit.py (SMILES vs geometry consistency).
- FK-Flow comparison: scripts/fkflow_import.py (run where Zenodo is reachable) then scripts/fkflow_compare.py.
- Slash commands: /second-opinion <ids | red | amber | mirror-repaired>, /new-rule "<rule>".
- Never use the reference TS (score_against_reference) to make a decision; it is for evaluation only.
- Data: 1,194 RDB7 test reactions with 8 GoFlow guesses each; rules are gated on 9,540 training references.
- Claims: keep geometry errors (D-MAE > 0.2 A) and handedness errors (wrong stereo) separate. The rules fix handedness
  (17.7% -> 2.3%), not geometry (14.1% -> 11.7%). 90% of wrong-handed guesses are whole mirror images: same energy and
  frequencies, i.e. the TS of the mirror-image reaction. Certified trust (conformal selection) is an average guarantee
  that assumes exchangeable data. Results are mostly retrospective, plus new TS optimisations:
  DFT wB97X/def2-SVP (scripts/dft_prospective.py, results/dft/): mirror pair 3509 0.0006 kcal/mol apart, both 1 imag;
  found reference TS 6/6 trusted vs 0/6 re-generate (5 of those still pass a frequency check).
  GFN2-xTB (scripts/xtb_prospective.py): mirror pairs 34/40 identical saddles; found reference TS 77% / 60% / 0% by verdict.
