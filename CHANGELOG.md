# Changelog

## 0.3.0 (2026-11-19)

- **Check a reaction from three `.xyz` files, no SMILES needed.** Without an atom-mapped SMILES, bonds, bond orders and
  stereocentres are read from the reactant and product geometries (`tsso/infer.py`, needs RDKit) and the inferred SMILES
  is returned for checking. On the 1,194 RDB7 test reactions this gives the same verdict as the dataset's own SMILES for
  99.7% and the same recommended guess for 98.7% (`scripts/validate_geometry_input.py`). Works in `tsso check`, Python
  (`second_opinion(None, ...)`), the MCP tool `triage_structures` and the installed app; the in-browser app still needs
  the SMILES. New `--charge` and `--multiplicity` options.
- **Batch mode.** `tsso batch <folder or CSV>` (and the MCP tool `triage_batch`) triages many reactions into one
  `verdicts.csv`, with the recommended geometry for each and ORCA inputs for the ones to check with DFT or re-generate.
  `batch_disagreement_rank` ranks model disagreement within your batch, for generators other than GoFlow.
- **One-line install, and agent set-up.** `pip install git+https://github.com/tammy357/ts-second-opinion` now works
  outside a clone: checking your own reactions needs no extra data, and the benchmark commands fetch the files they
  need on first use. `tsso mcp-config` prints settings for Claude Desktop and Cursor, `tsso skill --install` adds an
  agent skill to Claude Code, and `server.json` plus a release workflow prepare PyPI and MCP Registry publishing.
- **"Was this verdict right?"** Every result (command line, report card, app, MCP) carries a link to a pre-filled
  GitHub issue form. It shares the verdict and the rules broken, never your structures.
- **Clearer claims.** Certified trust is described everywhere as an average guarantee that assumes new reactions resemble
  the calibration set; the harder-split results note their overlap with the gate's training data; RESULTS.md now
  matches the recorded `/second-opinion` agent run.

## 0.2.0 (2026-10-22)

- OSI Hackathon 2026 release: gated LLM-written rules, mirror repair, certified trust, stereo-aware leaderboard, mirror
  test, RDB7 data audit, MCP server, marimo app (also in the browser), and new DFT and GFN2-xTB calculations.
