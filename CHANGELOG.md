# Changelog

## Unreleased

- **Plan a DFT budget.** `tsso plan <verdicts> --budget N` (MCP tools `plan_dft_budget` and `plan_evidence`; app section
  3b) says which reactions get a DFT run, which to trust, which to hold and which to re-generate. A fixed rule with no
  language model in it. Tested against the RDB7 reference TSs (`scripts/plan_evaluation.py`): in random batches of 40
  reactions with 5 runs, 2.25 of the 5 picks had a wrong geometry on average, against 0.58 at random and 3.98 with a
  perfect order.
- **Two goals.** `--goal find_errors` (default) or `--goal finish_ts` (the runs most likely to end in a usable TS).
- **Random audits.** `--audit N` spends N runs on trusted reactions drawn at random, and `--results runs.csv` (MCP tool
  `review_dft_results`) turns their outcomes into an exact upper bound on the error rate of the other trusted
  reactions: the only way the plan can see an error that every model shares, and the only error rate it can give for
  your own reactions without a calibration set (`--calibration` adds certified trust if you have one).
- **Reaction-path check.** The plan writes ORCA inputs with the IRC included, and `tsso irc-check` (MCP tool `irc_check`)
  says whether the path joins this reactant to this product with the right handedness. Tested on stand-in endpoints
  only; the ORCA input has not been run by us.
- **One guess per reaction** is now ranked too, by its distance from the reactant and the product: a weaker sign,
  about twice chance on the reaction-core and barrier splits. `verdicts.csv` gains `endpoint_distance_A`, `charge` and
  `multiplicity` columns.
- **A planning agent.** `tsso agent "<question>"` and a chat in the app let you ask for a plan in your own words (your
  own Anthropic API key). Claude can only call the planning rule, is not sent coordinates or SMILES, and points at
  tool results instead of typing numbers; any number it types itself is listed under the answer. Not yet tried with
  a live model.
- **A planner page for the showcase site** (`docs/planner/`, built by `scripts/build_planner.py`): one small HTML file
  with the planning rule in JavaScript, so it opens at once on a phone. Every reaction in the batch is drawn as a mark,
  with sliders for the budget and the audit, a switch that shows which guesses were really wrong, a place to mark runs
  right or wrong, ready answers to common questions, and the chat with Claude. A test keeps its rule equal to the
  Python one.
- RESULTS.md lists what the plan does not do and what is and is not new compared with published work.

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
- **Works where WebGL is off.** The 3D views on the showcase page, in report cards and in the app show a short note
  instead; before, the rest of the showcase page stayed empty in such browsers. The in-browser app also hides the
  charge and multiplicity inputs, which apply only without a SMILES.
- **Plain-language explainer** on the showcase page, and a roadmap of open problems in the README.

## 0.2.0 (2026-10-22)

- OSI Hackathon 2026 release: gated LLM-written rules, mirror repair, certified trust, stereo-aware leaderboard, mirror
  test, RDB7 data audit, MCP server, marimo app (also in the browser), and new DFT and GFN2-xTB calculations.
