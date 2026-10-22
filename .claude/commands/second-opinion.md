---
description: Triage machine-learned transition-state guesses and explain which to trust, repair or send to DFT
argument-hint: "[reaction ids, or 'red' / 'amber' / 'mirror-repaired']"
---
You are a careful computational chemist using the ts-second-opinion MCP tools.

Task: give a second opinion on the TS guesses for: $ARGUMENTS
(If that is a verdict name, call list_reactions with it and take the first 5 reactions.)

For each reaction:
1. Call triage_reaction (do not use score_against_reference for the decision).
2. Explain the verdict in two or three plain sentences: which rules failed, whether a mirror-image repair was applied
   and why that is safe (interatomic distances unchanged), and what the model disagreement says. If it was mirror-repaired,
   say that the unrepaired guess has exactly the same energy and frequencies as the repaired one (it is the TS of the
   mirror-image reaction), so a DFT frequency check would not have caught it; only the stereo rule did.
3. For amber or red verdicts, call write_orca_input and give the file path.
4. If Valency tools are available (names containing "valency" or "semantic_search_papers"), run one
   semantic_search_papers query with source=chemrxiv about the reaction type or the issue found (for example
   "stereochemistry errors in machine-learned transition states"), and cite at most two relevant preprints with links.
   Say clearly that these are related reading, not evidence for this specific reaction. Skip this if no such tools exist.

Finish with a short table: reaction id | verdict | action, then call score_against_reference=True for each reaction
and report, separately and labelled "evaluation only", whether the recommended guess matched the reference DFT TS.
Write the report to results/second_opinion_report.md.
