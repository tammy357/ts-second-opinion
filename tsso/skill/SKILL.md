---
name: ts-second-opinion
description: Check machine-learned transition-state (TS) guesses before spending DFT. Finds wrong-handed (mirror-image) and broken guesses, repairs whole mirror images, and says which guesses to trust, check with DFT or re-generate. Use when someone has TS guesses from a generator (GoFlow, React-OT, TSDiff or any other) as .xyz files, wants a batch of reactions triaged, wants a generator scored with a stereo-aware metric, or wants a plain-English chemistry rule turned into a tested check.
---

# TS Second Opinion

TS Second Opinion checks machine-learned transition-state guesses with plain-English chemistry rules that were turned
into code and kept only if they accept at least 99% of 9,540 real DFT transition states (RDB7). It repairs guesses
that are whole mirror images (same distances, opposite handedness), uses disagreement between several guesses to flag
uncertain ones, and returns a verdict per reaction:

- **green: trust.** Passes every rule and the guesses agree.
- **amber: check with DFT.** Passes the rules, but the guesses disagree a lot.
- **red: re-generate.** Still breaks a rule after repair.

## Run it

Install once: `pip install git+https://github.com/tammy357/ts-second-opinion` (or use the MCP server below).

One reaction, from three .xyz files (the reactant, the product and one or more TS guesses, all listing the same atoms
in the same order; multi-frame guess files are fine):

```bash
tsso check --reactant reactant.xyz --product product.xyz guesses.xyz --report report.html
tsso --json check --reactant reactant.xyz --product product.xyz guesses.xyz        # machine-readable
```

Without `--smiles`, bonds, bond orders and stereocentres are read from the reactant and product geometries (RDKit),
and the inferred atom-mapped SMILES is returned so a person can check it. Give `--charge` and `--multiplicity` for
charged or open-shell systems (defaults 0 and 1). With an atom-mapped SMILES (every atom, hydrogens too, mapped
1..N), pass `--smiles` instead.

Many reactions: a folder with one sub-folder per reaction (reactant.xyz, product.xyz, the guess files, optional
reaction.smi), or a CSV with columns id, reactant, product, guesses:

```bash
tsso batch my_reactions/ --out tsso_batch        # writes tsso_batch/verdicts.csv plus geometries and ORCA inputs
```

A fixed number of DFT runs to spend: `tsso plan tsso_batch --budget 5` (MCP tool plan_dft_budget) says which reactions
get a run (the ones whose guesses disagree most), which to trust, which to hold and which to re-generate, and writes
`dft_plan.csv` plus ORCA inputs that include the reaction path (IRC). It is a fixed rule: do not pick the reactions
yourself. Options: `--goal finish_ts` (the runs most likely to end in a usable TS; finds fewer errors) and `--audit N`
(N runs on trusted reactions drawn at random). After the runs: `tsso irc-check` per reaction (does the path join this
reactant to this product, with the right handedness?), then `tsso plan ... --results runs.csv`.

Score a generator's test-set predictions with and without handedness: `tsso evaluate preds.xyz refs.xyz reactions.smi`.

From Python: `from tsso import second_opinion; second_opinion(None, "reactant.xyz", "product.xyz", ["guesses.xyz"])`.

As MCP tools: `tsso mcp` (or `uvx --from git+https://github.com/tammy357/ts-second-opinion tsso mcp`). Tools include
triage_structures (one reaction; the SMILES is optional), triage_batch (a folder or CSV), evaluate_predictions,
gate_check and save_rule (turn a rule into a tested check).

## How to report results

1. Lead with the verdict and the reasons the tool gives. If a mirror-image repair was applied, explain that the
   unrepaired guess has the same energy and frequencies (it is the TS of the mirror-image reaction), so a DFT
   frequency check would not have caught it; only the stereo rule did.
2. For amber or red verdicts, point to the recommended geometry and the ORCA input.
3. Without a SMILES, show the inferred SMILES and suggest the person checks the bonds it implies.
4. Offer the feedback_url ("Was this verdict right?"): it opens a GitHub issue form and shares no structures.

## Keep the claims honest

- Never use a reference TS (score_against_reference) to make a decision; it is for evaluation only.
- Keep geometry errors (D-MAE > 0.2 Å) and handedness errors (a stereocentre the wrong way round) separate. The rules
  fix handedness (17.7% → 2.3% of final answers on RDB7), not geometry (14.1% → 11.7%).
- The rules were tested on small gas-phase H/C/N/O molecules (RDB7) with GoFlow guesses. For other chemistry, test
  them on the person's own reference TSs first (`tsso gate`).
- The amber cut-off was calibrated on GoFlow/RDB7. For another generator, use batch_disagreement_rank in batch output,
  or recalibrate.
- A DFT budget plan finds likely errors; it cannot see an error every guess shares, and it was tested on GoFlow/RDB7
  only (quote plan_evidence as benchmark results, not predictions). Reactions on hold are not trusted.
- An error rate for the trusted reactions comes only from random audits (review_dft_results) or from a calibration set
  of the person's own. Quote the audit bound with its confidence level, say that it is a bound and that few audits
  make it wide, and never derive a rate from the targeted runs: they were picked because they looked doubtful.
- Tell the person to confirm each optimised TS with an IRC calculation and irc_check: one imaginary frequency is not
  enough, and without the SMILES the check cannot see a mirror-image TS. The ORCA input with IRC has not been run by
  the authors; suggest trying one first.
- Certified trust (conformal selection) is an average guarantee that holds only if new reactions resemble the
  calibration set; say so whenever you quote it.
