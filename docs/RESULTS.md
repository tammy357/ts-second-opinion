# Detailed results

Everything behind the numbers in the [README](../README.md). The main results are retrospective: GoFlow's released predictions for the 1,194 RDB7 test reactions, with rules gated only on the 9,540 training reactions. New DFT transition-state optimisations (16 reactions) and GFN2-xTB ones (106 reactions) test the two central claims with new calculations; see [New DFT calculations](#new-calculations-optimising-the-guesses-with-dft) and [New xTB calculations](#new-calculations-optimising-the-guesses-with-gfn2-xtb). Every number is reproduced by the scripts listed in the README.

**Contents:** [Results on RDB7](#results-on-rdb7) · [Geometry-only input](#geometry-only-input-reading-the-reaction-from-its-structures) · [New DFT calculations](#new-calculations-optimising-the-guesses-with-dft) · [New xTB calculations](#new-calculations-optimising-the-guesses-with-gfn2-xtb) · [Is it only GoFlow?](#is-it-only-goflow) · [Closest prior work](#closest-prior-work-found-with-valency-mcp-bond) · [Partial flips](#can-partial-flips-be-repaired) · [Data audit](#data-audit-where-rdb7s-smiles-and-geometries-disagree) · [Leaderboard](#full-stereo-aware-leaderboard) · [Limitations](#limitations-in-detail)

## Results on RDB7

![final answers wrong, by kind](../results/fig_final_answer.png)
![certified trust](../results/fig_certify.png)

### What we found

Tested retrospectively on GoFlow's released predictions for the **1,194 RDB7 test reactions**
(8 GoFlow models trained with different seeds). Rules were written, tuned and gated only on the
**9,540 RDB7 training reactions**. We count two kinds of error against the reference DFT TS, because they
matter differently: **geometry** (interatomic-distance error, D-MAE, above 0.2 Å), which affects the barrier,
and **handedness** (a stereocentre the wrong way round), which usually leaves the energy unchanged but makes it
the transition state of a different reaction.

| | |
|---|---|
| **Handedness** | |
| GoFlow guesses with a stereocentre the wrong way round | **17.7%**, and 85% of them still count as accurate by D-MAE |
| …that are whole mirror images, with the same energy and frequencies as the right TS | **90%**; the other 10% flip only some stereocentres (a diastereomer, with a different energy) |
| …caught by two LLM-written stereo rules | **93%**, with 0.2% false alarms |
| Final answers with the wrong handedness: one model → after rules and mirror-image repair | **17.7% → 2.3%** |
| **Geometry** | |
| Final answers with a geometry error: one model → consensus of 8 → after the rules | 14.1% → 11.3% → 11.7% (the rules do not fix geometry) |
| Geometry errors among reactions marked trust / check with DFT / re-generate (80% / 20% / 0.5% of reactions) | **4.4% / 40% / 67%** |
| **Certified trust:** reactions trusted at a 5% target average geometry-error rate (assumes new reactions resemble the calibration set) | **80%** (actual 4.9%, averaged over 2,000 calibration splits) |
| **Both kinds together** | |
| Final answers wrong in either way: one model → TS Second Opinion | 29% → 13%, with no DFT |

Two honest negatives:

- **Shape rules (clashes, bond lengths, flatness) rarely fire** on GoFlow's guesses (under 1.2% each)
  and add nothing beyond model disagreement for geometry errors (AUROC 0.86 alone, 0.85 with the rules). GoFlow's
  shapes are physically sensible; its geometry errors are subtle, and disagreement is the better signal for them.
- The planarity rule never fires on these small RDB7 molecules (at most 7 heavy atoms). It is kept
  because it passed the gate, but it does no work here.

### Benchmark finding: how far the usual TS metrics miss mirror images

This blind spot is known. FK-Flow (Mark, Galustian, Kovar & Heid, 2025) notes that "the chirality error mode cannot be seen in
the mean absolute error of the inter-atomic distances (DMAE)", reports GoFlow chirality errors of 9.0% (persistent) and 27.8%
(thresholded) by its definitions, and cuts them to 0.7% and 5.8% by Feynman–Kac steering during sampling, which needs the
model and a g-xTB potential. Our handedness measure is defined differently (every checked stereocentre against the reference
TS), so the numbers are not directly comparable. What we add: where the blind spot sits in shared evaluation code, how much it
inflates reported success, and that most of these errors would pass a DFT frequency check too.


The RMSD code shared by OA-ReactDiff, React-OT and GoFlow is called with `ignore_chirality=True`
(it takes the better of the guess and its mirror image), and D-MAE compares interatomic distances, which a
mirror image keeps exactly. So a TS with the wrong handedness can score as perfect:

| | |
|---|---|
| Mirror image of the reference TS, for the 447 test reactions with stereocentres | D-MAE **0.000 Å**, RMSD with mirror allowed **0.000 Å** (proper-rotation RMSD: median 1.50 Å) |
| Wrong-handed GoFlow guesses that would count as accurate (D-MAE ≤ 0.2 Å) | **85%** (median D-MAE 0.085 Å, vs 0.081 Å for right-handed guesses) |
| GoFlow "success rate": D-MAE ≤ 0.2 Å only → also requiring correct stereo | **85.9% → 70.8%** |

### Why handedness errors slip past distance scores and even a DFT frequency check

Energy-based checks miss these errors too. 90% of GoFlow's wrong-handed guesses (1,514 of 1,690) flip every
stereocentre, so the guess is the mirror image of a correctly handed one. Mirror images have identical energies,
forces and vibrational frequencies, so a DFT TS optimisation and frequency check treats the guess exactly as it
treats the right one: it finds a genuine saddle point with one imaginary mode, but for the mirror-image reaction. An IRC from
it leads to the mirror images of the intended reactant and product, so it is caught only if the endpoints are compared with
their stereochemistry. The barrier is right; the stereochemistry is not. That matters whenever stereochemistry
does (stereoselectivity, chiral catalysts or enzymes, stereo-labelled datasets, mechanisms). The other 10% flip only some
stereocentres. They correspond to a diastereomer, whose energy does differ, and mirroring cannot fix them, so they are sent
for re-generation. A rule that compares handedness with the reaction as written catches both cheaply, before any DFT is spent.

We suggest reporting stereo correctness alongside D-MAE and RMSD, and judging confidence estimates with
risk–coverage curves (how often the answers you choose to trust are wrong):

![risk-coverage](../results/fig_risk_coverage.png)

| Share of trusted answers that are wrong | trust 50% | trust 80% | trust all | AURC |
|---|---|---|---|---|
| Geometry errors only: one GoFlow model, ranked by its disagreement with the other 7 | 2.5% | 6.2% | 14.1% | 0.039 |
| Geometry errors only: TS Second Opinion (rules + repair + disagreement) | **2.2%** | **4.6%** | **11.7%** | **0.032** |
| Geometry or handedness: one GoFlow model, ranked by disagreement | 18.9% | 22.6% | 29.2% | 0.205 |
| Geometry or handedness: TS Second Opinion | **3.9%** | **6.4%** | **13.5%** | **0.051** |

For geometry alone, most of the ranking power comes from model disagreement. The large gap in the last two rows comes
from handedness errors, which a model's disagreement cannot see.

### Certified trust: a stated error rate, not just a score

A verdict of "trust" is only useful if you know how often it is wrong. Conformal selection (Jin & Candès, *JMLR* 2023)
turns the triage into a guarantee. Calibrate on reactions whose answers you know (for example, from earlier DFT). For each new
reaction, compare its trust score (rules first, then model disagreement) with the calibration reactions that turned out
wrong, and trust the reactions picked by the Benjamini–Hochberg procedure. If new reactions are drawn like the calibration ones,
the expected share of trusted reactions that are wrong is at most the target, whatever the model or the sample size
(`tsso/certify.py`, `scripts/certify.py`, `tsso certify --target 0.05`, MCP tool `certify_trust`).

Tested by splitting the 1,194 test reactions at random into calibration and new halves 2,000 times:

| Target (average error rate among trusted reactions) | 3% | 5% | 10% |
|---|---|---|---|
| Geometry errors: reactions trusted / actual error rate | 60% / 2.9% | **80% / 4.9%** | 97% / 9.9% |
| Geometry or handedness: reactions trusted / actual error rate | 19% / 2.3% | 63% / 4.9% | 93% / 9.9% |

What the guarantee does and does not say:

- **It is an average.** In a single calibration draw the rate can be higher: for a 5% target, 95% of draws stayed at or below 7.4%.
- **The correction matters most at strict targets.** A naive cut-off (trust everything scoring above the point where
  the calibration error rate is under target) does about as well on average here, but it overshoots at strict targets
  (3% target → 4.1% actual, geometry or handedness) and comes with no guarantee.
- **It assumes nothing changed.** When we certified guesses from a different ensemble (5 GoFlow runs of four model sizes) using the
  8-seed calibration, the error rate stayed under target (2.5% at a 5% target) but only 53% of reactions were trusted;
  recalibrating on that ensemble restored 74% at 4.8%. A change in the other direction could break the guarantee, so
  recalibrate whenever the generator, the data or the rules change.
- **The trust score was designed with this test set in view** (see Limitations), which the theory does not allow; the
  repeated random splits are the practical check.

We did not find prior work applying conformal selection to machine-learned TS guesses, but our search was brief; conformal
methods have been used for generative molecular design ([ConfHit](https://arxiv.org/abs/2603.07371),
[SCoRE](https://arxiv.org/abs/2603.24704)).

### Where should the DFT go? Testing the budget plan

`tsso plan` (MCP tool `plan_dft_budget`, app section 3b) takes the verdicts for a set of reactions and the number of
DFT runs you can afford. The rule is fixed: re-generate what still breaks a rule; give the runs to the reactions whose
guesses disagree most; hold the other "check with DFT" reactions; trust the rest. It has no fitted parameters and never
sees a reference TS. Using model disagreement to decide where an expensive calculation goes is borrowed from active
learning ("query by committee", [Seung, Opper & Sompolinsky 1992](https://doi.org/10.1145/130385.130417); for learned
potentials, [Schran, Brezina & Marsalek 2020](https://doi.org/10.1063/5.0016004)). There it chooses training data;
here it chooses which answers to check.

**The test** (`scripts/plan_evaluation.py`, `results/plan_evaluation.json`). Afterwards, we used the reference TSs to
count how many picked reactions really had a wrong recommended guess (D-MAE above 0.2 Å), since those are the ones a
DFT run with a reaction-path check would expose. In batches the size a person brings (random draws from the 1,194 test
reactions, 5,000 draws each, 5 runs):

| Reactions in the batch | Picks with a wrong geometry (mean of 5) | Picked at random | Perfect order | Batches where no pick was wrong |
|---|---|---|---|---|
| 20 | 1.61 | 0.57 | 2.26 | 16% |
| 40 | 2.25 | 0.58 | 3.98 | 5% |
| 100 | 3.06 | 0.57 | 4.99 | 1% |

Over all 1,194 reactions at once (8 GoFlow seed models; 11.5% of the reactions that can get a run have a wrong
geometry):

| DFT runs | Picks with a wrong geometry | 95% interval for that share | Picked at random | Any "check with DFT" reaction | Perfect order |
|---|---|---|---|---|---|
| 5 | 5 (100%) | 60–100% | 0.57 | 1.98 | 5 |
| 10 | 9 (90%) | 70–100% | 1.14 | 3.97 | 10 |
| 25 | 21 (84%) | 64–100% | 2.86 | 9.92 | 25 |
| 50 | 32 (64%) | 52–82% | 5.72 | 19.83 | 50 |
| 100 | 52 (52%) | 42–63% | 11.45 | 39.66 | 100 |

"Any check-with-DFT reaction" is what you would get by working through that verdict in no particular order, so the
ranking within the verdict is what adds most at small budgets. Picking the least disputed first found 0 wrong
guesses at every budget. With a different ensemble (five GoFlow runs of four model sizes, same reactions) the first
10 picks held 6 wrong guesses and the first 50 held 33 (against 9 and 32 above), so the gain does not depend on
one set of seeds, though the top of the order is noisier.

**Five limits of the first version, and what was done about each.** All numbers are from
`results/plan_evaluation.json`.

*1. Errors that all the models share.* Model disagreement measures spread, not bias, so the trusted reactions keep a
4.4% geometry-error rate (6.2% counting handedness) whatever the budget. The same has been reported for learned
potentials ([Knowing when to trust machine-learned interatomic potentials, 2026](https://arxiv.org/abs/2605.00640)).
We looked for a reference-free signal that ranks the wrong guesses hidden among trusted reactions. Thirteen were
tried on the test reactions; the best, how far the guess sits from the nearer of the reactant and the product
(D-MAE), reached an AUROC of 0.746 among trusted reactions against 0.714 for disagreement. Adding it to the
order by rank found 3, 7 and 16 wrong guesses in the first 5, 10 and 25 picks against 5, 9 and 21 for
disagreement alone, so it was **not** adopted for the main order. What the plan offers instead is a **random audit**
(`--audit N`): N of the runs go to trusted reactions drawn by lot, and their outcomes give an exact upper bound
(hypergeometric, as in acceptance sampling) on how many of the other trusted reactions are wrong. It needs no
assumption about the models or the chemistry, only that the draw was random and that every audit is reported.

| Reactions in the batch | Random audits | Bound held (claimed: at least 95%) | Mean bound on the error rate | True error rate of the rest |
|---|---|---|---|---|
| 100 | 10 | 100.0% | 35.5% | 6.2% |
| 100 | 30 | 99.3% | 21.3% | 6.2% |
| 300 | 30 | 99.9% | 19.4% | 6.2% |
| 300 | 60 | 97.9% | 14.9% | 6.2% |

The bound held as often as it should. It is wide: about 34 audits that all come back right are needed to bound
100 trusted reactions at 10%, so with five runs an audit says almost nothing, and the plan says so. Every audit is
also a run not spent on a doubtful reaction. Choosing test points under a labelling budget, and the bias that
choosing them non-randomly introduces, is the subject of active testing
([Kossen, Farquhar, Gal & Rainforth 2021](https://arxiv.org/abs/2103.05331)); the audit here is the simplest
unbiased version of that idea, not a new method.

*2. Finding errors is not the same as finishing TSs.* A run started from a wrong guess rarely reaches the right TS:
in the GFN2-xTB tests, "check with DFT" guesses with a wrong geometry reached the reference TS in 1 of 9 runs, those with
a right geometry in 17 of 21. There is now a second goal, `--goal finish_ts`, which gives the runs to the *least*
disputed "check with DFT" reactions. The two goals pull in opposite directions, and the person has to choose:

| DFT runs | Wrong guesses found: find_errors | finish_ts | Picks with a right guess: find_errors | finish_ts |
|---|---|---|---|---|
| 5 | 5 | 2 | 0 | 3 |
| 10 | 9 | 3 | 1 | 7 |
| 25 | 21 | 6 | 4 | 19 |
| 50 | 33 | 14 | 17 | 36 |
| 100 | 53 | 31 | 47 | 69 |

"Right guess" counts geometry and handedness. It is a proxy for a run that finishes in the right TS, not a count of
finished runs: no DFT was run for this table.

*3. A frequency check is not enough.* Of the 62 xTB optimisations that ended with exactly one imaginary frequency, 21
were not the reference TS (5 trust, 11 check, 5 re-generate; some of this is xTB differing from DFT). The plan now
writes ORCA inputs that run the TS optimisation, the frequencies and the reaction path (IRC) in one job, and
`tsso irc-check` (MCP tool `irc_check`, `tsso/irc.py`) compares the two ends of the path with the reactant and the
product. Comparing the bonding at IRC ends is standard in automated workflows
(for example [this 2020 PeerJ Physical Chemistry workflow](https://peerj.com/articles/pchem-15/), which compares the adjacency matrices of the IRC ends). The addition here is handedness: a mirror-image TS
has the same bonds at both ends of its path, so a bonds-only check passes it, and at least one recent benchmark
removes stereochemistry before comparing the ends ([TSBench, 2026](https://arxiv.org/abs/2609.08503)). On stand-in
endpoints (the RDB7 reactant and product geometries, not IRC output) the check called 1194 of 1194 true pairs
right in either order, 1194 of 1194 paths that return to the reactant wrong, and 676 of 676 mirror-image pairs
wrong when given the SMILES, against 0 of 676 without it. With 0.03 Å of noise on every coordinate, 16 true pairs were
called wrong, so bonds read from distances are fragile near a threshold. **Neither the ORCA input nor the check has
been run on a real IRC calculation**: the input follows the ORCA manual, which allows IRC together with OptTS and
Freq, and real path ends are less tidy than optimised minima.

*4. One generator family, one dataset.* Still true, and no new DFT was run. Two things were added from the data in the
repository. First, fewer guesses (random subsets of the 8 seed models; mean and range over three subsets):

| Guesses per reaction | Wrong among the first 5 picks | first 10 | first 25 |
|---|---|---|---|
| 2 | 2.3 (1–3) | 5.0 (2–7) | 14.3 (11–16) |
| 3 | 3.3 (2–4) | 7.7 (5–9) | 17.7 (13–20) |
| 4 | 4.7 (4–5) | 8.7 (8–9) | 19.0 (15–22) |
| 6 | 4.7 (4–5) | 9.0 (7–10) | 21.0 (18–23) |
| 8 (all) | 5 | 9 | 21 |

Two guesses already help; four give most of the benefit. Second, the harder reaction-core and barrier splits have one
GoFlow model each, so there is no disagreement there. With one guess per reaction the plan now ranks by the endpoint
distance described above. On reactions outside the random-split test set, and with a different trained model, the
first picks held: reaction-core split, 4 of 10 (random 2.0), 11 of 25 (random 5.1), 24 of 50 (random 10.1), 40 of 100 (random 20.2); barrier split,
5 of 10 (random 2.5), 13 of 25 (random 6.1), 23 of 50 (random 12.3), 46 of 100 (random 24.5). That is about twice chance, far weaker than disagreement, with wide
intervals at 10 picks; it is a fallback, and the plan labels it as the weaker sign. The signal was chosen on the
random-split test reactions, so these two splits are its only fair test.

*5. An error rate for your own reactions.* Certified trust (conformal selection,
[Jin & Candès 2023](https://jmlr.org/papers/v24/22-1176.html)) holds only if new reactions are exchangeable with a
calibration set, and a visitor's reactions have none. There are now two honest routes. If you have reactions with
known answers from the same generator, `tsso plan ... --target 0.05 --calibration known.csv` certifies against them
(reactions that still break a rule are never certified). If you do not, the random audit above is the only source of
an error rate, and `tsso plan ... --results runs.csv` reports it after the runs. Without either, the plan states that
it gives no error rate.

**Is this new?** Partly, and less than first claimed. Using a learned confidence score to decide which generated TSs
go to DFT is published: OA-ReactDiff pairs its generator with a confidence model and sends 14% of reactions to DFT
([Duan, Du, Jia & Kulik 2023](https://arxiv.org/abs/2304.06174)). Ranking by committee disagreement is query by
committee ([Seung, Opper & Sompolinsky 1992](https://doi.org/10.1145/130385.130417)), random audits are acceptance
sampling, and bond-based IRC checks are routine. What we did not find in a short search (about ten queries, not a systematic review), and so describe as
"not found" rather than "first": (a) a budget plan that works after generation on any model's guesses with no
trained confidence model, tested for how many wrong guesses its picks contain; (b) an exact, assumption-free bound on
the error rate of the *accepted* TS guesses from a random audit inside the same budget; (c) a reaction-path check
that compares handedness as well as bonds; (d) an agent whose numbers are copied from tool results by reference. Each
is a small step, and (d) has not been tested with a live model.

**The planning agent** (`tsso agent`, the chat in app section 3b, `tsso/agent.py`) lets a person ask for a plan in their
own words. Three properties are enforced in code and covered by tests (`tests/test_plan.py`): Claude can only call the
planning rule and read-only tools, so it cannot choose the reactions; tool results are stripped of coordinates, SMILES
and ORCA inputs before Claude sees them (verdicts, rule names, atom labels such as "C3–H7" and disagreement values are
sent); and Claude is told not to type numbers but to point at a field of a tool result, for example
`{{1:counts.hold}}`, which the program replaces with the value. A pointer to a field that does not exist is shown as
such. Every number Claude types anyway is listed under the answer, split into "also present in a tool result" and
"found nowhere". The first group is a weak check: once the benchmark evidence has been read, about 31% of arbitrary
whole numbers below 100 match something by chance, which is why pointing replaced matching. Three limits remain. A
copied value can be the wrong field for its sentence. Whether Claude keeps to the pointer format has not been
measured, because no live run could be made where this was built (no API key). And no answer has been scored for
quality; the plan table is the record.

### New calculations: optimising the guesses with DFT

Everything above is retrospective. To test the two central claims at the level of theory that matters, we optimised TS
guesses to saddle points with DFT: ωB97X/def2-SVP in PySCF (density fitting; restricted or unrestricted Kohn-Sham from
the reaction's charge and multiplicity), Sella (order 1, internal coordinates, fmax 0.05 eV/Å, at most 80 steps), then
analytic Hessians and harmonic frequencies (`scripts/dft_prospective.py`; one JSON file per reaction in `results/dft/`,
summary in `results/dft_summary.json`). RDB7's reference TSs are ωB97X-D3/def2-TZVP, so ours is the same functional
family with a smaller basis and no dispersion. A guess "found the reference TS" if it reached a first-order saddle point
within 0.2 Å (D-MAE) of the RDB7 TS with the right handedness. Each optimisation took 3–45 minutes on 1 CPU core.

![new DFT calculations](../results/fig_dft.png)

**A mirror image passes a DFT frequency check.** Reaction 3509's best guess is a whole mirror image, which the rules
repair by reflection. We optimised both. Both converged to first-order saddle points: 654i and 656i cm⁻¹, energies
0.0006 kcal/mol apart, both 0.018 Å (D-MAE) from the reference. The unrepaired one kept the wrong handedness; the
repaired one is the reference TS. At DFT level the two are indistinguishable by energy, frequencies and distances; only
a stereo check tells them apart.

**The verdicts hold at DFT level.** We optimised the recommended guess for the first 6 ✓ trust reactions by id (not
mirror-repaired; those are covered above), the first three ! check with DFT reactions by id, and all 6 ✕ re-generate reactions in the test set.

| Verdict | Found the reference TS | First-order saddle point (passes a frequency check) | Wrong handedness after optimisation | Median D-MAE to the RDB7 TS |
|---|---|---|---|---|
| ✓ trust | 6 of 6 | 6 of 6 | 0 of 6 | 0.018 Å |
| ! check with DFT | 1 of 3 | 3 of 3 | 0 of 3 | 0.206 Å |
| ✕ re-generate | 0 of 6 | 5 of 6 | 5 of 6 | 0.258 Å |

Every trusted guess converged in 5–14 steps to the reference TS (D-MAE 0.009–0.026 Å). No re-generate guess
did, yet 5 of the 6 still converged to a first-order saddle point, so a DFT frequency check would have accepted them.
All of them kept a stereocentre the wrong way round: reaction 3923 (783i cm⁻¹, D-MAE 0.25 Å); reaction 3968 (609i cm⁻¹, D-MAE 0.27 Å); reaction 5144 (587i cm⁻¹, D-MAE 0.31 Å); reaction 5877 (674i cm⁻¹, D-MAE 0.18 Å); reaction 8502 (490i cm⁻¹, D-MAE 0.21 Å). The rules flagged all of them before any DFT was spent. The ! check with DFT guesses sit in between, as the verdict intends: all three became first-order saddle points with the
right handedness, one reached the reference TS (reaction 28, 0.037 Å) and two ended just outside the 0.2 Å cut-off
(reactions 72 and 24, 0.206 and 0.231 Å). These are the cases the verdict sends to DFT, and DFT is what settles them.

### New calculations: optimising the guesses with GFN2-xTB

Next, a cheaper level on a larger sample. We optimised TS guesses to
saddle points with GFN2-xTB (tblite) and Sella (order 1, internal coordinates, fmax 0.03 eV/Å, at most 400 steps), then
computed numerical frequencies (`scripts/xtb_prospective.py`; results in `results/xtb_prospective_*.json`). GFN2-xTB is
a semi-empirical level, not DFT, and the reference TSs are DFT, so these are supporting evidence, not a DFT validation.

**A mirror image passes a frequency check.** For the first 40 test reactions (by id) whose best guess was mirror-repaired,
we optimised both the unrepaired and the repaired guess. In 34 of 40, both reached a first-order saddle point (one
imaginary frequency) with the same energy (within 0.01 kcal/mol) and the same imaginary frequency (within 6 cm⁻¹). After
optimisation, the unrepaired guess had the wrong handedness in 40 of 40 and the repaired guess the right one in 40 of 40.
In the other 6, at least one optimisation stopped before converging or the two found different saddle points. So the
symmetry argument holds in practice: an optimiser and a frequency check accept the mirror image, and only a stereo check
catches it.

**The verdicts hold after optimisation.** For the first 30 reactions (by id) of each verdict (all 6 for re-generate;
mirror-repaired reactions excluded, as they are covered above), we optimised the recommended guess. A guess "found the
reference TS" if it reached a first-order saddle point within 0.2 Å (D-MAE) of the DFT TS with the right handedness.

| Verdict | Found the reference TS | First-order saddle point | Median D-MAE to the DFT TS |
|---|---|---|---|
| ✓ trust | 23 of 30 (77%) | 28 of 30 | 0.07 Å |
| ! check with DFT | 18 of 30 (60%) | 29 of 30 | 0.18 Å |
| ✕ re-generate | 0 of 6 (0%) | 5 of 6 | 0.32 Å |

A cheap optimisation does not rescue bad geometry: of the 13 sampled guesses with a geometry error before optimisation
(D-MAE > 0.2 Å), only 1 reached the reference TS after GFN2-xTB optimisation, so sending flagged geometry to DFT, not
to a quick xTB fix, is the right advice. Five of the six "re-generate" guesses still optimised to a saddle point with one imaginary frequency, so a frequency check
would have accepted them, but none was the reference TS and five had the wrong handedness. The samples are small and
deterministic (first by id), not random.

### Is it only GoFlow?

No. On Transition1x, a different dataset, the ReCurveflow study (Baek, Gim & Kang, arXiv:2608.20869, 2026; their Table 1,
native split) reports reflection-invariant RMSD against chirality-preserving RMSD for eight generators. Three lose
0.15–0.35 Å once mirror images are not allowed: GoFlow 0.33 → 0.61 Å, MolGen 0.27 → 0.62 Å and TSDiff 0.41 → 0.56 Å.
Their headline metric is still the reflection-invariant one. Neither that study nor the 2026 robustness study of React-OT
(Darouich et al., arXiv:2601.16469) checks, repairs or filters TS guesses for handedness after generation.

### Closest prior work, found with Valency MCP Bond

A literature scan with the Valency MCP Bond server (semantic search across ChemRxiv, arXiv and PubMed, 2023 onwards)
found stereo-aware *generation*: RitS (Nikitin, Anstine & Isayev, ChemRxiv 2026, doi:10.26434/chemrxiv.15001681)
encodes product chirality in the input graph so that the model generates the intended diastereomeric TS, and FK-Flow
steers GoFlow at inference time. Both act inside one model. It also found a 2026 review of TS evaluation practice (Sain,
ChemRxiv 2026, doi:10.26434/chemrxiv.15009215) asking what reported accuracies certify, and a validation-focused TS
generator (Kairos-GFM, bioRxiv 2026) showing that low RMSD does not guarantee the intended reaction channel. We found no
tool that checks, repairs and triages the handedness of any generator's TS guesses after generation, gives a stated
error rate for the ones to trust, or tests TS metrics for mirror blindness.

### Can partial flips be repaired?

Mostly not usefully. We tried a local "umbrella" inversion: reflect the wrong stereocentre and its smallest terminal
branch (usually a hydrogen) through the plane of its other neighbours, which inverts that centre and keeps every bond to
it exactly. Applied to single GoFlow guesses (9,552), it fixed the handedness of 40 of the 49 guesses it could edit and
broke none, cutting wrong-handed single guesses from 3.3% to 2.8%. But 41 of the 49 already had a geometry error
(D-MAE > 0.2 Å), and their mean D-MAE rose from 0.26 to 0.33 Å. In the 8-model pipeline only 5 final answers had a
detectable wrong centre. Partly flipped guesses are usually broken in other ways too, so the verdict stays "re-generate",
and the inversion is not part of the pipeline. The app offers it as an experimental button on partly flipped
guesses (`tsso/stereo_edit.py`), with the rules checking the result: on example 5877 the inverted centre comes out
right but the edit creates a clash, so the verdict stays "re-generate".

### Data audit: where RDB7's SMILES and geometries disagree

The checks that test TS guesses assume the atom-mapped SMILES and the reactant and product geometries describe the same
molecules. `tsso audit` tests that assumption directly: a structure is flagged when two atoms sit at single-bond distance
(below 1.10x the sum of their covalent radii) but the SMILES does not bond them, or when the SMILES bonds two atoms that
are more than 1.60x apart. In **51 of the 10,734 RDB7 reactions** used here (10 test, 41 training; 0.5%), and only ever in
the **product** geometry, the two disagree: for example, in reaction 6295 the SMILES puts H11 on C2, but in the
geometry it sits 1.08 A from C3 and 2.1 A from C2. Reaction 1093, one of the app's worked examples, is flagged (an N5-O6
bond at 1.45 A that the SMILES says is broken). The thresholds are deliberately conservative, so borderline cases (such
as a 1.69 A C-C contact) are not counted. The full list is in `results/rdb7_audit.csv`
(`python scripts/audit_rdb7.py`). The same check runs on every uploaded reaction and shows a warning when it fires. We
do not know whether these entries reflect the SMILES or the geometry being wrong; either way, a model trained or scored
on them sees an inconsistent target.

### Does it generalise?

All of these use predictions the GoFlow authors released; nothing was retrained.

| Setting | Stereo errors | Caught by the stereo rules | Wrong answers (geometry or handedness) before → after |
|---|---|---|---|
| 5 GoFlow runs of four model sizes (64–512 hidden features; random split), used as an ensemble | 15–18% per model | 92% (0.2% false alarms) | 31% → 14% (AUROC 0.70 → 0.91) |
| 1 GoFlow model, reaction-core split (harder) | 16% | 91% (0.1% false alarms) | 33% → 22% (rules + mirror only, no ensemble) |
| 1 GoFlow model, barrier-height split (harder) | 17% | 89% (0.1% false alarms) | 36% → 27% (rules + mirror only, no ensemble) |
| No ML: midpoint of aligned reactant and product | 1% | – | rules flag 43% of guesses; 76% of flagged vs 38% of unflagged are wrong |

So the stereo failure is not one unlucky model: every released GoFlow model size and split has it at about the same rate.
Caveat: the rules were gated on the random-split training set, which contains about 80% of the test reactions of
the two other splits. The gate only checks that a rule accepts real TSs, and the stereo rules have no tuned thresholds
beyond small signed-volume cut-offs, but these split results are not fully independent of the gate.

### Does it matter who writes the check?

For the same six rules we compared Claude's checks, written from the plain-English sentence, with drafts written directly
as Python during development (by the author with an AI coding assistant, before the rule-to-code pipeline existed).
Both versions passed the gate, and they perform about the same on GoFlow's test guesses (`scripts/ablation_authors.py`):

| Final answers after triage | Claude-authored checks | Directly written drafts |
|---|---|---|
| Geometry wrong | 11.7% | 11.6% |
| Handedness wrong | 2.3% | 1.4% |
| Product-stereo check: real TSs accepted / handedness errors caught | 99.8% / 87% | 99.1% / 93% |

So Claude's checks are not better than checks written directly as code. What the LLM adds is access: a chemist who does
not write Python gets a tested check from one sentence. What makes any check trustworthy, however it was written, is the gate.

### What did not work, and why that matters

- **Shape rules** (see above) rarely fire on GoFlow guesses and add nothing for geometry errors.
- **Claude proposing its own rules.** Asked for 8 new rules (`scripts/discover_rules.py`), Claude produced plausible ones
  (bond-angle ranges, torsions, near-linear hydrogen transfer, balanced progress of forming and breaking bonds). All 8 passed the gate, some after rewrites,
  but none improved triage on the held-out half of the test reactions, so none were added. The gate plus held-out
  selection is what stops plausible-sounding but useless rules from creeping in.
- **Cheap physics single points.** GFN2-xTB energies of the reference TSs track the DFT barriers well (Spearman 0.88),
  but on raw guesses the xTB "barrier" and RMS gradient flag bad geometries poorly (AUROC 0.63 and 0.65), adding them to
  disagreement makes ranking worse (0.79–0.81 vs 0.86), and choosing the lowest-energy guess is worse than the consensus.
  In a quick 5-reaction trial, xTB TS optimisation (Sella) moved 4 guesses further from the DFT TS, so xTB is not a good
  stand-in for the DFT check here (`scripts/xtb_check.py`).

### Geometry-only input: reading the reaction from its structures

Since version 0.3, a reaction can be checked from three `.xyz` files without an atom-mapped SMILES (`tsso/infer.py`).
Bonds come from interatomic distances (a pair is bonded below a cut-off of 1.13 to 1.30 times the sum of covalent radii,
depending on the elements, with a valence check that drops the relatively longest contact of an over-bonded atom, as in
strained cages), bond orders from RDKit's xyz2mol for closed-shell structures, and stereocentres from RDKit's assignment
from 3D. The cut-offs were chosen on RDB7's 10,734 reactant and product pairs, so this is not an independent test of
them; it measures how closely the inferred reaction reproduces the dataset's own SMILES.

| On the 1,194 RDB7 test reactions | Same as with the dataset's SMILES |
|---|---|
| Bonds, reactant / product | 100.0% / 98.6% |
| Stereocentres, reactant / product | 99.3% / 97.1% |
| Verdict | **99.7%** |
| Recommended guess (same model and same mirror repair) | **98.7%** |

| Final answers, against the reference TS (evaluation only) | With the SMILES | From geometry alone |
|---|---|---|
| Wrong handedness | 2.3% | 2.6% |
| Geometry wrong (D-MAE > 0.2 Å) | 11.7% | 11.6% |
| Trusted (green) | 79.7% | 79.6% |
| Wrong among trusted | 6.2% | 6.3% |

The outcome differs in 15 reactions (1 of them flagged by the data audit), and in every one the
inferred bonds or stereocentres differ from the dataset's SMILES. Bond orders were assigned for
87.4% of reactions; for the rest (some carbenes and zwitterions, whose RDB7 SMILES are
charge-separated) every bond is kept single, so the planarity rule has nothing to check there. Exact bond orders match
the dataset's SMILES for 79.2% of reactions, mostly because of how aromatic and
charge-separated structures are written. Reproduce with `python scripts/validate_geometry_input.py`.

## Full stereo-aware leaderboard

| Method | Model size | Test set | D-MAE (Å) | RMSD, mirror allowed (Å) | Success (D-MAE ≤ 0.2 Å) | Correct handedness | **Success with correct handedness** |
|---|---|---|---|---|---|---|---|
| GoFlow, 8 seeds (mean ± s.d.) | 256 | random split (1,194) | 0.110 ± 0.003 | 0.356 ± 0.006 | 85.9% ± 0.5% | 82.3% ± 0.7% | **70.8% ± 0.6%** |
| GoFlow, size run 1 | 64 | random split (1,194) | 0.163 | 0.494 | 72.2% | 85.0% | **61.9%** |
| GoFlow, size run 2 | 128 | random split (1,194) | 0.124 | 0.403 | 82.5% | 83.8% | **69.4%** |
| GoFlow, size run 3 | 256 | random split (1,194) | 0.102 | 0.335 | 88.4% | 82.4% | **72.9%** |
| GoFlow, size run 4 | 256 | random split (1,194) | 0.105 | 0.343 | 86.1% | 83.1% | **71.3%** |
| GoFlow, size run 5 | 512 | random split (1,194) | 0.114 | 0.357 | 85.0% | 83.1% | **70.8%** |
| GoFlow, 1 model | 256 | reaction-core split (1,106) | 0.138 | 0.425 | 79.2% | 83.8% | **67.3%** |
| GoFlow + TS Second Opinion (rules + mirror repair, no ensemble) | 256 | reaction-core split (1,106) | 0.138 | 0.425 | 79.2% | 97.6% | **77.8%** |
| GoFlow, 1 model | 256 | barrier-height split (1,192) | 0.149 | 0.432 | 74.7% | 83.5% | **64.0%** |
| GoFlow + TS Second Opinion (rules + mirror repair, no ensemble) | 256 | barrier-height split (1,192) | 0.149 | 0.432 | 74.7% | 96.4% | **73.2%** |
| No ML: midpoint of aligned reactant and product | – | random split (1,194) | 0.234 | 0.488 | 45.6% | 98.8% | **45.5%** |
| GoFlow 8 seeds + TS Second Opinion (final answer) | 256 | random split (1,194) | 0.099 | 0.314 | 88.3% | 97.7% | **86.5%** |
| GoFlow 5 size runs + TS Second Opinion (final answer) | – | random split (1,194) | 0.104 | 0.327 | 87.7% | 97.7% | **86.1%** |

## How it works, in detail

```
plain-English rule ──► Claude writes check(rxn, xyz) ──► gate: accept ≥ 99% of 9,540 real TSs,
        ▲                                                   never crash, and catch ≥ 1% of
        └──────── feedback with failing examples ◄──────── mirrored or noisy TSs  (else rewrite)

8 GoFlow guesses ──► run every accepted check ──► consensus among passing guesses
                 ──► if a stereo rule still fails, try the mirror image (same distances)
                 ──► verdict:  ✓ trust  ·  ! check with DFT (top 20% model disagreement)  ·  ✕ re-generate
                 ──► ORCA TS-optimisation input for flagged reactions
                 ──► optional: certified trust, calibrated on reactions with known answers
```

- **Rules → code** (`tsso/author.py`, `tsso/llm.py`): Claude via the Claude Code CLI (`claude -p`, works
  with a Pro/Max plan) or the Anthropic API. Checks only see what is known before any quantum chemistry:
  reactant and product geometries, bonds and stereo tags. The reference TS and the model predictions are
  hidden, and code with imports, file access or double-underscore names is rejected (`tsso/sandbox.py`).
- **The gate** catches both failure modes of LLM-written code: too strict (rejects real TSs) and too loose
  (never fires). In the final run, 5 of 6 rules passed first time and `planarity` needed two rewrites
  (in an earlier run, `reacting_bonds` also needed one). Logs: `results/authoring_log*.jsonl`.
- **Mirror-image repair**: if every guess has a stereocentre the wrong way round, the mirror image keeps
  every interatomic distance (so geometric accuracy is unchanged) and flips handedness. It is applied only if
  the mirrored guess then passes every rule; it was applied to 80 reactions.
- **Disagreement** is the mean D-MAE between the chosen guess and the other 7 models. The review cut-off is a
  budget (the most disputed 20% of reactions), set from the disagreement values alone, never from reference answers.

## The agent layer in detail

TS Second Opinion is also an MCP server, so Claude Code, Claude Desktop or any MCP client can use it as tools:
`triage_structures` (your own reaction and guesses, from any generator), `evaluate_predictions`, `certify_trust`,
`list_reactions`, `triage_reaction`, `list_rules`, `get_rule_source`, `gate_check`, `save_rule`, `write_orca_input`,
`benchmark_summary`, `triage_batch`, `plan_dft_budget`, `plan_evidence`, `review_dft_results` and `irc_check`.

```bash
pip install -e .            # installs the `tsso` command
claude                      # in this folder: .mcp.json starts `tsso mcp` automatically
> /new-rule "A hydrogen that is being transferred sits roughly on the line between donor and acceptor."
> /second-opinion 3509 1093 239
```

In `/new-rule` the agent writes the check, calls `gate_check`, reads which real transition states it wrongly rejected,
revises and re-tests, then saves it. In our recorded run (`results/agent_new_rule_transcript.md`) the first version
rejected 36% of real TSs (it forgot that 1,2- and 1,3-hydrogen shifts are bent); the second passed. The agent then
reported, honestly, that the rule never fires on GoFlow guesses and adds little. `/second-opinion` triages reactions,
explains each verdict, writes ORCA inputs for flagged ones, and, if the Valency MCP server is connected, cites related
ChemRxiv preprints as background reading. In the recorded run (`results/agent_second_opinion_transcript.md`,
22 October 2026) Valency was connected and supplied the literature steps; the triage steps called the Python functions
behind the MCP tools directly. Reaction 3509 is a good example: all eight GoFlow models put C2 the wrong way round, and the
mirror repair gives a guess with D-MAE 0.019 Å and correct stereochemistry.

The same functions are available as a CLI: `tsso check`, `tsso evaluate`, `tsso certify --target 0.05`,
`tsso triage 3509 --score`, `tsso list --verdict red`, `tsso rules`, `tsso gate "<rule>" check.py`, `tsso author "<rule>"`,
`tsso orca 3509` and `tsso summary` (add `--json` to any of them). An agent like [TSAgent](https://arxiv.org/abs/2605.14154), which runs DFT TS searches
but has no check on ML guesses beforehand and no stereo check, could call these tools before spending DFT.

## Where foundation models fit

- **Language models** write and revise the checks and explain verdicts; they do not judge geometry directly.
- **Interatomic-potential foundation models** are the natural next stage after triage: a 2026 benchmark found
  MACE-OMol25 locating TSs with a 96.6% success rate using fewer than four DFT gradients per reaction
  ([arXiv:2604.00405](https://arxiv.org/abs/2604.00405)). A pipeline would be: generator → TS Second Opinion → MLIP refinement
  (MACE-OMol25 or UMA) for amber cases → DFT only where needed. We did not run an MLIP here (it needs PyTorch, which was not
  available in our build environment).
- **Steering**: because Feynman-Kac steering accepts non-differentiable potentials, the same gated checks could steer
  GoFlow during generation once model weights are available, instead of filtering afterwards.

## Limitations in detail

- **Retrospective, apart from new DFT and GFN2-xTB optimisations.** The headline numbers mean close to RDB7's reference TS
  geometry, not that a TS optimisation from the guess converges. The DFT check is small (16 reactions: the first by id for trust and check, all six
  re-generate, not a random sample; ωB97X/def2-SVP rather than RDB7's ωB97X-D3/def2-TZVP) and
  confirms, rather than proves, the symmetry argument (an achiral Hamiltonian gives mirror images identical energies).
- **Whether a handedness error counts as wrong depends on what you need.** For a barrier height alone, a whole mirror image
  gives the right number; the geometry-only rows are the fair comparison for that use.
- **GoFlow's released predictions were aligned to the reference before saving.** GoFlow's evaluation code rotates each
  sample onto the reference TS and relabels symmetry-equivalent atoms (for example, the hydrogens of a methyl group) to
  match it, then saves the sample closest to the median. Our disagreement scores inherit this relabelling, which may make
  models look slightly more consistent than they would without a reference. Structures you triage yourself are not relabelled.
- **One dataset, one model family.** Small gas-phase molecules (H, C, N, O; at most 7 heavy atoms), GoFlow runs of
  four model sizes and two extra splits, plus a naive interpolation baseline. Other generators (React-OT, TsDiff) were not
  run. `tsso check` accepts any elements, but the rules were only gated on RDB7: test them on reference TSs from your own
  chemistry with `tsso gate` before relying on them there.
- **Eight models cost eight trainings.** Comparing against only 3 other models gives similar disagreement quality
  (AUROC 0.84 vs 0.86 for geometry errors), but a single model with many samples was not tested because
  the released files keep only one sample per model.
- **Some choices were made after looking at test results.** The product-stereocentre rule was added after we saw
  stereo errors at product stereocentres in the test predictions, and "wrong" uses a 0.2 Å threshold we chose
  (about twice GoFlow's average error). Every rule's code and threshold were fixed using training references only.
- **Some RDB7 reference TSs have inconsistent chirality.** FK-Flow reports RDB7 entries whose TS chirality disagrees with
  the reactant and product, and some with impossible chirality. The gate tolerates up to 1% rejected reference TSs, which
  leaves room for such entries, but a correct rule could still look slightly worse because of them.
- **Rules written by an LLM can be wrong.** The gate tests them on real transition states, but it cannot prove a rule
  is chemically complete. Every check is plain Python in `tsso/checks/` for a chemist to read.
- One of the 6 reactions marked "re-generate" (RDB7 1093) looks like a dataset inconsistency rather than a model error:
  its atom-mapped SMILES says the N5–O6 bond breaks, but the reactant, TS and product geometries all keep it near 1.4 Å.

## Related work in detail

- **GoFlow**: Galustian, Mark, Karwounopoulos, Kovar & Heid, *Digital Discovery* 2025, [code](https://github.com/heid-lab/goflow) (MIT).
  All model predictions used here are the test-set outputs the authors released in that repository.
- **FK-Flow**: Mark, Galustian, Kovar & Heid, [arXiv:2509.01543](https://arxiv.org/abs/2509.01543), [code](https://github.com/heid-lab/fkflow).
  It showed that many of GoFlow's errors on chiral reactions come from predicting the wrong chirality, and fixed them by Feynman-Kac steering of the model.
  TS Second Opinion catches and repairs the same failure after the fact, without access to the model, and lets
  chemists add new rules in plain English. The same checks could serve as steering potentials (FK steering does not need
  differentiable potentials), which is the natural next step once model weights are available.
- **RDB7**: Spiekermann, Pattanaik & Green, *Sci. Data* 2022 ([Zenodo](https://zenodo.org/records/13328872)); reference geometries and barriers.
- **TS-DFM** (Luo, Gu & Sun, [*Nat. Commun.* 2026](https://www.nature.com/articles/s41467-026-74101-0)) notes that distance
  geometry "cannot distinguish enantiomers" and anchors its reconstruction to the reactant and product, but reports no
  stereo-correctness metric, which supports reporting one.
- **TSAgent** (Madhavan et al., [arXiv:2605.14154](https://arxiv.org/abs/2605.14154)) automates DFT TS searches with LLM agents;
  it does not check ML guesses before DFT or check stereochemistry. TS Second Opinion is designed to sit in front of such agents.
- LLM-written reward and validation code has precedent in robotics (e.g. Eureka, Ma et al. 2023), and an LLM agent that writes
  reaction-classification rules and verifies each against a patent corpus was reported in 2026
  ([arXiv:2607.01061](https://arxiv.org/abs/2607.01061)). We did not find prior work applying the idea to 3D transition-state
  geometries, but our search was brief.
- **Conformal selection**: Jin & Candès, "Selection by prediction with conformal p-values", *JMLR* 2023.
- **MCP for chemistry** is not new: for example, El Agente Quntur ([arXiv:2602.04850](https://arxiv.org/abs/2602.04850)) is a
  quantum-chemistry research agent whose tools include IRC set-up from transition states, and several 2025 hackathon projects
  used MCP. To our knowledge, none triages machine-learned TS guesses before DFT.
- The idea grew out of a Q&A with Dr Esther Heid (TU Wien) at the AIchemy webinar on 21 October 2026,
  about LLM-written steering potentials and about predicting which TS guesses will fail before DFT. Any errors here are ours.
- 3D viewer: [3Dmol.js](https://3dmol.csb.pitt.edu/) (BSD-3, vendored in `tsso/assets`).
