# TS Second Opinion

<!-- mcp-name: io.github.tammy357/ts-second-opinion -->

[![tests](https://github.com/tammy357/ts-second-opinion/actions/workflows/tests.yml/badge.svg)](https://github.com/tammy357/ts-second-opinion/actions/workflows/tests.yml)
![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)
![MCP server](https://img.shields.io/badge/MCP-server-8A2BE2)
![CPU only](https://img.shields.io/badge/runs%20on-a%20laptop%20CPU-lightgrey)

**Check, repair and triage machine-learned transition-state (TS) guesses before you spend DFT.**
Works with any TS generator and any elements, from a web app, the command line, Python, or an AI agent over MCP.
Installed on your computer, all it needs is three `.xyz` files: the reactant, the product and your TS guesses (the
in-browser app also asks for an atom-mapped SMILES).

> **Open Scientific Intelligence (OSI) Hackathon 2026** · Visionary Award · solo entry · sponsor tracks: marimo, Valency MCP Bond<br>
> **[▶ Live app](https://tammy357.github.io/ts-second-opinion/app/)** (runs in your browser, nothing to install) ·
> **[Video (4 min)](https://youtu.be/zaH6y7H471k)** · [Showcase page](https://tammy357.github.io/ts-second-opinion/) ·
> [Detailed results](https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md)

## In 60 seconds

- **The problem.** About **1 in 6** GoFlow transition-state guesses has a stereocentre the wrong way round. Most are exact
  mirror images, so they pass D-MAE, mirror-allowed RMSD and **even a DFT frequency check**. It is not only GoFlow: an
  independent 2026 study sees the same blind spot in MolGen and TSDiff.
- **What it does.** Chemists write rules in plain English; **Claude turns each into code, kept only if it passes 9,540
  real DFT transition states**. The rules check every guess, mirror-repair whole mirror images, and give a verdict
  (✓ trust, ! check with DFT, ✕ re-generate) with a **stated error rate** for what to trust (conformal selection).
- **What it found.** Catches **93%** of wrong-handed guesses (0.2% false alarms); wrong-handed answers **17.7% → 2.3%**;
  **80%** of reactions trusted at a 5% target error rate (an average guarantee, if new reactions resemble the
  calibration set).
- **New DFT confirms it.** ωB97X/def2-SVP transition-state optimisations: **6 of 6 trusted guesses** reached the
  reference TS, **0 of 6 re-generate guesses** did (though 5 of them pass a DFT frequency check), and an unrepaired
  mirror image is a true DFT saddle point **0.0006 kcal/mol** from the right one. GFN2-xTB on 106 more reactions agrees.
- **Also new.** A data audit (51 RDB7 entries contradict their own SMILES), a mirror test for any TS metric, a
  stereo-aware leaderboard, and panels where anyone can drop in their own reactions, test sets and datasets.
- **New in 0.3.** Check a reaction from **three `.xyz` files, no SMILES needed** (same verdict as the SMILES route on
  99.7% of the test reactions); triage **whole folders** into one verdict table; set up an AI agent with one command; and a
  "Was this verdict right?" link on every result. See [CHANGELOG.md](https://github.com/tammy357/ts-second-opinion/blob/main/CHANGELOG.md).
- **Limits.** The DFT sample is small (16 reactions, smaller basis than RDB7's); one dataset, RDB7 (also tested on its
  harder splits, which partly overlap the data the rules were gated on); geometry errors are flagged for DFT, not fixed.

![a near-perfect TS guess with the wrong handedness](https://raw.githubusercontent.com/tammy357/ts-second-opinion/main/results/fig_hero_mirror.png)

## Why it matters

Models such as [GoFlow](https://github.com/heid-lab/goflow) predict a TS in milliseconds instead of hours of DFT, but they
cannot tell you when they are wrong. On the RDB7 benchmark, **about one in six GoFlow guesses has a stereocentre the wrong
way round**. Most are exact mirror images of the right TS, so they pass distance-based scores, chirality-blind RMSD and
even a DFT frequency check: they are transition states of the mirror-image reaction.
[FK-Flow](https://arxiv.org/abs/2509.01543) (from the GoFlow group) first measured these errors and fixed them by steering
GoFlow during generation. TS Second Opinion works afterwards instead, on any generator's output, and also tells you which
guesses to trust, at an error rate that holds on average when new reactions resemble the calibration set.
It is not only GoFlow: on another dataset (Transition1x), an independent 2026 study
([ReCurveflow](https://arxiv.org/abs/2608.20869)) found that MolGen, GoFlow and TSDiff all lose 0.15–0.35 Å of
RMSD once mirror images are no longer allowed.

| On 1,194 held-out RDB7 reactions (8 GoFlow models) | |
|---|---|
| Wrong-handed guesses caught by two LLM-written rules | **93%** (0.2% false alarms) |
| Final answers with the wrong handedness, before → after checking and mirror repair | **17.7% → 2.3%** |
| Geometry errors among guesses marked "trust" vs those sent to DFT | **4.4% vs 40%** |
| Reactions you can trust at a 5% target average geometry-error rate (conformal selection; assumes new reactions resemble the calibration set) | **80%** |
| "Success" (D-MAE ≤ 0.2 Å) that is also correctly handed: GoFlow → with TS Second Opinion | **71% → 86.5%** |
| **New DFT** (ωB97X/def2-SVP TS optimisation): guesses that reach the reference TS, ✓ trust / ! check / ✕ re-generate | **6 of 6 / 1 of 3 / 0 of 6** |
| **New GFN2-xTB** TS optimisation, same measure, ✓ trust / ! check / ✕ re-generate (106 reactions) | **77% / 60% / 0%** |

The rules repair handedness; they do not fix geometry errors (14.1% → 11.7%), which are flagged for DFT instead (a quick
xTB optimisation rescued only 1 of 13 flagged geometries). On RDB7's harder splits (reaction types GoFlow never saw in
training, and the barrier split), the rules still catch 90% and 89% of wrong-handed guesses at 0.1% false alarms. These
two numbers are not a fully independent test: the training references the rules were gated on include about 80% of
those splits' test reactions (details in [docs/RESULTS.md](https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md#does-it-generalise)).

**New DFT calculations** test the central claims directly. Optimised with ωB97X/def2-SVP (the functional family of RDB7's references)
and checked with analytic frequencies, the unrepaired mirror image of reaction 3509 is a first-order saddle point
0.0006 kcal/mol and 2 cm⁻¹ from the repaired one, so a DFT frequency check passes it; only the stereo check sees it.
All 6 trusted guesses reached the reference TS and 1 of 3 "check with DFT" guesses did; none of the 6 re-generate
guesses did, yet 5 of them still passed the frequency check, all with the wrong handedness. GFN2-xTB on 106 more reactions agrees (mirror pairs identical in 34 of 40).

![new DFT calculations](https://raw.githubusercontent.com/tammy357/ts-second-opinion/main/results/fig_dft.png)

Full results, caveats and every number's source: **[docs/RESULTS.md](https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md)**.

### How it differs from what exists

| | Sees a wrong-handed TS? | Works on any generator's existing output? | Repairs it? | Says which guesses to trust, at a stated error rate? | Audits datasets and metrics? |
|---|---|---|---|---|---|
| Usual scores (D-MAE, RMSD with the mirror allowed) | ✕ | ✓ | ✕ | ✕ | ✕ |
| DFT frequency check | ✕ (a mirror image passes) | ✓ | ✕ | ✕ | ✕ |
| [FK-Flow](https://arxiv.org/abs/2509.01543): steers GoFlow while it generates | ✓ | ✕ (inside one model) | prevents | — | — |
| [RitS](https://doi.org/10.26434/chemrxiv.15001681/v1): stereo-aware generator | ✓ | ✕ (its own model) | prevents | — | — |
| **TS Second Opinion** | ✓ | ✓ | ✓ (whole mirror images) | ✓ (conformal selection) | ✓ (`tsso audit`, mirror test) |

The approaches are complementary: a stereo-aware generator makes fewer errors, and TS Second Opinion catches what any
generator still gets wrong, before DFT is spent.

### Who can use it

- **Anyone running ML transition-state models** (catalysis, combustion and atmospheric chemistry, process chemistry,
  reaction-network exploration): check, repair and decide what to send to DFT.
- **Anyone building or benchmarking a 3D generative model of chiral molecules**: the mirror test tells you in one click
  whether your metric can see handedness, and `tsso evaluate` scores a model with and without it.
- **Anyone building a reaction dataset**: drop it into the app's "Audit your own reaction dataset" panel (or
  `tsso.service.audit_files`) to find entries whose geometry contradicts their own SMILES.
- **Anyone validating scientific ML outputs with LLMs**: rules written in plain English, turned into code by an LLM, and
  kept only if they pass real reference data, is a pattern that carries over to other structure predictions.

## Quickstart

**Use it now, nothing to install:** **[open the live app](https://tammy357.github.io/ts-second-opinion/app/)**, or
[![Open in molab](https://marimo.io/molab-shield.svg)](https://molab.marimo.io/github/tammy357/ts-second-opinion/blob/main/apps/demo.py/wasm).
Both run the marimo app in your browser with WebAssembly (Pyodide), so uploaded files never leave your computer. In the
browser, a small pure-Python SMILES reader stands in for RDKit; it gives identical results on every RDB7 reaction
(`tests/test_smiles_lite.py`), but there the atom-mapped SMILES is still needed. Installed on your computer, three
`.xyz` files are enough. To write new rules with Claude in the browser, paste your own Anthropic API key; it is sent
only to api.anthropic.com.

**Install** the library, the `tsso` command and the MCP server (Python 3.10+; no GPU and no PyTorch needed):

```bash
pip install git+https://github.com/tammy357/ts-second-opinion
```

Or work from a clone, which also has the app, the scripts and every result file:

```bash
git clone https://github.com/tammy357/ts-second-opinion && cd ts-second-opinion
pip install -e .
marimo run apps/demo.py             # opens the app in your browser
```

![the app: examples, verdict card and linked 3D views](https://raw.githubusercontent.com/tammy357/ts-second-opinion/main/results/demo.gif)

**What you can do in the live app**

| Section | Try this |
|---|---|
| 1 · Check a TS guess | Five worked examples, or upload your own reaction: verdict, reasons, linked 3D views, report card, ORCA input, a data check, and Claude's explanation |
| 2 · Rule Studio | Replay a recorded run where Claude's first rule was rejected and rewritten, or write a new rule with Claude (your API key) |
| 3 · How much can you trust? | Pick an error rate; conformal selection says what share of reactions you can trust |
| 4 · Stereo-aware leaderboard | Every released GoFlow model, plus **score your own model** (or try 200 GoFlow predictions) |
| 5 · Beyond one model and one dataset | How it differs from prior work, the field-wide evidence, **the mirror test** for any metric, **audit your own dataset**, the new DFT and xTB results, an experimental repair for partly flipped guesses |
| 6 · Results | The figures |

### Check one reaction from three `.xyz` files

```bash
cd examples/rdb7_3509
tsso check --reactant reactant.xyz --product product.xyz guesses.xyz --report report.html
```

```
Your reaction: ✓ TRUST  (trust: passes every rule and the models agree)
  SMILES read from the geometries: [C:1]([C@@:2]1([H:10])[C:3]([H:11])([H:12])[C@:4]1([C:5]#[C:6][H:14])[H:13])([H:7])([H:8])[H:9]>>...
  note: Bond orders could not be assigned for the product (no closed-shell bonding pattern fits, as for some carbenes
        and zwitterions), so the planarity rule has nothing to check there. Give an atom-mapped SMILES if you need it.
  recommended guess: model 8  (mirror-repaired)
  - 0 of 8 model guesses pass every rule.
  - Rejected guesses broke: chirality (8).
  - The best guess had a stereocentre the wrong way round; its mirror image fixes this without changing any
    interatomic distance. The unrepaired guess has the same energy and frequencies, so a DFT frequency check would
    not have caught it.
  - Models disagree most about H8-H14 (not a forming/breaking bond): spread 0.09 A.
  - Mean disagreement of the chosen guess with the other models: 0.024 A (review cut-off 0.109 A).
  Was this verdict right? Tell us (opens a GitHub form, shares no structures): https://github.com/tammy357/...
  report card: report.html
```

Give the reactant, the product and one or more TS guesses from any generator (several files, or one multi-frame
`.xyz`); several guesses (models, seeds or samples) add the model-disagreement signal. Every file must list the same
atoms in the same order. Without a SMILES, bonds, bond orders and stereocentres are read from the reactant and product
geometries (RDKit), and the inferred atom-mapped SMILES is printed so you can check it. Add `--charge` and
`--multiplicity` for charged or open-shell systems (defaults 0 and 1), and `--write-xyz best.xyz --orca best.inp` for
the recommended geometry and an ORCA OptTS + Freq input.

**How well does reading the geometry work?** On all 1,194 RDB7 test reactions it gives the same verdict as the
dataset's own SMILES for **99.7%** and the same recommended guess for **98.7%**. Final answers are wrong-handed **2.6%**
of the time from geometry alone, against 2.3% with the SMILES, and **6.3%** of trusted answers are wrong, against 6.2%
(checked against the reference TSs, for evaluation only). Every difference traces to a bond or a stereocentre that the
geometry reading sees differently from the SMILES. Bond orders could not be assigned for 12.6% of reactions (some
carbenes and zwitterions); there the planarity rule has nothing to check. Reproduce with
`python scripts/validate_geometry_input.py` (results in
[`results/geometry_input_validation.json`](https://github.com/tammy357/ts-second-opinion/blob/main/results/geometry_input_validation.json)).

**With an atom-mapped SMILES** (every atom mapped, hydrogens included; charge and multiplicity are read from it), add
`--smiles "<mapped reaction>"` or `--smiles reaction.smi`. Atoms are then listed by atom-map number by default; add
`--order smiles` if your files follow the SMILES order. See `examples/` for ready-made files.

### Check many reactions at once

```bash
tsso batch my_reactions/ --out tsso_batch
```

`my_reactions/` holds one sub-folder per reaction with `reactant.xyz`, `product.xyz`, the TS guess files and, if you
have it, `reaction.smi`. A CSV with columns `id, reactant, product, guesses` (and optional `smiles, charge,
multiplicity`) works too. You get `tsso_batch/verdicts.csv` with one row per reaction, the recommended geometry for each,
and ORCA inputs for the reactions to check with DFT or re-generate (`--reports` adds a report card per reaction). The
review cut-off for "check with DFT" was calibrated on GoFlow; for another generator, sort by
`batch_disagreement_rank` (the share of your batch with lower model disagreement) and review from the top.
[`examples/batch`](https://github.com/tammy357/ts-second-opinion/tree/main/examples/batch) has four reactions laid out this way (one of each verdict) and a matching CSV:
try `tsso batch examples/batch`.

### From Python

```python
from tsso import second_opinion
result = second_opinion(None, "reactant.xyz", "product.xyz", ["model_a.xyz", "model_b.xyz"], report="report.html")
print(result["verdict"], result["reasons"], result["smiles"])     # pass a mapped SMILES instead of None to use it
```

### From an AI agent

- **Claude Desktop or Cursor:** `tsso mcp-config` (or `--client cursor`) prints the settings that start the MCP server
  with [uv](https://docs.astral.sh/uv/), so nothing is installed by hand.
- **Claude Code:** `claude mcp add ts-second-opinion -- uvx --from git+https://github.com/tammy357/ts-second-opinion tsso mcp`,
  plus `tsso skill --install` for an agent skill that explains when and how to use the tool and how to report results
  without overclaiming. In a clone, run `claude` here, approve the server, and try `/second-opinion 3509` or
  `/new-rule "A transferred hydrogen sits roughly on the line between donor and acceptor."`.
- **Tools** include `triage_structures` (one reaction; the SMILES is optional), `triage_batch` (a folder or CSV),
  `evaluate_predictions`, `certify_trust`, `gate_check` and `save_rule`.

### Was this verdict right?

Every result, on the command line, in the report card, in the app and over MCP, carries a link that opens a pre-filled
GitHub issue form with the verdict, the rules broken and the version. It never includes your structures, and nothing is
sent until you submit the form. It is the quickest way to report a verdict that looks wrong (or right).

## More tools for your own data

To score a generator's test-set predictions with and without handedness (also in the live app, under the leaderboard):

```bash
tsso evaluate predictions.xyz references.xyz reactions.smi      # multi-frame xyz, one frame per reaction
```

Handedness stays at 82–85% correct for every released GoFlow model size (64 to 512 hidden features), so a bigger model
does not fix it. The full leaderboard is in [docs/RESULTS.md](https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md#full-stereo-aware-leaderboard); to add a
model, see [CONTRIBUTING.md](https://github.com/tammy357/ts-second-opinion/blob/main/CONTRIBUTING.md).

**Mirror test for metrics** (`tsso.mirror_test`, and in the app): scores each of 447 reference TSs against its own
mirror image. D-MAE and RMSD with the mirror allowed call all 447 pairs identical; RMSD with proper rotations sees
every one. Paste your own metric to test it.

**Data audit** (`tsso audit`): in 51 of the 10,734 RDB7 reactions used here (0.5%), the product geometry's bonding
contradicts the product SMILES, for example a hydrogen sitting on a different carbon. The list is in
[`results/rdb7_audit.csv`](https://github.com/tammy357/ts-second-opinion/blob/main/results/rdb7_audit.csv); every uploaded reaction gets the same check.

**Benchmark data** ([`benchmark/`](https://github.com/tammy357/ts-second-opinion/blob/main/benchmark/README.md)): handedness labels for all 9,552 GoFlow guesses on the 1,194 RDB7
test reactions, plus a mirror-image unit test for TS metrics (447 mirrored reference TSs that a handedness-blind metric
scores as perfect).

## How it works

```
plain-English rule ──► Claude writes check(rxn, xyz) ──► gate: accept ≥ 99% of 9,540 real DFT TSs, never crash,
        ▲                                                   and catch ≥ 1% of mirrored or distorted TSs
        └──────── feedback with the failing examples ◄──── (else rewrite)

TS guesses ──► run every accepted check ──► consensus among passing guesses
           ──► if a stereo rule still fails, try the mirror image (identical distances, opposite handedness)
           ──► verdict ✓ / ! / ✕ (! = most-disputed 20% by model disagreement) ──► report card + ORCA input
           ──► optional: certified trust, calibrated on reactions with known answers
```

### How LLMs are used

- **Rules → tested code.** A chemist writes a rule as one sentence; Claude writes the Python; the gate tests it on 9,540
  real DFT transition states and sends the failures back until it passes. Checks run in a sandbox that hides the
  reference TS. The app's Rule Studio replays real runs where the first attempt was rejected.
- **Agents.** An MCP server and two Claude Code commands let an agent triage reactions, certify, evaluate generators and
  write new rules (`/new-rule`, `/second-opinion`). With **Valency MCP Bond** connected, `/second-opinion` adds related
  reading from ChemRxiv, arXiv and PubMed to each report; see the recorded run in
  [results/second_opinion_report.md](https://github.com/tammy357/ts-second-opinion/blob/main/results/second_opinion_report.md). Valency also found the closest prior work:
  [RitS](https://doi.org/10.26434/chemrxiv.15001681/v1) (2026) makes TS *generation* stereo-aware; TS Second Opinion
  checks and repairs any generator's output afterwards.
- **Explanations.** Claude explains each verdict in plain English, from the facts the tool computed.
- **Honest comparison.** Claude's checks perform about the same as the same rules written directly as code. The value of
  the LLM is access (no coding needed); what makes any check trustworthy is the gate. Claude also proposed 8 rules of its
  own: all passed the gate, none helped on held-out data, so none were kept.

## OSI judging criteria at a glance

| Criterion | |
|---|---|
| **Potential for impact** | Clear users, usable today in a browser with their own data: chemists deciding which ML transition states to send to DFT (check, repair, trust at a stated error rate); TS-model developers, who get a stereo-aware score in one click; and dataset curators (audit your own dataset). New DFT confirms the verdicts: 6 of 6 trusted guesses reach the reference TS, 0 of 6 re-generate ones (GFN2-xTB on 106 reactions: 77% vs 0%). |
| **Innovativeness** | The obvious fixes are retraining or steering one model, or a DFT frequency check, which a mirror image passes. Instead, a model-agnostic second-opinion layer: chemists' rules in plain English become Claude-written code that must pass a data gate, a reflection repairs whole-mirror errors at no cost, and conformal selection states the error rate. |
| **Novelty** | Post-hoc handedness checking and repair for any generator (prior work, FK-Flow and RitS, acts inside one model); a stereo-aware leaderboard and a mirror test for any TS metric; an audit that found 51 RDB7 reactions whose product geometry contradicts its SMILES. Open avenues: other generators and datasets, larger DFT validation, and the gated rules as steering potentials during generation. |
| **Relevance** | A specific problem in chemistry: transition states gate kinetics, catalysis and reaction-network studies, and every wrong guess sent to DFT costs hours. |

## Roadmap: help wanted

These are the next steps, roughly in order. Issues and pull requests are very welcome; see
[CONTRIBUTING.md](https://github.com/tammy357/ts-second-opinion/blob/main/CONTRIBUTING.md).

- **A larger, random DFT sample** at RDB7's own level of theory (ωB97X-D3/def2-TZVP), to turn "6 of 6 trusted guesses
  reach the reference TS" into an estimate with an uncertainty.
- **More generators and datasets.** Transition1x, and generators such as React-OT, TSDiff and MolGen. If you have
  test-set predictions, `tsso evaluate` scores them, and a [leaderboard entry](https://github.com/tammy357/ts-second-opinion/blob/main/.github/ISSUE_TEMPLATE/leaderboard_entry.md)
  adds them here.
- **Reference TSs from other chemistry** (surfaces, organometallics, larger molecules), to re-gate the rules with
  `tsso gate` beyond small gas-phase H/C/N/O molecules.
- **A refinement stage** with a machine-learned interatomic potential (for example MACE-OMol25 or UMA) for the
  "check with DFT" cases, so DFT is spent only where it is needed.
- **Steering with the gated checks** during generation (Feynman-Kac steering accepts non-differentiable potentials).
- **Geometry-only input in the browser app**, which needs a pure-Python stereocentre reader in place of RDKit.

## Reproduce every number

```bash
pip install -e ".[dev]"
python -m pytest -q tests                # all tests
python scripts/evaluate.py               # -> results/metrics.json
python scripts/triage_all.py 0.2         # verdicts for all 1,194 reactions -> results/triage.csv
python scripts/benchmark.py              # mirror blindness, handedness split, risk-coverage -> results/benchmark.json
python scripts/certify.py                # certified trust -> results/certify.json
python scripts/leaderboard.py            # stereo-aware leaderboard -> results/leaderboard.md
python scripts/ablation_authors.py       # Claude-written vs directly written checks -> results/ablation_authors.json
python scripts/build_benchmark.py        # benchmark/ data files
python scripts/audit_rdb7.py             # RDB7 data audit -> results/rdb7_audit.csv
python scripts/validate_geometry_input.py   # geometry-only input vs the SMILES route -> results/geometry_input_validation.json
python scripts/xtb_prospective.py mirror 40; python scripts/xtb_prospective.py verdicts 30   # new xTB TS optimisations
pip install -e ".[dft]"; python scripts/dft_prospective.py mirror:3509 guess:97 guess:5877   # new DFT (about 15 min per job)
python scripts/dft_summary.py; python scripts/make_dft_figure.py   # -> results/dft_summary.json, results/fig_dft.png
python scripts/build_site.py             # the static showcase in docs/ (GitHub Pages)
python scripts/build_wasm.py             # the in-browser app in docs/app/
```

Writing new rules needs the `claude` CLI or an `ANTHROPIC_API_KEY` (`pip install -e ".[api]"`). Data files are derived
from GoFlow's released predictions and RDB7; `scripts/export_data.py` and `scripts/build_reference_set.py` rebuild them
from a clone of the GoFlow repository.

## Limitations

- **Mostly retrospective, with a small DFT check.** The main numbers compare guesses with RDB7's reference TSs. New DFT
  optimisations (16 reactions: the first by id for ✓ and !, all six ✕, not a random sample; ωB97X/def2-SVP, no
  dispersion, smaller basis than RDB7's ωB97X-D3/def2-TZVP) and GFN2-xTB (106 reactions) confirm the mirror-image
  argument and the verdicts; a larger, random DFT sample is the next step.
- **One dataset, one model family.** Small H/C/N/O molecules (RDB7) and GoFlow models, also tested on RDB7's harder
  splits (which partly overlap the gate's training data). `tsso check` accepts any elements, but the rules were gated
  on RDB7: test them on your own reference TSs with `tsso gate` first.
- **Reading bonds from geometry is a heuristic.** Distance cut-offs and RDKit's bond-order and stereocentre
  assignment were checked on RDB7 only. For unusual bonding (metals, hypervalent atoms, open-shell systems), check the
  inferred SMILES or give your own.
- **Mirror-image errors leave the barrier unchanged.** They matter whenever stereochemistry does; for barrier heights
  alone, use the geometry-only numbers.
- **Some choices were made after seeing test results**, and the guarantee is an average that assumes new reactions look
  like the calibration ones. Details: [docs/RESULTS.md](https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md#limitations-in-detail).

## Credits and citation

GoFlow (Galustian, Mark, Karwounopoulos, Kovar & Heid, *Digital Discovery* 2025; MIT) supplied every model prediction
used here. RDB7 (Spiekermann, Pattanaik & Green, *Sci. Data* 2022) supplied the reference structures. FK-Flow (Mark,
Galustian, Kovar & Heid, 2025) first measured GoFlow's chirality errors. The idea grew out of a Q&A with Dr Esther Heid
(TU Wien) at the AIchemy webinar on 21 October 2026. Conformal selection follows Jin & Candès (*JMLR* 2023). The 3D viewer
is [3Dmol.js](https://3dmol.csb.pitt.edu/) (BSD-3). More related work: [docs/RESULTS.md](https://github.com/tammy357/ts-second-opinion/blob/main/docs/RESULTS.md#related-work-in-detail).

If you use TS Second Opinion, please cite it ([CITATION.cff](https://github.com/tammy357/ts-second-opinion/blob/main/CITATION.cff); preprint:
[doi:10.26434/chemrxiv.15011459](https://doi.org/10.26434/chemrxiv.15011459/v1)) together with GoFlow and RDB7.
Contributions (new rules, new generators on the leaderboard) are welcome: [CONTRIBUTING.md](https://github.com/tammy357/ts-second-opinion/blob/main/CONTRIBUTING.md).

```
apps/demo.py         the app (marimo)            tsso/          library, CLI (tsso) and MCP server
examples/            ready-to-run inputs         tsso/checks/   the six gated, Claude-written rules
tsso/infer.py        bonds from geometry         tsso/batch.py  many reactions at once
tsso/skill/          the agent skill             server.json    MCP Registry entry
scripts/             experiments and figures     results/       metrics, figures, leaderboard, agent transcripts
docs/RESULTS.md      detailed results            data/          GoFlow predictions + RDB7 references (gzipped JSON)
docs/index.html      the showcase page           benchmark/     stereo-aware labels and the mirror-image unit test
docs/app/            the in-browser app (WASM)   results/dft/, results/xtb_prospective_*.json   the new DFT and xTB calculations
```

Code: MIT licence. Data files are derived from GoFlow (MIT) and RDB7; please cite both.
