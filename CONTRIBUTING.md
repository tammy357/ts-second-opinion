# Contributing

Thank you for helping. The two most useful contributions are **new rules** and **new generators on the leaderboard**.

## Set up

```bash
git clone https://github.com/tammy357/ts-second-opinion && cd ts-second-opinion
pip install -e ".[dev]"
python -m pytest -q tests
```

## Add a rule (no coding needed)

1. Write the rule as one plain-English sentence a chemist would check, e.g.
   *"A hydrogen that is being transferred sits roughly on the line between donor and acceptor."*
2. Let Claude write and test it: `tsso author "<your rule>"`, or `/new-rule "<your rule>"` in Claude Code
   (this repository's MCP server and commands load automatically). Writing it yourself is fine too:
   `tsso gate "<your rule>" my_check.py`.
3. A check is kept only if it passes the gate: it accepts at least 99% of 9,540 real DFT transition states,
   never crashes, and rejects at least 1% of mirrored or noisy ones. Accepted checks are saved in `tsso/checks_user/`.
4. Open a pull request with the file. We move a rule into `tsso/checks/` (used for triage) only when it also helps
   on held-out reactions, and we say in the README when it does not.

## Add a TS generator to the leaderboard

Export your model's test-set predictions and the reference TSs as multi-frame `.xyz` files (one frame per reaction,
same order, atoms listed by atom-map number) plus one atom-mapped reaction SMILES per line, then run

```bash
tsso evaluate predictions.xyz references.xyz reactions.smi --json > my_model.json
```

and open a pull request (or an issue using the "Leaderboard entry" template) with the numbers and a link to the model.

## Ground rules

- Never use the reference TS to make a decision; it is for scoring only.
- Report what did not work as well as what did.
- Keep checks as plain Python that a chemist can read: no imports, no file access (`tsso/sandbox.py` enforces this).
