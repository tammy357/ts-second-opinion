# Batch example

Four RDB7 test reactions with eight GoFlow guesses each, laid out the way `tsso batch` reads a folder: one sub-folder per
reaction, holding `reactant.xyz`, `product.xyz` and the guesses (here one multi-frame `guesses.xyz`). Only `rxn_3509` has
a `reaction.smi`; for the other three, bonds and stereocentres are read from the geometries.

    tsso batch examples/batch --out tsso_batch

gives one of each verdict: `rxn_3509` trust (after mirror repair), `rxn_1292` check with DFT, and `rxn_1093` and
`rxn_5877` re-generate. The results are in `tsso_batch/verdicts.csv`, with the recommended geometry for every reaction
and an ORCA input for the ones to check or re-generate.

`reactions.csv` shows the other input format, a table of file paths (relative to the CSV):

    tsso batch examples/batch/reactions.csv --out tsso_batch
