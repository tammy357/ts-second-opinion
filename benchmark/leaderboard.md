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
