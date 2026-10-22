RULE = "Bonds that form or break have lengths in between their reactant and product values."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.2%; rejects 1% of mirror images and 12% of noisy geometries.

def check(rxn, xyz):
    """Each forming/breaking bond has a TS length between its reactant and product lengths (with a tolerance)."""
    if not rxn.changed_bonds:
        return True, "rule not applicable"
    tol = 0.35
    for i, j in rxn.changed_bonds:
        dr = float(np.linalg.norm(rxn.reactant_xyz[i] - rxn.reactant_xyz[j]))
        dp = float(np.linalg.norm(rxn.product_xyz[i] - rxn.product_xyz[j]))
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        lo, hi = min(dr, dp), max(dr, dp)
        if d < lo - tol or d > hi + tol:
            return False, (f"{rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} is {d:.2f} A in the TS, "
                           f"outside reactant {dr:.2f} / product {dp:.2f} A (tolerance {tol} A)")
    return True, "all changing bonds lie between reactant and product lengths"
