RULE = "Bonds that survive the reaction keep a length close to their lengths in the reactant and product."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.6%; rejects 0% of mirror images and 99% of noisy geometries.

def check(rxn, xyz):
    """Bonds present in both reactant and product keep lengths close to their reactant/product values."""
    changed = set(tuple(sorted(b)) for b in rxn.changed_bonds)
    shared = [b for b in rxn.reactant_bonds if b in rxn.product_bonds and b not in changed]
    if not shared:
        return True, "rule not applicable"
    worst = None
    for (i, j) in shared:
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        dr = float(np.linalg.norm(rxn.reactant_xyz[i] - rxn.reactant_xyz[j]))
        dp = float(np.linalg.norm(rxn.product_xyz[i] - rxn.product_xyz[j]))
        lo, hi = min(dr, dp), max(dr, dp)
        # bonds at the reaction centre may change order and stretch more
        centre = i in rxn.reaction_centre or j in rxn.reaction_centre
        tol_lo = 0.30 if centre else 0.20
        tol_hi = 0.45 if centre else 0.30
        if d < lo - tol_lo:
            excess = lo - d
        elif d > hi + tol_hi:
            excess = d - hi
        else:
            continue
        if worst is None or excess > worst[0]:
            worst = (excess, i, j, d, lo, hi)
    if worst is not None:
        _, i, j, d, lo, hi = worst
        return False, (f"{rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} is {d:.2f} A "
                       f"(reactant/product range {lo:.2f}-{hi:.2f} A)")
    return True, "all surviving bonds are close to their reactant/product lengths"
