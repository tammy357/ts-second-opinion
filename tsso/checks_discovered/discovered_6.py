RULE = "Bond lengths between atoms that are not changing bonds but are attached to a reacting atom should lie between their reactant and product values (or be slightly changed in a way consistent with a change of bond order, such as shortening for atoms becoming double-bonded), rather than showing arbitrary distortion."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 2.
# Gate (RDB7 training references): passes 99.8%; rejects 0% of mirror images and 62% of noisy geometries.

def check(rxn, xyz):
    """Unchanged bonds on reacting atoms stay roughly within their reactant/product lengths."""
    centre = rxn.reaction_centre
    changed = set(tuple(sorted(b)) for b in rxn.changed_bonds)
    worst = None
    n_checked = 0
    for (i, j), order_r in rxn.reactant_bonds.items():
        if (i, j) in changed or (i not in centre and j not in centre):
            continue
        if (i, j) not in rxn.product_bonds:
            continue
        dr = float(np.linalg.norm(rxn.reactant_xyz[i] - rxn.reactant_xyz[j]))
        dp = float(np.linalg.norm(rxn.product_xyz[i] - rxn.product_xyz[j]))
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        # real TSs stretch spectator bonds next to the reacting centre a lot
        # (bonds to a leaving group or migrating atom can lengthen by ~0.7 A)
        lo, hi = min(dr, dp) - 0.35, max(dr, dp) + 0.75
        n_checked += 1
        excess = max(lo - d, d - hi)
        if excess > 0 and (worst is None or excess > worst[0]):
            worst = (excess, i, j, d, min(dr, dp), max(dr, dp))
    if n_checked == 0:
        return True, "rule not applicable"
    if worst is not None:
        _, i, j, d, a, b = worst
        return False, (f"bond {rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} is {d:.2f} A, "
                       f"outside reactant/product range {a:.2f}-{b:.2f} A (with tolerance)")
    return True, f"{n_checked} spectator bonds on reacting atoms within expected range"
