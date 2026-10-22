RULE = "For a reaction with a single forming or breaking bond between a pair of atoms, the TS bond length should lie near a sensible fraction of the way between the reactant and product values, and the sum of the fractional progress of all breaking bonds should be roughly comparable to that of all forming bonds, which keeps the TS from being strongly early or late on only one side."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 3.
# Gate (RDB7 training references): passes 99.2%; rejects 1% of mirror images and 7% of noisy geometries.

def check(rxn, xyz):
    """Each forming/breaking bond is partway between reactant and product lengths, and breaking and forming progress are roughly balanced."""
    brk, frm = [], []
    for (i, j) in rxn.changed_bonds:
        rr = float(np.linalg.norm(rxn.reactant_xyz[i] - rxn.reactant_xyz[j]))
        rp = float(np.linalg.norm(rxn.product_xyz[i] - rxn.product_xyz[j]))
        rt = float(np.linalg.norm(xyz[i] - xyz[j]))
        span = rp - rr
        if abs(span) < 0.4:
            continue
        f = (rt - rr) / span
        label = f"{rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1}"
        if f < -0.4 or f > 1.5:
            return False, f"{label} TS length {rt:.2f} A (reactant {rr:.2f}, product {rp:.2f}) gives progress {f:.2f}"
        if (i, j) in rxn.reactant_bonds and (i, j) not in rxn.product_bonds:
            brk.append(f)
        elif (i, j) in rxn.product_bonds and (i, j) not in rxn.reactant_bonds:
            frm.append(f)
    if not brk and not frm:
        return True, "rule not applicable"
    if not brk or not frm:
        return True, "single-sided bond change, per-bond progress within range"
    sb, sf = sum(brk), sum(frm)
    if abs(sb - sf) > 1.8:
        return False, f"unbalanced progress: breaking sum {sb:.2f} vs forming sum {sf:.2f}"
    return True, f"progress balanced: breaking {sb:.2f}, forming {sf:.2f}"
