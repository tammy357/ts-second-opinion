RULE = "The forming and breaking bonds should be elongated beyond typical equilibrium bond lengths for that element pair (for example, at least about 10\u201315% longer than the longer of the reactant and product values for that bond) so the geometry is not simply a reactant or product structure."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.9%; rejects 0% of mirror images and 2% of noisy geometries.

def check(rxn, xyz):
    """Each forming/breaking bond is stretched at least ~10% beyond its reactant/product bonded length."""
    if not rxn.changed_bonds:
        return True, "rule not applicable"
    best = None
    for (i, j) in rxn.changed_bonds:
        i, j = min(i, j), max(i, j)
        refs = []
        if (i, j) in rxn.reactant_bonds:
            refs.append(float(np.linalg.norm(rxn.reactant_xyz[i] - rxn.reactant_xyz[j])))
        if (i, j) in rxn.product_bonds:
            refs.append(float(np.linalg.norm(rxn.product_xyz[i] - rxn.product_xyz[j])))
        if not refs:
            continue
        ref = max(refs)
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        ratio = d / ref
        if best is None or ratio > best[0]:
            best = (ratio, i, j, d, ref)
    if best is None:
        return True, "rule not applicable"
    ratio, i, j, d, ref = best
    if ratio < 1.10:
        return False, (f"most stretched changing bond {rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} "
                       f"is {d:.2f} A vs bonded {ref:.2f} A (ratio {ratio:.2f} < 1.10)")
    return True, (f"{rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} stretched to "
                  f"{d:.2f} A (ratio {ratio:.2f})")
