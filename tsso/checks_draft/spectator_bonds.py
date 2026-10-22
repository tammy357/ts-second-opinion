RULE = ("Bonds that exist in both the reactant and the product stay within {slack:.0%} of the range "
        "spanned by their reactant and product lengths.")
PARAMS = {"slack": 0.30}

def check(rxn, xyz, slack=PARAMS["slack"]):
    bad = []
    for (i, j) in set(rxn.reactant_bonds) & set(rxn.product_bonds):
        r = np.linalg.norm(rxn.reactant_xyz[i] - rxn.reactant_xyz[j])
        p = np.linalg.norm(rxn.product_xyz[i] - rxn.product_xyz[j])
        t = np.linalg.norm(xyz[i] - xyz[j])
        lo, hi = min(r, p) * (1 - slack), max(r, p) * (1 + slack)
        if not (lo <= t <= hi):
            bad.append((abs(t - (lo if t < lo else hi)), i, j, t, lo, hi))
    if bad:
        _, i, j, t, lo, hi = max(bad)
        return False, f"bond {rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} is {t:.2f} A, outside {lo:.2f}-{hi:.2f} A"
    return True, "bonds that survive the reaction keep sensible lengths"
