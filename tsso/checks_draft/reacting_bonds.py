RULE = ("Every bond that forms or breaks has a transition-state length between its reactant and product "
        "lengths, allowing {slack:.0%} either side.")
PARAMS = {"slack": 0.20}

def check(rxn, xyz, slack=PARAMS["slack"]):
    bad = []
    for (i, j) in rxn.changed_bonds:
        r = np.linalg.norm(rxn.reactant_xyz[i] - rxn.reactant_xyz[j])
        p = np.linalg.norm(rxn.product_xyz[i] - rxn.product_xyz[j])
        t = np.linalg.norm(xyz[i] - xyz[j])
        lo, hi = min(r, p) * (1 - slack), max(r, p) * (1 + slack)
        if not (lo <= t <= hi):
            bad.append((abs(t - (lo if t < lo else hi)), i, j, t, lo, hi))
    if bad:
        _, i, j, t, lo, hi = max(bad)
        verb = "forming" if (i, j) in rxn.product_bonds else "breaking"
        return False, f"{verb} bond {rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} is {t:.2f} A, outside {lo:.2f}-{hi:.2f} A"
    return True, "forming and breaking bonds sit between reactant and product"
