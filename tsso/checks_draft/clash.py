RULE = ("Atoms that are not bonded in either the reactant or the product should not be closer than "
        "{heavy} Angstrom (two heavy atoms) or {light} Angstrom (if either is hydrogen).")
PARAMS = {"heavy": 1.6, "light": 1.1}

def check(rxn, xyz, heavy=PARAMS["heavy"], light=PARAMS["light"]):
    bonded = set(rxn.reactant_bonds) | set(rxn.product_bonds)
    worst = None
    for i in range(rxn.n):
        for j in range(i + 1, rxn.n):
            if (i, j) in bonded:
                continue
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            limit = light if 1 in (rxn.atoms[i], rxn.atoms[j]) else heavy
            if d < limit and (worst is None or d - limit < worst[0]):
                worst = (d - limit, i, j, d, limit)
    if worst:
        _, i, j, d, limit = worst
        return False, f"{rxn.symbols[i]}{i + 1} and {rxn.symbols[j]}{j + 1} are {d:.2f} A apart (limit {limit:.2f} A) but never bonded"
    return True, "no unbonded atoms are unrealistically close"
