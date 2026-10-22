RULE = ("Atoms with three bonded neighbours that keep a double or aromatic bond through the reaction "
        "stay flat: the atom sits within {tol} Angstrom of the plane of its neighbours.")
PARAMS = {"tol": 0.35}

def check(rxn, xyz, tol=PARAMS["tol"]):
    nb = rxn.neighbours('reactant')
    worst = None
    for a in range(rxn.n):
        if a in rxn.reaction_centre or len(nb[a]) != 3:
            continue
        keeps_pi = any(rxn.reactant_bonds.get(tuple(sorted((a, b))), 1) > 1 and
                       rxn.product_bonds.get(tuple(sorted((a, b))), 1) > 1 for b in nb[a])
        if not keeps_pi:
            continue
        p = xyz[nb[a]]
        normal = np.cross(p[1] - p[0], p[2] - p[0])
        norm = np.linalg.norm(normal)
        if norm < 1e-6:
            continue
        h = abs(np.dot(xyz[a] - p[0], normal / norm))
        if h > tol and (worst is None or h > worst[0]):
            worst = (h, a)
    if worst:
        h, a = worst
        return False, f"{rxn.symbols[a]}{a + 1} is {h:.2f} A out of the plane of its neighbours"
    return True, "unreacting double bonds stay flat"
