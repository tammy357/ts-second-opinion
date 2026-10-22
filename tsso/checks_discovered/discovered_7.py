RULE = "A reacting atom whose bond order or hybridization changes (such as a carbon going from sp3 to sp2) should show a degree of pyramidalization at the TS intermediate between its reactant and product values, measured by an improper torsion or the sum of bond angles around it."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.1%; rejects 1% of mirror images and 21% of noisy geometries.

def check(rxn, xyz):
    """Pyramidalization (sum of bond angles) of atoms changing hybridization lies between reactant and product values."""
    nr = rxn.neighbours('reactant')
    npd = rxn.neighbours('product')

    def order_sum(bonds, a):
        return sum(o for (i, j), o in bonds.items() if a in (i, j))

    def angle_sum(x, a, nb):
        tot = 0.0
        for k in range(3):
            for l in range(k + 1, 3):
                u = x[nb[k]] - x[a]
                v = x[nb[l]] - x[a]
                c = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
                tot += math.degrees(math.acos(max(-1.0, min(1.0, c))))
        return tot

    tol = 20.0
    checked = 0
    for a in sorted(rxn.reaction_centre):
        if rxn.atoms[a] == 1:
            continue
        r_nb = nr.get(a, [])
        p_nb = npd.get(a, [])
        d_ord = abs(order_sum(rxn.reactant_bonds, a) - order_sum(rxn.product_bonds, a))
        if len(r_nb) == len(p_nb) and d_ord < 0.4:
            continue
        common = [j for j in r_nb if j in p_nb]
        if len(common) != 3:
            continue
        sr = angle_sum(rxn.reactant_xyz, a, common)
        sp = angle_sum(rxn.product_xyz, a, common)
        if abs(sr - sp) < 10.0:
            continue
        st = angle_sum(xyz, a, common)
        checked += 1
        lo, hi = min(sr, sp) - tol, max(sr, sp) + tol
        if st < lo or st > hi:
            return False, (f"{rxn.symbols[a]}{a + 1} angle sum {st:.1f} deg at TS outside "
                           f"reactant {sr:.1f} / product {sp:.1f} range")
    if checked == 0:
        return True, "rule not applicable"
    return True, "pyramidalization of changing atoms is intermediate"
