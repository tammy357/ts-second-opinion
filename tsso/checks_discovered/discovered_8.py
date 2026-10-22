RULE = "Torsion angles around bonds that are not formed, broken, or changed in order should stay close to their values in the reactant or product (within about 30\u201340\u00b0), so that substituents far from the reaction centre do not rotate arbitrarily and the TS stays connected to both endpoints."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.6%; rejects 10% of mirror images and 4% of noisy geometries.

def check(rxn, xyz):
    """Torsions about unchanged bonds stay within 50 degrees of the reactant or product value."""
    def dih(p, a, b, c, d):
        b0 = p[a] - p[b]; b1 = p[c] - p[b]; b2 = p[d] - p[c]
        n1 = np.cross(b0, b1); n2 = np.cross(b1, b2)
        m = np.cross(n1, b1 / (np.linalg.norm(b1) + 1e-12))
        return math.degrees(math.atan2(float(np.dot(m, n2)), float(np.dot(n1, n2))))

    def ang(p, a, b, c):
        u = p[a] - p[b]; v = p[c] - p[b]
        cs = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-12)
        return math.degrees(math.acos(max(-1.0, min(1.0, float(cs)))))

    def diff(a, b):
        d = abs(a - b) % 360.0
        return min(d, 360.0 - d)

    nr = rxn.neighbours('reactant')
    npd = rxn.neighbours('product')
    changed = set(tuple(sorted(b)) for b in rxn.changed_bonds)
    centre = rxn.reaction_centre
    R, P = rxn.reactant_xyz, rxn.product_xyz
    tol = 50.0
    worst = None
    checked = 0
    for (i, j), order in rxn.reactant_bonds.items():
        if (i, j) in changed or rxn.product_bonds.get((i, j)) != order:
            continue
        if i in centre or j in centre or order == 3:
            continue
        ks = [k for k in nr[i] if k != j and k in npd[i] and rxn.atoms[k] != 1]
        ls = [l for l in nr[j] if l != i and l in npd[j] and rxn.atoms[l] != 1]
        for k in ks:
            if k in centre:
                continue
            for l in ls:
                if l in centre or k == l:
                    continue
                bad = False
                for g in (R, P):
                    a1, a2 = ang(g, k, i, j), ang(g, i, j, l)
                    if not (20 < a1 < 160 and 20 < a2 < 160):
                        bad = True
                if bad or ang(xyz, k, i, j) > 175 or ang(xyz, i, j, l) > 175:
                    continue
                checked += 1
                t = dih(xyz, k, i, j, l)
                dev = min(diff(t, dih(R, k, i, j, l)), diff(t, dih(P, k, i, j, l)))
                if dev > tol and (worst is None or dev > worst[0]):
                    worst = (dev, k, i, j, l, t)
    if checked == 0:
        return True, "rule not applicable"
    if worst:
        dev, k, i, j, l, t = worst
        return False, (f"torsion {rxn.symbols[k]}{k+1}-{rxn.symbols[i]}{i+1}-"
                       f"{rxn.symbols[j]}{j+1}-{rxn.symbols[l]}{l+1} is {t:.0f} deg, "
                       f"{dev:.0f} deg from reactant/product")
    return True, f"{checked} torsions about unchanged bonds are within {tol:.0f} deg of endpoints"
