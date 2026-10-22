RULE = "Atoms that are directly bonded to each other and not part of the reaction should keep their bond angles within a sensible range (for example, sp3 centres near 109\u00b0, sp2 near 120\u00b0, sp near 180\u00b0) relative to the reactant and product values, with deviations allowed to grow only for atoms next to the forming or breaking bonds."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.8%; rejects 0% of mirror images and 64% of noisy geometries.

def check(rxn, xyz):
    """Bond angles at atoms outside the reaction centre stay near reactant/product values."""
    centre = set(rxn.reaction_centre)
    nr = rxn.neighbours('reactant')
    npd = rxn.neighbours('product')
    near = set()
    for a in centre:
        near.update(nr.get(a, []))
        near.update(npd.get(a, []))

    def ang(p, c, i, j):
        u = p[i] - p[c]
        v = p[j] - p[c]
        cs = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v) + 1e-12)
        return math.degrees(math.acos(max(-1.0, min(1.0, float(cs)))))

    worst = None
    checked = 0
    for c in range(rxn.n):
        if c in centre:
            continue
        shared = sorted(set(nr.get(c, [])) & set(npd.get(c, [])))
        shared = [k for k in shared if k not in centre]
        for a in range(len(shared)):
            for b in range(a + 1, len(shared)):
                i, j = shared[a], shared[b]
                ref = [ang(rxn.reactant_xyz, c, i, j), ang(rxn.product_xyz, c, i, j)]
                cur = ang(xyz, c, i, j)
                dev = min(abs(cur - r) for r in ref)
                adjacent = c in near or i in near or j in near
                tol = 40.0 if adjacent else 25.0
                checked += 1
                excess = dev - tol
                if excess > 0 and (worst is None or excess > worst[0]):
                    worst = (excess, c, i, j, cur, ref, dev, tol)
    if checked == 0:
        return True, "rule not applicable"
    if worst is not None:
        _, c, i, j, cur, ref, dev, tol = worst
        return False, (f"angle {rxn.symbols[i]}{i + 1}-{rxn.symbols[c]}{c + 1}-{rxn.symbols[j]}{j + 1} "
                       f"is {cur:.0f} deg vs reactant {ref[0]:.0f}/product {ref[1]:.0f} "
                       f"(deviation {dev:.0f} > {tol:.0f})")
    return True, f"{checked} spectator bond angles are within tolerance"
