RULE = "A hydrogen that is being transferred between two heavy atoms sits roughly on the line between them at the transition state."
# Written by an agent via the ts-second-opinion MCP server on 2026-10-21.
# Gate (RDB7 training references): passes 99.9%; rejects 0% of mirror images and 1% of noisy geometries.

def check(rxn, xyz):
    """A transferred H is roughly collinear with donor and acceptor, unless a small ring TS forces it to bend."""
    nb_r = rxn.neighbours('reactant')
    nb_p = rxn.neighbours('product')
    checked = 0
    for h in range(rxn.n):
        if rxn.symbols[h] != 'H':
            continue
        heavy_r = [a for a in nb_r.get(h, []) if rxn.symbols[a] != 'H']
        heavy_p = [a for a in nb_p.get(h, []) if rxn.symbols[a] != 'H']
        if len(heavy_r) != 1 or len(heavy_p) != 1 or heavy_r[0] == heavy_p[0]:
            continue
        d, a = heavy_r[0], heavy_p[0]
        # 1,2- and 1,3-shifts (donor and acceptor bonded or sharing a neighbour) go through 3-/4-membered rings and must bend
        close = (a in nb_r.get(d, [])) or (a in nb_p.get(d, [])) \
            or (set(nb_r.get(d, [])) & set(nb_r.get(a, []))) or (set(nb_p.get(d, [])) & set(nb_p.get(a, [])))
        if close:
            continue
        v1 = xyz[d] - xyz[h]
        v2 = xyz[a] - xyz[h]
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            continue
        checked += 1
        ang = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2) / (n1 * n2))))))
        if ang < 90.0:
            return False, (f"H{h + 1} moves from {rxn.symbols[d]}{d + 1} to {rxn.symbols[a]}{a + 1} "
                           f"but the angle {rxn.symbols[d]}{d + 1}-H{h + 1}-{rxn.symbols[a]}{a + 1} "
                           f"is only {ang:.0f} degrees at the TS (expected >= 90)")
    if checked == 0:
        return True, "rule not applicable (no transferred hydrogen between distant heavy atoms)"
    return True, f"{checked} transferred hydrogen(s) lie roughly between donor and acceptor"
