RULE = "A transferring atom (such as H moving between two heavy atoms) should have its breaking and forming bond lengths both stretched, and the donor\u2013H\u2013acceptor angle should be fairly close to linear (roughly 150\u2013180\u00b0 for intermolecular-like transfers, and not below about 100\u00b0 in small rings)."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 3.
# Gate (RDB7 training references): passes 99.9%; rejects 0% of mirror images and 13% of noisy geometries.

def check(rxn, xyz):
    """Transferring H atoms have both bonds not compressed, stretched overall, and a donor-H-acceptor angle that fits the ring size."""
    rb = rxn.reactant_bonds
    rad = {1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66}
    nb = rxn.neighbours('reactant')
    checked = 0
    for h in range(rxn.n):
        if rxn.atoms[h] != 1:
            continue
        donors = [b[0] if b[1] == h else b[1] for b in rxn.changed_bonds if h in b and tuple(b) in rb]
        accs = [b[0] if b[1] == h else b[1] for b in rxn.changed_bonds if h in b and tuple(b) not in rb]
        if not donors or not accs:
            continue
        d, a = donors[0], accs[0]
        if d == a:
            continue
        checked += 1
        d1 = float(np.linalg.norm(xyz[h] - xyz[d]))
        d2 = float(np.linalg.norm(xyz[h] - xyz[a]))
        n1 = 0.31 + rad.get(int(rxn.atoms[d]), 0.7)
        n2 = 0.31 + rad.get(int(rxn.atoms[a]), 0.7)
        if d1 + d2 < 1.0 * (n1 + n2):
            return False, f"H{h + 1} bonds to {rxn.symbols[d]}{d + 1} ({d1:.2f} A) and {rxn.symbols[a]}{a + 1} ({d2:.2f} A) are not stretched"
        if d1 < 0.85 * n1 or d2 < 0.85 * n2:
            return False, f"H{h + 1} bond lengths {d1:.2f}/{d2:.2f} A to {rxn.symbols[d]}{d + 1}/{rxn.symbols[a]}{a + 1} are shorter than normal bonds"
        v1, v2 = xyz[d] - xyz[h], xyz[a] - xyz[h]
        c = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12))
        ang = math.degrees(math.acos(max(-1.0, min(1.0, c))))
        # shortest path length donor -> acceptor in reactant graph, skipping the transferring H
        dist, front, L = {d: 0}, [d], None
        while front and L is None:
            nxt = []
            for u in front:
                for w in nb[u]:
                    if w == h or w in dist:
                        continue
                    dist[w] = dist[u] + 1
                    if w == a:
                        L = dist[w]
                    nxt.append(w)
            front = nxt
        if L is None:
            limit = 90.0
        elif L == 1:
            limit = 35.0
        elif L == 2:
            limit = 45.0
        elif L == 3:
            limit = 62.0
        else:
            limit = 72.0
        if ang < limit:
            return False, f"{rxn.symbols[d]}{d + 1}-H{h + 1}-{rxn.symbols[a]}{a + 1} angle is {ang:.0f} deg (min {limit:.0f})"
    if checked == 0:
        return True, "rule not applicable"
    return True, "transferring H atoms are stretched and sufficiently linear"
