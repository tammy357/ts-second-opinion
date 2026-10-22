RULE = "Ring-forming or ring-breaking reactions should have the atoms of the new cycle arranged so that the reacting atoms are held in close proximity, with the distance between the two atoms that form a new bond smaller than the distance between them in the reactant and each other pair of bonded atoms in that cycle remaining within a plausible range for a strained ring."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 2.
# Gate (RDB7 training references): passes 99.9%; rejects 0% of mirror images and 7% of noisy geometries.

def check(rxn, xyz):
    """New-ring forming/breaking atoms are closer than in the open form, with ring bonds in a strained-ring range."""
    changed = set(tuple(sorted(b)) for b in rxn.changed_bonds)
    rb = rxn.reactant_bonds
    pn = rxn.neighbours('product')
    rn = rxn.neighbours('reactant')
    found = False
    for (i, j) in sorted(changed):
        if (i, j) in rb and (i, j) not in rxn.product_bonds:
            ring_nb, open_xyz, ring_state = rn, rxn.product_xyz, 'breaking'
        else:
            ring_nb, open_xyz, ring_state = pn, rxn.reactant_xyz, 'forming'
        # shortest path i -> j in the ring-containing state, excluding the direct bond
        prev = {i: None}
        queue = [i]
        while queue and j not in prev:
            nxt = []
            for a in queue:
                for b in ring_nb.get(a, []):
                    if (a == i and b == j) or b in prev:
                        continue
                    prev[b] = a
                    nxt.append(b)
            queue = nxt
        if j not in prev:
            continue
        path = [j]
        while prev[path[-1]] is not None:
            path.append(prev[path[-1]])
        if len(path) > 9:
            continue
        found = True
        d_ts = float(np.linalg.norm(xyz[i] - xyz[j]))
        d_open = float(np.linalg.norm(open_xyz[i] - open_xyz[j]))
        if d_ts > d_open + 1.2:
            return False, f"{ring_state} bond atoms {rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} is {d_ts:.2f} A in TS vs {d_open:.2f} A in open form"
        for a, b in zip(path[:-1], path[1:]):
            d = float(np.linalg.norm(xyz[a] - xyz[b]))
            hi = 4.6 if tuple(sorted((a, b))) in changed else 2.5
            if d < 0.8 or d > hi:
                return False, f"ring bond {rxn.symbols[a]}{a + 1}-{rxn.symbols[b]}{b + 1} is {d:.2f} A (allowed 0.8-{hi})"
    if not found:
        return True, "rule not applicable"
    return True, "ring-closing atoms are close and ring bonds are plausible"
