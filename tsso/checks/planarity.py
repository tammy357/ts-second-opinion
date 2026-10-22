RULE = "Double bonds and aromatic rings that do not take part in the reaction stay flat."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 3.
# Gate (RDB7 training references): passes 100.0%; rejects 0% of mirror images and 1% of noisy geometries.

def check(rxn, xyz):
    """Spectator C=C/C=N-type double bonds and aromatic rings stay planar."""
    rb, pb, centre, sy = rxn.reactant_bonds, rxn.product_bonds, rxn.reaction_centre, rxn.symbols
    nb, npd = rxn.neighbours('reactant'), rxn.neighbours('product')
    X = [[float(c) for c in row] for row in xyz]

    def sub(a, b): return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    def unit(v):
        m = math.sqrt(dot(v, v)) + 1e-12
        return [v[0] / m, v[1] / m, v[2] / m]
    def vec(a, b): return unit(sub(X[b], X[a]))

    def dihedral(a, b, c, d):
        b0, b1, b2 = sub(X[a], X[b]), vec(b, c), sub(X[d], X[c])
        v = [b0[k] - dot(b0, b1) * b1[k] for k in range(3)]
        w = [b2[k] - dot(b2, b1) * b1[k] for k in range(3)]
        return math.degrees(math.atan2(dot(cross(b1, v), w), dot(v, w)))

    found, worst, msg = False, 0.0, ""
    dbl = [k for k, q in rb.items() if q == 2]
    for (i, j) in dbl:
        if pb.get((i, j)) != 2 or i in centre or j in centre:
            continue
        if sum(1 for k in dbl if i in k) > 1 or sum(1 for k in dbl if j in k) > 1:
            continue  # cumulated double bonds
        si = [a for a in nb[i] if a != j and a in npd[i] and a not in centre]
        sj = [b for b in nb[j] if b != i and b in npd[j] and b not in centre]
        vij = vec(i, j)
        for a in si:
            for b in sj:
                if abs(dot(vec(i, a), vij)) > 0.97 or abs(dot(vec(j, b), vij)) > 0.97:
                    continue
                found = True
                d = abs(dihedral(a, i, j, b))
                dev = min(d, 180.0 - d)
                if dev > worst:
                    worst = dev
                    msg = (f"double bond {sy[i]}{i + 1}={sy[j]}{j + 1} twisted: dihedral "
                           f"{sy[a]}{a + 1}-{sy[b]}{b + 1} deviates {dev:.0f} deg from planar")
    if worst > 40.0:
        return False, msg

    parent = {}
    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for (i, j), o in rb.items():
        if o == 1.5 and pb.get((i, j)) == 1.5:
            parent[find(i)] = find(j)
    groups = {}
    for a in list(parent):
        groups.setdefault(find(a), []).append(a)
    for atoms in groups.values():
        if len(atoms) < 5 or any(a in centre for a in atoms):
            continue
        found = True
        cen = [sum(X[a][k] for a in atoms) / len(atoms) for k in range(3)]
        P = [sub(X[a], cen) for a in atoms]
        C = [[sum(p[r] * p[c] for p in P) for c in range(3)] for r in range(3)]
        tr = C[0][0] + C[1][1] + C[2][2]
        M = [[(tr if r == c else 0.0) - C[r][c] for c in range(3)] for r in range(3)]
        n = [0.577, 0.5, 0.64]
        for _ in range(100):
            n = unit([dot(M[r], n) for r in range(3)])
        devs = [abs(dot(p, n)) for p in P]
        k = devs.index(max(devs))
        if devs[k] > 0.3:
            a = atoms[k]
            return False, f"aromatic ring atom {sy[a]}{a + 1} is {devs[k]:.2f} A out of the ring plane"
    if not found:
        return True, "rule not applicable"
    return True, "spectator double bonds and aromatic rings are planar"
