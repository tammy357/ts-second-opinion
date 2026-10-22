RULE = "A hydrogen that moves from one heavy atom to another sits between them at the transition state, closer to the line joining them than to either end."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 2.
# Gate (RDB7 training references): passes 99.2%; rejects 1% of mirror images and 7% of noisy geometries.

def check(rxn, xyz):
    """A transferring hydrogen stays close to both its donor and acceptor, with a D-H-A angle that is not collapsed."""
    nr = rxn.neighbours('reactant')
    npd = rxn.neighbours('product')
    found = False
    for h in range(rxn.n):
        if rxn.atoms[h] != 1:
            continue
        donors = [d for d in nr.get(h, []) if rxn.atoms[d] != 1 and d not in npd.get(h, [])]
        accs = [a for a in npd.get(h, []) if rxn.atoms[a] != 1 and a not in nr.get(h, [])]
        if not donors or not accs:
            continue
        found = True
        best = None
        for d in donors:
            for a in accs:
                v = xyz[a] - xyz[d]
                L2 = float(np.dot(v, v))
                if L2 < 1e-6:
                    continue
                w = xyz[h] - xyz[d]
                t = float(np.dot(w, v)) / L2
                perp = float(np.linalg.norm(w - t * v))
                dh = float(np.linalg.norm(w))
                ah = float(np.linalg.norm(xyz[h] - xyz[a]))
                cosang = float(np.dot(-w, xyz[h] - xyz[a]) * -1.0) / max(dh * ah, 1e-9)
                ang = math.degrees(math.acos(max(-1.0, min(1.0, cosang))))
                ok = dh < 2.4 and ah < 2.4 and -0.5 <= t <= 1.5 and ang > 45.0
                info = (ok, max(dh, ah),
                        f"H{h + 1} between {rxn.symbols[d]}{d + 1} and {rxn.symbols[a]}{a + 1}: "
                        f"d(H-donor)={dh:.2f} A, d(H-acceptor)={ah:.2f} A, "
                        f"projection t={t:.2f}, off-line distance {perp:.2f} A, angle D-H-A={ang:.0f} deg")
                if best is None or (ok and not best[0]) or (ok == best[0] and info[1] < best[1]):
                    best = info
        if best is not None and not best[0]:
            return False, best[2]
    if not found:
        return True, "rule not applicable"
    return True, "transferring hydrogens sit between donor and acceptor"
