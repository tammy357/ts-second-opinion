RULE = "Stereocentres that are not part of the reaction keep the handedness they have in the reactant."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 100.0%; rejects 25% of mirror images and 0% of noisy geometries.

def check(rxn, xyz):
    """Stereocentres outside the reaction centre keep their reactant handedness (sign of signed volume)."""
    nb = rxn.neighbours('reactant')
    checked = 0

    def vol(c, ns, X):
        if len(ns) >= 4:
            o = X[ns[3]]
            a, b, d = X[ns[0]] - o, X[ns[1]] - o, X[ns[2]] - o
        else:
            o = X[c]
            a, b, d = X[ns[0]] - o, X[ns[1]] - o, X[ns[2]] - o
        return float(np.dot(a, np.cross(b, d)))

    for c in range(rxn.n):
        if rxn.reactant_chiral_tags[c] not in (1, 2) or c in rxn.reaction_centre:
            continue
        ns = list(nb.get(c, []))
        if len(ns) < 3:
            continue
        ns = ns[:4]
        v_r = vol(c, ns, rxn.reactant_xyz)
        if abs(v_r) < 0.3:
            continue
        checked += 1
        v_t = vol(c, ns, xyz)
        if v_t * v_r < 0 and abs(v_t) > 0.4:
            return False, (f"stereocentre {rxn.symbols[c]}{c + 1} inverted: signed volume "
                           f"{v_r:.2f} in reactant vs {v_t:.2f} in TS")
    if checked == 0:
        return True, "rule not applicable"
    return True, f"{checked} spectator stereocentre(s) keep their handedness"
