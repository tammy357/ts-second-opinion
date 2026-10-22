RULE = "Atoms that are stereocentres in the product already have the product's handedness at the transition state."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.8%; rejects 33% of mirror images and 2% of noisy geometries.

def check(rxn, xyz):
    """Product stereocentres have the same handedness (signed volume sign) in the TS as in the product."""
    nb = rxn.neighbours('product')
    centres = [i for i in range(rxn.n) if rxn.product_chiral_tags[i] in (1, 2)]

    def vol(c, pos, ns):
        v = []
        for k in ns[:3]:
            d = pos[k] - pos[c]
            v.append(d / (np.linalg.norm(d) + 1e-9))
        return float(np.dot(v[0], np.cross(v[1], v[2])))

    checked = 0
    for c in centres:
        ns = nb.get(c, [])
        if len(ns) < 3:
            continue
        ref = vol(c, rxn.product_xyz, ns)
        if abs(ref) < 0.1:
            continue
        ts = vol(c, xyz, ns)
        checked += 1
        # centres in the reaction centre may legitimately be near-planar at the TS
        thr = 0.3 if c in rxn.reaction_centre else 0.15
        if ts * ref < 0 and abs(ts) > thr:
            return False, (f"{rxn.symbols[c]}{c + 1} has inverted handedness vs product "
                           f"(signed volume {ts:+.2f} at TS, {ref:+.2f} in product)")
    if checked == 0:
        return True, "rule not applicable"
    return True, f"{checked} product stereocentre(s) have product handedness"
