RULE = "Atoms that are never bonded to each other, in the reactant or the product, are not squashed unrealistically close together."
# Written by Claude Code CLI default via claude on 2026-10-21, attempt 1.
# Gate (RDB7 training references): passes 99.3%; rejects 1% of mirror images and 24% of noisy geometries.

def check(rxn, xyz):
    """Atom pairs never bonded in reactant or product are not closer than a minimum contact distance."""
    bonded = set(rxn.reactant_bonds) | set(rxn.product_bonds)
    n = rxn.n
    if n < 2:
        return True, "rule not applicable"
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    heavy = rxn.atoms != 1
    worst = None
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in bonded:
                continue
            if heavy[i] and heavy[j]:
                lim = 1.55
            elif heavy[i] or heavy[j]:
                lim = 1.15
            else:
                lim = 1.05
            ratio = d[i, j] / lim
            if ratio < 1.0 and (worst is None or ratio < worst[0]):
                worst = (ratio, i, j, lim)
    if worst is not None:
        _, i, j, lim = worst
        return False, (f"non-bonded {rxn.symbols[i]}{i + 1}-{rxn.symbols[j]}{j + 1} "
                       f"is {d[i, j]:.2f} A apart (minimum {lim:.2f} A)")
    return True, "no non-bonded pair is unrealistically close"
