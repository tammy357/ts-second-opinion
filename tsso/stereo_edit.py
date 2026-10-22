"""Experimental: repair a partly flipped TS guess by inverting only its wrong stereocentres.

Mirroring fixes a guess whose every stereocentre is inverted (the TS of the mirror-image reaction). When only some are
inverted (a diastereomer), mirroring would break the correct ones. Here each wrong centre gets a local "umbrella"
inversion: the centre and its smallest terminal branch (usually one hydrogen) are reflected through the plane of its
other neighbours. That inverts the centre and keeps every bond to it, and inside the moved branch, exactly.

What we found on 9,552 single GoFlow guesses (docs/RESULTS.md, "Can partial flips be repaired?"): it fixed the
handedness of 40 of the 49 guesses it could edit and broke none, but 41 of the 49 already had geometry errors. So it is
offered as a diagnostic in the app, not used by the pipeline.
"""
import numpy as np

from .stereo import stereocentre_status


def _branch(rxn, c, n):
    """Atoms reachable from n without passing through c (reactant and product bonds); None if n is in a ring with c."""
    adj = {i: set() for i in range(rxn.n)}
    for a, b in list(rxn.reactant_bonds) + list(rxn.product_bonds):
        adj[a].add(b); adj[b].add(a)
    seen, stack = {n}, [n]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v == c and u != n:
                return None
            if v != c and v not in seen:
                seen.add(v); stack.append(v)
    return seen


def invert_centre(rxn, xyz, c, neighbours):
    """Umbrella inversion at atom c. Returns (new xyz, moved atoms) or (None, None) if every branch is in a ring."""
    x = np.array(xyz, float)
    if len(neighbours) >= 4:
        options = [(len(b), int(rxn.atoms[n]), n, b) for n in neighbours for b in [_branch(rxn, c, n)] if b is not None]
        if not options:
            return None, None
        _, _, n, moved = min(options, key=lambda o: (o[0], o[1]))
        plane = [k for k in neighbours if k != n][:3]
        moved = set(moved) | {c}
    else:
        plane, moved = list(neighbours[:3]), {c}
    p0, p1, p2 = x[plane[0]], x[plane[1]], x[plane[2]]
    normal = np.cross(p1 - p0, p2 - p0)
    normal /= np.linalg.norm(normal) + 1e-12
    for k in moved:
        x[k] = x[k] - 2 * np.dot(x[k] - p0, normal) * normal
    return x, sorted(moved)


def edit_stereo(rxn, xyz):
    """Invert every stereocentre the stereo rules call wrong. Returns (new xyz or None, [labels of edited centres])."""
    nbr, nbp = rxn.neighbours('reactant'), rxn.neighbours('product')
    x, edited = np.array(xyz, float), []
    for s in stereocentre_status(rxn, xyz):
        if s['ok'] is False:
            ns = nbr[s['atom']][:4] if s['compared_with'] == 'reactant' else nbp[s['atom']]
            y, _ = invert_centre(rxn, x, s['atom'], ns)
            if y is None:
                return None, edited
            x = y
            edited.append(s['label'])
    return x, edited
