RULE = "Stereocentres that are not part of the reaction keep the handedness they have in the reactant."

def check(rxn, xyz):
    nb = rxn.neighbours('reactant')
    flipped = []
    for a in range(rxn.n):
        if rxn.reactant_chiral_tags[a] not in (1, 2) or a in rxn.reaction_centre or len(nb[a]) < 3:
            continue
        ref = rxn.reactant_xyz[nb[a][:3]] - rxn.reactant_xyz[a]
        new = xyz[nb[a][:3]] - xyz[a]
        s_ref = np.dot(ref[0], np.cross(ref[1], ref[2]))
        s_new = np.dot(new[0], np.cross(new[1], new[2]))
        if np.sign(s_ref) != np.sign(s_new):
            flipped.append(a)
    if flipped:
        names = ", ".join(f"{rxn.symbols[a]}{a + 1}" for a in flipped)
        return False, f"stereocentre {names} has the opposite handedness to the reactant"
    return True, "spectator stereocentres keep their handedness"
