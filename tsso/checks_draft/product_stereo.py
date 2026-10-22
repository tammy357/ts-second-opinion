RULE = ("Atoms that are stereocentres in the product already have the product's handedness at the "
        "transition state, judged using their neighbours in the product.")

def check(rxn, xyz):
    nb = rxn.neighbours('product')
    flipped = []
    for a in range(rxn.n):
        if rxn.product_chiral_tags[a] not in (1, 2) or len(nb[a]) < 3:
            continue
        ref = rxn.product_xyz[nb[a][:3]] - rxn.product_xyz[a]
        new = xyz[nb[a][:3]] - xyz[a]
        if np.sign(np.dot(ref[0], np.cross(ref[1], ref[2]))) != np.sign(np.dot(new[0], np.cross(new[1], new[2]))):
            flipped.append(a)
    if flipped:
        names = ", ".join(f"{rxn.symbols[a]}{a + 1}" for a in flipped)
        return False, f"stereocentre {names} has the opposite handedness to the product"
    return True, "product stereocentres already have the product's handedness"
