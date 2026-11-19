"""After a DFT run: does the optimised TS join this reactant to this product, with the right handedness?

A TS optimisation that ends with one imaginary frequency has found a saddle point, not necessarily the one for this
reaction. The usual test is to follow the reaction path downhill in both directions (an intrinsic reaction
coordinate, IRC, calculation) and compare the bonding at its two ends with the reactant and the product; automated
workflows have done this with adjacency matrices for years (for example peerj.com/articles/pchem-15, 2020). `irc_check` does that, and adds the comparison those checks leave out: handedness. A mirror-image TS
has the same bonding at both ends of its path, so a bonding-only check passes it; here each stereocentre at each
end is compared with the reactant or product it should be.

    tsso irc-check --reactant r.xyz --product p.xyz --forward job_IRC_F_trj.xyz --backward job_IRC_B_trj.xyz \
                   --smiles reaction.smi

`ORCA_IRC_TEMPLATE` is an ORCA input that runs the TS optimisation, the frequencies and the IRC in one job (the
ORCA manual allows the IRC keyword together with OptTS and Freq). It has not been run by us: check it on one
reaction before queueing many.
"""
import numpy as np

from .infer import bonded_pairs
from .data import ATOMIC_NUMBER

ORCA_IRC_TEMPLATE = """! wB97X D4 def2-TZVP OptTS FREQ IRC
# TS optimisation, frequencies and reaction path (IRC) for {rid}; method follows GoFlow's tsopt.inp
# Guess: {source}
# Afterwards: tsso irc-check --reactant reactant.xyz --product product.xyz \\
#             --forward <this job>_IRC_F_trj.xyz --backward <this job>_IRC_B_trj.xyz

%maxcore 1500
%pal
  nprocs 2
end

%scf
  Convergence Tight
  MaxIter 500
end

%geom
  MaxIter 500
end

%irc
  MaxIter 40
end

* xyzfile {charge} {mult} {xyzname}
"""


def _last_frame(x):
    from .service import _frames
    frames = _frames(x)
    if not frames:
        raise ValueError('an xyz input is empty')
    return frames[-1]


def _signed_volume(xyz, centre, nbrs):
    v = xyz[nbrs[:3]] - xyz[centre]
    return float(np.dot(v[0], np.cross(v[1], v[2])))


def _handedness(tags, bonds, n, reference_xyz, end_xyz):
    """(centres compared, labels of the centres whose handedness differs) for one end of the path."""
    nb = {i: [] for i in range(n)}
    for i, j in bonds:
        nb[i].append(j); nb[j].append(i)
    compared, wrong = 0, []
    for i in range(n):
        if tags[i] in (1, 2) and len(nb[i]) >= 3:
            a, b = _signed_volume(reference_xyz, i, sorted(nb[i])), _signed_volume(end_xyz, i, sorted(nb[i]))
            if abs(a) < 0.05 or abs(b) < 0.05:          # nearly flat: no handedness to compare
                continue
            compared += 1
            if np.sign(a) != np.sign(b):
                wrong.append(i)
    return compared, wrong


def irc_check(reactant, product, forward, backward, mapped_reaction_smiles=None, order='map'):
    """Compare the two ends of a reaction path with the reactant and the product.

    reactant, product: xyz text or .xyz paths. forward, backward: the two IRC branches (xyz text or paths; for a
    trajectory the last frame is used). Every structure must list the same atoms in the same order.
    mapped_reaction_smiles: optional, and needed for the handedness check, because it says which atoms are
    stereocentres. Without it only the bonding is compared and 'handedness' is 'not checked'.

    Returns {'outcome': 'right' | 'wrong', ...}: 'right' means one end has the reactant's bonds, the other the
    product's, and (if checked) every stereocentre at each end has the handedness of the structure it matches.
    Bonds are read from distances (tsso.infer.bonded_pairs), so a path that stopped before its bonds finished forming
    can be reported as not matching: optimise the ends, or run more IRC steps, before calling a guess wrong.
    """
    try:
        (rs, rx), (ps, px) = _last_frame(reactant), _last_frame(product)
        (fs, fx), (bs, bx) = _last_frame(forward), _last_frame(backward)
    except (ValueError, OSError) as e:
        return {'error': f'could not read the xyz input: {e}'}
    if rs is None:
        return {'error': 'give the reactant as xyz text or an .xyz file, so the elements are known'}
    for name, (s, x) in (('product', (ps, px)), ('forward end', (fs, fx)), ('backward end', (bs, bx))):
        if len(x) != len(rx) or (s is not None and s != rs):
            return {'error': f'the {name} does not list the same atoms in the same order as the reactant'}
    atoms = [ATOMIC_NUMBER[s] for s in rs]
    bR, bP = set(bonded_pairs(atoms, rx)), set(bonded_pairs(atoms, px))
    bF, bB = set(bonded_pairs(atoms, fx)), set(bonded_pairs(atoms, bx))
    lab = lambda b: f'{rs[b[0]]}{b[0] + 1}-{rs[b[1]]}{b[1] + 1}'
    if bR == bP:
        return {'outcome': 'wrong', 'connects': False, 'handedness': 'not checked',
                'why': 'the reactant and the product have the same bonds, so bonding cannot tell the two ends apart '
                       '(a conformational change or a stereo-inversion); compare the structures by eye'}
    if bF == bR and bB == bP:
        ends = {'forward': 'reactant', 'backward': 'product'}
    elif bF == bP and bB == bR:
        ends = {'forward': 'product', 'backward': 'reactant'}
    else:
        def nearest(b):
            dr, dp = b ^ bR, b ^ bP
            name, diff = ('reactant', dr) if len(dr) <= len(dp) else ('product', dp)
            return {'nearest': name, 'bonds_that_differ': sorted(lab(x) for x in diff)}
        same = bF == bB
        return {'outcome': 'wrong', 'connects': False, 'handedness': 'not checked',
                'forward_end': nearest(bF), 'backward_end': nearest(bB),
                'why': ('both ends of the path have the same bonds, so this saddle point does not join the reactant to '
                        'the product' if same else
                        'the ends of the path do not have the bonds of this reactant and this product')}
    out = {'connects': True, 'ends': ends}
    if not (mapped_reaction_smiles or '').strip():
        out.update(outcome='right', handedness='not checked',
                   why='the path joins the reactant to the product; handedness was not checked (give the atom-mapped '
                       'SMILES to check it: a mirror-image TS has the same bonds at both ends)')
        return out
    from .service import build_reaction
    rxn, err = build_reaction(mapped_reaction_smiles, reactant, product, [forward, backward], order)
    if err:
        return {'error': err}
    endR, endP = (rxn.preds[-2], rxn.preds[-1]) if ends['forward'] == 'reactant' else (rxn.preds[-1], rxn.preds[-2])
    nR, wR = _handedness(rxn.reactant_chiral_tags, rxn.reactant_bonds, rxn.n, rxn.reactant_xyz, endR)
    nP, wP = _handedness(rxn.product_chiral_tags, rxn.product_bonds, rxn.n, rxn.product_xyz, endP)
    sym = rxn.symbols
    wrong = sorted({f'{sym[i]}{i + 1}' for i in wR + wP})
    if nR + nP == 0:
        out.update(outcome='right', handedness='no stereocentres', why='the path joins the reactant to the product')
    elif wrong:
        out.update(outcome='wrong', handedness='wrong', wrong_centres=wrong,
                   why=f"the path joins the right bonds, but {', '.join(wrong)} {'has' if len(wrong) == 1 else 'have'} the "
                       f"opposite handedness: this is the TS of a stereoisomer of the reaction (a mirror image or a "
                       f"diastereomer)")
    else:
        out.update(outcome='right', handedness='right', stereocentres_compared=nR + nP,
                   why='the path joins the reactant to the product with the right handedness')
    return out
