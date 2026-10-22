"""A small pure-Python reader for atom-mapped SMILES, used when RDKit is not installed (e.g. the in-browser app).

It reads exactly what tsso needs: atoms, map numbers, bonds and bond orders, tetrahedral stereo tags, formal charges,
and radical electrons computed the way RDKit's AssignRadicals does for bracket atoms. tests/test_smiles_lite.py
checks that it gives the same answers as RDKit on every RDB7 reaction and on charged, radical and aromatic cases.
"""

# symbol:valence list:outer electrons, taken from RDKit's periodic table
_TABLE = """H:1:1 He:0:2 Li:1,-1:1 Be:2:2 B:3:3 C:4:4 N:3:5 O:2:6 F:1:7 Ne:0:8 Na:1,-1:1 Mg:2,-1:2 Al:3:3 Si:4:4 P:3,5:5
S:2,4,6:6 Cl:1:7 Ar:0:8 K:1,-1:1 Ca:2,-1:2 Sc:-1:3 Ti:-1:4 V:-1:5 Cr:-1:6 Mn:-1:7 Fe:-1:8 Co:-1:9 Ni:-1:10
Cu:-1:11 Zn:-1:2 Ga:3:3 Ge:4:4 As:3,5:5 Se:2,4,6:6 Br:1:7 Kr:0:8 Rb:1,-1:1 Sr:2,-1:2 Y:-1:3 Zr:-1:4 Nb:-1:5
Mo:-1:6 Tc:-1:7 Ru:-1:8 Rh:-1:9 Pd:-1:10 Ag:-1:11 Cd:-1:2 In:3:3 Sn:2,4:4 Sb:3,5:5 Te:2,4,6:6 I:1,3,5:7
Xe:0,2,4,6:8 Cs:1:1 Ba:2,-1:2 La:-1:3 Ce:-1:4 Pr:-1:3 Nd:-1:4 Pm:-1:5 Sm:-1:6 Eu:-1:7 Gd:-1:8 Tb:-1:9 Dy:-1:10
Ho:-1:11 Er:-1:12 Tm:-1:13 Yb:-1:14 Lu:-1:15 Hf:-1:4 Ta:-1:5 W:-1:6 Re:-1:7 Os:-1:8 Ir:-1:9 Pt:-1:10 Au:-1:11
Hg:-1:2 Tl:-1:3 Pb:2,4:4 Bi:3,5:5 Po:2,4,6:6 At:1,3,5:7 Rn:0:8 Fr:1:1 Ra:2,-1:2 Ac:-1:3 Th:-1:4 Pa:-1:3 U:-1:4
Np:-1:5 Pu:-1:6 Am:-1:7 Cm:-1:8 Bk:-1:9 Cf:-1:10 Es:-1:11 Fm:-1:12 Md:-1:13 No:-1:14 Lr:-1:15 Rf:-1:2 Db:-1:2
Sg:-1:2 Bh:-1:2 Hs:-1:2 Mt:-1:2 Ds:-1:2 Rg:-1:2 Cn:-1:2 Nh:-1:2 Fl:-1:2 Mc:-1:2 Lv:-1:2 Ts:-1:2 Og:-1:2"""
Z, VALENCES, OUTER = {}, {}, {}
for _z, _e in enumerate(_TABLE.split(), start=1):
    _s, _v, _o = _e.split(':')
    Z[_s], VALENCES[_z], OUTER[_z] = _z, [int(v) for v in _v.split(',')], int(_o)
AROMATIC = {'b': 'B', 'c': 'C', 'n': 'N', 'o': 'O', 'p': 'P', 's': 'S', 'se': 'Se', 'as': 'As', 'te': 'Te'}
ORGANIC = ['Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I', 'b', 'c', 'n', 'o', 's', 'p', '*']
BOND = {'-': 'SINGLE', '=': 'DOUBLE', '#': 'TRIPLE', '$': 'QUADRUPLE', ':': 'AROMATIC', '/': 'SINGLE', '\\': 'SINGLE'}
CONTRIB = {'SINGLE': 1.0, 'DOUBLE': 2.0, 'TRIPLE': 3.0, 'QUADRUPLE': 4.0, 'AROMATIC': 1.5}
CW, CCW = 1, 2          # RDKit's CHI_TETRAHEDRAL_CW / _CCW


class Atom:
    __slots__ = ('z', 'aromatic', 'map', 'chiral', 'h', 'charge', 'bracket', 'radicals', 'written', 'created')

    def __init__(self, z, aromatic=False, bracket=False):
        self.z, self.aromatic, self.bracket = z, aromatic, bracket
        self.map = self.chiral = self.h = self.charge = self.radicals = 0
        self.written, self.created = [], []    # neighbours in SMILES order, and in the order bonds are made


class Mol:
    def __init__(self, atoms, bonds):
        self.atoms, self.bonds = atoms, bonds      # bonds: (i, j, 'SINGLE' | 'DOUBLE' | ...)


def _bad(smi):
    return ValueError(f'could not parse SMILES: {smi}')


def _bracket(text, smi):
    i = 0
    while i < len(text) and text[i].isdigit():           # isotope (ignored)
        i += 1
    if text[i:i + 1] == '*':
        a, i = Atom(0, bracket=True), i + 1
    elif len(text) >= i + 2 and text[i:i + 2] in AROMATIC:
        a, i = Atom(Z[AROMATIC[text[i:i + 2]]], True, True), i + 2
    elif text[i:i + 1] in AROMATIC:
        a, i = Atom(Z[AROMATIC[text[i]]], True, True), i + 1
    elif len(text) >= i + 2 and text[i:i + 2] in Z:
        a, i = Atom(Z[text[i:i + 2]], bracket=True), i + 2
    elif text[i:i + 1] in Z:
        a, i = Atom(Z[text[i]], bracket=True), i + 1
    else:
        raise _bad(smi)
    if text[i:i + 1] == '@':
        if text[i:i + 2] == '@@':
            a.chiral, i = CW, i + 2
        elif text[i + 1:i + 3] in ('TH', 'AL', 'SP', 'TB', 'OH'):
            cls, i = text[i + 1:i + 3], i + 3
            j = i
            while j < len(text) and text[j].isdigit():
                j += 1
            num, i = int(text[i:j] or 0), j
            a.chiral = {1: CCW, 2: CW}.get(num, 0) if cls == 'TH' else 3
        else:
            a.chiral, i = CCW, i + 1
    if text[i:i + 1] == 'H':
        j = i + 1
        while j < len(text) and text[j].isdigit():
            j += 1
        a.h, i = int(text[i + 1:j] or 1), j
    if text[i:i + 1] in ('+', '-'):
        sign = 1 if text[i] == '+' else -1
        j = i + 1
        while j < len(text) and text[j].isdigit():
            j += 1
        if j > i + 1:
            a.charge, i = sign * int(text[i + 1:j]), j
        else:
            while j < len(text) and text[j] == text[i]:
                j += 1
            a.charge, i = sign * (j - i), j
    if text[i:i + 1] == ':':
        j = i + 1
        while j < len(text) and text[j].isdigit():
            j += 1
        if j == i + 1:
            raise _bad(smi)
        a.map, i = int(text[i + 1:j]), j
    if i != len(text):
        raise _bad(smi)
    return a


def mol_from_smiles(smi):
    """Parse one side of a reaction (fragments separated by '.'), keeping explicit hydrogens as atoms."""
    s = smi.strip().split()[0] if smi.strip() else ''
    atoms, bonds, stack, rings = [], [], [], {}
    prev, pending, i = None, None, 0

    ring_bonds = []

    def add_bond(x, y, sym, ring=0):
        if x == y or any({x, y} == {b[0], b[1]} for b in bonds + [r[2] for r in ring_bonds]):
            raise _bad(smi)
        kind = BOND[sym] if sym and sym not in '/\\' else (     # '/' and '\\' act like no symbol, as in RDKit
            'AROMATIC' if atoms[x].aromatic and atoms[y].aromatic else 'SINGLE')
        if ring:        # like RDKit, ring closures are added last, in order of ring number
            ring_bonds.append((ring, len(ring_bonds), (x, y, kind)))
        else:
            bonds.append((x, y, kind))

    while i < len(s):
        c = s[i]
        if c == '(':
            if prev is None or s[i - 1] == '(':
                raise _bad(smi)
            stack.append(prev); i += 1
        elif c == ')':
            if not stack or pending or s[i - 1] == '(':
                raise _bad(smi)
            prev = stack.pop(); i += 1
        elif c in BOND:
            if pending or prev is None:
                raise _bad(smi)
            pending = c; i += 1
        elif c == '.':
            if pending or stack:
                raise _bad(smi)
            prev = None; i += 1
        elif c.isdigit() or c == '%':
            if prev is None:
                raise _bad(smi)
            if c == '%':
                if s[i + 1:i + 2] == '(':
                    j = s.find(')', i)
                    if j < 0 or not s[i + 2:j].isdigit():
                        raise _bad(smi)
                    num, i = int(s[i + 2:j]), j + 1
                elif s[i + 1:i + 3].isdigit() and len(s[i + 1:i + 3]) == 2:
                    num, i = int(s[i + 1:i + 3]), i + 3
                else:
                    raise _bad(smi)
            else:
                num, i = int(c), i + 1
            if num in rings:
                other, sym = rings.pop(num)
                if pending and sym and BOND[pending] != BOND[sym]:
                    raise _bad(smi)
                add_bond(other, prev, pending or sym, ring=num)
                slot = atoms[other].written.index(('ring', num))
                atoms[other].written[slot] = prev
                atoms[prev].written.append(other)
            else:
                rings[num] = (prev, pending)
                atoms[prev].written.append(('ring', num))
            pending = None
        else:
            if c == '[':
                j = s.find(']', i)
                if j < 0:
                    raise _bad(smi)
                a, i = _bracket(s[i + 1:j], smi), j + 1
            else:
                tok = next((t for t in ORGANIC if s.startswith(t, i)), None)
                if tok is None:
                    raise _bad(smi)
                a = Atom(0) if tok == '*' else Atom(Z[AROMATIC.get(tok, tok)], tok.islower())
                i += len(tok)
            atoms.append(a)
            if prev is not None:
                add_bond(prev, len(atoms) - 1, pending)
                atoms[prev].written.append(len(atoms) - 1)
                a.written.append(prev)
            prev, pending = len(atoms) - 1, None
    if rings or stack or pending:
        raise _bad(smi)
    bonds += [b for _, _, b in sorted(ring_bonds)]
    for x, y, _ in bonds:
        atoms[x].created.append(y)
        atoms[y].created.append(x)
    for a in atoms:            # RDKit's CW/CCW refer to the order bonds were made: flip on an odd permutation
        if a.chiral in (CW, CCW) and sorted(a.written) == sorted(a.created) and _odd(a.written, a.created):
            a.chiral = CW + CCW - a.chiral
    return Mol(atoms, bonds)


def _odd(order, target):
    perm, swaps = [target.index(x) for x in order], 0
    for i in range(len(perm)):
        while perm[i] != i:
            j = perm[i]
            perm[i], perm[j] = perm[j], perm[i]
            swaps += 1
    return swaps % 2 == 1


def assign_radicals(m):
    """Radical electrons for bracket atoms (no implicit hydrogens), following RDKit's AssignRadicals (2026 releases):
    a charged atom is treated like its isoelectronic neighbour, e.g. C+ like B and N- like O."""
    total = [0.0] * len(m.atoms)
    for x, y, kind in m.bonds:
        total[x] += CONTRIB[kind]
        total[y] += CONTRIB[kind]
    for k, a in enumerate(m.atoms):
        a.radicals = 0
        if not a.bracket or a.z == 0:
            continue
        valence, chg = int(total[k] + a.h + 0.1), a.charge
        n_eff, base = OUTER[a.z] - chg, (2 if a.z in (1, 2) else 8)
        if VALENCES[a.z] == [-1]:                      # metals: only a bare atom with an odd electron count
            rad = n_eff % 2 if valence == 0 and n_eff > 0 else 0
        else:
            if n_eff >= base // 2:                     # electron-rich: complete the octet (or duet)
                rad = base - n_eff - valence
            else:                                      # electron-poor: unpaired electrons left after bonding
                rad = n_eff - valence if valence <= n_eff else base - n_eff - valence
            if rad < 0:                                # hypervalent: use the next allowed valence
                rad = 0
                for v in (VALENCES[a.z] if len(VALENCES[a.z]) > 1 else []):
                    if v - valence + chg >= 0:
                        rad = v - valence + chg
                        break
            if n_eff >= base // 2 and n_eff - valence >= 0:
                rad = min(rad, n_eff - valence)
        a.radicals = max(rad, 0)
    return m
