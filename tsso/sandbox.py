"""Run rule checks (hand- or LLM-written Python) against a restricted view of a reaction.

A check is Python source defining:

    def check(rxn, xyz):
        return passed: bool, detail: str

`rxn` exposes only what is known before any quantum chemistry is run: atoms, symbols,
reactant/product geometries and bonds, the reaction centre and reactant/product stereo tags.
The reference transition state and the model predictions are hidden, so a check cannot
peek at the answer. `xyz` is the (N, 3) candidate transition state in Angstrom.
"""
import math
import re
import numpy as np

ALLOWED_BUILTINS = {k: __builtins__[k] if isinstance(__builtins__, dict) else getattr(__builtins__, k)
                    for k in ['abs', 'all', 'any', 'bool', 'dict', 'enumerate', 'float', 'int', 'len', 'list',
                              'max', 'min', 'range', 'round', 'set', 'sorted', 'str', 'sum', 'tuple', 'zip',
                              'isinstance', 'ValueError', 'Exception', 'True', 'False', 'None', 'reversed',
                              'frozenset', 'map', 'filter']}
FORBIDDEN = re.compile(r'(__|\bimport\b|\bopen\s*\(|\beval\s*\(|\bexec\s*\(|\bcompile\s*\(|\bglobals\b|'
                       r'\blocals\b|\bgetattr\b|\bsetattr\b|\bos\.|\bsys\.|subprocess|ts_ref|preds)')


class RxnView:
    """Read-only reaction information a check may use (no reference TS, no predictions)."""

    def __init__(self, rxn):
        self.atoms = rxn.atoms.copy()
        self.symbols = list(rxn.symbols)
        self.n = rxn.n
        self.reactant_xyz = rxn.reactant_xyz.copy()
        self.product_xyz = rxn.product_xyz.copy()
        self.reactant_bonds = dict(rxn.reactant_bonds)
        self.product_bonds = dict(rxn.product_bonds)
        self.changed_bonds = list(rxn.changed_bonds)
        self.reaction_centre = set(rxn.reaction_centre)
        self.reactant_chiral_tags = rxn.reactant_chiral_tags.copy()
        self.product_chiral_tags = rxn.product_chiral_tags.copy()
        self._nb = {'reactant': rxn.neighbours('reactant'), 'product': rxn.neighbours('product')}
        self.smiles = rxn.smiles

    def neighbours(self, which='reactant'):
        return {k: list(v) for k, v in self._nb[which].items()}


def compile_check(source):
    """Static screen + compile. Returns the check function or raises ValueError."""
    code_only = re.sub(r'(?s)""".*?"""|\'\'\'.*?\'\'\'|#[^\n]*', '', source)
    bad = FORBIDDEN.search(code_only)
    if bad:
        raise ValueError(f'disallowed construct in check: {bad.group(0)!r}')
    ns = {'np': np, 'math': math, '__builtins__': ALLOWED_BUILTINS}
    exec(compile(source, '<check>', 'exec'), ns)
    if 'check' not in ns or not callable(ns['check']):
        raise ValueError('source must define check(rxn, xyz)')
    return ns['check']


def run_check(fn, rxn, xyz):
    try:
        out = fn(RxnView(rxn), np.array(xyz, dtype=float))
        ok, detail = (out if isinstance(out, tuple) else (bool(out), ''))
        return bool(ok), str(detail)
    except Exception as e:  # a crashing check counts as a failure of the check, not the geometry
        return None, f'check crashed: {type(e).__name__}: {e}'
