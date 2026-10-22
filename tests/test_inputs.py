"""Uploads that used to slip through or confuse the reader (found by stress-testing the app's upload path)."""
from pathlib import Path
from tsso.service import build_reaction, amber_cutoff
from tsso.triage import load_checks, triage

EX = Path(__file__).resolve().parent.parent / 'examples' / 'rdb7_3509'
SMI = (EX / 'reaction.smi').read_text().strip()
R, P, G = ((EX / f).read_text() for f in ('reactant.xyz', 'product.xyz', 'guesses.xyz'))


def verdict(r, p, g):
    rxn, err = build_reaction(SMI, r, p, g)
    return err if err else triage(rxn, load_checks(), amber_cutoff()).verdict


def test_byte_order_mark_and_windows_line_endings_are_fine():
    assert verdict('﻿' + R, P.replace('\n', '\r\n'), ['﻿' + G]) == verdict(R, P, [G])


def test_non_finite_coordinates_are_refused():
    lines = G.splitlines()
    lines[2] = lines[2].split()[0] + ' nan 0.0 0.0'
    assert 'finite' in verdict(R, P, ['\n'.join(lines)])


def test_atom_count_line_must_match():
    assert 'coordinate lines follow' in verdict('20\n' + R.split('\n', 1)[1], P, [G])


def test_a_copy_of_the_reactant_or_product_is_not_trusted():
    assert verdict(R, P, [R]) == 'red' and verdict(R, P, [P]) == 'red'
