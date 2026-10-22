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


def test_audit_flags_smiles_geometry_contradictions():
    from tsso.audit import audit_reaction
    from tsso.data import load_reactions
    by = {r.id: r for r in load_reactions()}
    assert not audit_reaction(by[3509])
    assert any(i['atoms'] == 'N5-O6' and i['side'] == 'product' for i in audit_reaction(by[1093]))


def test_uploads_warn_when_the_product_geometry_contradicts_the_smiles():
    from tsso.service import second_opinion
    good = second_opinion(SMI, R, P, [G])
    swapped = second_opinion(SMI, R, R, [G])          # the reactant geometry given as the product
    assert good['input_warnings'] == [] and swapped['input_warnings']


def test_local_inversion_fixes_a_partly_flipped_guess_and_keeps_its_bonds():
    import numpy as np
    from tsso.data import load_reactions
    from tsso.metrics import handedness_correct
    from tsso.stereo_edit import edit_stereo
    r = {x.id: x for x in load_reactions()}[5877]
    t = triage(r, load_checks(), amber_cutoff())
    assert not handedness_correct(r, t.xyz, r.ts_ref_xyz)
    y, edited = edit_stereo(r, t.xyz)
    assert edited == ['C2'] and handedness_correct(r, y, r.ts_ref_xyz)
    for a, b in r.product_bonds:
        if 1 in (a, b):
            assert abs(np.linalg.norm(y[a] - y[b]) - np.linalg.norm(t.xyz[a] - t.xyz[b])) < 1e-9


def test_mirror_test_tells_blind_metrics_from_seeing_ones():
    from tsso.mirror_test import run_mirror_test, BUILT_IN, compile_metric, load_pairs
    import pytest
    pairs = load_pairs()
    assert not run_mirror_test(BUILT_IN['dmae'], pairs)['sees_handedness']
    assert not run_mirror_test(BUILT_IN['rmsd_mirror_allowed'], pairs)['sees_handedness']
    assert run_mirror_test(BUILT_IN['rmsd_proper_rotations'], pairs)['sees_handedness']
    with pytest.raises(ValueError):
        compile_metric('import os\ndef metric(a, b):\n    return 0.0\n')
