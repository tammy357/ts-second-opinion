"""Fast sanity tests: python -m pytest -q"""
import numpy as np
import pytest
from tsso.data import load_reactions
from tsso.sandbox import compile_check, run_check
from tsso.gate import load_check_sources
from tsso.triage import load_checks, triage

R = load_reactions()


def test_data_shapes():
    assert len(R) == 1194
    r = R[0]
    assert r.preds.shape == (8, r.n, 3) and r.ts_ref_xyz.shape == (r.n, 3)


def test_sandbox_blocks_imports_and_answer_peeking():
    for bad in ['import os\ndef check(rxn, xyz):\n    return True, ""',
                'def check(rxn, xyz):\n    return bool(rxn.ts_ref_xyz is None), ""',
                'def check(rxn, xyz):\n    return open("x").read(), ""']:
        with pytest.raises(ValueError):
            compile_check(bad)


def test_all_checks_accept_reference_ts_mostly():
    for name, src in load_check_sources().items():
        fn = compile_check(src)
        ok = sum(run_check(fn, r, r.ts_ref_xyz)[0] is True for r in R[:300])
        assert ok >= 0.97 * 300, name


def test_mirror_image_flagged_by_stereo_rule():
    fn = compile_check(load_check_sources()['chirality'])
    chiral = [r for r in R if any(r.reactant_chiral_tags[i] in (1, 2) and i not in r.reaction_centre for i in range(r.n))]
    flagged = sum(run_check(fn, r, r.ts_ref_xyz * np.array([1, 1, -1]))[0] is False for r in chiral[:100])
    assert flagged > 50


def test_triage_runs():
    checks = load_checks()
    t = triage(R[0], checks, amber_cutoff=0.11)
    assert t.verdict in ('green', 'amber', 'red') and t.xyz.shape == R[0].ts_ref_xyz.shape


def test_service_and_mcp_tools():
    from tsso import service
    t = service.triage_reaction(3509, score_against_reference=True)
    assert t['mirror_repaired'] and t['score_against_reference']['stereo_correct']
    r = R[0]
    xyz = lambda X: '\n'.join([str(r.n), ''] + [f'{s} {a} {b} {c}' for s, (a, b, c) in zip(r.symbols, X)])
    out = service.triage_structures(r.smiles, xyz(r.reactant_xyz), xyz(r.product_xyz), [xyz(p) for p in r.preds[:3]])
    assert out['verdict'] in ('green', 'amber', 'red')
    one = service.triage_structures(r.smiles, xyz(r.reactant_xyz), xyz(r.product_xyz), [xyz(r.preds[0])])
    assert one['disagreement_A'] is None
    import asyncio
    from tsso.mcp_server import mcp
    names = {t.name for t in asyncio.run(mcp.list_tools())}
    assert {'triage_reaction', 'gate_check', 'save_rule', 'triage_structures', 'certify_trust', 'evaluate_predictions'} <= names


def test_certified_trust_controls_error_rate_on_synthetic_data():
    from tsso.certify import certify
    rng = np.random.default_rng(1)
    errs, shares = [], []
    for _ in range(400):
        s = rng.normal(size=600)
        wrong = rng.random(600) < 1 / (1 + np.exp(3 * s))     # higher score -> less likely wrong
        sel = certify(s[:300], wrong[:300], s[300:], 0.05)
        errs.append(wrong[300:][sel].mean() if sel.any() else 0.0); shares.append(sel.mean())
    assert np.mean(errs) <= 0.05 + 0.01 and np.mean(shares) > 0.2


def test_certify_service():
    from tsso import service
    out = service.certify_trust(0.05, 'geometry')
    assert out['trusted'] > 0 and 0 < out['share_trusted'] <= 1
    assert out['over_2000_random_splits']['mean_error_among_trusted'] <= 0.05


def test_any_element_and_plain_xyz_files(tmp_path):
    """SN2: Cl- + CH3Br -> CH3Cl + Br-. Toy geometries; checks that non-HCNO reactions run end to end."""
    from tsso import second_opinion
    smi = '[Cl-:1].[C:2]([H:3])([H:4])([H:5])[Br:6]>>[Cl:1][C:2]([H:3])([H:4])[H:5].[Br-:6]'
    H = [[0.36, 1.03, 0.0], [0.36, -0.51, 0.89], [0.36, -0.51, -0.89]]
    def frame(cl, br, hx):
        rows = [('Cl', [cl, 0, 0]), ('C', [0, 0, 0])] + [('H', [hx, y, z]) for _, y, z in H] + [('Br', [br, 0, 0])]
        return '6\n\n' + '\n'.join(f'{e} {a:.3f} {b:.3f} {c:.3f}' for e, (a, b, c) in rows)
    (tmp_path / 'r.xyz').write_text(frame(-3.2, 1.94, -0.36))
    (tmp_path / 'p.xyz').write_text(frame(-1.78, 3.4, 0.36))
    (tmp_path / 'g.xyz').write_text('\n'.join(frame(-2.3 + d, 2.5 - d, 0.0) for d in (0.0, 0.03, -0.02)))
    out = second_opinion(smi, tmp_path / 'r.xyz', tmp_path / 'p.xyz', tmp_path / 'g.xyz')
    assert 'error' not in out, out
    assert out['verdict'] in ('green', 'amber', 'red') and out['n_guesses'] == 3
    bad = second_opinion(smi.replace('[H:5]', '[H]'), tmp_path / 'r.xyz', tmp_path / 'p.xyz', tmp_path / 'g.xyz')
    assert 'atom-map' in bad['error']


def test_stereo_aware_metrics_see_mirror_images():
    from tsso.metrics import evaluate
    chiral = [r for r in R if r.id == 3509]
    mirror = [r.ts_ref_xyz * np.array([1, 1, -1]) for r in chiral]
    rep = evaluate(mirror, [r.ts_ref_xyz for r in chiral], chiral)
    assert rep['mean_dmae_A'] < 1e-9 and rep['mean_rmsd_mirror_allowed_A'] < 1e-6 and rep['mean_rmsd_A'] > 0.5
    assert rep['success_dmae_only'] == 1.0 and rep['success_dmae_and_handedness'] == 0.0


def test_evaluate_files_cli_example():
    from tsso import service
    out = service.evaluate_files('examples/evaluate/predictions.xyz', 'examples/evaluate/references.xyz',
                                 'examples/evaluate/reactions.smi')
    assert out['n'] == 200 and out['success_dmae_and_handedness'] < out['success_dmae_only']


def test_report_card(tmp_path):
    from tsso import second_opinion
    d = 'examples/rdb7_3509/'
    out = second_opinion(open(d + 'reaction.smi').read().strip(), d + 'reactant.xyz', d + 'product.xyz', d + 'guesses.xyz',
                         report=tmp_path / 'r.html')
    page = (tmp_path / 'r.html').read_text()
    assert out['report'] and 'Trust' in page and 'inverted' in page and '$3Dmol' in page
