"""Tests for 0.3: geometry-only input, batch mode, feedback links, agent setup and release metadata."""
import asyncio, csv, json, re
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

import numpy as np

import tsso
from tsso import service
from tsso.chem import parse_reaction
from tsso.data import load_reactions
from tsso.infer import bonded_pairs, reaction_smiles_from_geometry

ROOT = Path(__file__).resolve().parent.parent
EX = ROOT / 'examples' / 'rdb7_3509'
R = {r.id: r for r in load_reactions()}


def centres(tags):
    return sorted(np.nonzero(np.isin(np.asarray(tags), (1, 2)))[0].tolist())


def bonds(rows):
    return {(a, b) for a, b, _ in rows}


def test_geometry_route_matches_the_smiles_route_on_the_example():
    files = (EX / 'reactant.xyz', EX / 'product.xyz', [EX / 'guesses.xyz'])
    given = service.second_opinion((EX / 'reaction.smi').read_text().split()[0], *files)
    read = service.second_opinion(None, *files)
    assert read['smiles_source'].startswith('read from') and given['smiles_source'] == 'given'
    for key in ('verdict', 'chosen_guess', 'mirror_repaired'):
        assert read[key] == given[key], key
    assert read['verdict'] == 'green' and read['mirror_repaired']
    assert read['inference_notes'] and (read['charge'], read['multiplicity']) == (0, 1)


def test_inferred_bonds_and_stereocentres_match_the_dataset_on_a_sample():
    ids = sorted(R)[::30]
    same = np.zeros(3, int)
    for i in ids:
        r = R[i]
        atoms, chiral, rb, pb, pchiral = parse_reaction(
            reaction_smiles_from_geometry(r.symbols, r.reactant_xyz, r.product_xyz)['smiles'], with_product_chirality=True)
        assert list(atoms) == list(r.atoms)
        same += [bonds(rb) == set(r.reactant_bonds), bonds(pb) == set(r.product_bonds),
                 centres(chiral) == centres(r.reactant_chiral_tags)]
    assert same[0] >= len(ids) - 1 and same[1] >= len(ids) - 5 and same[2] >= len(ids) - 2, same


def test_a_strained_cage_is_not_over_bonded():
    r = R[5877]          # a transannular contact once gave a carbon five neighbours and hid its stereocentres
    for X in (r.reactant_xyz, r.product_xyz):
        counts = np.bincount(np.array(bonded_pairs(r.atoms, X)).ravel(), minlength=r.n)
        assert all(counts[i] <= 4 for i in range(r.n) if r.atoms[i] == 6)
    atoms, chiral, rb, pb, pchiral = parse_reaction(
        reaction_smiles_from_geometry(r.symbols, r.reactant_xyz, r.product_xyz)['smiles'], with_product_chirality=True)
    assert centres(pchiral) == centres(r.product_chiral_tags) == [1, 2, 4]


def test_inference_notes_open_shell_and_unchanged_structures():
    r = R[3509]
    open_shell = reaction_smiles_from_geometry(r.symbols, r.reactant_xyz, r.product_xyz, multiplicity=2)
    assert not open_shell['bond_orders'] and any('open-shell' in n for n in open_shell['notes'])
    same = reaction_smiles_from_geometry(r.symbols, r.reactant_xyz, r.reactant_xyz)
    assert any('No bond forms or breaks' in n for n in same['notes'])


def test_without_a_smiles_every_file_must_list_the_same_atoms():
    lines = (EX / 'product.xyz').read_text().splitlines()
    body = lines[2:]
    swapped = '\n'.join(lines[:2] + [body[-1]] + body[1:-1] + [body[0]])     # first and last atoms swapped
    res = service.second_opinion(None, (EX / 'reactant.xyz').read_text(), swapped, [(EX / 'guesses.xyz').read_text()])
    assert 'error' in res and 'same atoms in the same order' in res['error']
    res = service.second_opinion(None, np.zeros((14, 3)), np.zeros((14, 3)), [np.zeros((14, 3))])
    assert 'error' in res and 'elements are known' in res['error']


def test_batch_folder_and_csv(tmp_path):
    from tsso.batch import triage_batch
    out = triage_batch(ROOT / 'examples' / 'batch', tmp_path / 'folder')
    assert out['reactions'] == 4 and out['verdicts'] == {'green': 1, 'amber': 1, 'red': 2}
    rows = {r['id']: r for r in csv.DictReader(open(tmp_path / 'folder' / 'verdicts.csv'))}
    assert rows['rxn_3509']['verdict'] == 'green' and rows['rxn_3509']['smiles_source'] == 'given'
    assert rows['rxn_1292']['smiles_source'].startswith('read from')
    assert (tmp_path / 'folder' / 'rxn_1292' / 'orca_optts.inp').is_file()
    assert not (tmp_path / 'folder' / 'rxn_3509' / 'orca_optts.inp').exists()
    out = triage_batch(ROOT / 'examples' / 'batch' / 'reactions.csv', tmp_path / 'csv')
    assert out['reactions'] == 2 and out['checked'] == 2


def test_batch_reports_reactions_it_cannot_check(tmp_path):
    from tsso.batch import triage_batch
    (tmp_path / 'in' / 'half').mkdir(parents=True)
    (tmp_path / 'in' / 'half' / 'reactant.xyz').write_text((EX / 'reactant.xyz').read_text())
    out = triage_batch(tmp_path / 'in', tmp_path / 'out')
    assert out['checked'] == 0 and 'missing product' in out['error_rows'][0]['error']
    assert 'error' in triage_batch(tmp_path / 'nothing_here', tmp_path / 'out')


def test_command_line_check_without_a_smiles(capsys):
    from tsso.cli import main
    main(['check', '--reactant', str(EX / 'reactant.xyz'), '--product', str(EX / 'product.xyz'), str(EX / 'guesses.xyz')])
    text = capsys.readouterr().out
    assert 'TRUST' in text and 'SMILES read from the geometries' in text and 'Was this verdict right?' in text


def test_feedback_link_shares_the_verdict_but_no_structures():
    res = service.second_opinion(None, EX / 'reactant.xyz', EX / 'product.xyz', [EX / 'guesses.xyz'])
    url = res['feedback_url']
    query = parse_qs(urlparse(url).query)
    assert url.startswith('https://github.com/tammy357/ts-second-opinion/issues/new?')
    assert query['template'] == ['verdict_feedback.yml'] and query['verdict'] == ['green']
    assert query['version'] == [tsso.__version__]
    text = unquote(url)
    first_coordinate = (EX / 'reactant.xyz').read_text().splitlines()[2].split()[1]
    assert first_coordinate not in text and '[C' not in text and res['smiles'] not in text
    form = (ROOT / '.github' / 'ISSUE_TEMPLATE' / 'verdict_feedback.yml').read_text()
    fields = set(re.findall(r'^\s+id: (\w+)', form, re.M))
    assert set(query) - {'template', 'title'} <= fields


def test_report_card_shows_the_inferred_smiles_and_the_feedback_link(tmp_path):
    path = tmp_path / 'report.html'
    service.second_opinion(None, EX / 'reactant.xyz', EX / 'product.xyz', [EX / 'guesses.xyz'], report=str(path))
    html = path.read_text()
    assert 'SMILES read from the geometries' in html and 'Was this verdict right?' in html
    assert 'verdict_feedback.yml' in html


def test_agent_skill_and_mcp_settings(tmp_path):
    from tsso import agents
    packaged = agents.SKILL.read_text()
    assert packaged.startswith('---\nname: ts-second-opinion\n')
    assert (ROOT / '.claude' / 'skills' / 'ts-second-opinion' / 'SKILL.md').read_text() == packaged
    assert agents.install_skill(tmp_path).read_text() == packaged
    assert agents.mcp_settings('claude-code').startswith('claude mcp add ts-second-opinion -- uvx --from git+')
    for client in ('claude-desktop', 'cursor'):
        server = json.loads(agents.mcp_settings(client))['mcpServers']['ts-second-opinion']
        assert server['command'] == 'uvx' and server['args'][-2:] == ['tsso', 'mcp']


def test_data_paths_in_a_clone():
    from tsso import paths
    assert paths.IN_CLONE
    assert paths.data_file('goflow8_rdb7_test.json.gz') == ROOT / 'data' / 'goflow8_rdb7_test.json.gz'
    assert paths.output_dir() == ROOT / 'results'


def test_mcp_tools_include_batch_and_make_the_smiles_optional():
    from tsso.mcp_server import mcp
    tools = {t.name: t for t in asyncio.run(mcp.list_tools())}
    assert 'triage_batch' in tools
    schema = tools['triage_structures'].inputSchema
    assert 'mapped_reaction_smiles' in schema['properties']
    assert 'mapped_reaction_smiles' not in schema.get('required', [])


def test_release_metadata_agrees():
    version = tsso.__version__
    pyproject = (ROOT / 'pyproject.toml').read_text()
    assert re.search(r'^version = "([^"]+)"', pyproject, re.M).group(1) == version
    server = json.loads((ROOT / 'server.json').read_text())
    assert server['version'] == version and server['packages'][0]['version'] == version
    assert server['packages'][0]['identifier'] == re.search(r'^name = "([^"]+)"', pyproject, re.M).group(1)
    assert len(server['description']) <= 100
    assert f"mcp-name: {server['name']}" in (ROOT / 'README.md').read_text()
    assert re.search(r'^version: ([\d.]+)', (ROOT / 'CITATION.cff').read_text(), re.M).group(1) == version
    assert f'## {version} ' in (ROOT / 'CHANGELOG.md').read_text()
