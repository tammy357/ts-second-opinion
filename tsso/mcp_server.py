"""TS Second Opinion as an MCP server, so any agent (Claude Code, Claude Desktop, ...) can use it as tools.

Run:  tsso mcp            (after `pip install -e .`)
      python -m tsso.mcp_server
"""
from mcp.server.fastmcp import FastMCP
from . import service

mcp = FastMCP('ts-second-opinion', instructions=(
    'Checks machine-learned transition-state (TS) guesses before any DFT is spent. Rules are plain-English chemistry '
    'checks turned into Python and kept only if they accept >= 99% of 9,540 real DFT TSs (the gate). '
    'Typical flow: list_reactions -> triage_reaction -> explain; to add a rule: write a check, call gate_check, '
    'revise from its feedback, then save_rule. Never use score_against_reference to make decisions; it is for evaluation.'))


@mcp.tool()
def list_reactions(verdict: str | None = None, mirror_repaired: bool | None = None, limit: int = 20, offset: int = 0) -> dict:
    """List reactions in the RDB7 test set (1,194 reactions, 8 GoFlow guesses each) with their triage verdict.
    verdict: 'green' (trust), 'amber' (check with DFT) or 'red' (re-generate). mirror_repaired: filter repaired ones."""
    return service.list_reactions(verdict, mirror_repaired, limit, offset)


@mcp.tool()
def triage_reaction(reaction_id: int, score_against_reference: bool = False) -> dict:
    """Triage one RDB7 test reaction: run every gated rule on all 8 GoFlow guesses, pick the consensus among guesses
    that pass, repair a wrong-handed stereocentre by mirroring if that makes it pass, and return a verdict with
    reasons and the recommended TS geometry (xyz). score_against_reference adds D-MAE and stereo correctness against
    the reference DFT TS (evaluation only)."""
    return service.triage_reaction(reaction_id, score_against_reference)


@mcp.tool()
def triage_structures(mapped_reaction_smiles: str, reactant_xyz: str, product_xyz: str, guesses_xyz: list[str],
                      order: str = 'map') -> dict:
    """Triage your own reaction, from any TS generator and with any elements. Give an atom-mapped reaction SMILES
    with every atom (hydrogens too) mapped 1..N, the reactant and product geometries and one or more TS guesses, as
    xyz text (multi-frame xyz allowed). order='map' if atoms are listed by atom-map number, 'smiles' if listed in the
    order they appear in the SMILES. Several guesses enable the model-disagreement signal. Charge and multiplicity
    are read from the SMILES."""
    return service.triage_structures(mapped_reaction_smiles, reactant_xyz, product_xyz, guesses_xyz, order)


@mcp.tool()
def list_rules() -> dict:
    """List the plain-English rules, their gate results, and whether each is used in triage."""
    return service.list_rules()


@mcp.tool()
def get_rule_source(name: str) -> dict:
    """Show the Python code of a rule, so you can read it or use it as an example."""
    return service.get_rule_source(name)


@mcp.tool()
def gate_check(rule_text: str, python_source: str) -> dict:
    """Test a candidate rule check against 9,540 real DFT transition states (RDB7 training set) and distorted ones.
    python_source must define check(rxn, xyz) -> (passed: bool, detail: str); numpy is available as np and math as
    math; no imports. rxn offers: n, atoms, symbols, reactant_xyz, product_xyz, reactant_bonds, product_bonds
    ({(i, j): order}), changed_bonds, reaction_centre, neighbours('reactant'|'product'), reactant_chiral_tags,
    product_chiral_tags. Accepted if it passes >= 99% of real TSs, never crashes, and rejects >= 1% of mirrored or
    noisy TSs. Returns failing examples and feedback to revise the code."""
    return service.gate_check(rule_text, python_source)


@mcp.tool()
def save_rule(name: str, rule_text: str, python_source: str, use_in_triage: bool = False) -> dict:
    """Save a rule check that passes the gate (re-tested on save). By default it goes to tsso/checks_user and is not
    used in triage; set use_in_triage=True only when a person has asked for it."""
    return service.save_rule(name, rule_text, python_source, use_in_triage)


@mcp.tool()
def write_orca_input(reaction_id: int) -> dict:
    """Write an ORCA TS-optimisation input (wB97X-D4/def2-TZVP OptTS Freq, as in GoFlow's tsopt.inp) for the
    recommended guess of a test reaction. Returns the file paths and input text."""
    return service.write_orca_input(reaction_id)


@mcp.tool()
def audit_rdb7_dataset() -> dict:
    """Data-quality audit of RDB7: reactions whose reactant or product geometry contradicts the atom-mapped SMILES
    (atoms at single-bond distance that the SMILES does not bond, or a SMILES bond whose atoms are far apart).
    Uploaded reactions get the same check automatically (input_warnings in triage_structures)."""
    return service.audit_rdb7(limit=60)


@mcp.tool()
def evaluate_predictions(predictions_xyz: str, references_xyz: str, reaction_smiles: list[str], order: str = 'map') -> dict:
    """Stereo-aware evaluation of any TS generator: multi-frame xyz text of predictions and of reference TSs (one
    frame per reaction, same order) plus the atom-mapped reaction SMILES. Returns D-MAE, RMSD with and without the
    mirror image allowed, success rate, share with correct handedness, and success with correct handedness."""
    return service.evaluate_files(predictions_xyz, references_xyz, reaction_smiles, order)


@mcp.tool()
def certify_trust(target_error_rate: float = 0.05, count: str = 'geometry', seed: int = 0) -> dict:
    """Certified trust (conformal selection): calibrate on a random half of the RDB7 test reactions, then return
    which of the other half can be trusted so that, on average, at most target_error_rate of the trusted ones are
    wrong. count: 'geometry' (D-MAE > 0.2 A; what matters for barriers) or 'geometry_or_handedness'. The guarantee
    assumes new reactions are drawn like the calibration ones; recalibrate when the generator or data change."""
    return service.certify_trust(target_error_rate, count, seed)


@mcp.tool()
def benchmark_summary() -> dict:
    """Headline results, split into geometry errors and handedness errors: rates before and after, how many
    wrong-handed guesses the rules catch and how many are whole mirror images (same energy, so even DFT would not
    flag them), certified-trust numbers, and risk-coverage curves."""
    return service.benchmark_summary()


def main():
    mcp.run()


if __name__ == '__main__':
    main()
