"""TS Second Opinion: check, repair and triage machine-learned transition-state guesses before spending DFT.

    from tsso import second_opinion
    result = second_opinion(None, "reactant.xyz", "product.xyz", "guesses.xyz")   # or an atom-mapped SMILES first
    print(result["verdict"], result["reasons"])
"""
__version__ = "0.3.0"


def second_opinion(*args, **kwargs):
    """Triage TS guesses for any reaction from any generator. See tsso.service.second_opinion for arguments."""
    from .service import second_opinion as _so
    return _so(*args, **kwargs)
