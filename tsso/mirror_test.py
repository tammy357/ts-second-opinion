"""The mirror test: does a TS metric see handedness?

Every RDB7 reference TS with a stereocentre (447 in the test set) is compared with its own mirror image, which has the
wrong handedness and is the TS of a different (mirror-image) reaction. A metric that scores these pairs as identical
cannot see handedness errors. Files: benchmark/mirror_test_references.xyz (mirrored) and
benchmark/mirror_test_true_references.xyz (the originals, same order).

    from tsso.mirror_test import run_mirror_test, BUILT_IN
    run_mirror_test(BUILT_IN['dmae'])
"""
import math, re
from pathlib import Path
import numpy as np

from .geometry import dmae
from .metrics import rmsd
from .sandbox import ALLOWED_BUILTINS, FORBIDDEN

ROOT = Path(__file__).resolve().parent.parent
IDENTICAL_BELOW = 0.01      # a score below this (Angstrom-like units) counts as "scored as the same structure"
BUILT_IN = {
    'dmae': lambda pred, ref: dmae(pred, ref),
    'rmsd_mirror_allowed': lambda pred, ref: rmsd(pred, ref, True),
    'rmsd_proper_rotations': lambda pred, ref: rmsd(pred, ref, False),
}
TEMPLATE = '''def metric(pred, ref):
    """pred, ref: (N, 3) arrays in Angstrom, atoms in the same order. Return a score; lower = closer."""
    d = lambda x: np.linalg.norm(x[:, None] - x[None], axis=-1)
    return float(np.abs(d(pred) - d(ref)).mean())
'''


def compile_metric(source):
    """Static screen + compile a user metric (same restrictions as gated checks). Returns metric(pred, ref)."""
    code_only = re.sub(r'(?s)""".*?"""|\'\'\'.*?\'\'\'|#[^\n]*', '', source)
    bad = FORBIDDEN.search(code_only)
    if bad:
        raise ValueError(f'disallowed construct: {bad.group(0)!r}')
    ns = {'np': np, 'math': math, '__builtins__': ALLOWED_BUILTINS}
    exec(compile(source, '<metric>', 'exec'), ns)
    if not callable(ns.get('metric')):
        raise ValueError('the code must define metric(pred, ref)')
    return ns['metric']


def load_pairs():
    from .service import read_xyz_frames
    bench = ROOT / 'benchmark'
    mirrored = [f[1] for f in read_xyz_frames((bench / 'mirror_test_references.xyz').read_text())]
    true = [f[1] for f in read_xyz_frames((bench / 'mirror_test_true_references.xyz').read_text())]
    return mirrored, true


def run_mirror_test(metric, pairs=None):
    mirrored, true = pairs or load_pairs()
    scores = np.array([float(metric(m, t)) for m, t in zip(mirrored, true)])
    same = int((scores < IDENTICAL_BELOW).sum())
    return {'n': len(scores), 'scored_as_identical': same, 'share_scored_as_identical': same / len(scores),
            'median_score': float(np.median(scores)), 'sees_handedness': same / len(scores) < 0.5}
