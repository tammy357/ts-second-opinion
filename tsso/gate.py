"""The gate: a check is trusted only if it accepts (almost) every reference DFT transition state."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from .sandbox import compile_check, run_check

CHECK_DIR = Path(__file__).resolve().parent / 'checks'
GATE_PASS_RATE = 0.99


@dataclass
class GateResult:
    name: str
    pass_rate: float
    crash_rate: float
    failures: list          # [(reaction id, detail)]
    accepted: bool

    def summary(self):
        verdict = 'ACCEPTED' if self.accepted else 'REJECTED'
        return (f'{self.name}: {verdict} - passes {100 * self.pass_rate:.1f}% of reference TSs'
                f' (crashes {100 * self.crash_rate:.1f}%)')


def gate(name, source, reactions, threshold=GATE_PASS_RATE, kwargs=None):
    fn = compile_check(source)
    if kwargs:
        base = fn
        fn = lambda rxn, xyz: base(rxn, xyz, **kwargs)
    ok = crash = 0
    failures = []
    for r in reactions:
        passed, detail = run_check(fn, r, r.ts_ref_xyz)
        if passed is None:
            crash += 1; failures.append((r.id, detail))
        elif passed:
            ok += 1
        else:
            failures.append((r.id, detail))
    n = len(reactions)
    res = GateResult(name, ok / n, crash / n, failures, ok / n >= threshold and crash == 0)
    return res, fn


def load_check_sources():
    return {p.stem: p.read_text() for p in sorted(CHECK_DIR.glob('*.py')) if not p.stem.startswith('_')}
