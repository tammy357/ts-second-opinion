"""Ask an LLM (Claude) to turn a plain-English rule into a Python check.

Backends, tried in this order by `backend="auto"`:
  api    - Anthropic API, if the `anthropic` package and ANTHROPIC_API_KEY are available
  claude - the Claude Code CLI (`claude -p`), which works with a Claude Pro/Max plan
  cache  - checks already saved in tsso/checks (offline demo; only for rules seen before)
"""
import os, re, shutil, subprocess
from pathlib import Path

CHECK_DIR = Path(__file__).resolve().parent / 'checks'
API_MODEL = os.environ.get('TSSO_MODEL', 'claude-sonnet-5-5')
CLI_MODEL = os.environ.get('TSSO_CLAUDE_MODEL')          # e.g. "sonnet"; default = CLI default

SYSTEM = """You write small, careful Python validators for chemistry. Reply with exactly one ```python code block and nothing else."""

API_DOC = '''Write a function

    def check(rxn, xyz):
        """<one line saying what is checked>"""
        ...
        return passed, detail      # passed: bool, detail: short human-readable str

that tests whether a candidate transition-state (TS) geometry obeys the RULE below.

Inputs
- xyz: numpy array (N, 3), candidate TS coordinates in Angstrom.
- rxn: read-only reaction info (the true TS is NOT available):
    rxn.n                    number of atoms
    rxn.atoms                numpy int array of atomic numbers (H=1, C=6, N=7, O=8)
    rxn.symbols              list of element symbols
    rxn.reactant_xyz, rxn.product_xyz   numpy (N, 3) optimised reactant / product geometries
    rxn.reactant_bonds, rxn.product_bonds   dict {(i, j): bond order} with i < j; order 1, 2, 3 or 1.5 (aromatic)
    rxn.changed_bonds        list of (i, j) bonds that form or break
    rxn.reaction_centre      set of atom indices in forming/breaking bonds
    rxn.neighbours(which)    dict {atom: sorted list of bonded atoms}, which = 'reactant' or 'product'
    rxn.reactant_chiral_tags, rxn.product_chiral_tags   numpy int arrays; 1 or 2 = tetrahedral stereocentre, 0 = none
                             (RDKit CW/CCW tags refer to SMILES neighbour order, so compare handedness by geometry,
                             e.g. the sign of a signed volume, rather than by the tag value)
- Atom i is shown to people as symbol + str(i + 1), e.g. "C3".

Rules for the code
- `np` (numpy) and `math` are already available. Do not import anything. No file, network or system access.
- No double-underscore names. Keep it under 60 lines and fast (it runs on ~10,000 reactions).
- Do not wrap the function in try/except: errors must surface so they can be fixed.
- Choose sensible numeric thresholds yourself. Real TSs stretch forming/breaking bonds a lot.
- The check is tested on ~10,000 real DFT transition states (it must accept at least 99%) and on distorted ones
  (mirror images and 0.3 Angstrom random noise), where it should catch clear violations of its rule.
- If the rule does not apply to a reaction (e.g. no stereocentres), return (True, "rule not applicable").
- The detail string for a failure should name the atoms involved and the measured value.

Example of the expected style (for a different rule: "every hydrogen stays close to some heavy atom"):

```python
def check(rxn, xyz):
    """Every hydrogen has a heavy atom within 1.6 Angstrom."""
    for h in range(rxn.n):
        if rxn.atoms[h] != 1:
            continue
        heavy = [j for j in range(rxn.n) if rxn.atoms[j] != 1]
        if not heavy:
            return True, "rule not applicable"
        d = min(float(np.linalg.norm(xyz[h] - xyz[j])) for j in heavy)
        if d > 1.6:
            return False, f"H{h + 1} is {d:.2f} A from the nearest heavy atom"
    return True, "every hydrogen is near a heavy atom"
```
'''


def build_prompt(rule, feedback=None, previous=None):
    p = f'{API_DOC}\nRULE: {rule}\n'
    if previous and feedback:
        p += f'\nYour previous attempt was:\n```python\n{previous}\n```\n\nIt was tested on real reference transition states and {feedback}\n'
    return p


def extract_code(text):
    m = re.search(r'```(?:python)?\s*\n(.*?)```', text, re.S)
    return (m.group(1) if m else text).strip() + '\n'


def available_backends():
    out = []
    if os.environ.get('ANTHROPIC_API_KEY'):
        try:
            import anthropic  # noqa: F401
            out.append('api')
        except ImportError:
            pass
    if shutil.which('claude'):
        out.append('claude')
    out.append('cache')
    return out


def ask(prompt, backend='auto', timeout=240):
    """Return (code, backend_used, model_label)."""
    if backend == 'auto':
        backend = available_backends()[0]
    if backend == 'api':
        import anthropic
        client = anthropic.Anthropic()
        msg = client.messages.create(model=API_MODEL, max_tokens=2000, system=SYSTEM,
                                     messages=[{'role': 'user', 'content': prompt}])
        return extract_code(msg.content[0].text), 'api', API_MODEL
    if backend == 'claude':
        cmd = ['claude', '-p'] + (['--model', CLI_MODEL] if CLI_MODEL else [])
        res = subprocess.run(cmd, input=SYSTEM + '\n\n' + prompt, capture_output=True, text=True, timeout=timeout)
        if res.returncode != 0:
            raise RuntimeError(f'claude CLI failed: {res.stderr.strip()[:300]}')
        return extract_code(res.stdout), 'claude', CLI_MODEL or 'Claude Code CLI default'
    if backend == 'cache':
        raise LookupError('cache backend has no new code to offer; use lookup_cached()')
    raise ValueError(f'unknown backend {backend}')


API_URL = 'https://api.anthropic.com/v1/messages'


def request_body(prompt, model=None):
    return {'model': model or API_MODEL, 'max_tokens': 2000, 'system': SYSTEM, 'messages': [{'role': 'user', 'content': prompt}]}


def request_headers(key):
    # the last header lets a web page call the API directly; the key stays in the visitor's own browser
    return {'x-api-key': key, 'anthropic-version': '2023-06-01', 'content-type': 'application/json',
            'anthropic-dangerous-direct-browser-access': 'true'}


async def ask_async(prompt, key, fetch, model=None):
    """Ask Claude over plain HTTPS with the caller's own API key (for the in-browser app, where fetch is
    pyodide.http.pyfetch). Returns (code, backend_used, model_label)."""
    import json
    resp = await fetch(API_URL, method='POST', headers=request_headers(key), body=json.dumps(request_body(prompt, model)))
    data = json.loads(await resp.string())
    if not resp.ok:
        msg = (data.get('error') or {}).get('message', str(data)[:200])
        raise RuntimeError(f'the Anthropic API said: {msg}')
    text = ''.join(b.get('text', '') for b in data.get('content', []))
    return extract_code(text), 'api (browser)', data.get('model', model or API_MODEL)


def normalise(rule):
    return re.sub(r'\s+', ' ', rule.strip().lower())


def lookup_cached(rule):
    for p in sorted(CHECK_DIR.glob('*.py')):
        src = p.read_text()
        m = re.search(r'^RULE = (.+)$', src, re.M)
        if m:
            try:
                text = eval(m.group(1), {'__builtins__': {}})  # a plain string literal written by save_check
            except Exception:
                continue
            if normalise(text) == normalise(rule):
                return src, p.stem
    return None, None
