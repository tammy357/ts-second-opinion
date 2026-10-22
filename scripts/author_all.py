"""Author every rule in rules.md with the LLM pipeline (tsso.author)."""
import re, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.author import author

backend = sys.argv[1] if len(sys.argv) > 1 else 'auto'
only = set(sys.argv[2:])
for line in Path('rules.md').read_text().splitlines():
    m = re.match(r'- (\w+): (.+)', line)
    if not m or (only and m.group(1) not in only):
        continue
    name, rule = m.groups()
    t = time.time()
    print(f'\n[{name}] {rule}')
    out = author(rule, name=name, backend=backend)   # saves into tsso/checks (used by triage)
    print(f'  {"ACCEPTED" if out["accepted"] else "REJECTED"} after {out["attempts"]} attempt(s), {time.time() - t:.0f}s')
