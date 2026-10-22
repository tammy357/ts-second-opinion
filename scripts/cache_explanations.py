"""Pre-generate Claude's plain-English explanations for the demo's example reactions (so the hosted demo needs no key).

Usage: python scripts/cache_explanations.py  -> results/explanations.json   (needs the claude CLI)
"""
import json, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tsso.data import load_reactions
from tsso.stereo import stereocentre_status
from tsso.triage import load_checks, triage

EXAMPLES = [3509, 5877, 1292, 2577, 1093]
cut = json.load(open('results/triage_summary.json'))['amber_cutoff_A']


def prompt(r, t):
    st = stereocentre_status(r, t.xyz)
    stereo = '; '.join(f"{s['label']} vs {s['compared_with']}: {'right' if s['ok'] else 'inverted' if s['ok'] is False else 'too flat to judge'}"
                       for s in st) or 'none'
    return ('You are helping a computational chemist decide whether to trust a machine-learned transition-state (TS) guess. '
            'In 3 short sentences, plain English, no headings, no bullet points, explain the verdict and what to do next. '
            'Use only these facts; do not invent numbers, and do not suggest a repair the facts do not support.\n'
            f'Reaction (atom-mapped SMILES): {r.smiles}\nVerdict: {t.verdict} (green = trust, amber = check with DFT, red = re-generate)\n'
            f'Stereocentres of the recommended guess: {stereo}\n' + '\n'.join(t.reasons))


R = {r.id: r for r in load_reactions()}; checks = load_checks(); out = {}
for rid in EXAMPLES:
    r = R[rid]; t = triage(r, checks, cut)
    res = subprocess.run(['claude', '-p'], input=prompt(r, t), capture_output=True, text=True, timeout=300)
    out[str(rid)] = {'verdict': t.verdict, 'explanation': res.stdout.strip()}
    print(rid, t.verdict, res.stdout.strip()[:200].replace('\n', ' '))
json.dump({'model': 'Claude via the Claude Code CLI', 'examples': out}, open('results/explanations.json', 'w'), indent=2)
