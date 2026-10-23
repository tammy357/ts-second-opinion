"""Figure for the new DFT calculations (results/dft/*.json -> results/fig_dft.png). Same style as make_figures.py."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.edgecolor': GRID,
                     'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})

rows = [json.loads(p.read_text()) for p in sorted((ROOT / 'results' / 'dft').glob('*.json'))]
points = []                                            # (row, dmae, right-handed, n_imaginary, label)
for r in rows:
    if 'recommended' in r and 'error' not in r['recommended']:
        g = r['recommended']
        points.append((0 if r['verdict'] == 'green' else 1 if r['verdict'] == 'amber' else 2,
                       g['dmae_to_reference_ts_A'], g['handedness_correct'], g['n_imaginary'], r['reaction_id']))
    if 'unrepaired' in r and 'error' not in r['unrepaired'] and 'error' not in r.get('repaired', {'error': 1}):
        mirror = r
has_amber = any(p[0] == 1 for p in points)
n = {k: sum(p[0] == k for p in points) for k in (0, 1, 2)}
ylab = [f'✓ trust ({n[0]})'] + ([f'! check ({n[1]})'] if has_amber else []) + [f'✕ re-generate ({n[2]})']
ypos = {0: 0, 1: 1, 2: 2 if has_amber else 1}

fig, ax = plt.subplots(figsize=(7.6, 4.4), facecolor=SURF)
ax.set_facecolor(SURF)
ax.axvspan(0.005, 0.2, color='#eef4fb', zorder=0)
ax.axvline(0.2, color=INK2, lw=1, ls=(0, (3, 3)), zorder=1)
ax.text(0.19, -0.62, 'D-MAE ≤ 0.2 Å', ha='right', va='center', color=INK2, fontsize=9)
for k, d, ok, nimag, rid in points:
    y = ypos[k] + (0.14 if ok else -0.14)
    ax.scatter(d, y, s=90, marker='o' if ok else 'X', color=BLUE if ok else ORANGE, edgecolor=SURF, lw=1.5, zorder=3)
    if nimag != 1:
        ax.annotate(f'{nimag} imaginary frequencies', (d, y), xytext=(10, -11), textcoords='offset points',
                    ha='right', va='top', color=INK2, fontsize=8.5)
# the 3509 mirror pair
ym = max(ypos.values()) + 1.1
u, p = mirror['unrepaired'], mirror['repaired']
ax.scatter(p['dmae_to_reference_ts_A'], ym + 0.14, s=90, marker='o', color=BLUE, edgecolor=SURF, lw=1.5, zorder=3)
ax.scatter(u['dmae_to_reference_ts_A'], ym - 0.14, s=90, marker='X', color=ORANGE, edgecolor=SURF, lw=1.5, zorder=3)
ax.text(0.03, ym, f"repaired vs unrepaired: {abs(u['energy_kcal'] - p['energy_kcal']):.4f} kcal/mol apart,\n"
        f"{u['imaginary_cm1']:.0f}i vs {p['imaginary_cm1']:.0f}i cm⁻¹: a frequency check passes both",
        va='center', color=INK, fontsize=9,
        bbox=dict(fc=SURF, ec=GRID, boxstyle='round,pad=0.4'))
ax.set_xscale('log'); ax.set_xlim(0.005, 1.0)
ax.set_xticks([0.01, 0.03, 0.1, 0.2, 0.5, 1.0]); ax.set_xticklabels(['0.01', '0.03', '0.1', '0.2', '0.5', '1'])
ax.minorticks_off()
ax.set_yticks(list(range(len(ylab))) + [ym]); ax.set_yticklabels(ylab + [f"mirror pair\n(reaction {mirror['reaction_id']})"])
ax.set_ylim(ym + 0.6, -0.9)
ax.set_xlabel('D-MAE to the RDB7 reference TS after DFT optimisation (Å)')
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.tick_params(length=0)
ax.xaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
fig.text(0.02, 0.965, 'New DFT: the verdicts hold, and a frequency check passes mirror images',
         color=INK, fontsize=12.5, fontweight='bold', va='top')
fig.text(0.02, 0.905, 'ωB97X/def2-SVP (PySCF) transition-state optimisation with Sella from each recommended guess, then\n'
         'analytic frequencies. One imaginary frequency unless marked.', color=INK2, fontsize=9, va='top')
ax.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=BLUE, ms=9), plt.Line2D([], [], marker='X', ls='', color=ORANGE, ms=9)],
          labels=['right handedness', 'wrong handedness'], frameon=False, loc='upper right', labelcolor=INK2, fontsize=9)
fig.subplots_adjust(left=0.2, right=0.97, top=0.78, bottom=0.15); fig.savefig(ROOT / 'results' / 'fig_dft.png', dpi=200, facecolor=SURF); plt.close(fig)
print('wrote results/fig_dft.png')
