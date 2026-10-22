"""Make the final-answer figure from results/metrics.json (light theme, validated palette)."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

m = json.load(open('results/metrics.json'))
SURF, INK, INK2, GRID = '#fcfcfb', '#0b0b0b', '#52514e', '#e4e3df'
BLUE, ORANGE = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.edgecolor': GRID,
                     'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})


def rounded_bar(ax, x, h, w, color):
    ax.add_patch(FancyBboxPatch((x - w / 2, 0), w, h, boxstyle='round,pad=0,rounding_size=0.012',
                                mutation_aspect=1 / 30, fc=color, ec=SURF, lw=2))


def style(ax, title, sub):
    ax.set_facecolor(SURF)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    ax.tick_params(length=0)
    ax.set_title(title, loc='left', color=INK, fontsize=13, fontweight='bold', pad=40)
    ax.text(0, 1.03, sub, transform=ax.transAxes, color=INK2, fontsize=9.5, va='bottom')


# Figure: how often the final answer is wrong, split into geometry and handedness errors
fa = m['final_answer']
stages = ['single_model_avg', 'consensus', 'check_aware_consensus', 'check_aware_plus_mirror']
labels = ['one GoFlow\nmodel', 'consensus\nof 8 models', '+ rules pick\na passing guess', '+ mirror-image\nrepair']
geo = [100 * fa[k]['geometry_bad_rate'] for k in stages]
hand = [100 * fa[k]['stereo_bad_rate'] for k in stages]
fig, ax = plt.subplots(figsize=(7.4, 4.8), facecolor=SURF)
w = 0.34
for i in range(4):
    rounded_bar(ax, i - w / 2 - 0.01, geo[i], w, BLUE)
    rounded_bar(ax, i + w / 2 + 0.01, hand[i], w, ORANGE)
    ax.text(i - w / 2 - 0.01, geo[i] + 0.6, f'{geo[i]:.0f}%', ha='center', color=INK2, fontsize=10)
    ax.text(i + w / 2 + 0.01, hand[i] + 0.6, f'{hand[i]:.0f}%', ha='center', color=INK, fontsize=10,
            fontweight='bold' if i in (0, 3) else 'normal')
ax.set_xlim(-0.6, 3.6); ax.set_ylim(0, 22); ax.set_yticks(range(0, 21, 5))
ax.set_xticks(range(4)); ax.set_xticklabels(labels)
ax.set_ylabel('final answers that are wrong (%)')
style(ax, 'Rules fix handedness; geometry errors are flagged, not fixed',
      'Geometry: D-MAE > 0.2 Å. Handedness: a stereocentre the wrong way round.\n1,194 RDB7 test reactions, 8 GoFlow models, no DFT.')
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=BLUE), plt.Rectangle((0, 0), 1, 1, fc=ORANGE)],
          labels=['geometry wrong', 'handedness wrong'], frameon=False, loc='upper right', labelcolor=INK2)
fig.tight_layout(); fig.savefig('results/fig_final_answer.png', dpi=200, facecolor=SURF); plt.close(fig)
print('wrote results/fig_final_answer.png')
