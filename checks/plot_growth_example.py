"""Draw a fixed, illustrative realization of the tree/occupancy coupling.

This figure is a schematic example, not a sampled or fitted data set.
Run directly, or through verify_and_plot.py.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT = Path(__file__).resolve().parents[1]


def plot_growth_example():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'pdf.fonttype': 42, 'savefig.bbox': 'tight'})
    blue, orange, ink = '#176b89', '#b75c23', '#30383d'
    fig, (tree, draws) = plt.subplots(
        1, 2, figsize=(6.35, 1.85), gridspec_kw={'width_ratios': [1, 1.6]})
    tree.set(xlim=(0, 3.0), ylim=(-0.2, 2.4), aspect='equal')
    tree.axis('off')
    tree.text(0.04, 2.27, '(a) A tree with four vertices', fontsize=9, weight='bold')
    pos = {1: (1.25, 1.84), 2: (1.25, 1.15),
           3: (0.58, 0.46), 4: (1.93, 0.46)}
    for parent, child in [(1, 2), (2, 3), (2, 4)]:
        a, b = pos[parent], pos[child]
        tree.plot([a[0], b[0]], [a[1], b[1]],
                  color=blue if child == 4 else ink,
                  lw=2.2 if child == 4 else 1.2, zorder=1)
    for i, (x, y) in pos.items():
        leaf = i in (3, 4)
        tree.add_patch(Circle((x, y), .17, facecolor='#f7dfc9' if leaf else '#edf0f2',
                              edgecolor=orange if leaf else ink, lw=1.2, zorder=2))
        tree.text(x, y, str(i), ha='center', va='center', fontsize=9, zorder=3)
    tree.text(2.06, 1.00, 'new\nedge', color=blue, fontsize=8, va='center')
    tree.text(1.25, -0.10, r'Leaves: $3,4$;  $I_4=2$, $L_4=2$',
              ha='center', fontsize=8.5)

    draws.set(xlim=(0, 4.8), ylim=(-0.2, 2.4))
    draws.axis('off')
    draws.text(0, 2.27, '(b) One unavailable draw is replaced', fontsize=9, weight='bold')
    xs = (1.90, 2.53, 3.16)
    draws.text(0, 1.90, r'Arrival time $t$', fontsize=8.5)
    for x, t in zip(xs, (2, 3, 4)):
        draws.text(x, 1.90, str(t), ha='center', fontsize=8.5)
    draws.text(0, 1.38, r'Zipf draw $X_t$', va='center', fontsize=8.5)
    draws.text(0, .55, r'Parent $P_t$', va='center', fontsize=8.5)
    for row, values in ((1.38, (1, 2, 6)), (.55, (1, 2, 2))):
        for j, (x, value) in enumerate(zip(xs, values)):
            color = (orange if row > 1 else blue) if j == 2 else ink
            draws.text(x, row, str(value), ha='center', va='center', color=color,
                       bbox={'boxstyle': 'round,pad=.22', 'fc': 'white', 'ec': color,
                             'lw': 1.1}, fontsize=9)
    draws.annotate('', xy=(xs[2], .78), xytext=(xs[2], 1.14),
                   arrowprops={'arrowstyle': '->', 'color': blue, 'lw': 1.2})
    draws.text(3.55, 1.38, r'$K=3$', va='center', fontsize=9)
    draws.text(3.55, .55, r'$I=2$', va='center', fontsize=9)
    draws.text(3.45, .97, r'$6\geq4$: unavailable', fontsize=7.6, va='center', color=orange)
    draws.text(0, -.10, r'Replacement $2\in\{1,2,3\}$ is already a parent.', fontsize=8.5)
    fig.subplots_adjust(left=.005, right=.995, bottom=.08, top=.98, wspace=.06)
    (ROOT/'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'figures'/'growth_example.pdf')
    fig.savefig(ROOT/'figures'/'growth_example.png', dpi=200)
    plt.close(fig)


if __name__ == '__main__':
    plot_growth_example()
