import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import LogLocator, LogFormatterMathtext

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['mathtext.fontset'] = 'stix'

DARK_BLUE  = '#003087'
SKY_BLUE   = '#56B4E9'
ORANGE     = '#E69F00'
VERMILLION = '#D55E00'
GREEN      = '#009E73'
PINK       = '#CC79A7'

FONT_LABEL = 24
FONT_TICK  = 18
FONT_LEG   = 14

FIT_COLOR = 'k'


def load_dat_blocks(filepath):
    blocks = []
    current_rows = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                if current_rows:
                    blocks.append(np.array(current_rows, dtype=float))
                    current_rows = []
                continue
            if line.startswith('#'):
                continue
            try:
                current_rows.append([float(v) for v in line.split()])
            except ValueError:
                pass
    if current_rows:
        blocks.append(np.array(current_rows, dtype=float))
    return blocks


# ── x-axis definitions ───────────────────────────────────────────────────────

def compute_x_Reb(block):
    Re, B, uh = block[:, 0], block[:, 2], block[:, 29]
    return (uh**3) * Re / B

def compute_xerr_Reb(block):
    Re, B, uh, uh_err = block[:, 0], block[:, 2], block[:, 29], block[:, 30]
    return uh_err * 3 * (uh**2) * Re / B

def compute_x_ReG(block):
    return block[:, 17] / block[:, 2]

def compute_xerr_ReG(block):
    return block[:, 18] / block[:, 2]


# ── panel plotter ─────────────────────────────────────────────────────────────

def plot_panel(ax, blocks, datasets, x_fn, x_err_fn,
               xlabel, show_ylabel, x_power, legend_label=None, x_scale=1.0):
    """
    x_power = -0.25  for both panels → plots VTurb = 1 - (x/x_scale)^{-1/4}
    """
    for idx, label, color, marker in datasets:
        block = blocks[idx]
        x    = x_fn(block)
        xerr = x_err_fn(block)
        y    = block[:, 19]    # VTurb
        yerr = block[:, 20]    # VTurb_err
        mask = (np.isfinite(x) & np.isfinite(xerr) &
                np.isfinite(y) & np.isfinite(yerr) & (x > 0))
        ax.errorbar(x[mask], y[mask], xerr=xerr[mask], yerr=yerr[mask],
                    fmt=marker, color=color, markersize=8,
                    capsize=3, linestyle='none', label=label,
                    elinewidth=1.2, capthick=1.2)

    ax.set_xscale('log')
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel(xlabel, fontsize=FONT_LABEL)
    ax.tick_params(labelsize=FONT_TICK, which='both')

    if show_ylabel:
        ax.set_ylabel(r'$V_\mathrm{Turb}$', fontsize=FONT_LABEL)

    # ── theory curve: VTurb = 1 - x^x_power ─────────────────────────────────
    xs = []
    for idx, *_ in datasets:
        block = blocks[idx]
        x = x_fn(block)
        mask = np.isfinite(x) & (x > 0)
        xs.append(x[mask])
    xs = np.concatenate(xs)

    if len(xs) >= 2:
        x_line = np.logspace(np.log10(xs.min()), np.log10(xs.max()), 300)
        y_line = 1.0 - (x_line / x_scale) ** x_power
        h, = ax.plot(x_line, y_line, '--', color=FIT_COLOR, linewidth=2.0, zorder=5)
        lbl = legend_label if legend_label is not None else r'$1 - Re_G/Re_B^{*\,1/4}$'
        ax.legend([h], [lbl],
                  loc='upper left', fontsize=FONT_LEG,
                  frameon=True, framealpha=0.85)


def main():
    blocks = load_dat_blocks('tavg_bflux.dat')
    datasets = [
        (0, 'Steady (300,30)',   ORANGE,    'o'),
        (1, 'Steady (600,30)',   SKY_BLUE,  'o'),
        (2, 'Steady (600,60)',   DARK_BLUE, 'o'),
        (3, 'Steady (1000,10)',  VERMILLION,'o'),
        (4, 'Steady (1000,100)', GREEN,     'o'),
        (5, 'Steady (600,600)', PINK,      'o'),
        (6, 'Stoch (600,60)',    DARK_BLUE, 'D'),
        (7, 'Stoch (1000,100)', GREEN,     'D'),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    fig.subplots_adjust(wspace=0.08, bottom=0.30)

    print("Plotting theory curves (VTurb = 1 - x^{-1/4}):")

    plot_panel(axes[0], blocks, datasets,
               compute_x_Reb, compute_xerr_Reb,
               xlabel=r'$Re_B^*$',
               show_ylabel=True,
               x_power=-0.25,
               legend_label=r'$1 - (Re_B^*/33)^{-1/4}$',
               x_scale=33)

    plot_panel(axes[1], blocks, datasets,
               compute_x_ReG, compute_xerr_ReG,
               xlabel=r'$Re_G$',
               show_ylabel=False,
               x_power=-0.25,
               legend_label=r'$1 - Re_G^{-1/4}$')

    # force every integer power of 10 to show on the ReG axis
    axes[1].xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=20))
    axes[1].xaxis.set_major_formatter(LogFormatterMathtext())
    axes[1].tick_params(labelsize=FONT_TICK, which='both')

    # shared legend below both panels
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='lower center', ncol=4,
               fontsize=FONT_LEG,
               bbox_to_anchor=(0.5, 0.02),
               frameon=True)

    fig.savefig('VTurb_combined.pdf', bbox_inches='tight')
    fig.savefig('VTurb_combined.png', dpi=150, bbox_inches='tight')
    print("Saved VTurb_combined.pdf and VTurb_combined.png")


if __name__ == '__main__':
    main()
