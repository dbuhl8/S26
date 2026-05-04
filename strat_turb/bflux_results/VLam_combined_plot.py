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


def fit_power_law(blocks, datasets, x_fn):
    """Fit VLam = C * x^(-1/4) by least-squares in log-log space."""
    xs, ys = [], []
    for idx, *_ in datasets:
        block = blocks[idx]
        x = x_fn(block)
        vlam = 1.0 - block[:, 19]
        mask = np.isfinite(x) & np.isfinite(vlam) & (x > 0) & (vlam > 0)
        xs.append(x[mask])
        ys.append(vlam[mask])
    xs = np.concatenate(xs)
    ys = np.concatenate(ys)
    if len(xs) < 2:
        return None
    # log(vlam) = log(C) - 0.25*log(x)  =>  log(C) = mean(log(vlam) + 0.25*log(x))
    C = np.exp(np.mean(np.log(ys) + 0.25 * np.log(xs)))
    return C, xs.min(), xs.max()


def plot_panel(ax, blocks, datasets, x_fn, x_err_fn, xlabel, x_sym, show_ylabel):
    all_y = []
    for idx, label, color, marker in datasets:
        block = blocks[idx]
        x    = x_fn(block)
        xerr = x_err_fn(block)
        vturb = block[:, 19]
        verr  = block[:, 20]
        y    = 1.0 - vturb       # VLam
        yerr = verr               # d(1 - VTurb) = dVTurb
        mask = (np.isfinite(x) & np.isfinite(xerr) &
                np.isfinite(y) & np.isfinite(yerr) & (x > 0) & (y > 0))
        ax.errorbar(x[mask], y[mask], xerr=xerr[mask], yerr=yerr[mask],
                    fmt=marker, color=color, markersize=8,
                    capsize=3, linestyle='none', label=label,
                    elinewidth=1.2, capthick=1.2)
        all_y.append(y[mask] - yerr[mask])   # include error bar extents
        all_y.append(y[mask] + yerr[mask])

    # set y limits with a small log-space margin around the data
    all_y = np.concatenate(all_y)
    all_y = all_y[all_y > 0]
    log_pad = 0.08   # padding in decades
    y_lo = 10 ** (np.log10(all_y.min()) - log_pad)
    y_hi = 10 ** (np.log10(all_y.max()) + log_pad)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel(xlabel, fontsize=FONT_LABEL)
    ax.tick_params(labelsize=FONT_TICK, which='both')

    if show_ylabel:
        ax.set_ylabel(r'$V_\mathrm{Lam}$', fontsize=FONT_LABEL)

    result = fit_power_law(blocks, datasets, x_fn)
    if result is not None:
        C, x_min, x_max = result
        x_line = np.logspace(np.log10(x_min), np.log10(x_max), 300)
        h, = ax.plot(x_line, C * x_line**(-0.25),
                     '--', color='k', linewidth=2.0, zorder=5)
        fit_label = fr'$\propto ({x_sym})^{{-1/4}}$'
        ax.legend(handles=[h], labels=[fit_label],
                  loc='upper right', fontsize=FONT_LEG,
                  frameon=True, framealpha=0.85)
        print(f"  {xlabel} fit: VLam = {C:.4f} * {x_sym}^(-1/4)")


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

    plot_panel(axes[0], blocks, datasets,
               compute_x_Reb, compute_xerr_Reb,
               xlabel=r'$Re_B^*$', x_sym=r'Re_B^*',
               show_ylabel=True)

    plot_panel(axes[1], blocks, datasets,
               compute_x_ReG, compute_xerr_ReG,
               xlabel=r'$Re_G$', x_sym=r'Re_G',
               show_ylabel=False)

    axes[1].xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks=20))
    axes[1].xaxis.set_major_formatter(LogFormatterMathtext())
    axes[1].tick_params(labelsize=FONT_TICK, which='both')

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='lower center', ncol=4,
               fontsize=FONT_LEG,
               bbox_to_anchor=(0.5, 0.02),
               frameon=True)

    fig.savefig('VLam_combined.pdf', bbox_inches='tight')
    fig.savefig('VLam_combined.png', dpi=150, bbox_inches='tight')
    print("Saved VLam_combined.pdf and VLam_combined.png")


if __name__ == '__main__':
    main()
