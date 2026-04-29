import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

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

# VTurb range used to fit the kick-up power law — adjust as needed
FIT_VTURB_MIN = 0.01
FIT_VTURB_MAX = 0.80


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
    return block[:, 17] / block[:, 2]          # mdisp / B

def compute_xerr_ReG(block):
    return block[:, 18] / block[:, 2]          # mdisp_err / B


# ── kick-up power-law fit ────────────────────────────────────────────────────

def fit_kickup(blocks, datasets, x_fn):
    """Fit VTurb = a * x^b in log-log space, restricted to the transition region."""
    xs, ys = [], []
    for idx, *_ in datasets:
        block = blocks[idx]
        x = x_fn(block)
        y = block[:, 19]   # VTurb
        mask = (np.isfinite(x) & np.isfinite(y) & (x > 0)
                & (y >= FIT_VTURB_MIN) & (y <= FIT_VTURB_MAX))
        xs.append(x[mask])
        ys.append(y[mask])
    xs = np.concatenate(xs)
    ys = np.concatenate(ys)
    slope, intercept = np.polyfit(np.log10(xs), np.log10(ys), 1)
    return 10**intercept, slope, xs.min(), xs.max()


# ── panel plotter ─────────────────────────────────────────────────────────────

def plot_panel(ax, blocks, datasets, x_fn, x_err_fn,
               xlabel, x_sym, show_ylabel, show_legend):

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
    if show_legend:
        ax.legend(fontsize=FONT_LEG, loc='upper left')

    # ── kick-up power law ────────────────────────────────────────────────────
    a, b, x_lo, x_hi = fit_kickup(blocks, datasets, x_fn)
    x_line = np.logspace(np.log10(x_lo), np.log10(x_hi), 300)
    y_line = a * x_line**b
    ax.plot(x_line, y_line, 'k--', linewidth=2.0, zorder=5)

    # Annotate with the fitted slope; wrap x_sym in () so nested ^ is valid LaTeX
    slope_num = int(round(1.0 / b)) if abs(1.0/b - round(1.0/b)) < 0.07 else None
    if slope_num is not None:
        label_str = fr'$\propto ({x_sym})^{{1/{slope_num}}}$'
    else:
        label_str = fr'$\propto ({x_sym})^{{{b:.2f}}}$'

    ax.text(0.97, 0.05, label_str,
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=FONT_LEG + 1,
            bbox=dict(facecolor='white', alpha=0.85,
                      edgecolor='lightgray', boxstyle='round,pad=0.3'))

    print(f"  {xlabel}: a = {a:.4e},  slope = {b:.4f}")


def main():
    blocks = load_dat_blocks('tavg_bflux.dat')
    datasets = [
        (0, 'Steady (300,30)',   ORANGE,    'o'),
        (1, 'Steady (600,30)',   SKY_BLUE,  'o'),
        (2, 'Steady (600,60)',   DARK_BLUE, 'o'),
        (3, 'Steady (1000,10)',  VERMILLION,'o'),
        (4, 'Steady (1000,100)', GREEN,     'o'),
        (5, 'Stoch (600,60)',    DARK_BLUE, '^'),
        (6, 'Stoch (1000,100)', GREEN,     '^'),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    fig.subplots_adjust(wspace=0.08)

    print("Kick-up power-law fits "
          f"(VTurb in [{FIT_VTURB_MIN}, {FIT_VTURB_MAX}]):")

    plot_panel(axes[0], blocks, datasets,
               compute_x_Reb, compute_xerr_Reb,
               xlabel=r'$Re_B^*$', x_sym=r'Re_B^*',
               show_ylabel=True, show_legend=True)

    plot_panel(axes[1], blocks, datasets,
               compute_x_ReG, compute_xerr_ReG,
               xlabel=r'$Re_G$', x_sym=r'Re_G',
               show_ylabel=False, show_legend=False)

    fig.savefig('VTurb_combined.pdf', bbox_inches='tight')
    fig.savefig('VTurb_combined.png', dpi=150, bbox_inches='tight')
    print("Saved VTurb_combined.pdf and VTurb_combined.png")


if __name__ == '__main__':
    main()
