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

# VTurb thresholds for the two fit regions
FIT_VTURB_MIN = 0.01
FIT_VTURB_MID = 0.65   # boundary between kick-up and taper regions

FIT_COLOR_KICKUP = 'k'
FIT_COLOR_TAPER  = '#AA3377'   # muted rose (Paul Tol qualitative palette)


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


# ── logarithmic fit over a VTurb sub-range ───────────────────────────────────

def fit_region(blocks, datasets, x_fn, vturb_lo, vturb_hi):
    """Fit VTurb = a*ln(x) + b for points with VTurb in [vturb_lo, vturb_hi].
    Returns (a, b, x_min, x_max) or None if fewer than 2 points qualify."""
    xs, ys = [], []
    for idx, *_ in datasets:
        block = blocks[idx]
        x = x_fn(block)
        y = block[:, 19]
        mask = (np.isfinite(x) & np.isfinite(y) & (x > 0)
                & (y >= vturb_lo) & (y <= vturb_hi))
        xs.append(x[mask])
        ys.append(y[mask])
    xs = np.concatenate(xs)
    ys = np.concatenate(ys)
    if len(xs) < 2:
        return None
    a, b = np.polyfit(np.log(xs), ys, 1)
    return a, b, xs.min(), xs.max()


# ── panel plotter ─────────────────────────────────────────────────────────────

def plot_panel(ax, blocks, datasets, x_fn, x_err_fn,
               xlabel, x_sym, show_ylabel, kickup_vmin=None):

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

    # ── kick-up fit ──────────────────────────────────────────────────────────
    vmin_k = kickup_vmin if kickup_vmin is not None else FIT_VTURB_MIN
    res_k = fit_region(blocks, datasets, x_fn, vmin_k, FIT_VTURB_MID)
    res_t = fit_region(blocks, datasets, x_fn, FIT_VTURB_MID, 1.05)

    # Determine non-overlapping x extents
    if res_k is not None and res_t is not None:
        a_k, b_k, x_lo_k, x_hi_k = res_k
        a_t, b_t, x_lo_t, x_hi_t = res_t
        if x_lo_t >= x_hi_k:
            x_end_k, x_start_t = x_hi_k, x_lo_t
        else:
            # Ranges overlap — split at the log-midpoint
            x_split = 10 ** ((np.log10(x_lo_t) + np.log10(x_hi_k)) / 2)
            x_end_k, x_start_t = x_split, x_split
    elif res_k is not None:
        a_k, b_k, x_lo_k, x_hi_k = res_k
        x_end_k = x_hi_k
    elif res_t is not None:
        a_t, b_t, x_lo_t, x_hi_t = res_t
        x_start_t = x_lo_t

    fit_handles, fit_labels = [], []

    if res_k is not None:
        x_line = np.logspace(np.log10(x_lo_k), np.log10(x_end_k), 300)
        h, = ax.plot(x_line, a_k * np.log(x_line) + b_k,
                     '--', color=FIT_COLOR_KICKUP, linewidth=2.0, zorder=5,
                     label='_nolegend_')
        fit_handles.append(h)
        fit_labels.append(fr'$\propto {a_k:.2f}\ln({x_sym})$')
        print(f"  kick-up {xlabel}: a = {a_k:.4f},  b = {b_k:.4f}")

    if res_t is not None:
        x_line = np.logspace(np.log10(x_start_t), np.log10(x_hi_t), 300)
        h, = ax.plot(x_line, a_t * np.log(x_line) + b_t,
                     '--', color=FIT_COLOR_TAPER, linewidth=2.0, zorder=5,
                     label='_nolegend_')
        fit_handles.append(h)
        fit_labels.append(fr'$\propto {a_t:.2f}\ln({x_sym})$')
        print(f"  taper   {xlabel}: a = {a_t:.4f},  b = {b_t:.4f}")

    if res_k is not None and res_t is not None:
        ReG_crit = np.exp((b_t - b_k) / (a_k - a_t))
        print(f"  critical {xlabel}: {xlabel} = {ReG_crit:.4g}  (kick-up / taper intersection)")

    if fit_handles:
        ax.legend(handles=fit_handles, labels=fit_labels,
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
        (5, 'Stoch (600,60)',    DARK_BLUE, '^'),
        (6, 'Stoch (1000,100)', GREEN,     '^'),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    fig.subplots_adjust(wspace=0.08, bottom=0.30)

    print(f"Logarithmic fits (kick-up: VTurb in [{FIT_VTURB_MIN}, {FIT_VTURB_MID}],"
          f" taper: VTurb in [{FIT_VTURB_MID}, 1.05]):")

    plot_panel(axes[0], blocks, datasets,
               compute_x_Reb, compute_xerr_Reb,
               xlabel=r'$Re_B^*$', x_sym=r'Re_B^*',
               show_ylabel=True, kickup_vmin=0.05)

    plot_panel(axes[1], blocks, datasets,
               compute_x_ReG, compute_xerr_ReG,
               xlabel=r'$Re_G$', x_sym=r'Re_G',
               show_ylabel=False)

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
