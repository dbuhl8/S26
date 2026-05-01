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
FONT_LEG   = 14   # panel annotation text boxes
FONT_KEY   = 17   # figure-level dataset legend


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


def valid_mask(x, xerr, y, yerr):
    return (np.isfinite(x) & np.isfinite(xerr) & np.isfinite(y) & np.isfinite(yerr)
            & (x > 0) & (y > 0))


def fit_amplitude(x, y, slope):
    """Fix slope, fit amplitude in log-log space."""
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    return np.exp(np.mean(np.log(y[mask]) - slope * np.log(x[mask])))


def fit_powerlaw(x, y):
    """Free power-law fit (both slope and amplitude) in log-log space."""
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    slope, logA = np.polyfit(np.log10(x[mask]), np.log10(y[mask]), 1)
    return 10**logA, slope


# ── Fr* = sqrt(B) / uh ──────────────────────────────────────────────────────

def compute_x_Fr(block):
    return np.sqrt(block[:, 2]) / block[:, 29]


def compute_xerr_Fr(block):
    B, uh, uh_err = block[:, 2], block[:, 29], block[:, 30]
    return np.sqrt(B) * uh_err / uh**2


# ── Fr_turb = mdisp * uh / (Re * sqrt(B) * U_tot) ───────────────────────────

def _FrT(block):
    Re, B = block[:, 0], block[:, 2]
    ux, uy, uz = block[:, 5], block[:, 7], block[:, 9]
    mdisp, uh = block[:, 17], block[:, 29]
    U_tot = np.sqrt(ux**2 + uy**2 + uz**2)
    return mdisp * uh / (Re * np.sqrt(B) * U_tot)


def compute_x_FrT(block):
    return 1.0 / _FrT(block)


def compute_xerr_FrT(block):
    ux, ux_err = block[:, 5],  block[:, 6]
    uy, uy_err = block[:, 7],  block[:, 8]
    uz, uz_err = block[:, 9],  block[:, 10]
    mdisp, mdisp_err = block[:, 17], block[:, 18]
    uh,    uh_err    = block[:, 29], block[:, 30]
    U2  = ux**2 + uy**2 + uz**2
    FrT = _FrT(block)
    err = FrT * np.sqrt(
        (mdisp_err / mdisp)**2 + (uh_err / uh)**2 +
        (ux * ux_err / U2)**2  + (uy * uy_err / U2)**2 + (uz * uz_err / U2)**2
    )
    return err / FrT**2          # σ(1/FrT)


# ── Print fit comparison table ───────────────────────────────────────────────

def print_fit_comparison(fit_results, label):
    row_names = ['Total (free fit)', 'Turb  (fixed slope)', 'Lam   (fixed slope)']
    print(f"\n{'='*68}")
    print(f"  Fit comparison — {label}")
    print(f"{'='*68}")
    print(f"  {'':22s}  {'wrms':>22s}  {'brms':>22s}")
    print(f"  {'':22s}  {'amplitude':>11s} {'slope':>9s}  {'amplitude':>11s} {'slope':>9s}")
    print(f"  {'-'*64}")
    for row, name in enumerate(row_names):
        w = fit_results.get((row, 0))
        b = fit_results.get((row, 1))
        ws = f"  {w[0]:.4e}  {w[1]:+.3f}" if w else f"  {'N/A':>11s}  {'N/A':>7s}"
        bs = f"  {b[0]:.4e}  {b[1]:+.3f}" if b else f"  {'N/A':>11s}  {'N/A':>7s}"
        print(f"  {name:22s}{ws}{bs}")
    print(f"{'='*68}\n")


# ── Figure builder ───────────────────────────────────────────────────────────

def make_combined_figure(blocks, datasets, wrms_panels, brms_panels,
                         compute_x_fn, compute_xerr_fn, xlabel, x_sym, filename):
    """
    3 rows × 2 cols:
      col 0 = wrms   col 1 = brms
      row 0 = total  row 1 = turb  row 2 = lam

    panel spec: (ycol, yerrcol, ylabel, slope, slope_label)
      slope = float  → fix slope, fit amplitude only
      slope = 'free' → fit both slope and amplitude
      slope = None   → no fit line
    """
    fig, axes = plt.subplots(3, 2, figsize=(13, 16), sharex=True,
                             constrained_layout=False)
    fig.subplots_adjust(hspace=0.08, wspace=0.32,
                        left=0.10, right=0.97, top=0.91, bottom=0.13)

    legend_handles = []
    legend_labels  = []
    fit_results = {}   # (row, col) -> (amplitude, slope)

    for row in range(3):
        for col, panels in enumerate([wrms_panels, brms_panels]):
            ycol, yerrcol, ylabel, slope, slope_label = panels[row]
            ax = axes[row, col]

            all_x, all_y = [], []
            all_y_lo, all_y_hi = [], []
            for blk_idx, label, color, marker in datasets:
                block = blocks[blk_idx]
                x    = compute_x_fn(block)
                xerr = compute_xerr_fn(block)
                y    = block[:, ycol]
                yerr = block[:, yerrcol]
                mask = valid_mask(x, xerr, y, yerr)
                if not mask.any():
                    continue
                cont = ax.errorbar(
                    x[mask], y[mask], xerr=xerr[mask], yerr=yerr[mask],
                    fmt=marker, color=color, markersize=7,
                    capsize=3, linestyle='none', label=label,
                    elinewidth=1.2, capthick=1.2
                )
                if row == 0 and col == 0:
                    legend_handles.append(cont)
                    legend_labels.append(label)
                all_x.append(x[mask])
                all_y.append(y[mask])
                all_y_lo.append(y[mask] - yerr[mask])
                all_y_hi.append(y[mask] + yerr[mask])

            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.tick_params(labelsize=FONT_TICK, which='both')
            ax.set_ylabel(ylabel, fontsize=FONT_LABEL)

            if slope is not None and all_x:
                x_cat = np.concatenate(all_x)
                y_cat = np.concatenate(all_y)

                if slope == 'free':
                    a, use_slope = fit_powerlaw(x_cat, y_cat)
                    annot = fr'$\propto ({x_sym})^{{{use_slope:.2f}}}$'
                else:
                    a = fit_amplitude(x_cat, y_cat, slope)
                    use_slope = slope
                    annot = slope_label

                fit_results[(row, col)] = (a, use_slope)
                x_line = np.logspace(np.log10(x_cat.min()),
                                     np.log10(x_cat.max()), 300)
                ax.plot(x_line, a * x_line**use_slope,
                        color='black', linestyle='--', linewidth=2.0)
                ax.text(0.97, 0.97, annot,
                        transform=ax.transAxes, ha='right', va='top',
                        fontsize=FONT_LEG,
                        bbox=dict(facecolor='white', alpha=0.85,
                                  edgecolor='lightgray', boxstyle='round,pad=0.3'))

            if all_y_lo and all_y_hi:
                y_lo_cat = np.concatenate(all_y_lo)
                y_hi_cat = np.concatenate(all_y_hi)
                pos_lo = y_lo_cat[y_lo_cat > 0]
                log_pad = 0.15
                log_ymin = np.log10(pos_lo.min() if len(pos_lo) else y_hi_cat.min()) - log_pad
                log_ymax = np.log10(y_hi_cat.max()) + log_pad
                ax.set_ylim(10**log_ymin, 10**log_ymax)

            if row == 2:
                ax.set_xlabel(xlabel, fontsize=FONT_LABEL)
            else:
                ax.tick_params(labelbottom=False)

    axes[0, 0].set_title(r'$w_\mathrm{rms}$', fontsize=FONT_LABEL + 2, pad=6)
    axes[0, 1].set_title(r'$b_\mathrm{rms}$', fontsize=FONT_LABEL + 2, pad=6)

    fig.legend(legend_handles, legend_labels,
               loc='lower center', ncol=4,
               bbox_to_anchor=(0.5, 0.01),
               fontsize=FONT_KEY,
               framealpha=0.95, edgecolor='gray',
               handlelength=1.5, handletextpad=0.5, columnspacing=1.0)

    fig.savefig(filename + '.pdf', bbox_inches='tight')
    fig.savefig(filename + '.png', dpi=150, bbox_inches='tight')
    print(f"Saved {filename}.pdf and {filename}.png")

    print_fit_comparison(fit_results, filename)


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

    # ── Fr* panels ──────────────────────────────────────────────────────────
    wrms_Fr = [
        (9,  10, r'$w_\mathrm{rms}$',          'free', None),
        (21, 22, r'$w_\mathrm{rms,\,turb}$',   -0.5,  r'$\propto (Fr^*)^{-1/2}$'),
        (23, 24, r'$w_\mathrm{rms,\,lam}$',    -1.0,  r'$\propto (Fr^*)^{-1}$'),
    ]
    brms_Fr = [
        (13, 14, r'$b_\mathrm{rms}$',          'free', None),
        (25, 26, r'$b_\mathrm{rms,\,turb}$',   -1.5,  r'$\propto (Fr^*)^{-3/2}$'),
        (27, 28, r'$b_\mathrm{rms,\,lam}$',    -1.0,  r'$\propto (Fr^*)^{-1}$'),
    ]

    make_combined_figure(
        blocks, datasets, wrms_Fr, brms_Fr,
        compute_x_Fr, compute_xerr_Fr,
        r'$(Fr^*)^{-1}$', r'Fr^*{}^{-1}',
        'brms_wrms_combined_Fr'
    )

    # ── Fr_turb panels ───────────────────────────────────────────────────────
    wrms_FrT = [
        (9,  10, r'$w_\mathrm{rms}$',          'free', None),
        (21, 22, r'$w_\mathrm{rms,\,turb}$',   -0.5,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-1/2}$'),
        (23, 24, r'$w_\mathrm{rms,\,lam}$',    -1.0,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-1}$'),
    ]
    brms_FrT = [
        (13, 14, r'$b_\mathrm{rms}$',          'free', None),
        (25, 26, r'$b_\mathrm{rms,\,turb}$',   -1.5,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-3/2}$'),
        (27, 28, r'$b_\mathrm{rms,\,lam}$',    -1.0,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-1}$'),
    ]

    make_combined_figure(
        blocks, datasets, wrms_FrT, brms_FrT,
        compute_x_FrT, compute_xerr_FrT,
        r'$Fr_\mathrm{turb}^{-1}$', r'Fr_\mathrm{turb}^{-1}',
        'brms_wrms_combined_turbFr'
    )


if __name__ == '__main__':
    main()
