import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['mathtext.fontset'] = 'stix'

# Wong (2011) colorblind-safe 8-color palette
DARK_BLUE  = '#003087'
SKY_BLUE   = '#56B4E9'
ORANGE     = '#E69F00'
VERMILLION = '#D55E00'
GREEN      = '#009E73'
PINK       = '#CC79A7'


def load_dat_blocks(filepath):
    """Parse tavg_bflux.dat into a list of numpy arrays, one per index block."""
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


def compute_FrT(block):
    """Fr_turb = mdisp * uh_rms / (Re * sqrt(B) * sqrt(ux^2 + uy^2 + uz^2))"""
    Re    = block[:, 0]   # col 1
    B     = block[:, 2]   # col 3
    ux    = block[:, 5]   # col 6
    uy    = block[:, 7]   # col 8
    uz    = block[:, 9]   # col 10
    mdisp = block[:, 17]  # col 18
    uh    = block[:, 29]  # col 30
    U_tot = np.sqrt(ux**2 + uy**2 + uz**2)
    return mdisp * uh / (Re * np.sqrt(B) * U_tot)


def compute_FrT_err(block):
    """σ_FrT = FrT * sqrt((σ_mdisp/mdisp)^2 + (σ_uh/uh)^2
                         + (ux*σ_ux/U²)^2 + (uy*σ_uy/U²)^2 + (uz*σ_uz/U²)^2)"""
    ux        = block[:, 5]   # col 6
    ux_err    = block[:, 6]   # col 7
    uy        = block[:, 7]   # col 8
    uy_err    = block[:, 8]   # col 9
    uz        = block[:, 9]   # col 10
    uz_err    = block[:, 10]  # col 11
    mdisp     = block[:, 17]  # col 18
    mdisp_err = block[:, 18]  # col 19
    uh        = block[:, 29]  # col 30
    uh_err    = block[:, 30]  # col 31
    U_tot2    = ux**2 + uy**2 + uz**2
    FrT = compute_FrT(block)
    return FrT * np.sqrt(
        (mdisp_err / mdisp)**2 +
        (uh_err    / uh)**2 +
        (ux * ux_err / U_tot2)**2 +
        (uy * uy_err / U_tot2)**2 +
        (uz * uz_err / U_tot2)**2
    )


def compute_x(block):
    """Inverse turbulent Froude number: 1 / Fr_turb"""
    return 1.0 / compute_FrT(block)


def compute_xerr(block):
    """σ(1/FrT) = σ_FrT / FrT^2"""
    FrT     = compute_FrT(block)
    FrT_err = compute_FrT_err(block)
    return FrT_err / FrT**2


def valid_mask(x, xerr, y, yerr):
    return (np.isfinite(x) & np.isfinite(xerr) & np.isfinite(y) & np.isfinite(yerr)
            & (x > 0) & (y > 0))


def fit_amplitude(x, y, slope):
    """Fit amplitude a for y = a * x^slope in log space."""
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    return np.exp(np.mean(np.log(y[mask]) - slope * np.log(x[mask])))


def plot_panel(ax, blocks, datasets, ycol, yerrcol):
    """Plot all datasets on ax; return concatenated valid (x, y) for fitting."""
    all_x, all_y = [], []
    for idx, label, color, marker in datasets:
        block = blocks[idx]
        x    = compute_x(block)
        xerr = compute_xerr(block)
        y    = block[:, ycol]
        yerr = block[:, yerrcol]
        mask = valid_mask(x, xerr, y, yerr)
        if mask.any():
            ax.errorbar(x[mask], y[mask], xerr=xerr[mask], yerr=yerr[mask],
                        fmt=marker, color=color, markersize=6,
                        capsize=3, linestyle='none', label=label)
            all_x.append(x[mask])
            all_y.append(y[mask])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(labelsize=15)
    ax.set_xlabel(r'$Fr_\mathrm{turb}^{-1}$', fontsize=20)
    return (np.concatenate(all_x) if all_x else np.array([]),
            np.concatenate(all_y) if all_y else np.array([]))


def add_power_law(ax, x_all, y_all, slope, line_label):
    a = fit_amplitude(x_all, y_all, slope)
    x_line = np.logspace(np.log10(x_all.min()), np.log10(x_all.max()), 300)
    ax.plot(x_line, a * x_line**slope,
            color='black', linestyle='--', linewidth=1.5, label=line_label)
    ax.legend(fontsize=12, loc='best')


def make_figure(blocks, datasets, panels, filename):
    """panels: list of (ycol, yerrcol, ylabel, slope_or_None, slope_label_or_None)"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for ax, (ycol, yerrcol, ylabel, slope, slope_label) in zip(axes, panels):
        x_all, y_all = plot_panel(ax, blocks, datasets, ycol, yerrcol)
        ax.set_ylabel(ylabel, fontsize=20)
        if slope is not None and len(x_all) > 0:
            add_power_law(ax, x_all, y_all, slope, slope_label)

    # dataset legend on first panel only
    axes[0].legend(fontsize=11, loc='best')

    fig.tight_layout()
    fig.savefig(filename + '.pdf')
    fig.savefig(filename + '.png', dpi=150)
    print(f"Saved {filename}.pdf and {filename}.png")


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

    # (ycol, yerrcol, ylabel, slope, slope_label)
    brms_panels = [
        (13, 14, r'$b_\mathrm{rms}$',        None,  None),
        (25, 26, r'$b_\mathrm{rms,\,turb}$', -1.5,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-3/2}$'),
        (27, 28, r'$b_\mathrm{rms,\,lam}$',  -1.0,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-1}$'),
    ]

    wrms_panels = [
        (9,  10, r'$w_\mathrm{rms}$',        None,  None),
        (21, 22, r'$w_\mathrm{rms,\,turb}$', -0.5,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-1/2}$'),
        (23, 24, r'$w_\mathrm{rms,\,lam}$',  -1.0,  r'$\propto (Fr_\mathrm{turb}^{-1})^{-1}$'),
    ]

    make_figure(blocks, datasets, brms_panels, 'brms_vs_turbFr')
    make_figure(blocks, datasets, wrms_panels, 'wrms_vs_turbFr')


if __name__ == '__main__':
    main()
