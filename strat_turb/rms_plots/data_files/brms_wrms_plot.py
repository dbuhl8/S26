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


def compute_x(block):
    # (Fr*)^{-1} = sqrt(B) / uh_rms
    B  = block[:, 2]   # col 3
    uh = block[:, 29]  # col 30
    return np.sqrt(B) / uh


def compute_xerr(block):
    # xerr = sqrt(B) * uh_err / uh^2
    B      = block[:, 2]   # col 3
    uh     = block[:, 29]  # col 30
    uh_err = block[:, 30]  # col 31
    return np.sqrt(B) * uh_err / uh**2


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
    ax.set_xlabel(r'$(Fr^*)^{-1}$', fontsize=20)
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
        (13, 14, r'$b_\mathrm{rms}$',           None,   None),
        (25, 26, r'$b_\mathrm{rms,\,turb}$',    -1.5,   r'$\propto (Fr^*)^{-3/2}$'),
        (27, 28, r'$b_\mathrm{rms,\,lam}$',     -1.0,   r'$\propto (Fr^*)^{-1}$'),
    ]

    wrms_panels = [
        (9,  10, r'$w_\mathrm{rms}$',           None,   None),
        (21, 22, r'$w_\mathrm{rms,\,turb}$',    -0.5,   r'$\propto (Fr^*)^{-1/2}$'),
        (23, 24, r'$w_\mathrm{rms,\,lam}$',     -1.0,   r'$\propto (Fr^*)^{-1}$'),
    ]

    make_figure(blocks, datasets, brms_panels, 'brms_vs_Fr')
    make_figure(blocks, datasets, wrms_panels, 'wrms_vs_Fr')


if __name__ == '__main__':
    main()
