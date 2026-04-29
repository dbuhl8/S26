import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

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


def compute_Fr_eff_inv(block):
    """(Fr*)^{-1} = sqrt(B) / uh_rms"""
    B  = block[:, 2]   # col 3
    uh = block[:, 29]  # col 30
    return np.sqrt(B) / uh


def recon_brms(X, c_turb, c_lam):
    """b_pred = sqrt(VTurb*(c_turb * FrT_inv^(3/2))^2 + VLam*(c_lam * FrT_inv)^2)"""
    FrT_inv, VTurb = X
    VLam = 1.0 - VTurb
    return np.sqrt(VTurb * (c_turb * FrT_inv**1.5)**2
                 + VLam  * (c_lam  * FrT_inv)**2)


def recon_wrms(X, c_turb, c_lam):
    """w_pred = sqrt(VTurb*(c_turb * FrT_inv^(1/2))^2 + VLam*(c_lam * FrT_inv)^2)"""
    FrT_inv, VTurb = X
    VLam = 1.0 - VTurb
    return np.sqrt(VTurb * (c_turb * FrT_inv**0.5)**2
                 + VLam  * (c_lam  * FrT_inv)**2)


def valid_mask(FrT_inv, VTurb, y):
    return (np.isfinite(FrT_inv) & np.isfinite(VTurb) & np.isfinite(y)
            & (FrT_inv > 0) & (y > 0) & (VTurb >= 0) & (VTurb <= 1))


def collect_xy(blocks, datasets, ycol):
    """Pool (Fr*^{-1}, y) pairs across all datasets; discard non-finite/non-positive."""
    xs, ys = [], []
    for idx, *_ in datasets:
        block = blocks[idx]
        x = compute_Fr_eff_inv(block)
        y = block[:, ycol]
        m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        xs.append(x[m]); ys.append(y[m])
    return np.concatenate(xs), np.concatenate(ys)


def fit_amplitude(x, y, slope):
    """Fit amplitude a for  y = a * x^slope  in log space (slope fixed)."""
    return np.exp(np.mean(np.log(y) - slope * np.log(x)))


def make_parity_panel(ax, blocks, datasets,
                      ycol_total, ycol_turb, slope_turb, ycol_lam, slope_lam,
                      recon_func, actual_label, pred_label, quantity_name):
    # Fit each coefficient from its own decomposed component with fixed slope.
    # c_turb from turb_component = c_turb * Fr^slope_turb
    # c_lam  from lam_component  = c_lam  * Fr^slope_lam
    x_t, y_t = collect_xy(blocks, datasets, ycol_turb)
    x_l, y_l = collect_xy(blocks, datasets, ycol_lam)
    c_turb = fit_amplitude(x_t, y_t, slope_turb)
    c_lam  = fit_amplitude(x_l, y_l, slope_lam)
    print(f"{quantity_name}: c_turb = {c_turb:.4e},  c_lam = {c_lam:.4e}")

    all_actual, all_pred = [], []
    has_open = False

    for idx, label, color, marker in datasets:
        block   = blocks[idx]
        FrT_inv = compute_Fr_eff_inv(block)
        VTurb   = block[:, 19]
        y       = block[:, ycol_total]
        ReG     = block[:, 17] / block[:, 2]   # mdisp / B
        m = valid_mask(FrT_inv, VTurb, y)
        if not m.any():
            continue

        y_pred  = recon_func((FrT_inv[m], VTurb[m]), c_turb, c_lam)
        ReG_m   = ReG[m]
        hi      = ReG_m >= 1   # filled
        lo      = ReG_m <  1   # open

        # Re_G >= 1: filled markers (carry the dataset label)
        if hi.any():
            ax.scatter(y_pred[hi], y[m][hi], color=color, marker=marker,
                       s=60, zorder=3, label=label)

        # Re_G < 1: open (hollow) markers, same color/shape, no extra label
        if lo.any():
            ax.scatter(y_pred[lo], y[m][lo],
                       facecolors='none', edgecolors=color, marker=marker,
                       s=60, zorder=3, linewidths=1.5)
            has_open = True

        all_actual.append(y[m])
        all_pred.append(y_pred)

    # 1:1 reference line over the full data range
    combined = np.concatenate(all_actual + all_pred)
    vmin, vmax = combined.min() * 0.8, combined.max() * 1.2
    ref = np.logspace(np.log10(vmin), np.log10(vmax), 200)
    ax.plot(ref, ref, 'k--', linewidth=1.5, label='1:1', zorder=2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(pred_label,   fontsize=20)
    ax.set_ylabel(actual_label, fontsize=20)
    ax.tick_params(labelsize=15)
    ax.set_aspect('equal', adjustable='datalim')

    return c_turb, c_lam, has_open


def main():
    blocks   = load_dat_blocks('tavg_bflux.dat')
    datasets = [
        (0, 'Steady (300,30)',   ORANGE,    'o'),
        (1, 'Steady (600,30)',   SKY_BLUE,  'o'),
        (2, 'Steady (600,60)',   DARK_BLUE, 'o'),
        (3, 'Steady (1000,10)',  VERMILLION,'o'),
        (4, 'Steady (1000,100)', GREEN,     'o'),
        (5, 'Stoch (600,60)',    DARK_BLUE, '^'),
        (6, 'Stoch (1000,100)', GREEN,     '^'),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    _, _, has_open_b = make_parity_panel(
        axes[0], blocks, datasets,
        ycol_total=13, ycol_turb=25, slope_turb=1.5,
                       ycol_lam=27,  slope_lam=1.0,
        recon_func=recon_brms,
        actual_label=r'$b_\mathrm{rms}$',
        pred_label=r'$b_\mathrm{pred}$',
        quantity_name='brms',
    )
    make_parity_panel(
        axes[1], blocks, datasets,
        ycol_total=9,  ycol_turb=21, slope_turb=0.5,
                       ycol_lam=23,  slope_lam=1.0,
        recon_func=recon_wrms,
        actual_label=r'$w_\mathrm{rms}$',
        pred_label=r'$w_\mathrm{pred}$',
        quantity_name='wrms',
    )

    # legend on first panel only; append open-marker proxy if needed
    handles, labels = axes[0].get_legend_handles_labels()
    if has_open_b:
        proxy = Line2D([0], [0], marker='o', color='dimgray',
                       markerfacecolor='none', markeredgewidth=1.5,
                       markersize=8, linestyle='none',
                       label=r'open: $Re_G < 1$')
        handles.append(proxy)
        labels.append(r'open: $Re_G < 1$')
    axes[0].legend(handles, labels, fontsize=11, loc='upper left')

    fig.tight_layout()
    fig.savefig('reconstruction_plot.pdf')
    fig.savefig('reconstruction_plot.png', dpi=150)
    print("Saved reconstruction_plot.pdf and reconstruction_plot.png")


if __name__ == '__main__':
    main()
