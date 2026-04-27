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


def compute_Fr_eff_inv_err(block):
    """σ((Fr*)^{-1}) = sqrt(B) * uh_err / uh^2"""
    B      = block[:, 2]   # col 3
    uh     = block[:, 29]  # col 30
    uh_err = block[:, 30]  # col 31
    return np.sqrt(B) * uh_err / uh**2


def compute_bflux(block):
    """Buoyancy flux = tdisp / Pe"""
    return block[:, 15] / block[:, 1]   # col 16 / col 2


def compute_bflux_err(block):
    """σ = tdisp_err / Pe"""
    return block[:, 16] / block[:, 1]   # col 17 / col 2


def fit_amplitude(x, y, slope):
    """Fit amplitude a for y = a * x^slope in log space (same method as brms_wrms_plot.py)."""
    mask = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
    return np.exp(np.mean(np.log(y[mask]) - slope * np.log(x[mask])))


def collect_component(blocks, datasets, ycol, slope):
    """Pool all valid (Fr_eff_inv, y) and return the fitted amplitude."""
    all_x, all_y = [], []
    for idx, *_ in datasets:
        block = blocks[idx]
        x    = compute_Fr_eff_inv(block)
        xerr = compute_Fr_eff_inv_err(block)
        y    = block[:, ycol]
        # Use the same valid mask as brms_wrms_plot.py (requires finite x, xerr, y > 0)
        mask = (np.isfinite(x) & np.isfinite(xerr) & np.isfinite(y)
                & (x > 0) & (y > 0))
        all_x.append(x[mask])
        all_y.append(y[mask])
    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)
    return fit_amplitude(all_x, all_y, slope)


def valid_mask_bflux(Fr_eff_inv, VTurb, y):
    return (np.isfinite(Fr_eff_inv) & np.isfinite(VTurb) & np.isfinite(y)
            & (Fr_eff_inv > 0) & (y > 0) & (VTurb >= 0) & (VTurb <= 1))


def predict_bflux(Fr_eff_inv, VTurb, C_turb, C_lam):
    """flux_pred = VTurb*C_turb*(Fr*_inv)^2 + VLam*C_lam*(Fr*_inv)^2"""
    VLam = 1.0 - VTurb
    return VTurb * C_turb * Fr_eff_inv**(-2) + VLam * C_lam * Fr_eff_inv**(-2)


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

    # Fit amplitudes from the individual component plots, matching brms_wrms_plot.py exactly
    # Columns (0-based): turb_brms=25, lam_brms=27, turb_wrms=21, lam_wrms=23
    a_b_turb = collect_component(blocks, datasets, ycol=25, slope=-1.5)
    a_b_lam  = collect_component(blocks, datasets, ycol=27, slope=-1.0)
    a_w_turb = collect_component(blocks, datasets, ycol=21, slope=-0.5)
    a_w_lam  = collect_component(blocks, datasets, ycol=23, slope=-1.0)

    C_turb = a_w_turb * a_b_turb
    C_lam  = a_w_lam  * a_b_lam

    print(f"a_b_turb = {a_b_turb:.4e},  a_b_lam = {a_b_lam:.4e}")
    print(f"a_w_turb = {a_w_turb:.4e},  a_w_lam = {a_w_lam:.4e}")
    print(f"C_turb   = {C_turb:.4e},    C_lam   = {C_lam:.4e}")

    fig, ax = plt.subplots(figsize=(7, 6))

    all_pred, all_actual = [], []

    for idx, label, color, marker in datasets:
        block      = blocks[idx]
        Fr_eff_inv = compute_Fr_eff_inv(block)
        VTurb      = block[:, 19]       # vturb (col 20)
        y          = compute_bflux(block)
        y_err      = compute_bflux_err(block)
        m = valid_mask_bflux(Fr_eff_inv, VTurb, y)
        if not m.any():
            continue

        y_pred = predict_bflux(Fr_eff_inv[m], VTurb[m], C_turb, C_lam)

        ax.errorbar(y[m], y_pred, yerr=y_err[m],
                    fmt=marker, color=color, markersize=7,
                    capsize=3, linestyle='none', label=label, zorder=3)

        all_pred.append(y_pred)
        all_actual.append(y[m])

    # 1:1 reference line over the union of pred and actual ranges
    combined = np.concatenate(all_pred + all_actual)
    combined = combined[combined > 0]
    ref = np.logspace(np.log10(combined.min() * 0.8),
                      np.log10(combined.max() * 1.2), 200)
    ax.plot(ref, ref, 'k--', linewidth=1.5, label='1:1', zorder=2)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\langle wb \rangle$', fontsize=20)
    ax.set_ylabel(r'$\langle wb \rangle_\mathrm{pred}$', fontsize=20)
    ax.legend(fontsize=12, loc='upper left')
    ax.tick_params(labelsize=15)
    ax.set_aspect('equal', adjustable='datalim')

    fig.tight_layout()
    fig.savefig('bflux_reconstruction.pdf')
    fig.savefig('bflux_reconstruction.png', dpi=150)
    print("Saved bflux_reconstruction.pdf and bflux_reconstruction.png")


if __name__ == '__main__':
    main()
