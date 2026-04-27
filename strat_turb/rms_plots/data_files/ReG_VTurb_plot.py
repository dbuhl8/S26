import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times New Roman']
mpl.rcParams['mathtext.fontset'] = 'stix'  # matches Times for math

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
    # Re_G = mdisp / B
    mdisp = block[:, 17]  # col 18
    B     = block[:, 2]   # col 3
    return mdisp / B


def compute_xerr(block):
    # B is a fixed simulation parameter, so xerr = mdisp_err / B
    mdisp_err = block[:, 18]  # col 19
    B         = block[:, 2]   # col 3
    return mdisp_err / B


def valid_mask(x, xerr, y, yerr):
    return np.isfinite(x) & np.isfinite(xerr) & np.isfinite(y) & np.isfinite(yerr) & (x > 0)


def main():
    blocks = load_dat_blocks('tavg_bflux.dat')

    # (block_index, label, color, marker)
    # Steady: circle ('o'), Stochastic: triangle ('^')
    # Matching (Re,Pe) pairs share a color across steady/stoch runs
    datasets = [
        (0, 'Steady (300,30)',   ORANGE,    'o'),
        (1, 'Steady (600,30)',   SKY_BLUE,  'o'),
        (2, 'Steady (600,60)',   DARK_BLUE, 'o'),
        (3, 'Steady (1000,10)',  VERMILLION,'o'),
        (4, 'Steady (1000,100)', GREEN,     'o'),
        (5, 'Stoch (600,60)',    DARK_BLUE, '^'),
        (6, 'Stoch (1000,100)', GREEN,     '^'),
    ]

    fig, ax = plt.subplots(figsize=(9, 6))

    for idx, label, color, marker in datasets:
        block = blocks[idx]
        x    = compute_x(block)
        xerr = compute_xerr(block)
        y    = block[:, 19]   # vturb     (col 20)
        yerr = block[:, 20]   # vturb_err (col 21)

        mask = valid_mask(x, xerr, y, yerr)
        ax.errorbar(x[mask], y[mask], xerr=xerr[mask], yerr=yerr[mask],
                    fmt=marker, color=color, markersize=8,
                    capsize=3, linestyle='none', label=label)

    ax.set_xscale('log')
    ax.set_ylim(-0.05, 1)
    ax.set_xlabel(r'$Re_G$', fontsize=28)
    ax.set_ylabel(r'$V_\mathrm{Turb}$', fontsize=28)
    ax.legend(loc='upper left', fontsize=16)
    ax.tick_params(labelsize=22)

    plt.tight_layout()
    plt.savefig('ReG_VTurb_colorblind.pdf')
    plt.savefig('ReG_VTurb_colorblind.png', dpi=150)
    print("Saved ReG_VTurb_colorblind.pdf and ReG_VTurb_colorblind.png")


if __name__ == '__main__':
    main()
