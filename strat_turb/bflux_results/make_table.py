"""Generate sim_table.tex — APJ-style deluxetable* of simulation diagnostics."""
import numpy as np

def load_dat_blocks(filepath):
    blocks, current_rows = [], []
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


blocks = load_dat_blocks('tavg_bflux.dat')

# (block_idx, display_type, nominal_Re, nominal_Pe)
SIM_META = [
    (0, 'Steady', 300,   30),
    (1, 'Steady', 600,   30),
    (2, 'Steady', 600,   60),
    (3, 'Steady', 1000,  10),
    (4, 'Steady', 1000, 100),
    (5, 'Stoch.', 600,   60),
    (6, 'Stoch.', 1000, 100),
]

def fmt(val, dec):
    if not np.isfinite(val):
        return r'\ldots'
    return f'{val:.{dec}f}'

lines = []

lines += [
    r'\begin{deluxetable*}{llccrrrrrr}',
    r'\tablecaption{Time-averaged diagnostics for all simulations.'
    r' \label{tab:sims}}',
    r'\tablewidth{0pt}',
    r'\tabletypesize{\small}',
    r'\tablehead{',
    r'  \colhead{Type} &',
    r'  \colhead{$(Re,\,Pe)$} &',
    r'  \colhead{$Pr$} &',
    r'  \colhead{$Fr^{*}$} &',
    r'  \colhead{$u_{h,\mathrm{rms}}$} &',
    r'  \colhead{$w_\mathrm{rms}$} &',
    r'  \colhead{$b_\mathrm{rms}$} &',
    r'  \colhead{$V_\mathrm{Turb}$} &',
    r'  \colhead{$\mathcal{M}_\mathrm{disp}$} &',
    r'  \colhead{$\mathcal{T}_\mathrm{disp}$}',
    r'}',
    r'\startdata',
]

first_group = True
for idx, sim_type, Re_nom, Pe_nom in SIM_META:
    block = blocks[idx]
    valid_rows = []
    for row in block:
        Re, Pe, B = row[0], row[1], row[2]
        uh,  uz   = row[29], row[9]
        brms      = row[13]
        vt        = row[19]
        mdisp     = row[17]
        tdisp     = row[15]
        if not (np.isfinite(Re) and Re > 0 and
                np.isfinite(uh) and uh > 0 and
                np.isfinite(B)  and B  > 0):
            continue
        Fr = uh / np.sqrt(B)
        Pr = Pe / Re
        valid_rows.append((Pr, Fr, uh, uz, brms, vt, mdisp, tdisp))

    if not valid_rows:
        continue

    if not first_group:
        lines.append(r'\tableline')
    first_group = False

    run_label = f'$({Re_nom},{Pe_nom})$'
    for i, (Pr, Fr, uh, uz, brms, vt, mdisp, tdisp) in enumerate(valid_rows):
        type_col = sim_type if i == 0 else ''
        run_col  = run_label  if i == 0 else ''
        lines.append(
            f'  {type_col} & {run_col} & '
            f'{fmt(Pr,3)} & {fmt(Fr,3)} & '
            f'{fmt(uh,4)} & {fmt(uz,4)} & '
            f'{fmt(brms,4)} & {fmt(vt,3)} & '
            f'{fmt(mdisp,1)} & {fmt(tdisp,3)} \\\\'
        )

lines += [
    r'\enddata',
    r'\tablecomments{$Fr^* \equiv u_{h,\mathrm{rms}}/\sqrt{B}$ where $B$ is the'
    r' dimensionless buoyancy parameter. $Pr = Pe/Re$.'
    r' $\mathcal{M}_\mathrm{disp}$ and $\mathcal{T}_\mathrm{disp}$ are the'
    r' momentum and thermal displacement diagnostics (see text).'
    r' Rows with $V_\mathrm{Turb} = 0$ are laminar states.}',
    r'\end{deluxetable*}',
]

tex = '\n'.join(lines) + '\n'

with open('sim_table.tex', 'w') as f:
    f.write(tex)

print("Written sim_table.tex")
print(f"Total data rows: {sum(1 for l in lines if l.strip().endswith('\\\\'))}")
