import dbuhlMod as db
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from netCDF4 import MFDataset

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'stix'
plt.rcParams["font.size"] = 16

def make_norm(maxnum):
    return colors.Normalize(vmin=-maxnum, vmax=maxnum)

ptstp = -1

# open simdat files
fn = 'simdat*.cdf'
cdf_file = MFDataset(fn)

print(cdf_file.variables)

# grid
x = np.array(cdf_file.variables["x"])
y = np.array(cdf_file.variables["y"])
z = np.array(cdf_file.variables["z"])
t = np.array(cdf_file.variables["t"][:])
Nx = len(x)
Ny = len(y)
Nz = len(z)
Nt = len(t)
print('Nt={}, Nx={}, Ny={}, Nz={}'.format(Nt, Nx, Ny, Nz))

gx = 4 * np.pi
gy = 4 * np.pi
gz = np.pi
dx = gx / Nx
dy = gy / Ny
dz = gz / Nz
dt = t[1] - t[0]

# load velocity fields — indexed [t, z, y, x]
ux = np.array(cdf_file.variables["ux"][:])
uy = np.array(cdf_file.variables["uy"][:])
uz = np.array(cdf_file.variables["uz"][:])

# vorticity components (sign convention matches db.FD6X/FD6Y: returns -deriv)
# wz = -d(uy)/dx + d(ux)/dy  (i.e. -omega_z in standard notation)
# wy = -d(ux)/dz + d(uz)/dx
# wx = -d(uz)/dy + d(uy)/dz
wz = db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
wy = db.FD6Z(ux, Nz, dz)    - db.FD6X(uz, Nx, dx)
wx = db.FD6Y(uz, Ny, dy) - db.FD6Z(uy, Nz, dz)

del ux, uy

# colour limits
wzmax = np.abs(wz).max()
wymax = np.abs(wy).max()
wxmax = np.abs(wx).max()
uzmax = np.abs(uz).max() / 2

cmap = 'RdBu_r'

# face slices:
#   xy top face  (z = top, index -1) : field[ptstp, -1, :, :]
#   xz front face (y = 0)            : field[ptstp, :,  0, :]
#   yz side face  (x = 0)            : field[ptstp, :,  :, 0]

fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 1])

# --- row 1: vorticity ---
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

norm_wz = colors.Normalize(vmin=-wzmax, vmax=wzmax)
norm_wy = colors.Normalize(vmin=-wymax, vmax=wymax)
norm_wx = colors.Normalize(vmin=-wxmax, vmax=wxmax)

im_wz = ax1.imshow(wz[ptstp, -1, :, :], norm=norm_wz,
                   cmap=cmap, origin='lower', aspect='auto')
fig.colorbar(im_wz, ax=ax1, fraction=0.046, pad=0.04)
ax1.set_title(r'$\omega_z$ (top, $z=0$)')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$', rotation=0, labelpad=10)
ax1.set_xticks([])
ax1.set_yticks([])

im_wy = ax2.imshow(wy[ptstp, :, 0, :], norm=norm_wy,
                   cmap=cmap, origin='lower', aspect='auto')
fig.colorbar(im_wy, ax=ax2, fraction=0.046, pad=0.04)
ax2.set_title(r'$\omega_y$ (front, $y=0$)')
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$z$', rotation=0, labelpad=10)
ax2.set_xticks([])
ax2.set_yticks([])

im_wx = ax3.imshow(wx[ptstp, :, :, 0], norm=norm_wx,
                   cmap=cmap, origin='lower', aspect='auto')
fig.colorbar(im_wx, ax=ax3, fraction=0.046, pad=0.04)
ax3.set_title(r'$\omega_x$ (side, $x=0$)')
ax3.set_xlabel(r'$y$')
ax3.set_ylabel(r'$z$', rotation=0, labelpad=10)
ax3.set_xticks([])
ax3.set_yticks([])

# --- row 2: vertical velocity ---
ax4 = fig.add_subplot(gs[1, 0])
ax5 = fig.add_subplot(gs[1, 1])
ax6 = fig.add_subplot(gs[1, 2])

norm_uz = colors.Normalize(vmin=-uzmax, vmax=uzmax)

im_uz1 = ax4.imshow(uz[ptstp, -1, :, :], norm=norm_uz,
                    cmap=cmap, origin='lower', aspect='auto')
fig.colorbar(im_uz1, ax=ax4, fraction=0.046, pad=0.04)
ax4.set_title(r'$u_z$ (top, $z=0$)')
ax4.set_xlabel(r'$x$')
ax4.set_ylabel(r'$y$', rotation=0, labelpad=10)
ax4.set_xticks([])
ax4.set_yticks([])

im_uz2 = ax5.imshow(uz[ptstp, :, 0, :], norm=norm_uz,
                    cmap=cmap, origin='lower', aspect='auto')
fig.colorbar(im_uz2, ax=ax5, fraction=0.046, pad=0.04)
ax5.set_title(r'$u_z$ (front, $y=0$)')
ax5.set_xlabel(r'$x$')
ax5.set_ylabel(r'$z$', rotation=0, labelpad=10)
ax5.set_xticks([])
ax5.set_yticks([])

im_uz3 = ax6.imshow(uz[ptstp, :, :, 0], norm=norm_uz,
                    cmap=cmap, origin='lower', aspect='auto')
fig.colorbar(im_uz3, ax=ax6, fraction=0.046, pad=0.04)
ax6.set_title(r'$u_z$ (side, $x=0$)')
ax6.set_xlabel(r'$y$')
ax6.set_ylabel(r'$z$', rotation=0, labelpad=10)
ax6.set_xticks([])
ax6.set_yticks([])

fig.tight_layout()
plt.savefig('simdat_vort_pseudoplot.png', dpi=800, bbox_inches='tight')
plt.savefig('simdat_vort_pseudoplot.pdf', bbox_inches='tight')
