import dbuhlMod as db
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import MFDataset

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'stix'
plt.rcParams["font.size"] = 16

ptstp = -1

# ── open simdat files ────────────────────────────────────────────────────────
fn = 'simdat*.cdf'
cdf_file = MFDataset(fn)

print(cdf_file.variables)

# grid
x = np.array(cdf_file.variables["x"])
y = np.array(cdf_file.variables["y"])
z = np.array(cdf_file.variables["z"])
t = np.array(cdf_file.variables["t"][:])
Nx, Ny, Nz = len(x), len(y), len(z)
Nt = len(t)
print('Nt={}, Nx={}, Ny={}, Nz={}'.format(Nt, Nx, Ny, Nz))

gx, gy, gz = 4*np.pi, 4*np.pi, np.pi
dx, dy, dz = gx/Nx, gy/Ny, gz/Nz

# ── load fields ──────────────────────────────────────────────────────────────
ux = np.array(cdf_file.variables["ux"][:])
uy = np.array(cdf_file.variables["uy"][:])
uz = np.array(cdf_file.variables["uz"][:])

# vorticity components (sign convention matches db.FD6X/FD6Y: returns -deriv)
wz = db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
wy = db.FD6Z(ux, Nz, dz) - db.FD6X(uz, Nx, dx)
wx = db.FD6Y(uz, Ny, dy) - db.FD6Z(uy, Nz, dz)

u_rms = np.sqrt(np.mean(ux**2 + uy**2 + uz**2))
del ux, uy

# colour limits
# vorticity clim: U_rms / gz gives the physically motivated shear scale;
# vort_clim_factor lets you brighten (<1) or darken (>1) the vorticity colours
vort_clim_factor = 1.0
wmax  = vort_clim_factor * u_rms / gz
uzmax = np.abs(uz).max() / 2

cmap_obj = plt.get_cmap('RdBu_r')

# ── helpers ──────────────────────────────────────────────────────────────────
def face_rgba(data, vmax):
    norm = colors.Normalize(vmin=-vmax, vmax=vmax)
    return cmap_obj(norm(data))

def clean_ax(ax):
    ax.set_axis_off()

def plot_bbox(ax, top_data, front_data, side_data, vmax_top, vmax_front, vmax_side):
    """
    Render three bounding-box faces on a 3-D axes (full resolution).

    Faces:
      top   – XY plane at z = gz  (physical top of domain, z index -1)
      front – XZ plane at y = 0   (index  0 along y)
      side  – YZ plane at x = 0   (index  0 along x)
    """
    # z_wall extends to gz so the wall faces are flush with the top face
    z_wall = np.linspace(0, gz, Nz)

    # front face ── XZ plane at y = 0
    Xf, Zf = np.meshgrid(x, z_wall)
    Yf = np.zeros_like(Xf)
    ax.plot_surface(Xf, Yf, Zf,
                    facecolors=face_rgba(front_data, vmax_front),
                    shade=False, rcount=Nz, ccount=Nx)

    # side face ─── YZ plane at x = 0
    Ys, Zs = np.meshgrid(y, z_wall)
    Xs = np.zeros_like(Ys)
    ax.plot_surface(Xs, Ys, Zs,
                    facecolors=face_rgba(side_data, vmax_side),
                    shade=False, rcount=Nz, ccount=Ny)

    # top face ─── XY plane at z = gz (plotted last so painter's algorithm keeps it on top)
    Xt, Yt = np.meshgrid(x, y)
    Zt = gz * np.ones_like(Xt)
    ax.plot_surface(Xt, Yt, Zt,
                    facecolors=face_rgba(top_data, vmax_top),
                    shade=False, rcount=Ny, ccount=Nx)

    ax.set_xlim(0, gx)
    ax.set_ylim(0, gy)
    ax.set_zlim(0, gz)
    ax.set_box_aspect((4, 4, 1))   # match domain proportions: 4π × 4π × π
    ax.view_init(elev=35, azim=225)
    clean_ax(ax)

# ── single figure: vorticity (left) and uz (right) ───────────────────────────
fig = plt.figure(figsize=(18, 8))
fig.patch.set_facecolor('white')

ax_v = fig.add_subplot(121, projection='3d')
ax_u = fig.add_subplot(122, projection='3d')

# vorticity bounding box — all faces share the same colour scale
plot_bbox(ax_v,
          wz[ptstp, -1, :, :],   # top  (z = gz): ω_z
          wy[ptstp, :,  0, :],   # front (y = 0): ω_y
          wx[ptstp, :,  :, 0],   # side  (x = 0): ω_x
          wmax, wmax, wmax)
ax_v.set_title('Vorticity', pad=10)

# uz bounding box
plot_bbox(ax_u,
          uz[ptstp, -1, :, :],
          uz[ptstp, :,  0, :],
          uz[ptstp, :,  :, 0],
          uzmax, uzmax, uzmax)
ax_u.set_title(r'$u_z$', pad=10)

# ── colorbars ─────────────────────────────────────────────────────────────────
# vorticity: single shared colorbar (between the two subplots)
cax_v = fig.add_axes([0.50, 0.20, 0.012, 0.55])
sm_v = plt.cm.ScalarMappable(cmap=cmap_obj,
                              norm=colors.Normalize(vmin=-wmax, vmax=wmax))
fig.colorbar(sm_v, cax=cax_v, label=r'$\omega$')

# uz: single colorbar to the right of the right subplot
cax_u = fig.add_axes([0.92, 0.20, 0.012, 0.55])
sm_uz = plt.cm.ScalarMappable(cmap=cmap_obj,
                               norm=colors.Normalize(vmin=-uzmax, vmax=uzmax))
fig.colorbar(sm_uz, cax=cax_u, label=r'$u_z$')

plt.savefig('simdat_3dbbox.png', dpi=800, bbox_inches='tight')
plt.savefig('simdat_3dbbox.pdf', bbox_inches='tight')
