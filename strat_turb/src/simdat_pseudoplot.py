import dbuhlMod as db
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap
from netCDF4 import MFDataset

# Example CVD-safe hex colors (Blue to Orange)
clors = ["#0072B2", "#E69F00"] 
custom_cmap = LinearSegmentedColormap.from_list("cb_friendly", clors)

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"]
plt.rcParams["mathtext.fontset"] = 'stix'

FONT_LABEL = 24
FONT_TICK  = 18
FONT_LEG   = 14

num_contours = 500

# ── open simdat files ────────────────────────────────────────────────────────
fn = 'simdat*.cdf'
cdf_file = MFDataset(fn)

print(cdf_file.variables)

# grid
x = np.array(cdf_file.variables["x"])
y = np.array(cdf_file.variables["y"])
z = np.array(cdf_file.variables["z"])
t_last = float(cdf_file.variables["t"][-1])
Nx, Ny, Nz = len(x), len(y), len(z)
print('t={:.4f}, Nx={}, Ny={}, Nz={}'.format(t_last, Nx, Ny, Nz))

gx = float(cdf_file.variables["Gammax"][()])
gy = float(cdf_file.variables["Gammay"][()])
gz = float(cdf_file.variables["Gammaz"][()])
Fr = 1./np.sqrt(float(cdf_file.variables["B_therm"][()]))
dx, dy, dz = gx/Nx, gy/Ny, gz/Nz

z_phys = np.linspace(0, gz, Nz)

# ── load fields (last timestep only) ─────────────────────────────────────────
ux = np.array(cdf_file.variables["ux"][-1:])
uy = np.array(cdf_file.variables["uy"][-1:])
uz = np.array(cdf_file.variables["uz"][-1:])

u_rms = np.sqrt(np.mean(ux**2 + uy**2 + uz**2))
uh_rms = np.sqrt(np.mean(ux**2 + uy**2))

# ── vorticity on each face (slice first, differentiate in-plane only) ─────────
# top face  (XY at z = gz, index -1): ω_z = ∂uy/∂x − ∂ux/∂y   → (Ny, Nx)
wz_top = (db.FD6X(uy[:, 0:1, :, :], Nx, dx)
        - db.FD6Y(ux[:, 0:1, :, :], Ny, dy))[0, 0]
#print('Shape of wz_top: ', np.shape(wz_top))
#print('Mean of wz_top: ', np.mean(wz_top))
#print('STD of wz_top: ', np.std(wz_top))

# front face (XZ at y = 0,  index  0): ω_y = ∂ux/∂z − ∂uz/∂x  → (Nz, Nx)
uy_front = np.concatenate((uy[:,:,0:5,:],uy[:,:,Ny-5:,:]),axis=2)
ux_front = np.concatenate((ux[:,:,0:5,:],ux[:,:,Ny-5:,:]),axis=2)
uy_side = np.concatenate((uy[:,:,:,0:5],uy[:,:,:,Nx-5:]),axis=3)
ux_side = np.concatenate((ux[:,:,:,0:5],ux[:,:,:,Nx-5:]),axis=3)

del ux
del uy

wy_front = (db.FD6X(uy_front, Nx, dx)
          - db.FD6Y(ux_front, 8, dy))[0, :, 0]

# side face  (YZ at x = 0,  index  0): ω_x = ∂uz/∂y − ∂uy/∂z  → (Nz, Ny)
wx_side = (db.FD6X(uy_side, 8, dx)
         - db.FD6Y(ux_side, Ny, dy))[0, :, :, 0]

# ── uz face slices ────────────────────────────────────────────────────────────
uz_top   = uz[0, -1, :, :]
uz_front = uz[0, :,  0, :]
uz_side  = uz[0, :,  :, 0]
uzmax = np.abs(uz).max() / 2
del uz

# colour limits
wmax = max(np.abs(wz_top).max(), np.abs(wy_front).max(), np.abs(wx_side).max())
wmax = np.sqrt(uh_rms/Fr)

# ── 2D meshgrids for each face ───────────────────────────────────────────────
X_xy, Y_xy = np.meshgrid(x, y)        # (Ny, Nx)
X_xz, Z_xz = np.meshgrid(x, z_phys)  # (Nz, Nx)
Y_yz, Z_yz = np.meshgrid(y, z_phys)  # (Nz, Ny)

# ── bounding box edges ───────────────────────────────────────────────────────
vertices = np.array([
    [0,   0,   0 ],  # i0: bottom front left
    [gx,  0,   0 ],  # i1: bottom front right
    [gx,  gy,  0 ],  # i2: bottom back right (hidden from default view)
    [0,   gy,  0 ],  # i3: bottom back left
    [0,   0,   gz],  # i4: top front left
    [gx,  0,   gz],  # i5: top front right
    [gx,  gy,  gz],  # i6: top back right
    [0,   gy,  gz],  # i7: top back left
])
edges = [
    [3, 0], [0, 1], [1, 5], [5, 6],
    [6, 7], [7, 4], [4, 5], [7, 3], [4, 0]
]

def draw_bbox(ax):
    for e in edges:
        v1, v2 = vertices[e[0]], vertices[e[1]]
        ax.plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]],
                color='black', linewidth=1, zorder=3)

# ── box panel plotter ─────────────────────────────────────────────────────────
def plot_bbox(ax, top_data, front_data, side_data, vmax, cmap):
    """
    Render three boundary faces of the domain using contourf + zdir.

    top   – XY plane at z = gz  (z index -1): contourf zdir='z'
    front – XZ plane at y = 0   (y index  0): contourf zdir='y'
    side  – YZ plane at x = 0   (x index  0): contourf zdir='x'
    """
    ax.computed_zorder = False
    kw = dict(levels=num_contours, vmin=-vmax, vmax=vmax, cmap=cmap, zorder=1)

    # top face: XY at z = gz
    ax.contourf(X_xy, Y_xy, top_data,
                zdir='z', offset=gz, **kw)

    # front face: XZ at y = 0  (args: X, DATA, Z)
    ax.contourf(X_xz, front_data, Z_xz,
                zdir='y', offset=0, **kw)

    # side face: YZ at x = 0  (args: DATA, Y, Z)
    ax.contourf(side_data, Y_yz, Z_yz,
                zdir='x', offset=0, **kw)

    draw_bbox(ax)

    ax.set_xlim(0, gx)
    ax.set_ylim(0, gy)
    ax.set_zlim(0, gz)
    ax.set_box_aspect((gx/gz, gy/gz, 1), zoom=0.95)
    ax.view_init(25, 240, 0)
    ax.set_xlabel(r'$x$', labelpad=-5,  fontsize=FONT_LABEL)
    ax.set_ylabel(r'$y$', labelpad=-15, fontsize=FONT_LABEL)
    ax.set_zlabel(r'$z$', labelpad=-15, fontsize=FONT_LABEL)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

# ── figure ────────────────────────────────────────────────────────────────────
CMAP_VORT = 'RdBu_r'
CMAP_UZ   = 'PuOr'

fig = plt.figure(figsize=(16, 7))
fig.patch.set_facecolor('white')
fig.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.12, wspace=0.05)

ax_v = fig.add_subplot(121, projection='3d')
ax_u = fig.add_subplot(122, projection='3d')

# vorticity box: ω_z on top, ω_y on front, ω_x on side
plot_bbox(ax_v,
          wz_top, wy_front, wx_side,
          wmax, CMAP_VORT)

# uz box
plot_bbox(ax_u,
          uz_top, uz_front, uz_side,
          uzmax, CMAP_UZ)

# ── colorbars (horizontal, below each subplot) ────────────────────────────────
cax_v = fig.add_axes([0.06, 0.06, 0.38, 0.03])
sm_v = plt.cm.ScalarMappable(cmap=CMAP_VORT,
                              norm=colors.Normalize(vmin=-wmax, vmax=wmax))
cb_v = fig.colorbar(sm_v, cax=cax_v, orientation='horizontal')
cb_v.ax.tick_params(labelsize=FONT_TICK)
cb_v.set_label(r'$\omega_z$', fontsize=FONT_LABEL)

cax_u = fig.add_axes([0.56, 0.06, 0.38, 0.03])
sm_uz = plt.cm.ScalarMappable(cmap=CMAP_UZ,
                               norm=colors.Normalize(vmin=-uzmax, vmax=uzmax))
cb_u = fig.colorbar(sm_uz, cax=cax_u, orientation='horizontal')
cb_u.ax.tick_params(labelsize=FONT_TICK)
cb_u.set_label(r'$w$', fontsize=FONT_LABEL)

fig.tight_layout()

plt.savefig('simdat_3dbbox.png', dpi=800, bbox_inches='tight')
plt.savefig('simdat_3dbbox.pdf', bbox_inches='tight')
