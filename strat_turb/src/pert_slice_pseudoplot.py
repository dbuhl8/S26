import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from netCDF4 import MFDataset 

# this file is meant to work for a simdat.cdf file. 
# different configurations exsit for slice files which are 2D in nature. 
def FDX(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dx_field[:,:,0] = (field[:,:,1]-field[:,:,0])/dx
    dx_field[:,:,nx-1] = (field[:,:,nx-1]-field[:,:,nx-2])/dx
    # centered finite difference formula (inside)
    for i in range(nx-2):
        dx_field[:,:,i+1] = (0.5/dx)*(field[:,:,i+2]-field[:,:,i])
    return dx_field

def FDY(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative # with respect to x_comp
    dy_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    dy_field[:,0,:] = (field[:,1,:]-field[:,0,:])/dy
    dy_field[:,ny-1,:] = (field[:,ny-1,:]-field[:,ny-2,:])/dy
    # centered finite difference formula (inside)
    for i in range(ny-2):
        dy_field[:,i+1,:] = (0.5/dy)*(field[:,i+2,:]-field[:,i,:])
    return dy_field

def make_norm(maxnum):
    return colors.Normalize(vmin=-maxnum,vmax=maxnum)

def sclmap(cmap, maxnum):
    return mpl.cm.ScalarMappable(cmap=cmap,norm=make_norm(maxnum))

# opening slice data files
fnxy = 'XYSLICE*.cdf'
fnxz = 'XZSLICE*.cdf'
fnyz = 'YZSLICE*.cdf'
cdf_filexy = MFDataset(fnxy)
cdf_filexz = MFDataset(fnxz)
cdf_fileyz = MFDataset(fnyz)
dtype = np.float32
cmap = 'RdYlBu_r'
num_contours = 100 # enforces smoothness

x = np.array(cdf_filexy.variables["x"])
y = np.array(cdf_filexy.variables["y"])
z = np.array(cdf_filexz.variables["z"])
t = np.array(cdf_filexy.variables["t"][:])
Nx = len(x) # using MFDataset (these might be extra long)
Ny = len(y)
Nz = len(z)
Nt = len(t)
gx = 4*np.pi
gy = 2*np.pi
gz = np.pi
dx = gx/Nx
dy = gy/Ny
dz = gz/Nz
dt = t[1]-t[0]

#these arrays are indexed by [t,z,y,x]
ux_xy =  np.array(cdf_filexy.variables["ux"][:])
ux_xz =  np.array(cdf_filexz.variables["ux"][:])
ux_yz =  np.array(cdf_fileyz.variables["ux"][:])

uy_xy =  np.array(cdf_filexy.variables["uy"][:])
uy_xz =  np.array(cdf_filexz.variables["uy"][:])
uy_yz =  np.array(cdf_fileyz.variables["uy"][:])

uz_xy =  np.array(cdf_filexy.variables["uz"][:])
uz_xz =  np.array(cdf_filexz.variables["uz"][:])
uz_yz =  np.array(cdf_fileyz.variables["uz"][:])

temp_xy =  np.array(cdf_filexy.variables["Temp"][:])
temp_xz =  np.array(cdf_filexz.variables["Temp"][:])
temp_yz =  np.array(cdf_fileyz.variables["Temp"][:])

# Finding max values for the colorbar
uxmax = np.max([ux_xz.max(), ux_yz.max(), ux_xy.max()])
uymax = np.max([uy_xz.max(), uy_yz.max(), ux_xy.max()])
uzmax = np.max([uz_xz.max(), uz_yz.max(), uz_xy.max()])
tempmax = np.max([temp_xz.max(), temp_yz.max(), temp_xy.max()])

#print('Maximum; ux: ', uxmax, ', uy: ',\
    #uymax, ', uz: ', uzmax, ', temp: ',tempmax)

kwux = {
    'vmin': -uxmax,
    'vmax': uxmax,
    #'norm': make_norm(uxmax),
    'cmap': cmap
}

kwuy = {
    'vmin': -uymax,
    'vmax': uymax,
    #'norm': make_norm(uymax),
    'cmap': cmap
}

kwuz = {
    'vmin': -uzmax,
    'vmax': uzmax,
    #'norm': make_norm(tempmax),
    'cmap': cmap
}

kwtemp = {
    'vmin': -tempmax,
    'vmax': tempmax,
    #'norm': make_norm(uzmax),
    'cmap': cmap
}

# plotting with matplotlib
XX, YY, ZZ = np.meshgrid(x,y, -z)

"""
    # can loop over this to create a movie through the vertical domain
    pc1 = ax[0,0].pcolor(XX[:,:],
        YY[:,:], ux[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-uxmax,vmax=uxmax), cmap='seismic')
    fig.colorbar(pc1, ax=ax[0,0])
    ax[0,0].set_title(r'$u_x$')

    pc2 = ax[0,1].pcolor(XX[:,:],
        YY[:,:], uy[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-uymax,vmax=uymax), cmap='seismic')
    fig.colorbar(pc2, ax=ax[0,1])
    ax[0,1].set_title(r'$u_y$')

    pc3 = ax[1,0].pcolor(XX[:,:],
        YY[:,:], wz[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-wzmax,vmax=wzmax), cmap='seismic')
    fig.colorbar(pc3, ax=ax[1,0])
    ax[1,0].set_title(r'$\omega_z$')

    pc4 = ax[1,1].pcolor(XX[:,:],
        YY[:,:], uz[-1,0,:,:].T,
        norm=colors.Normalize(vmin=-uzmax,vmax=uzmax), cmap='seismic')
    fig.colorbar(pc4, ax=ax[1,1])
    ax[1,1].set_title(r'$u_z$')
    ------------------------------------------------------------------------------  
    # plot using imshow (much faster than pcolor for evenly spaced discretization)
    # axis style is worse however
    pc1 = ax[0,0].imshow(ux[0,:,:].T,
        norm=colors.Normalize(vmin=-uxmax,vmax=uxmax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc1, ax=ax[0,0])
    ax[0,0].set_title(r'$u_x$')

    pc2 = ax[0,1].imshow(uy[0,:,:].T,
        norm=colors.Normalize(vmin=-uymax,vmax=uymax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc2, ax=ax[0,1])
    ax[0,1].set_title(r'$u_y$')

    pc3 = ax[1,0].imshow(wz[0,:,:].T,
        norm=colors.Normalize(vmin=-wzmax,vmax=wzmax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc3, ax=ax[1,0])
    ax[1,0].set_title(r'$\omega_z$')

    pc4 = ax[1,1].imshow(uz[0,:,:].T,
        norm=colors.Normalize(vmin=-uzmax,vmax=uzmax), cmap='seismic',
        origin='lower')
    fig.colorbar(pc4, ax=ax[1,1])
    ax[1,1].set_title(r'$u_z$')
    -----------------------------------------------------------------------------
    need to implement this for each subplot
    import matplotlib.pyplot as plt
    import numpy as np

    # Define dimensions
    Nx, Ny, Nz = 100, 300, 500
    X, Y, Z = np.meshgrid(np.arange(Nx), np.arange(Ny), -np.arange(Nz))

    # Create fake data
    data = (((X+100)**2 + (Y-20)**2 + 2*Z)/1000+1)

    kw = {
        'vmin': data.min(),
        'vmax': data.max(),
        'levels': np.linspace(data.min(), data.max(), 10),
    }

    # Create a figure with 3D ax
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111, projection='3d')

    # Plot contour surfaces
    _ = ax.contourf(
        X[:, :, 0], Y[:, :, 0], data[:, :, 0],
        zdir='z', offset=0, **kw
    )
    _ = ax.contourf(
        X[0, :, :], data[0, :, :], Z[0, :, :],
        zdir='y', offset=0, **kw
    )
    C = ax.contourf(
        data[:, -1, :], Y[:, -1, :], Z[:, -1, :],
        zdir='x', offset=X.max(), **kw
    )
    # --


    # Set limits of the plot from coord limits
    xmin, xmax = X.min(), X.max()
    ymin, ymax = Y.min(), Y.max()
    zmin, zmax = Z.min(), Z.max()
    ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

    # Plot edges
    edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
    ax.plot([xmax, xmax], [ymin, ymax], 0, **edges_kw)
    ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
    ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

    # Set labels and zticks
    ax.set(
        xlabel='X [km]',
        ylabel='Y [km]',
        zlabel='Z [m]',
        zticks=[0, -150, -300, -450],
    )

    # Set zoom and angle view
    ax.view_init(40, -30, 0)
    ax.set_box_aspect(None, zoom=0.9)

    # Colorbar
    fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Name [units]')

    # Show Figure
    plt.show()
"""

def plot_frame(frame):
    # makes a box plot for each subplot
    fig, ax = plt.subplots(2, 2,subplot_kw=dict(projection='3d'))

    # ux box plot
    ux_top = ax[0,0].contourf(XX[:,:,0],YY[:,:,0],ux_xy[frame,:,:], levels=num_contours,
        zdir='z', offset=0, **kwux)
    ux_xzside = ax[0,0].contourf(XX[0,:,:],ux_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwux)
    ux_yzside = ax[0,0].contourf(ux_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwux)
    ax[0,0].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
    ax[0,0].set_title(r"$u_x$")
    ax[0,0].view_init(40, 240, 0)
    ax[0,0].set_box_aspect((4, 4, 1), zoom=0.9)
    fig.colorbar(sclmap(cmap,uxmax), ax=ax[0,0])

    # uy box plot
    uy_top = ax[0,1].contourf(XX[:,:,0],YY[:,:,0],uy_xy[frame,:,:], levels=num_contours,
        zdir='z', offset=0, **kwuy)
    uy_xzside = ax[0,1].contourf(XX[0,:,:],uy_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuy)
    uy_yzside = ax[0,1].contourf(uy_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuy)
    ax[0,1].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
    ax[0,1].set_title(r"$u_y$")
    ax[0,1].view_init(40, 240, 0)
    ax[0,1].set_box_aspect((4, 4, 1), zoom=0.9)
    #fig.colorbar(uy_top, ax=ax[0,1], norm=make_norm(uymax), fraction=0.02, pad=0.1)
    fig.colorbar(sclmap(cmap,uymax), ax=ax[0,1])


    # temp box plot
    temp_top = ax[1,0].contourf(XX[:,:,0],YY[:,:,0],temp_xy[frame,:,:],
        levels=num_contours, zdir='z', offset=0, **kwtemp)
    temp_xzside = ax[1,0].contourf(XX[0,:,:],temp_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwtemp)
    temp_yzside = ax[1,0].contourf(temp_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwtemp)
    ax[1,0].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
    ax[1,0].set_title(r"$T'$")
    ax[1,0].view_init(40, 240, 0)
    ax[1,0].set_box_aspect((4, 4, 1), zoom=0.9)
    #fig.colorbar(temp_top, ax=ax[1,0], norm=make_norm(tempmax), fraction=0.02, pad=0.1)
    fig.colorbar(sclmap(cmap,tempmax), ax=ax[1,0])


    # uz box plot
    uz_top = ax[1,1].contourf(XX[:,:,0],YY[:,:,0],uz_xy[frame,:,:], levels=num_contours,
        zdir='z', offset=0, **kwuz)
    uz_xzside = ax[1,1].contourf(XX[0,:,:],uz_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuz)
    uz_yzside = ax[1,1].contourf(uz_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuz)
    ax[1,1].set(xlim=[XX.min(),XX.max()],ylim=[YY.min(), YY.max()], zlim=[ZZ.min(), ZZ.max()])
    ax[1,1].set_title(r"$u_z$")
    ax[1,1].view_init(40, 240, 0)
    ax[1,1].set_box_aspect((4, 4, 1), zoom=0.9)
    #fig.colorbar(uz_top, ax=ax[1,1], norm=make_norm(uzmax), fraction=0.02, pad=0.1)
    fig.colorbar(sclmap(cmap,uzmax), ax=ax[1,1])

    # Define the bounding box coordinates
    bbox_min_x = 0
    bbox_max_x = 4*np.pi
    bbox_min_y = 0
    bbox_max_y = 4*np.pi
    bbox_min_z = -np.pi
    bbox_max_z = 0

    # Draw the bounding box
    # Define the vertices of the bounding box
    vertices = [
        [bbox_min_x, bbox_min_y, bbox_min_z],#i0: bottom front center
        [bbox_max_x, bbox_min_y, bbox_min_z],#i1: bottom right
        [bbox_max_x, bbox_max_y, bbox_min_z],#i2: bottom back center(not visible)
        [bbox_min_x, bbox_max_y, bbox_min_z],#i3: bottom left
        [bbox_min_x, bbox_min_y, bbox_max_z],#i4: top front center
        [bbox_max_x, bbox_min_y, bbox_max_z],#i5: top right
        [bbox_max_x, bbox_max_y, bbox_max_z],#i6: top back center
        [bbox_min_x, bbox_max_y, bbox_max_z] #i7: top left
    ]

    # Define the edges of the bounding box
    edges = [
        [3,0], # bottom left to bottom front center
        [0,1], # bottom front center to bottom right
        [1,5], # bottom right to top right
        [5,6], # top right to top back center
        [6,7], # top back center to top left
        [7,4], # top left to top front center
        [4,5], # top front center to top right
        [7,3], # top left to bottom left
        [4,0] # top front center to bottom front center
    ]

    # Plot the edges
    for edge in edges:
        v1 = vertices[edge[0]]
        v2 = vertices[edge[1]]
        ax[0,0].plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]],
            color='black',linewidth=1,zorder=3)
        ax[0,1].plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]],
            color='black',linewidth=1,zorder=3)
        ax[1,0].plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]],
            color='black',linewidth=1,zorder=3)
        ax[1,1].plot([v1[0], v2[0]], [v1[1], v2[1]], [v1[2], v2[2]],
            color='black',linewidth=1,zorder=3)

    # Additional Plot Formating
    # makes room for the slider
    fig.subplots_adjust(0, .05, .90, .90,
    .05,.05)
    for axis_set in ax:
        for axis in axis_set:
            axis.set_xticks([])
            axis.set_yticks([])
            axis.set_zticks([])
    # horizaontally oriented slider
    taxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
    staxis = Slider(taxis, 'Time', t[0], t[-1], valinit=t[frame], valstep=dt)
    fig.tight_layout()
    fig.savefig('rot_Re1000_Om1_B100_pseudoplot_frame{:0{}d}.png'.format(frame,5),dpi=800)
    print('Done with frame: ', frame)
    plt.close('all')
    plt.clf()
    plt.cla()
    #gc.collect()
    return True


# movie down vertical extent of domain
def update_frame(frame):
    ux_top = ax[0,0].contourf(XX[:,:,0],YY[:,:,0],
        ux_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwux)
    ux_xzside = ax[0,0].contourf(XX[0,:,:],ux_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwux)
    ux_yzside = ax[0,0].contourf(ux_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwux)
    uy_top = ax[0,1].contourf(XX[:,:,0],YY[:,:,0],
        uy_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwuy)
    uy_xzside = ax[0,1].contourf(XX[0,:,:],uy_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuy)
    uy_yzside = ax[0,1].contourf(uy_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuy)
    temp_top = ax[1,0].contourf(XX[:,:,0],YY[:,:,0],
        temp_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwtemp)
    temp_xzside = ax[1,0].contourf(XX[0,:,:],temp_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwtemp)
    temp_yzside = ax[1,0].contourf(temp_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwtemp)
    uz_top = ax[1,1].contourf(XX[:,:,0],YY[:,:,0],
        uz_xy[frame,:,:], levels=num_contours, zdir='z', offset=0, **kwuz)
    uz_xzside = ax[1,1].contourf(XX[0,:,:],uz_xz[frame,:,:].T,
        ZZ[0,:,:], levels=num_contours, zdir='y', offset=0, **kwuz)
    uz_yzside = ax[1,1].contourf(uz_yz[frame,:,:].T,YY[:,-1,:],
        ZZ[:,-1,:], levels=num_contours,zdir='x', offset=0, **kwuz)
    staxis.set_val(t[frame]) 
    print('Done with frame: ', frame)
    return (ux_top, ux_xzside, ux_yzside, uy_top, uy_xzside, uy_yzside,
        temp_top, temp_xzside, temp_yzside, uz_top, uz_xzside, uz_yzside)


for i in range(Nt):
    gate = plot_frame(i) 
    
#ani = animation.FuncAnimation(fig=fig,
    #func=update_frame,frames=Nt,interval=100,blit=True)
#ani.save('psuedoplot.gif')
#plt.show()

