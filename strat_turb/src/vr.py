import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import dbuhlMod as db
from netCDF4 import MFDataset

fn = 'simdat2.cdf'
cdf_file = MFDataset(fn)

x = np.array(cdf_file.variables["x"])
y = np.array(cdf_file.variables["y"])
z = np.array(cdf_file.variables["z"])
t = np.array(cdf_file.variables["t"][:])
Nx = len(x)
Ny = len(y)
Nz = len(z)
Nt = len(t)
print('Number of Timesteps: ', Nt)
del x
del y
del z
gx = cdf_file.variables["Gammax"][0]
gy = cdf_file.variables["Gammay"][0]
gz = cdf_file.variables["Gammaz"][0]
dx = gx/Nx
dy = gy/Ny
dz = gz/Nz
Re = 1./cdf_file.variables["D_visc"][0]
Pe = 1./cdf_file.variables["D_therm"][0]
B = cdf_file.variables["B_therm"][0]
Fr = 1./np.sqrt(B)
print("Parameters: Re - ", Re, ", Pe - ", Pe, ', B - ', B)

ptstp = -1
ux =  np.array(cdf_file.variables["ux"][:])
ux_inst =  ux[ptstp,:,:,:]
del ux
uy =  np.array(cdf_file.variables["uy"][:])
uy_inst =  uy[ptstp,:,:,:]
del uy
uz =  np.array(cdf_file.variables["uz"][:])
uz_inst =  uz[ptstp,:,:,:]
del uz
temp =  np.array(cdf_file.variables["Temp"][:])
temp_inst =  temp[ptstp,:,:,:]
del temp

cdf_file.close()

dTdz = db.iFD6Z(temp_inst,Nz,dz)
tdisp = db.iFD6Z(temp_inst,Nz,dz)**2 + db.iFD6Y(temp_inst,Ny,dy)**2 +\
    db.iFD6X(temp_inst,Nx,dx)**2
mdisp = db.iFD6Z(ux_inst,Nz,dz)**2 + db.iFD6Y(ux_inst,Ny,dy)**2 +\
    db.iFD6X(ux_inst,Nx,dx)**2 + db.iFD6Z(uy_inst,Nz,dz)**2 + db.iFD6Y(uy_inst,Ny,dy)**2 +\
    db.iFD6X(uy_inst,Nx,dx)**2 +db.iFD6Z(uz_inst,Nz,dz)**2 + db.iFD6Y(uz_inst,Ny,dy)**2 +\
    db.iFD6X(uz_inst,Nx,dx)**2
Ri = (1 + dTdz)/((Fr**2)*np.minimum((db.iFD6Z(ux_inst,Nz,dz)**2 + db.iFD6Z(uy_inst,Nz,dz)**2),1e-4)) 
eta = (B*tdisp/Pe)/(B*tdisp/Pe + mdisp/Re)
vortz = db.iFD6X(uy_inst,Nx,dx) - db.iFD6Y(ux_inst,Ny,dy)

idx_lam = np.where(vortz**2 <= 1./Fr)
idx_turb = np.where(vortz**2 > 1./Fr)

# computing max
temp_max = np.abs(temp_inst).max()
dtemp_max = np.abs(dTdz).max()
tdisp_max = tdisp.max()
mdisp_max = mdisp.max()
vortz_max = 2*np.abs(ux_inst).max()

print('Avg TDisp: ',tdisp.mean())
print('Avg MDisp: ',mdisp.mean())
print('Avg Eta (local): ',eta.mean())
print('Avg Eta (global): ', (B*tdisp.mean()/Pe)/\
    (B*tdisp.mean()/Pe + mdisp.mean()/Re))
print('Avg Lam Eta (local): ', eta[idx_lam].mean(), ' +/- ', eta[idx_lam].std())
print('Avg Lam Eta (global): ', (B*tdisp[idx_lam].mean()/Pe)/\
    (B*tdisp[idx_lam].mean()/Pe + mdisp[idx_lam].mean()/Re))
print('Avg Lam TDisp: ',tdisp[idx_lam].mean())
print('Avg Lam MDisp: ',mdisp[idx_lam].mean())
print('Avg Turb Eta (local): ', eta[idx_turb].mean(), ' +/- ', eta[idx_turb].std())
print('Avg Turb Eta (global): ', (B*tdisp[idx_turb].mean()/Pe)/\
    (B*tdisp[idx_turb].mean()/Pe + mdisp[idx_turb].mean()/Re))
print('Avg Turb TDisp: ',tdisp[idx_turb].mean())
print('Avg Turb MDisp: ',mdisp[idx_turb].mean())
print('VLam: ', len(idx_lam[0])/(Nx*Ny*Nz))
print('VTurb: ', len(idx_turb[0])/(Nx*Ny*Nz))

# transposing 
vortz = np.transpose(vortz, axes=(2,1,0))
dTdz = np.transpose(dTdz, axes=(2,1,0))
Ri = np.transpose(Ri, axes=(2,1,0))
ux_inst = np.transpose(ux_inst,axes=(2,1,0))
temp_inst = np.transpose(temp_inst,axes=(2,1,0))
eta = np.transpose(eta,axes=(2,1,0))
mdisp = np.transpose(mdisp,axes=(2,1,0))
tdisp = np.transpose(tdisp,axes=(2,1,0))

pl = pv.Plotter(border=False)

#ux plot
#opacity = [.3,0,0,0,.3]
#pl.add_volume(ux_inst,clim=[-vortz_max/2,vortz_max/2],cmap='RdYlBu_r',opacity=opacity)
#print('UX max :', vortz_max/2)

#vort plot
#opacity = [.3,0,0,0,.3]
#pl.add_volume(vortz,clim=[-vortz_max,vortz_max],cmap='RdYlBu_r',opacity=opacity)
#print('Vortz max :', vortz_max)

# temp plot
#opacity = [.3,0,0.05,0,.3]
#pl.add_volume(temp_inst,clim=[-temp_max,temp_max],cmap='RdYlBu_r',opacity=opacity)

#Ri Plot
#opacity = [.1,.1,.3,0,0]
#pl.add_volume(Ri,clim=[-0.75,1.25],cmap='RdYlBu_r',opacity=opacity)

# mdisp Plot
#opacity = [0,0,.1,.3,.3]
#pl.add_volume(mdisp/Re,clim=[0,mdisp_max/(2*Re)],cmap='viridis_r',opacity=opacity)
#print('mdisp max: ', mdisp_max/(2*Re))

# tdisp Plot
opacity = [0,0,.1,.3,.3]
pl.add_volume(tdisp,clim=[0,tdisp_max/2],cmap='viridis_r',opacity=opacity)
#print('tdisp max: ', tdisp_max/2)

# eta Plot
#opacity = [.01,.01,.01,0,.4]
#pl.add_volume(eta,clim=[0,1],cmap='RdYlBu_r',opacity=opacity)

pl.add_bounding_box()
#pl.view_xz()
pl.set_viewup([0,1,0])
pl.show()

