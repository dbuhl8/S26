import dbuhlMod as db
import numpy as np
import os
import fnmatch
import collections
from netCDF4 import MFDataset
from netCDF4 import Dataset

# extracted data files
steady_tavg_file = "steady_with_new_simdats/steady_tavg_eta.dat"
stoch_tavg_file = "stochastic/stoch_tavg_eta.dat"

# read in tavg files (long process) store data in arrays
# Stoch (B, Re, Pe) = (10, 600, 60) ~ sim 1a
# vs 
# Steady (B, Re, Pe) = (30, 600, 60) ~ sim 1b
#
# and
# 
# Stoch (B, Re, Pe) = (180, 600, 60) ~ sim 2a
# vs 
# Steady (B, Re, Pe) = (400, 600, 60) ~ sim 2b

print("Begin gathering time averaged data")
# getting steady data
steady_data = open(steady_tavg_file, 'r')
need_more_data = 0 
for line in steady_data.readlines():
    if "#" in line:
        continue
    elif len(line.split()) == 0:
        continue
    elif need_more_data == 2:
        break
    data = [float(x) for x in line.split()]
    Re = data[0]
    B = data[1]
    Pe = data[3]
    if  ((Re == 600 and Pe == 60) and (B == 30)):
        sim1b_Re = 600.
        sim1b_Fr = 1./np.sqrt(30.)
        sim1b_Pe = 60.
        sim1b_uhrms = data[7]
        sim1b_wrms = data[11]
        sim1b_mdisp = data[15]
        need_more_data += 1
    elif  ((Re == 600 and Pe == 60) and (B == 400)):
        sim2b_Re = 600.
        sim2b_Fr = 1./np.sqrt(400.)
        sim2b_Pe = 60.
        sim2b_uhrms = data[7]
        sim2b_wrms = data[11]
        sim2b_mdisp = data[15]
        need_more_data += 1
steady_data.close()

stoch_data = open(stoch_tavg_file, 'r')
need_more_data = 0 
for line in stoch_data.readlines():
    if "#" in line:
        continue
    elif len(line.split()) == 0:
        continue
    elif need_more_data == 2:
        break
    data = [float(x) for x in line.split()]
    Re = data[0]
    B = data[1]
    Pe = data[3]
    if  ((Re == 600 and Pe == 60) and (B == 10)):
        sim1a_Re = 600.
        sim1a_Fr = 1./np.sqrt(30.)
        sim1a_Pe = 60.
        sim1a_uhrms = data[7]
        sim1a_wrms = data[11]
        sim1a_mdisp = data[15]
        need_more_data += 1
    elif  ((Re == 600 and Pe == 60) and (B == 180)):
        sim2a_Re = 600.
        sim2a_Fr = 1./np.sqrt(400.)
        sim2a_Pe = 60.
        sim2a_uhrms = data[7]
        sim2a_wrms = data[11]
        sim2a_mdisp = data[15]
        need_more_data += 1

stoch_data.close()

print("Finished gathering time averaged data")

# open desired spectra files
sim1a_fn = 'stochastic/nonrotating/B10Re600Pe60/XY_SPEC6'
sim1b_fn = 'steady_with_new_simdats/horizontal-shear/Re600_Pe60_B30/XY_SPEC15'
sim2a_fn = 'stochastic/nonrotating/B180Re600Pe60/XY_SPEC6'
sim2b_fn = 'steady_with_new_simdats/horizontal-shear/Re600_Pe60_B400/XY_SPEC16'

# get number of modes and timesteps from spectra files 
print("Begin gathering spectra info data")
counter = 0
stoch_filler = 0 
stoch_num_modes = 0 
sim1a_kx = np.array([])
sim1a_ky = np.array([])
file = open(sim1a_fn,'r')
for line in file.readlines():
    if "#" in line:
        stoch_filler += 1
        continue
    elif len(line.split()) == 0:
        counter += 1
        continue
    elif counter == 2:
        break
    stoch_num_modes += 1
    data = [float(x) for x in line.split()]
    sim1a_kx = np.append(sim1a_kx,data[0])
    sim1a_ky = np.append(sim1a_ky,data[1])
sim1a_nt = 5
#sim1a_nt = (len(file.readlines())-stoch_filler)/\
        #(stoch_num_modes+1)
file.close()
sim1a_kx = np.unique(sim1a_kx)
sim1a_ky = np.unique(sim1a_ky)

file = open(sim2a_fn,'r')
sim2a_nt = (len(file.readlines())-stoch_filler)/\
        (stoch_num_modes+1)
file.close()
sim2a_kx = sim1a_kx
sim2a_ky = sim1a_ky


# need to make four 3-dimensional arrays for energy [time, kx, ky]
sim1a_nkx = len(sim1a_kx)
sim1a_nky = len(sim1a_ky)
sim1a_Ex = np.zeros(sim1a_nt,sim1a_nkx,sim1a_nky)
sim1a_Ey = np.zeros_like(sim1a_Ex)
sim1a_Ey = np.zeros_like(sim1a_Ex)
sim1a_Ez = np.zeros_like(sim1a_Ex)
sim1a_Et = np.zeros_like(sim1a_Ex)
sim1a_E = np.zeros_like(sim1a_Ex)
sim1a_Kh = np.zeros(sim1a_nkx,sim1a_nky)
sim2a_nkx = sim1a_nkx
sim2a_nky = sim1a_nky
sim2a_Ex = np.zeros(sim2a_nt,sim2a_nkx,sim2a_nky)
sim2a_Ey = np.zeros_like(sim2a_Ex)
sim2a_Ey = np.zeros_like(sim2a_Ex)
sim2a_Ez = np.zeros_like(sim2a_Ex)
sim2a_Et = np.zeros_like(sim2a_Ex)
sim2a_E = np.zeros_like(sim2a_Ex)
sim2a_Kh = np.zeros(sim2a_nkx,sim2a_nky)
"""
counter = 0
steady_num_modes = 0
steady_filler = 0 
file = open(sim1b_fn, 'r')
for line in file.readlines():
    if "#" in line:
        steady_filler += 1
        continue
    elif len(line.split()) == 0:
        counter += 1
        continue
    elif counter == 2:
        break
    steady_num_modes += 1
    data = [float(x) for x in line.split()]
    sim1b_kx = np.append(sim1b_kx,data[0])
    sim1b_ky = np.append(sim1b_ky,data[1])
sim1b_nt = (len(file.readlines())-steady_filler)/\
        (steady_num_modes+1)
sim1b_kx = np.unique(sim1b_kx)
sim1b_ky = np.unique(sim1b_ky)
file.close()
file = open(sim2b_fn, 'r')
sim1b_nt = (len(file.readlines())-steady_filler)/\
        (steady_num_modes+1)
sim2b_kx = sim1b_kx
sim2b_ky = sim1b_ky
file.close()

sim1b_nkx = len(sim1b_kx)
sim1b_nky = len(sim1b_ky)
sim1b_Ex = np.zeros(sim1b_nt,sim1b_nkx,sim1b_nky)
sim1b_Ey = np.zeros_like(sim1b_Ex)
sim1b_Ey = np.zeros_like(sim1b_Ex)
sim1b_Ez = np.zeros_like(sim1b_Ex)
sim1b_Et = np.zeros_like(sim1b_Ex)
sim1b_E = np.zeros_like(sim1b_Ex)
sim1b_Kh = np.zeros(sim1b_nkx,sim1b_nky)
sim2b_nkx = sim1b_nkx
sim2b_nky = sim1b_nky
sim2b_Ex = np.zeros(sim2b_nt,sim2b_nkx,sim2b_nky)
sim2b_Ey = np.zeros_like(sim2b_Ex)
sim2b_Ey = np.zeros_like(sim2b_Ex)
sim2b_Ez = np.zeros_like(sim2b_Ex)
sim2b_Et = np.zeros_like(sim2b_Ex)
sim2b_E = np.zeros_like(sim2b_Ex)
sim2b_Kh = np.zeros(sim2b_nkx,sim2b_nky)
"""


# read spectra
i = -1
print("Starting sim1a")
with open(sim1a_fn,'r') as file:
    num_lines = (stoch_num_modes+1)*5
    lines = collections.deque(num_lines,file)
    for line in lines:
        if "#" in line:
            continue
        elif len(line.split()) == 0:
            i += 1
            continue
        data = [float(x) for x in line.split()]
        kx = data[0]
        ky = data[1]
        indx = np.where(sim1a_kx == kx)[0]
        indy = np.where(sim1a_ky == ky)[0]
        sim1a_Ex[i,indx,indy] = data[2] 
        sim1a_Ey[i,indx,indy] = data[3] 
        sim1a_Ez[i,indx,indy] = data[4] 
        sim1a_E[i,indx,indy] = data[5] 
        sim1a_Et[i,indx,indy] = data[10] 
        sim1a_Kh[indx,indy] = np.sqrt(kx**2 + ky**2)
print("Finished sim1a")
"""
i = -1
print("Starting sim1b")
file = open(sim1b_fn,'r')
for line in file.readlines():
    if "#" in line:
        continue
    elif len(line.split()) == 0:
        i += 1
        continue
    data = [float(x) for x in line.split()]
    kx = data[0]
    ky = data[1]
    indx = np.where(sim1b_kx == kx)[0]
    indy = np.where(sim1b_ky == ky)[0]
    sim1b_Ex[i,indx,indy] = data[2] 
    sim1b_Ey[i,indx,indy] = data[3] 
    sim1b_Ez[i,indx,indy] = data[4] 
    sim1b_E[i,indx,indy] = data[5] 
    sim1b_Et[i,indx,indy] = data[10] 
    sim1b_Kh[indx,indy] = np.sqrt(kx**2 + ky**2)
file.close()
print("Finished sim1b")

i = -1
print("Starting sim2a")
file = open(sim2a_fn,'r')
for line in file.readlines():
    if "#" in line:
        continue
    elif len(line.split()) == 0:
        i += 1
        continue
    data = [float(x) for x in line.split()]
    kx = data[0]
    ky = data[1]
    indx = np.where(sim2a_kx == kx)[0]
    indy = np.where(sim2a_ky == ky)[0]
    sim2a_Ex[i,indx,indy] = data[2] 
    sim2a_Ey[i,indx,indy] = data[3] 
    sim2a_Ez[i,indx,indy] = data[4] 
    sim2a_E[i,indx,indy] = data[5] 
    sim2a_Et[i,indx,indy] = data[10] 
    sim2a_Kh[indx,indy] = np.sqrt(kx**2 + ky**2)
file.close()
print("Finished sim2a")

i = -1
print("Starting sim2b")
file = open(sim2b_fn,'r')
for line in file.readlines():
    if "#" in line:
        continue
    elif len(line.split()) == 0:
        i += 1
        continue
    data = [float(x) for x in line.split()]
    kx = data[0]
    ky = data[1]
    indx = np.where(sim2b_kx == kx)[0]
    indy = np.where(sim2b_ky == ky)[0]
    sim2b_Ex[i,indx,indy] = data[2] 
    sim2b_Ey[i,indx,indy] = data[3] 
    sim2b_Ez[i,indx,indy] = data[4] 
    sim2b_E[i,indx,indy] = data[5] 
    sim2b_Et[i,indx,indy] = data[10] 
    sim2b_Kh[indx,indy] = np.sqrt(kx**2 + ky**2)
file.close()
print("Finished sim2b")
"""

# this function returns each field as a function of unique kh
# repeat kh fields are averaged together (i.e. this presumes the fields are
# horizontally isotropic)
def unique_khavg(kh_array,nt,*arg):
    kh_flat = np.flatten(kh_array)
    kh_unique, idx, old_idx, counts = np.unique(kh_flat,return_index=True,\
        return_inverse=True,return_counts=True)
    nkh = len(kh_unique)
    ret_args = []
    for i,field in enumerate(arg):
        new_arg = np.zeros((nt,nkh))
        temp_arg = np.zeros((nt,len(old_idx)))
        for j in range(nt):
            temp_arg[j,:] = arg[i][j,:,:].flatten()
        for j,oidx in enumerate(old_idx):
            new_arg[j,:] += temp_arg[oidx]
        new_arg /= counts
        ret_arg.append(new_arg) 
    return [kh_array] + ret_arg

# time average the energy arrays
def tavg(nt,nt_avg,*arg):
    for i,field in enumerate(arg):
        arg[i] = np.mean(field[nt-nt_avg:nt,:],axis=0)
    return arg

kh,Ex,Ey = unique_khavg(sim1a_kh,sim1a_Ex,sim1a_Ey)
Ex,Ey = tavg(sim1a_nt,5,Ex,Ey)

# plot using data from extracted data files
fig, ax = plt.subplots()
font = {'family': 'Times New Roman',
                        'size'   : 20}
ax.xlabel(r'$k_{h}$',**font)
ax.ylabel(r'$E$',**font)
ax.plot(kh,Ex+Ey,'k-',label=r'$E_{u_h}$')

fig.tight_layout()
plt.savefig('test_spectra.png',dpi=800)





