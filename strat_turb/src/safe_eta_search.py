import numpy as np
import dbuhlMod as db
import os
import fnmatch
from netCDF4 import MFDataset
from netCDF4 import Dataset


testing = True
stochastic = True  # True for stochastic forcing sims, False for steady sims

if stochastic:
    if testing:
        # only run this code on a single simulation set to ensure accuracy of code
        Re600_Pe60=['nonrotating/B300Re600Pe60/']
        Re600_Pe60_bounds=[[40,80]]

        simulations = [Re600_Pe60]
        bounds = [Re600_Pe60_bounds]
    else:
        Re600_Pe60=['B300Re600Pe60/','B180Re600Pe60/','B100Re600Pe60/','B30Re600Pe60/',\
                    'B10Re600Pe60/','B3Re600Pe60/','B1Re600Pe60/']
        Re1000_Pe100 = ['B100Re1000Pe100/', 'B10Re1000Pe100/']

        Re600_Pe60_bounds=[[40,80],[40,90],[45,160],\
            [35,80],[50,250],[55,110],[45,125]]
        Re1000_Pe100_bounds=[[40,70],[70,110]]

        simulations = [Re600_Pe60, Re1000_Pe100]
        bounds = [Re600_Pe60_bounds, Re1000_Pe100_bounds]

    io_filename = 'safe_nonrotating_stoch_eta.dat'
    tavg_filename = 'stoch_tavg_eta.dat'

else:  # steady forcing sims
    if testing:
        # only run this code on a single simulation set to ensure accuracy of code
        Re600_Pe60=['horizontal-shear/Re600_Pe60_B100/']
        Re600_Pe60_bounds=[[150,300]]

        simulations = [Re600_Pe60]
        bounds = [Re600_Pe60_bounds]
    else:
        Re600_Pe60=['horizontal-shear/Re600_Pe60_B100/',\
                    'horizontal-shear/Re600_Pe60_B0.1/',\
                    'horizontal-shear/Re600_Pe60_B1/',\
                    'horizontal-shear/Re600_Pe60_B10/',\
                    'horizontal-shear/Re600_Pe60_B1000/',\
                    'horizontal-shear/Re600_Pe60_B30/',\
                    'horizontal-shear/Re600_Pe60_B3000/',\
                    'horizontal-shear/Re600_Pe60_B400/',\
                    'horizontal-shear/Re600_Pe60_B6000/']
        Re600_Pe60_bounds=[[150,300],[275,400],[180,320],\
            [220,450],[500,900],[375,560],[650,900],\
            [225,460],[40,180]]

        Re600_Pe30=['horizontal-shear/Re600_Pe30_B10/',\
                    'horizontal-shear/Re600_Pe30_B100/',\
                    'horizontal-shear/Re600_Pe30_B30/']
        Re600_Pe30_bounds=[[260,300],[300,400],[400,480]]

        Re1000_Pe100=['horizontal-shear/Re1000_Pe100_B10/',\
                      'horizontal-shear/Re1000_Pe100_B100/',\
                      'horizontal-shear/Re1000_Pe100_B3/',\
                      'horizontal-shear/Re1000_Pe100_B30/',\
                      'horizontal-shear/Re1000_Pe100_B300/']
        Re1000_Pe100_bounds=[[35,65],[91,92],[85,86],[63,64],[91,92]]

        Re300_Pe30=['horizontal-shear/Re300_Pe30_B0.01/',\
                    'horizontal-shear/Re300_Pe30_B0.1/',\
                    'horizontal-shear/Re300_Pe30_B1/',\
                    'horizontal-shear/Re300_Pe30_B10/',\
                    'horizontal-shear/Re300_Pe30_B100/',\
                    'horizontal-shear/Re300_Pe30_B1000/',\
                    'horizontal-shear/Re300_Pe30_B10000/',\
                    'horizontal-shear/Re300_Pe30_B30/',\
                    'horizontal-shear/Re300_Pe30_B300/']
        Re300_Pe30_bounds=[[370,460],[260,460],[140,280],[60,160],\
            [300,1200],[800,1300],[1050,1450],[50,300],[300,900]]

        Re1000_Pe10=['horizontal-shear/Re1000_Pe10_B10/',\
                     'horizontal-shear/Re1000_Pe10_B100/',\
                     'horizontal-shear/Re1000_Pe10_B1000/',\
                     'horizontal-shear/Re1000_Pe10_B3/']
        Re1000_Pe10_bounds=[[35,40],[35,70],[94,95],[62,63]]

        simulations = [Re600_Pe60, Re600_Pe30, Re1000_Pe100, Re300_Pe30,\
            Re1000_Pe10]
        bounds = [Re600_Pe60_bounds, Re600_Pe30_bounds, Re1000_Pe100_bounds,\
            Re300_Pe30_bounds, Re1000_Pe10_bounds]

    io_filename = 'safe_nonrotating_steady_eta.dat'
    tavg_filename = 'steady_tavg_eta.dat'

# preparing output file
header_string = "#" + "{:<s}    "*36
tavg_header_string = "#" + "{:<s}    "*71
fmt_str = "{:.4e}    "*36
tavg_fmt_str = "{:.4e}    "*71
# 20 columns
header_string = header_string.format('Re','B','Pr', 'Pe', 'BPe', 't', 'uh_rms',\
    'vortz_rms', 'wrms',\
    'tdisp', 'mdisp', 'eta (local)', 'eta (global)', \
    'lam wrms', 'lam tdisp', 'lam mdisp',\
    'lam eta (local)', 'lam eta (global)', \
    'turb wrms', 'turb tdisp', 'turb mdisp',\
    'turb eta (local)', 'turb eta (global)', \
    'vlam', 'vturb',\
    'lam wrms wght', 'turb wrms wght',\
    'lam wrms eff wght', 'turb wrms eff wght',\
    'trms', 'lam trms', 'turb trms',\
    'lam trms wght', 'turb trms wght',\
    'lam trms eff wght', 'turb trms eff wght')
tavg_header_string = tavg_header_string.format('Re','B','Pr', 'Pe', 'BPe',\
    'lb','ub','uh_rms','uh_err','vortz_rms', 'vortz err','wrms','wrms err',\
    'tdisp','tdisp_err', 'mdisp','mdisp_err',\
    'eta (local)', 'eta (local) err', 'eta (global)','eta (global) err', \
    'lam wrms', 'lam wrms err',\
    'lam tdisp','lam tdisp err', 'lam mdisp','lam mdisp err',\
    'lam eta (local)', 'lam eta (local) err',\
    'lam eta (global)','lam eta (global) err',\
    'turb wrms', 'turb wrms err',\
    'turb tdisp','turb tdisp err', 'turb mdisp', 'turb mdisp err',\
    'turb eta (local)', 'turb eta (local) err', \
    'turb eta (global)', 'turb eta (global) err', \
    'vlam avg','vlam err', 'vturb avg', 'vturb err',\
    'lam wrms wght', 'lam wrms wght err', \
    'turb wrms wght', 'turb wrms wght err',\
    'lam wrms eff wght', 'lam wrms eff wght err',\
    'turb wrms eff wght', 'turb wrms eff wght err',\
    'trms', 'trms err',\
    'lam trms', 'lam trms err',\
    'turb trms', 'turb trms err',\
    'lam trms wght', 'lam trms wght err',\
    'turb trms wght', 'turb trms wght err',\
    'lam trms eff wght', 'lam trms eff wght err',\
    'turb trms eff wght', 'turb trms eff wght err')

io_file = open(io_filename,'w')
tavg_file = open(tavg_filename,'w')
io_file.write(header_string)
tavg_file.write(tavg_header_string)
io_file.write('\n')
tavg_file.write('\n')

index_counter = 0


for m, sim_set in enumerate(simulations):
    # loops over each parameter set (i.e. Re = 600, Pe = 60 is one iteration of
    # the loop)
    for n, dir_path in enumerate(sim_set):
        # loops over each B value in that parameter suite
        # the following 3 lines finds the file path of each individual simdat
        # file
        directory = os.listdir(dir_path)
        simdat_files = fnmatch.filter(directory,'simdat*.cdf')
        simdat_files = [dir_path+x for x in simdat_files]
        B = 0
        Re = 0
        Pe = 0
        t_tot = np.array([])
        uh_rms = np.array([])
        vortz_rms = np.array([])

        avg_wrms = np.array([])
        avg_lam_wrms = np.array([])
        avg_turb_wrms = np.array([])

        avg_lam_wrms_wght = np.array([])
        avg_lam_wrms_eff_wght = np.array([])
        avg_turb_wrms_wght = np.array([])
        avg_turb_wrms_eff_wght = np.array([])

        avg_trms = np.array([])
        avg_lam_trms = np.array([])
        avg_turb_trms = np.array([])

        avg_lam_trms_wght = np.array([])
        avg_lam_trms_eff_wght = np.array([])
        avg_turb_trms_wght = np.array([])
        avg_turb_trms_eff_wght = np.array([])

        avg_tdisp = np.array([])
        avg_lam_tdisp = np.array([])
        avg_turb_tdisp = np.array([])

        avg_mdisp = np.array([])
        avg_lam_mdisp = np.array([])
        avg_turb_mdisp = np.array([])

        avg_local_eta = np.array([])
        avg_local_lam_eta = np.array([])
        avg_local_turb_eta = np.array([])

        avg_global_eta = np.array([])
        avg_global_lam_eta = np.array([])
        avg_global_turb_eta = np.array([])

        vturb = np.array([])
        vlam = np.array([])
        np_tot = 0
        

        for j, fn in enumerate(simdat_files):
            # inject extraction routine here
            cdf_file = Dataset(fn,exclude='Chem')

            x = np.array(cdf_file.variables["x"])
            y = np.array(cdf_file.variables["y"])
            z = np.array(cdf_file.variables["z"])
            t = np.array(cdf_file.variables["t"])
            Nx = len(x)
            Ny = len(y)
            Nz = len(z)
            np_tot = Nx*Ny*Nz
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

            ux = np.array(cdf_file.variables["ux"])
            uy = np.array(cdf_file.variables["uy"])
            uz = np.array(cdf_file.variables["uz"])
            temp = np.array(cdf_file.variables["Temp"])

            cdf_file.close()

            #dTdz = db.FD6Z(temp_inst,Nz,dz)
            tdisp = db.FD6Z(temp,Nz,dz)**2 + db.FD6Y(temp,Ny,dy)**2 +\
                db.FD6X(temp,Nx,dx)**2
            trms = temp
            del temp

            mdisp = db.FD6Z(ux,Nz,dz)**2 + db.FD6Y(ux,Ny,dy)**2 +\
                db.FD6X(ux,Nx,dx)**2 + db.FD6Z(uy,Nz,dz)**2 +\
                db.FD6Y(uy,Ny,dy)**2 + db.FD6X(uy,Nx,dx)**2 +\
                db.FD6Z(uz,Nz,dz)**2 + db.FD6Y(uz,Ny,dy)**2 +\
                db.FD6X(uz,Nx,dx)**2
            #Ri = (1 + dTdz)/((Fr**2)*np.minimum((db.iFD6Z(ux_inst,Nz,dz)**2 +\
            #   db.iFD6Z(uy_inst,Nz,dz)**2),1e-8)) 
            eta = (B*tdisp/Pe)/(B*tdisp/Pe + mdisp/Re)
            vortz = db.FD6X(uy,Nx,dx) - db.FD6Y(ux,Ny,dy)
            uh = (ux**2+uy**2)**0.5
            wrms = uz
            del ux
            del uy
            del uz

            for i in range(Nt):
                # recording timestep and average uh
                t_tot = np.append(t_tot, t[i])
                uh_rms = np.append(uh_rms,db.rms(uh[i,:,:,:]))
                vortz_rms = np.append(vortz_rms,db.rms(vortz[i,:,:,:]))

                # eta = (thermal dissipation) / (total energy dissipation)
                # eta = (B/Pe)*tdisp / ((B/Pe)*tdisp + (1/Re)*mdisp)

                # identifying turbulent and laminar regions
                # threshold: vortz_eff = vortz/uh_rms, vortz**2 <= uh_rms/Fr
                idx_lam = np.where(vortz[i,:,:,:]**2 <= uh_rms[-1]/Fr)
                idx_turb = np.where(vortz[i,:,:,:]**2 > uh_rms[-1]/Fr)


                # computing local eta
                avg_local_eta = np.append(avg_local_eta,db.rms(eta[i,:,:,:]))
                avg_local_lam_eta = np.append(avg_local_lam_eta,\
                    db.rms(eta[i,idx_lam[0],idx_lam[1],idx_lam[2]]))
                avg_local_turb_eta = np.append(avg_local_turb_eta,\
                    db.rms(eta[i,idx_turb[0],idx_turb[1],idx_turb[2]]))

                # computing wrms (and weighted version)
                eps = 1e-8 # to protect against nans
                #vortz_rms = db.rms(vortz[i,:,:,:])
                vortz_inv_rms =\
                    np.sqrt(np.sum(1./np.maximum(vortz[i,:,:,:]**2,eps))/np_tot)

                avg_wrms = np.append(avg_wrms,db.rms(wrms[i,:,:,:]))
                avg_lam_wrms = np.append(avg_lam_wrms,\
                    db.rms(wrms[i,idx_lam[0],idx_lam[1],idx_lam[2]]))
                avg_turb_wrms = np.append(avg_turb_wrms,\
                    db.rms(wrms[i,idx_turb[0],idx_turb[1],idx_turb[2]]))

                avg_lam_wrms_wght = np.append(avg_lam_wrms_wght,\
                    np.sqrt(np.sum(wrms[i,:,:,:]**2/\
                    np.maximum(vortz[i,:,:,:]**2,eps))/np_tot)/vortz_inv_rms)
                avg_turb_wrms_wght = np.append(avg_turb_wrms_wght,\
                    db.rms(wrms[i,:,:,:]*vortz[i,:,:,:])/vortz_rms[-1])
                avg_lam_wrms_eff_wght = np.append(avg_lam_wrms_eff_wght,\
                    avg_lam_wrms_wght[-1]/uh_rms[-1])
                avg_turb_wrms_eff_wght = np.append(avg_turb_wrms_eff_wght,\
                    avg_turb_wrms_wght[-1]/uh_rms[-1])

                # computing trms (temperature rms)
                avg_trms = np.append(avg_trms,db.rms(trms[i,:,:,:]))
                avg_lam_trms = np.append(avg_lam_trms,\
                    db.rms(trms[i,idx_lam[0],idx_lam[1],idx_lam[2]]))
                avg_turb_trms = np.append(avg_turb_trms,\
                    db.rms(trms[i,idx_turb[0],idx_turb[1],idx_turb[2]]))

                avg_lam_trms_wght = np.append(avg_lam_trms_wght,\
                    np.sqrt(np.sum(trms[i,:,:,:]**2/\
                    np.maximum(vortz[i,:,:,:]**2,eps))/np_tot)/vortz_inv_rms)
                avg_turb_trms_wght = np.append(avg_turb_trms_wght,\
                    db.rms(trms[i,:,:,:]*vortz[i,:,:,:])/vortz_rms[-1])
                avg_lam_trms_eff_wght = np.append(avg_lam_trms_eff_wght,\
                    avg_lam_trms_wght[-1]/uh_rms[-1])
                avg_turb_trms_eff_wght = np.append(avg_turb_trms_eff_wght,\
                    avg_turb_trms_wght[-1]/uh_rms[-1])

                # computing avg tdisp and mdisp, NOTE that TDISP and MDisp do
                # not need to be sqrt_ed, so we use the mean() function
                avg_tdisp = np.append(avg_tdisp,db.mean(tdisp[i,:,:,:]))
                avg_lam_tdisp = np.append(avg_lam_tdisp,\
                    db.mean(tdisp[i,idx_lam[0],idx_lam[1],idx_lam[2]]))
                avg_turb_tdisp = np.append(avg_turb_tdisp,\
                    db.mean(tdisp[i,idx_turb[0],idx_turb[1],idx_turb[2]]))

                avg_mdisp = np.append(avg_mdisp,db.mean(mdisp[i,:,:,:]))
                avg_lam_mdisp = np.append(avg_lam_mdisp,\
                    db.mean(mdisp[i,idx_lam[0],idx_lam[1],idx_lam[2]]))
                avg_turb_mdisp = np.append(avg_turb_mdisp,\
                    db.mean(mdisp[i,idx_turb[0],idx_turb[1],idx_turb[2]]))

                # computing global eta
                avg_global_eta = np.append(avg_global_eta,\
                    (B*avg_tdisp[-1]/Pe)/\
                    (B*avg_tdisp[-1]/Pe + avg_mdisp[-1]/Re))
                avg_global_lam_eta = np.append(avg_global_lam_eta,\
                    (B*avg_lam_tdisp[-1]/Pe)/(B*avg_lam_tdisp[-1]/Pe +\
                    avg_lam_mdisp[-1]/Re))
                avg_global_turb_eta = np.append(avg_global_turb_eta,\
                    (B*avg_turb_tdisp[-1]/Pe)/\
                    (B*avg_turb_tdisp[-1]/Pe + avg_turb_mdisp[-1]/Re))

                vlam = np.append(vlam, len(idx_lam[0])/np_tot)
                vturb = np.append(vturb, len(idx_turb[0])/np_tot)
                            
            del uh
            del wrms
            del trms
            del eta
            del tdisp
            del mdisp
            del vortz
        
        # remove repeated timesteps and sort chronologically
        t, indices = np.unique(t_tot,return_index=True)
        del t_tot

        uh_rms = uh_rms[indices]
        vortz_rms = vortz_rms[indices]

        avg_wrms = avg_wrms[indices]
        avg_lam_wrms = avg_lam_wrms[indices]
        avg_turb_wrms = avg_turb_wrms[indices]

        avg_lam_wrms_wght = avg_lam_wrms_wght[indices]
        avg_lam_wrms_eff_wght = avg_lam_wrms_eff_wght[indices]
        avg_turb_wrms_wght = avg_turb_wrms_wght[indices]
        avg_turb_wrms_eff_wght = avg_turb_wrms_eff_wght[indices]

        avg_trms = avg_trms[indices]
        avg_lam_trms = avg_lam_trms[indices]
        avg_turb_trms = avg_turb_trms[indices]

        avg_lam_trms_wght = avg_lam_trms_wght[indices]
        avg_lam_trms_eff_wght = avg_lam_trms_eff_wght[indices]
        avg_turb_trms_wght = avg_turb_trms_wght[indices]
        avg_turb_trms_eff_wght = avg_turb_trms_eff_wght[indices]

        avg_tdisp = avg_tdisp[indices]
        avg_lam_tdisp = avg_lam_tdisp[indices]
        avg_turb_tdisp = avg_turb_tdisp[indices]

        avg_mdisp = avg_mdisp[indices]
        avg_lam_mdisp = avg_lam_mdisp[indices]
        avg_turb_mdisp = avg_turb_mdisp[indices]

        avg_local_eta = avg_local_eta[indices]
        avg_local_lam_eta = avg_local_lam_eta[indices]
        avg_local_turb_eta = avg_local_turb_eta[indices]

        avg_global_eta = avg_global_eta[indices]
        avg_global_lam_eta = avg_global_lam_eta[indices]
        avg_global_turb_eta = avg_global_turb_eta[indices]

        vlam = vlam[indices]
        vturb = vturb[indices]

        Nt = len(t)
        io_file.write('# Index {:03d}: Re = {:6.2f}, Pe = {:6.2f}, B = {:6.2f}'\
            .format(index_counter,Re,Pe,B))
        tavg_file.write('# Index {:03d}: Re = {:6.2f}, Pe = {:6.2f}'\
            .format(index_counter,Re,Pe))
        io_file.write('\n')
        tavg_file.write('\n')
        index_counter += 1
        # write to output file
        for i in range(Nt):
            io_file.write(fmt_str.format(Re,B,Pe/Re,Pe,B*Pe,\
            t[i],uh_rms[i],vortz_rms[i],avg_wrms[i],\
            avg_tdisp[i],avg_mdisp[i],avg_local_eta[i],avg_global_eta[i],\
            avg_lam_wrms[i],avg_lam_tdisp[i],avg_lam_mdisp[i],\
            avg_local_lam_eta[i],avg_global_lam_eta[i],\
            avg_turb_wrms[i],avg_turb_tdisp[i],avg_turb_mdisp[i],\
            avg_local_turb_eta[i],avg_global_turb_eta[i],\
            vlam[i],vturb[i],\
            avg_lam_wrms_wght[i],avg_turb_wrms_wght[i],\
            avg_lam_wrms_eff_wght[i],avg_turb_wrms_eff_wght[i],\
            avg_trms[i],\
            avg_lam_trms[i],avg_turb_trms[i],\
            avg_lam_trms_wght[i],avg_turb_trms_wght[i],\
            avg_lam_trms_eff_wght[i],avg_turb_trms_eff_wght[i],\
            ))
            io_file.write('\n')

        # do temporal averages and write to file
        lb, ub = bounds[m][n]
        tidx = np.where((lb <= t) & (t <= ub))

        uh_rms, uh_err = db.tavg(uh_rms, tidx)
        vortz_rms, vortz_err = db.tavg(vortz_rms, tidx)

        avg_wrms, err_wrms = db.tavg(avg_wrms, tidx)
        avg_lam_wrms, err_lam_wrms = db.discounted_tavg(avg_lam_wrms,\
            vlam, tidx)
        avg_turb_wrms, err_turb_wrms = db.discounted_tavg(avg_turb_wrms,\
            vturb, tidx)

        avg_lam_wrms_wght, err_lam_wrms_wght = db.tavg(avg_lam_wrms_wght, tidx)
        avg_turb_wrms_wght, err_turb_wrms_wght = db.tavg(avg_turb_wrms_wght, tidx)
        avg_lam_wrms_eff_wght, err_lam_wrms_eff_wght = \
            db.tavg(avg_lam_wrms_eff_wght, tidx)
        avg_turb_wrms_eff_wght, err_turb_wrms_eff_wght = \
            db.tavg(avg_turb_wrms_eff_wght, tidx)

        # TRMS
        avg_trms, err_trms = db.tavg(avg_trms, tidx)
        avg_lam_trms, err_lam_trms = db.discounted_tavg(avg_lam_trms, vlam, tidx)
        avg_turb_trms, err_turb_trms = db.discounted_tavg(avg_turb_trms, vturb, tidx)

        avg_lam_trms_wght, err_lam_trms_wght = db.tavg(avg_lam_trms_wght, tidx)
        avg_turb_trms_wght, err_turb_trms_wght = db.tavg(avg_turb_trms_wght, tidx)
        avg_lam_trms_eff_wght, err_lam_trms_eff_wght = \
            db.tavg(avg_lam_trms_eff_wght, tidx)
        avg_turb_trms_eff_wght, err_turb_trms_eff_wght = \
            db.tavg(avg_turb_trms_eff_wght, tidx)

        # TDISP
        avg_tdisp, err_tdisp = db.tavg(avg_tdisp, tidx)
        avg_lam_tdisp, err_lam_tdisp = db.discounted_tavg(avg_lam_tdisp, vlam, tidx)
        avg_turb_tdisp, err_turb_tdisp = db.discounted_tavg(avg_turb_tdisp, vturb, tidx)

        # MDISP
        avg_mdisp, err_mdisp = db.tavg(avg_mdisp, tidx)
        avg_lam_mdisp, err_lam_mdisp = db.discounted_tavg(avg_lam_mdisp, vlam, tidx)
        avg_turb_mdisp, err_turb_mdisp = db.discounted_tavg(avg_turb_mdisp, vturb, tidx)

        # Locally computed ETA
        avg_local_eta, err_local_eta = db.tavg(avg_local_eta, tidx)
        avg_local_lam_eta, err_local_lam_eta = db.discounted_tavg(avg_local_lam_eta,\
            vlam, tidx)
        avg_local_turb_eta, err_local_turb_eta = db.discounted_tavg(avg_local_turb_eta,\
            vturb, tidx)

        # Globally Computed ETA
        avg_global_eta, err_global_eta = db.tavg(avg_global_eta, tidx)
        avg_global_lam_eta, err_global_lam_eta = db.discounted_tavg(avg_global_lam_eta,\
            vlam, tidx)
        avg_global_turb_eta, err_global_turb_eta = db.discounted_tavg(avg_global_turb_eta,\
            vturb, tidx)

        # Volume Fractions
        vlam_avg, vlam_err = db.tavg(vlam, tidx)
        vturb_avg, vturb_err = db.tavg(vturb, tidx)


        tavg_file.write(tavg_fmt_str.format(Re,B,Pe/Re,Pe,B*Pe,lb,ub,\
            uh_rms,uh_err,vortz_rms,vortz_err,avg_wrms,err_wrms,\
            avg_tdisp,err_tdisp,avg_mdisp,err_mdisp,\
            avg_local_eta,err_local_eta,avg_global_eta,err_global_eta,\
            avg_lam_wrms,err_lam_wrms,\
            avg_lam_tdisp,err_lam_tdisp,avg_lam_mdisp,err_lam_mdisp,\
            avg_local_lam_eta,err_local_lam_eta,\
            avg_global_lam_eta,err_global_lam_eta,\
            avg_turb_wrms,err_turb_wrms,\
            avg_turb_tdisp,err_turb_tdisp,avg_turb_mdisp,err_turb_mdisp,\
            avg_local_turb_eta,err_local_turb_eta,\
            avg_global_turb_eta,err_global_turb_eta,\
            vlam_avg,vlam_err,vturb_avg,vturb_err,\
            avg_lam_wrms_wght,err_lam_wrms_wght,\
            avg_turb_wrms_wght,err_turb_wrms_wght,\
            avg_lam_wrms_eff_wght,err_lam_wrms_eff_wght,\
            avg_turb_wrms_eff_wght,err_turb_wrms_eff_wght,\
            avg_trms,err_trms,\
            avg_lam_trms,err_lam_trms,\
            avg_turb_trms,err_turb_trms,\
            avg_lam_trms_wght,err_lam_trms_wght,\
            avg_turb_trms_wght,err_turb_trms_wght,\
            avg_lam_trms_eff_wght,err_lam_trms_eff_wght,\
            avg_turb_trms_eff_wght,err_turb_trms_eff_wght))
        tavg_file.write('\n')
        io_file.write('\n\n\n')
    io_file.write('\n\n\n')
    tavg_file.write('\n\n\n')
    


