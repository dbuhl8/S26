import numpy as np
import dbuhlMod as db
import os
import fnmatch
from netCDF4 import MFDataset
from netCDF4 import Dataset


testing = False

if testing:
    # only run this code on a single simulation set to ensure accuracy of code
    Re600_Pe60=['nonrotating/B300Re600Pe60/']
    Re600_Pe60_bounds=[[40,80]]

    simulations = [Re600_Pe60]
    bounds = [Re600_Pe60_bounds]
else: 
    Re600_Pe60=['nonrotating/B300Re600Pe60/','nonrotating/B180Re600Pe60/',\
        'nonrotating/B100Re600Pe60/','nonrotating/B30Re600Pe60/',\
        'nonrotating/B10Re600Pe60/','nonrotating/B3Re600Pe60/',\
        'nonrotating/B1Re600Pe60/']
    Re1000_Pe100 = ['nonrotating/B100Re1000Pe100/', 'nonrotating/B10Re1000Pe100/']

    Re600_Pe60_bounds=[[40,80],[40,90],[45,160],\
        [35,80],[50,250],[55,110],[45,125]]
    Re1000_Pe100_bounds=[[40,70],[70,110]]

    simulations = [Re600_Pe60, Re1000_Pe100]
    bounds = [Re600_Pe60_bounds, Re1000_Pe100_bounds]

# preparing output file
header_string = "#" + "{:<s}    "*53
tavg_header_string = "#" + "{:<s}    "*101
fmt_str = "{:.4e}    "*53
tavg_fmt_str = "{:.4e}    "*101
# 20 columns
header_string = header_string.format('Re','B','Pr', 'Pe', 'BPe', 't',\
    'uh_rms','w_rms','temp_rms','vortz_rms',\
    'tdisp', 'mdisp', 'eta', 'TKE_disp', 'turb_dens_flux', 'turb_flux_coef',\
    'turb_Fr', 'grad_Ri', 'buoyancy_Re',\
    'LH', 'LZ', 'LT', 'LO', 'LE', 'LK',\
    'vlam', 'vturb')

io_file = open('length_scale_extraction.dat','w')
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
        w_rms = np.array([])

        avg_tdisp = np.array([])
        avg_mdisp = np.array([])
        avg_eta = np.array([])

        avg_TKE_disp = np.array([])
        avg_buoy_var_disp = np.array([])
        avg_turb_buoy_flux = np.array([])
        turb_flux_coef = np.array([])

        LH = np.array([])
        LZ = np.array([])
        LT = np.array([])
        LO = np.array([])
        LE = np.array([])
        LK = np.array([])

        ReB = np.array([])
        FrT = np.array([])
        RiG = np.array([])

        vlam = np.array([])
        vturb = np.array([])

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

            pert_uz = uz - np.mean(uz,axis=(2,3),keepdims=True)

            #dTdz = db.FD6Z(temp_inst,Nz,dz)
            tdisp = db.FD6Z(temp,Nz,dz)**2 + db.FD6Y(temp,Ny,dy)**2 +\
                db.FD6X(temp,Nx,dx)**2

            buoy_var_disp = tdisp/Pe
            turb_buoy_flux = pert_uz*temp
            del temp

            pert_ux = ux - np.mean(ux,axis=(2,3),keepdims=True)
            pert_uy = uy - np.mean(uy,axis=(2,3),keepdims=True)

            TKE = 0.5*(pert_ux**2 + pert_uy**2 + pert_uz**2)

            ux_x = db.FD6X(pert_ux,Nx,dx)
            ux_y = db.FD6Y(pert_ux,Ny,dy)
            ux_z = db.FD6Z(pert_ux,Nz,dz)
            del pert_ux
            uy_x = db.FD6X(pert_uy,Nx,dx)
            uy_y = db.FD6Y(pert_uy,Nx,dy)
            uy_z = db.FD6Z(pert_uy,Nx,dz)
            del pert_uy
            uz_x = db.FD6X(pert_uz,Nx,dx)
            uz_y = db.FD6Y(pert_uz,Ny,dy)
            uz_z = db.FD6Z(pert_uz,Nz,dz)
            del pert_uz
            TKE_disp = 2*(ux_x**2 + 0.5*ux_y**2 + 0.5*ux_z**2 \
                + 0.5*uy_x**2 + uy_y**2 + 0.5*uy_z**2\
                + 0.5*uz_x**2 + 0.5*uz_y**2 + uz_z**2\
                + ux_y*uy_x + ux_z*uz_x + uy_z*uz_y)/Re
            del ux_x 
            del ux_y
            del ux_z
            del uy_x
            del uy_y
            del uy_z
            del uz_x
            del uz_y
            del uz_z

            mdisp = db.FD6Z(ux,Nz,dz)**2 + db.FD6Y(ux,Ny,dy)**2 +\
                db.FD6X(ux,Nx,dx)**2 + db.FD6Z(uy,Nz,dz)**2 +\
                db.FD6Y(uy,Ny,dy)**2 + db.FD6X(uy,Nx,dx)**2 +\
                db.FD6Z(uz,Nz,dz)**2 + db.FD6Y(uz,Ny,dy)**2 +\
                db.FD6X(uz,Nx,dx)**2
            uh = (ux**2+uy**2)**0.5
            vortz = db.FD6X(uy,Nx,dx) - db.FD6Y(ux,Ny,dy)
            del ux
            del uy

            for i in range(Nt):
                # recording timestep and average uh
                t_tot = np.append(t_tot, t[i])
                uh_rms = np.append(uh_rms,db.rms(uh[i,:,:,:]))
                vortz_rms = np.append(vortz_rms,db.rms(vortz[i,:,:,:]))
                # computing wrms (and weighted version)
                w_rms = np.append(avg_wrms,db.rms(uz[i,:,:,:]))

                                # identifying turbulent and laminar regions
                # with rescaling to the effective Froude and vortz
                # effective vortz: vortz_eff = vortz/uh_rms
                # vortz**2 <= uh/Fr
                idx_lam = np.where(vortz[i,:,:,:]**2 <=\
                    uh_rms[-1]/Fr)
                idx_turb = np.where(vortz[i,:,:,:]**2 >\
                    uh_rms[-1]/Fr)

            
                # computing avg tdisp and mdisp, NOTE that TDISP and MDisp do
                # not need to be sqrt_ed, so we use the mean() function
                avg_tdisp = np.append(avg_tdisp,db.mean(tdisp[i,:,:,:]))

                avg_mdisp = np.append(avg_mdisp,db.mean(mdisp[i,:,:,:]))

                avg_TKE = np.append(avg_TKE,db.mean(TKE[i,:,:,:]))

                avg_TKE_disp = np.append(avg_TKE_disp,db.mean(TKE_disp[i,:,:,:])
                avg_buoy_var_disp = np.append(avg_buoy_var_disp,\
                    db.mean(buoy_var_disp[i,:,:,:]))
                avg_turb_buoy_flux = np.append(avg_turb_buoy_flux,\
                    db.mean(turb_buoy_flux[i,:,:,:]))
                turb_flux_coef = np.append(turb_flux_coef,\
                    avg_turb_buoy_flux[-1]/avg_TKE_disp[-1])

                # eta = (thermal dissipation) / (total energy dissipation)
                # eta = (B/Pe)*tdisp / ((B/Pe)*tdisp + (1/Re)*mdisp)
                eta = np.append(eta,\
                    (B*avg_tdisp[-1]/Pe)/\
                    (B*avg_tdisp[-1]/Pe + avg_mdisp[-1]/Re))

                ReB = np.append(ReB,avg_TKE_disp[-1]*Fr**2/Re)    
                FrT = np.append(FrT,avg_TKE_disp[-1]*Fr/avg_TKE[-1])

                LH = 

                vlam = np.append(vlam, len(idx_lam[0])/np_tot)
                vturb = np.append(vturb, len(idx_turb[0])/np_tot)

            del uh
            del uz
            del tdisp
            del mdisp
            del vortz
        
        # remove repeated timesteps and sort chronologically
        t, indices = np.unique(t_tot,return_index=True)
        del t_tot

        uh_rms = uh_rms[indices]
        vortz_rms = vortz_rms[indices]
        buoy_rms = buoy_rms[indices]
        w_rms = wrms[indices]

        tdisp = tdisp[indices]
        mdisp = mdisp[indices]
        eta = eta[indices]

        TKE_disp = TKE_disp[indices]
        buoy_var_disp = buoy_var_disp[indices]
        turb_dens_flux = turb_dens_flux[indices]
        turb_flux_coef = turb_flux_coef[indices]

        ReB = 
        FrT = 
        RiG = 

        LH = LH[indices]
        LZ = LZ[indices]
        LT = LT[indices]
        LO = LO[indices]
        LE = LE[indices]
        LK = LK[indices]

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
            avg_lam_Fr_wrms[i],avg_lam_Fr_tdisp[i],avg_lam_Fr_mdisp[i],\
            avg_local_lam_Fr_eta[i],avg_global_lam_Fr_eta[i],\
            avg_turb_Fr_wrms[i],avg_turb_Fr_tdisp[i],avg_turb_Fr_mdisp[i],\
            avg_local_turb_Fr_eta[i],avg_global_turb_Fr_eta[i],\
            avg_lam_Fr_vortz_wrms[i],\
            avg_lam_Fr_vortz_tdisp[i],avg_lam_Fr_vortz_mdisp[i],\
            avg_local_lam_Fr_vortz_eta[i],avg_global_lam_Fr_vortz_eta[i],\
            avg_turb_Fr_vortz_wrms[i],\
            avg_turb_Fr_vortz_tdisp[i],avg_turb_Fr_vortz_mdisp[i],\
            avg_local_turb_Fr_vortz_eta[i],avg_global_turb_Fr_vortz_eta[i],\
            vlam[i],vturb[i],\
            vlam_Fr[i],vturb_Fr[i],\
            vlam_Fr_vortz[i],vturb_Fr_vortz[i],\
            avg_lam_wrms_wght[i], avg_turb_wrms_wght[i],\
            avg_lam_wrms_eff_wght[i], avg_turb_wrms_eff_wght[i],\
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
        avg_lam_Fr_wrms, err_lam_Fr_wrms = db.discounted_tavg(avg_lam_Fr_wrms,\
            vlam_Fr, tidx)
        avg_lam_Fr_vortz_wrms, err_lam_Fr_vortz_wrms = db.discounted_tavg(\
            avg_lam_Fr_vortz_wrms, vlam_Fr_vortz, tidx)
        avg_turb_wrms, err_turb_wrms = db.discounted_tavg(avg_turb_wrms,\
            vturb, tidx)
        avg_turb_Fr_wrms, err_turb_Fr_wrms = db.discounted_tavg(\
            avg_turb_Fr_wrms, vturb_Fr, tidx)
        avg_turb_Fr_vortz_wrms, err_turb_Fr_vortz_wrms = db.discounted_tavg(\
            avg_turb_Fr_vortz_wrms, vturb_Fr_vortz, tidx)

        avg_lam_wrms_wght, err_lam_wrms_wght = db.tavg(avg_lam_wrms_wght,\
            tidx)
        avg_turb_wrms_wght, err_turb_wrms_wght = db.tavg(avg_turb_wrms_wght,\
            tidx)
        avg_lam_wrms_eff_wght, err_lam_wrms_eff_wght = \
            db.tavg(avg_lam_wrms_eff_wght, tidx)
        avg_turb_wrms_eff_wght, err_turb_wrms_eff_wght = \
            db.tavg(avg_turb_wrms_eff_wght, tidx)

        # TDISP
        avg_tdisp, err_tdisp = db.tavg(avg_tdisp, tidx)
        avg_lam_tdisp, err_lam_tdisp = db.discounted_tavg(avg_lam_tdisp,\
            vlam, tidx)
        avg_lam_Fr_tdisp, err_lam_Fr_tdisp = db.discounted_tavg(avg_lam_Fr_tdisp,\
            vlam_Fr, tidx)
        avg_lam_Fr_vortz_tdisp, err_lam_Fr_vortz_tdisp = db.discounted_tavg(\
            avg_lam_Fr_vortz_tdisp, vlam_Fr_vortz, tidx)
        avg_turb_tdisp, err_turb_tdisp = db.discounted_tavg(avg_turb_tdisp,\
            vturb, tidx)
        avg_turb_Fr_tdisp, err_turb_Fr_tdisp = db.discounted_tavg(\
            avg_turb_Fr_tdisp, vturb_Fr, tidx)
        avg_turb_Fr_vortz_tdisp, err_turb_Fr_vortz_tdisp = db.discounted_tavg(\
            avg_turb_Fr_vortz_tdisp, vturb_Fr_vortz, tidx)

        # MDISP
        avg_mdisp, err_mdisp = db.tavg(avg_mdisp, tidx)
        avg_lam_mdisp, err_lam_mdisp = db.discounted_tavg(avg_lam_mdisp,\
            vlam, tidx)
        avg_lam_Fr_mdisp, err_lam_Fr_mdisp = db.discounted_tavg(avg_lam_Fr_mdisp,\
            vlam_Fr, tidx)
        avg_lam_Fr_vortz_mdisp, err_lam_Fr_vortz_mdisp = db.discounted_tavg(\
            avg_lam_Fr_vortz_mdisp, vlam_Fr_vortz, tidx)
        avg_turb_mdisp, err_turb_mdisp = db.discounted_tavg(avg_turb_mdisp,\
            vturb, tidx)
        avg_turb_Fr_mdisp, err_turb_Fr_mdisp = db.discounted_tavg(\
            avg_turb_Fr_mdisp, vturb_Fr, tidx)
        avg_turb_Fr_vortz_mdisp, err_turb_Fr_vortz_mdisp = db.discounted_tavg(\
            avg_turb_Fr_vortz_mdisp, vturb_Fr_vortz, tidx)

        # Locally computed ETA
        avg_local_eta, err_local_eta = db.tavg(avg_local_eta, tidx)
        avg_local_lam_eta, err_local_lam_eta = db.discounted_tavg(avg_local_lam_eta,\
            vlam, tidx)
        avg_local_lam_Fr_eta, err_local_lam_Fr_eta = db.discounted_tavg(avg_local_lam_Fr_eta,\
            vlam_Fr, tidx)
        avg_local_lam_Fr_vortz_eta, err_local_lam_Fr_vortz_eta = db.discounted_tavg(\
            avg_local_lam_Fr_vortz_eta, vlam_Fr_vortz, tidx)
        avg_local_turb_eta, err_local_turb_eta = db.discounted_tavg(avg_local_turb_eta,\
            vturb, tidx)
        avg_local_turb_Fr_eta, err_local_turb_Fr_eta = db.discounted_tavg(\
            avg_local_turb_Fr_eta, vturb_Fr, tidx)
        avg_local_turb_Fr_vortz_eta, err_local_turb_Fr_vortz_eta = db.discounted_tavg(\
            avg_local_turb_Fr_vortz_eta, vturb_Fr_vortz, tidx)

        # Globally Computed ETA
        avg_global_eta, err_global_eta = db.tavg(avg_global_eta, tidx)
        avg_global_lam_eta, err_global_lam_eta = db.discounted_tavg(avg_global_lam_eta,\
            vlam, tidx)
        avg_global_lam_Fr_eta, err_global_lam_Fr_eta = db.discounted_tavg(avg_global_lam_Fr_eta,\
            vlam_Fr, tidx)
        avg_global_lam_Fr_vortz_eta, err_global_lam_Fr_vortz_eta = db.discounted_tavg(\
            avg_global_lam_Fr_vortz_eta, vlam_Fr_vortz, tidx)
        avg_global_turb_eta, err_global_turb_eta = db.discounted_tavg(avg_global_turb_eta,\
            vturb, tidx)
        avg_global_turb_Fr_eta, err_global_turb_Fr_eta = db.discounted_tavg(\
            avg_global_turb_Fr_eta, vturb_Fr, tidx)
        avg_global_turb_Fr_vortz_eta, err_global_turb_Fr_vortz_eta = db.discounted_tavg(\
            avg_global_turb_Fr_vortz_eta, vturb_Fr_vortz, tidx)

        # Volume Fractions
        vlam_avg, vlam_err = db.tavg(vlam,tidx)
        vlam_Fr_avg, vlam_Fr_err = db.tavg(vlam_Fr,tidx)
        vlam_Fr_vortz_avg, vlam_Fr_vortz_err = db.tavg(vlam_Fr_vortz,tidx)
        vturb_avg, vturb_err = db.tavg(vturb,tidx)
        vturb_Fr_avg, vturb_Fr_err = db.tavg(vturb_Fr,tidx)
        vturb_Fr_vortz_avg, vturb_Fr_vortz_err = db.tavg(vturb_Fr_vortz,tidx)

        """
        print('1:',Re,'2:',B,'3:',Pe/Re,'4:',Pe,'5:',B*Pe,'6:',lb,'7:',ub,'8:',\
            uh_rms,'9:',uh_err,'10:',vortz_rms,'11:',vortz_err,'12:',avg_wrms,'13:',err_wrms,'14:',\
            avg_tdisp,'15:',err_tdisp,'16:',avg_mdisp,'17:',err_mdisp,'18:',\
            avg_local_eta,'19:',err_local_eta,'20:',avg_global_eta,'21:',err_global_eta,'22:',\
            avg_lam_wrms,'23:',err_lam_wrms,'24:',\
            avg_lam_tdisp,'25:',err_lam_tdisp,'26:',avg_lam_mdisp,'27:',err_lam_mdisp,'28:',\
            avg_local_lam_eta,'29:',err_local_lam_eta,'30:',\
            avg_global_lam_eta,'31:',err_global_lam_eta,'32:',\
            avg_turb_wrms,'33:',err_turb_wrms,'34:',\
            avg_turb_tdisp,'35:',err_turb_tdisp,'36:',avg_turb_mdisp,'37:',err_turb_mdisp,'38:',\
            avg_local_turb_eta,'39:',err_local_turb_eta,'40:',\
            avg_global_turb_eta,'41:',err_global_turb_eta,'42:',\
            avg_lam_Fr_wrms,'42:',err_lam_Fr_wrms,'44:',\
            avg_lam_Fr_tdisp,'45:',err_lam_Fr_tdisp,'46:',\
            avg_lam_Fr_mdisp,'47:',err_lam_Fr_tdisp,'48:',\
            avg_local_lam_Fr_eta,'49:',err_local_lam_Fr_eta,'50:',\
            avg_global_lam_Fr_eta,'51:',err_global_lam_Fr_eta,'152:',\
            avg_turb_Fr_wrms,'53:',err_turb_Fr_wrms,'54:',\
            avg_turb_Fr_tdisp,'55:',err_turb_Fr_tdisp,'56:',\
            avg_turb_Fr_mdisp,'57:',err_turb_Fr_mdisp,'58:',\
            avg_local_turb_Fr_eta,'59:',err_local_turb_Fr_eta,'60:',\
            avg_global_turb_Fr_eta,'61:',err_global_turb_Fr_eta,'62:',\
            avg_lam_Fr_vortz_wrms,'63:', err_lam_Fr_vortz_wrms,'64:',\
            avg_lam_Fr_vortz_tdisp,'65:',err_lam_Fr_vortz_tdisp,'66:',\
            avg_lam_Fr_vortz_mdisp,'67:',err_lam_Fr_vortz_mdisp,'68:',\
            avg_local_lam_Fr_vortz_eta,'69:',err_local_lam_Fr_vortz_eta,'70:',\
            avg_global_lam_Fr_vortz_eta,'71:',err_global_lam_Fr_vortz_eta,'72:',\
            avg_turb_Fr_vortz_wrms,'73:', err_turb_Fr_vortz_wrms,'74:',\
            avg_turb_Fr_vortz_tdisp,'75:',err_turb_Fr_vortz_tdisp,'76:',\
            avg_turb_Fr_vortz_mdisp,'77:',err_turb_Fr_vortz_mdisp,'78:',\
            avg_local_turb_Fr_vortz_eta,'79:',err_local_turb_Fr_vortz_eta,'80:',\
            avg_global_turb_Fr_vortz_eta,'81:',err_global_turb_Fr_vortz_eta,'82:',\
            vlam_avg,'83:',vlam_err,'84:',vturb_avg,'85:',vturb_err,'86:',\
            vlam_Fr_avg,'87:', vlam_Fr_err,'88:',vturb_Fr_avg,'89:',vturb_Fr_err,'90:',\
            vlam_Fr_vortz_avg,'91:', vlam_Fr_vortz_err,'92:',\
            vturb_Fr_vortz_avg,'93:',vturb_Fr_vortz_err,'94:',\
            avg_lam_wrms_wght,'95:',err_lam_wrms_wght,'96:',\
            avg_turb_wrms_wght,'97:',err_turb_wrms_wght,'98:',\
            avg_lam_wrms_eff_wght,'99:',err_lam_wrms_eff_wght,'100:',\
            avg_turb_wrms_eff_wght,'101:',err_turb_wrms_eff_wght)
        """

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
            avg_lam_Fr_wrms,err_lam_Fr_wrms,\
            avg_lam_Fr_tdisp,err_lam_Fr_tdisp,\
            avg_lam_Fr_mdisp,err_lam_Fr_tdisp,\
            avg_local_lam_Fr_eta,err_local_lam_Fr_eta,\
            avg_global_lam_Fr_eta,err_global_lam_Fr_eta,\
            avg_turb_Fr_wrms,err_turb_Fr_wrms,\
            avg_turb_Fr_tdisp,err_turb_Fr_tdisp,\
            avg_turb_Fr_mdisp,err_turb_Fr_mdisp,\
            avg_local_turb_Fr_eta,err_local_turb_Fr_eta,\
            avg_global_turb_Fr_eta,err_global_turb_Fr_eta,\
            avg_lam_Fr_vortz_wrms, err_lam_Fr_vortz_wrms,\
            avg_lam_Fr_vortz_tdisp,err_lam_Fr_vortz_tdisp,\
            avg_lam_Fr_vortz_mdisp,err_lam_Fr_vortz_mdisp,\
            avg_local_lam_Fr_vortz_eta,err_local_lam_Fr_vortz_eta,\
            avg_global_lam_Fr_vortz_eta,err_global_lam_Fr_vortz_eta,\
            avg_turb_Fr_vortz_wrms, err_turb_Fr_vortz_wrms,\
            avg_turb_Fr_vortz_tdisp,err_turb_Fr_vortz_tdisp,\
            avg_turb_Fr_vortz_mdisp,err_turb_Fr_vortz_mdisp,\
            avg_local_turb_Fr_vortz_eta,err_local_turb_Fr_vortz_eta,\
            avg_global_turb_Fr_vortz_eta,err_global_turb_Fr_vortz_eta,\
            vlam_avg,vlam_err,vturb_avg,vturb_err,\
            vlam_Fr_avg, vlam_Fr_err,vturb_Fr_avg,vturb_Fr_err,\
            vlam_Fr_vortz_avg, vlam_Fr_vortz_err,\
            vturb_Fr_vortz_avg,vturb_Fr_vortz_err,\
            avg_lam_wrms_wght,err_lam_wrms_wght,\
            avg_turb_wrms_wght,err_turb_wrms_wght,\
            avg_lam_wrms_eff_wght,err_lam_wrms_eff_wght,\
            avg_turb_wrms_eff_wght,err_turb_wrms_eff_wght))
        tavg_file.write('\n')
        io_file.write('\n\n\n')
    io_file.write('\n\n\n')
    tavg_file.write('\n\n\n')
    


