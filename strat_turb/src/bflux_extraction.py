import numpy as np
import dbuhlMod as db
import os
import fnmatch
from netCDF4 import Dataset


testing = False

# ---------------------------------------------------------------------------
# Stochastic forcing simulations
# ---------------------------------------------------------------------------
if testing:
    stoch_Re600_Pe60 = ['nonrotating/B300Re600Pe60/']
    stoch_Re600_Pe60_bounds = [[40, 80]]
    stoch_Re600_Pe60_params = [(600, 60, 300)]

    stoch_simulations = [stoch_Re600_Pe60]
    stoch_bounds = [stoch_Re600_Pe60_bounds]
    stoch_params = [stoch_Re600_Pe60_params]
else:
    stoch_Re600_Pe60 = ['nonrotating/B1Re600Pe60/',   'nonrotating/B3Re600Pe60/',
                        'nonrotating/B30Re600Pe60/',  'nonrotating/B100Re600Pe60/',
                        'nonrotating/B180Re600Pe60/', 'nonrotating/B300Re600Pe60/']
    stoch_Re1000_Pe100 = ['nonrotating/B10Re1000Pe100/',
                          'nonrotating/B100Re1000Pe100/']

    stoch_Re600_Pe60_bounds = [[45, 125], [55, 110], [35, 80],
                               [45, 160], [40, 90],  [40, 80]]
    stoch_Re1000_Pe100_bounds = [[70, 110], [40, 70]]

    stoch_Re600_Pe60_params  = [(600, 60, 1),   (600, 60, 3),   (600, 60, 30),
                                (600, 60, 100), (600, 60, 180), (600, 60, 300)]
    stoch_Re1000_Pe100_params = [(1000, 100, 10), (1000, 100, 100)]

    stoch_Re300_Pe30 = ['nonrotating/B10Re300Pe30/']
    stoch_Re300_Pe30_bounds = [[50, 250]]
    stoch_Re300_Pe30_params = [(300, 30, 10)]

    stoch_simulations = [stoch_Re600_Pe60, stoch_Re1000_Pe100, stoch_Re300_Pe30]
    stoch_bounds = [stoch_Re600_Pe60_bounds, stoch_Re1000_Pe100_bounds, stoch_Re300_Pe30_bounds]
    stoch_params = [stoch_Re600_Pe60_params, stoch_Re1000_Pe100_params, stoch_Re300_Pe30_params]

# ---------------------------------------------------------------------------
# Steady forcing simulations
# ---------------------------------------------------------------------------
if testing:
    steady_Re600_Pe60 = ['horizontal-shear/Re600_Pe60_B100/']
    steady_Re600_Pe60_bounds = [[150, 300]]
    steady_Re600_Pe60_params = [(600, 60, 100)]

    steady_simulations = [steady_Re600_Pe60]
    steady_bounds = [steady_Re600_Pe60_bounds]
    steady_params = [steady_Re600_Pe60_params]
else:
    steady_Re600_Pe60 = ['horizontal-shear/Re600_Pe60_B0.1/',
                         'horizontal-shear/Re600_Pe60_B1/',
                         'horizontal-shear/Re600_Pe60_B10/',
                         'horizontal-shear/Re600_Pe60_B30/',
                         'horizontal-shear/Re600_Pe60_B100/',
                         'horizontal-shear/Re600_Pe60_B400/',
                         'horizontal-shear/Re600_Pe60_B1000/',
                         'horizontal-shear/Re600_Pe60_B3000/',
                         'horizontal-shear/Re600_Pe60_B6000/']
    steady_Re600_Pe60_bounds = [[275, 400], [180, 320], [220, 450],
                                [375, 560], [150, 300], [225, 460],
                                [500, 900], [650, 900], [40, 180]]
    steady_Re600_Pe60_params = [(600, 60, 0.1),  (600, 60, 1),    (600, 60, 10),
                                (600, 60, 30),   (600, 60, 100),  (600, 60, 400),
                                (600, 60, 1000), (600, 60, 3000), (600, 60, 6000)]

    steady_Re600_Pe30 = ['horizontal-shear/Re600_Pe30_B10/',
                         'horizontal-shear/Re600_Pe30_B30/',
                         'horizontal-shear/Re600_Pe30_B100/']
    steady_Re600_Pe30_bounds = [[260, 300], [400, 480], [300, 400]]
    steady_Re600_Pe30_params = [(600, 30, 10), (600, 30, 30), (600, 30, 100)]

    steady_Re1000_Pe100 = ['horizontal-shear/Re1000_Pe100_B3/',
                           'horizontal-shear/Re1000_Pe100_B10/',
                           'horizontal-shear/Re1000_Pe100_B30/',
                           'horizontal-shear/Re1000_Pe100_B100/',
                           'horizontal-shear/Re1000_Pe100_B300/']
    steady_Re1000_Pe100_bounds = [[65, 90], [30, 70], [50, 90], [90, 110], [109, 120]]
    steady_Re1000_Pe100_params = [(1000, 100, 3),  (1000, 100, 10), (1000, 100, 30),
                                  (1000, 100, 100), (1000, 100, 300)]

    steady_Re300_Pe30 = ['horizontal-shear/Re300_Pe30_B0.01/',
                         'horizontal-shear/Re300_Pe30_B0.1/',
                         'horizontal-shear/Re300_Pe30_B1/',
                         'horizontal-shear/Re300_Pe30_B10/',
                         'horizontal-shear/Re300_Pe30_B30/',
                         'horizontal-shear/Re300_Pe30_B100/',
                         'horizontal-shear/Re300_Pe30_B300/',
                         'horizontal-shear/Re300_Pe30_B1000/',
                         'horizontal-shear/Re300_Pe30_B10000/']
    steady_Re300_Pe30_bounds = [[370, 460], [260, 460], [140, 280], [60, 160],
                                [50, 300],  [300, 1200], [300, 900],
                                [800, 1300], [1050, 1450]]
    steady_Re300_Pe30_params = [(300, 30, 0.01), (300, 30, 0.1),   (300, 30, 1),
                                (300, 30, 10),   (300, 30, 30),    (300, 30, 100),
                                (300, 30, 300),  (300, 30, 1000),  (300, 30, 10000)]

    steady_Re1000_Pe10 = ['horizontal-shear/Re1000_Pe10_B3/',
                          'horizontal-shear/Re1000_Pe10_B100/',
                          'horizontal-shear/Re1000_Pe10_B1000/']
                          #'horizontal-shear/Re1000_Pe10_B10/',\
    steady_Re1000_Pe10_bounds = [[35, 40], [35, 70], [94, 95], [62, 63]]
    steady_Re1000_Pe10_params = [(1000, 10, 3), (1000, 10, 100), (1000, 10, 1000)]

    steady_Re600_Pe600 = ['horizontal-shear/Re600_Pe600_B100/']
    steady_Re600_Pe600_bounds = [[100, 125]]
    steady_Re600_Pe600_params = [(600, 600, 100)]

    steady_simulations = [steady_Re600_Pe60, steady_Re600_Pe30,
                          steady_Re1000_Pe100, steady_Re300_Pe30,
                          steady_Re1000_Pe10, steady_Re600_Pe600]
    steady_bounds = [steady_Re600_Pe60_bounds, steady_Re600_Pe30_bounds,
                     steady_Re1000_Pe100_bounds, steady_Re300_Pe30_bounds,
                     steady_Re1000_Pe10_bounds, steady_Re600_Pe600_bounds]
    steady_params = [steady_Re600_Pe60_params, steady_Re600_Pe30_params,
                     steady_Re1000_Pe100_params, steady_Re300_Pe30_params,
                     steady_Re1000_Pe10_params, steady_Re600_Pe600_params]

# ---------------------------------------------------------------------------
# Combine into a single list in the desired index order:
#   Index group 0: Re 300  Pe 30  steady
#   Index group 1: Re 600  Pe 30  steady
#   Index group 2: Re 600  Pe 60  steady
#   Index group 3: Re 1000 Pe 10  steady
#   Index group 4: Re 1000 Pe 100 steady
#   Index group 5: Re 600  Pe 600 steady  (Pr=1)
#   Index group 6: Re 600  Pe 60  stoch
#   Index group 7: Re 1000 Pe 100 stoch
#   Index group 8: Re 300  Pe 30  stoch
# (each group spans multiple indices, one per B value)
# ---------------------------------------------------------------------------
simulations = [steady_Re300_Pe30, steady_Re600_Pe30, steady_Re600_Pe60,
               steady_Re1000_Pe10, steady_Re1000_Pe100, steady_Re600_Pe600,
               stoch_Re600_Pe60, stoch_Re1000_Pe100, stoch_Re300_Pe30]
bounds = [steady_Re300_Pe30_bounds, steady_Re600_Pe30_bounds, steady_Re600_Pe60_bounds,
          steady_Re1000_Pe10_bounds, steady_Re1000_Pe100_bounds, steady_Re600_Pe600_bounds,
          stoch_Re600_Pe60_bounds, stoch_Re1000_Pe100_bounds, stoch_Re300_Pe30_bounds]
params = [steady_Re300_Pe30_params, steady_Re600_Pe30_params, steady_Re600_Pe60_params,
          steady_Re1000_Pe10_params, steady_Re1000_Pe100_params, steady_Re600_Pe600_params,
          stoch_Re600_Pe60_params, stoch_Re1000_Pe100_params, stoch_Re300_Pe30_params]

# flat dict: dir_path -> (Re, Pe, B)
sim_params = {path: p
              for sim_set, param_set in zip(simulations, params)
              for path, p in zip(sim_set, param_set)}
# Column indices (0-based) for quantities extracted from OUT* files.
# Both formats share: col 1 = t, col 5 = TEMPrms (brms).
# Steady  OUT format: (i7, 35E20.7)
# Stoch   OUT format: (i7, 50E20.7)
OUT_COLS = {
    'steady': {
        't':        1,
        'uxrms':   27,
        'uyrms':   28,
        'uzrms':   29,
        'vortzrms':32,
        'brms':     5,
        'tdisp':   33,
        'mdisp':   35,
    },
    'stoch': {
        't':        1,
        'uxrms':   35,
        'uyrms':   36,
        'uzrms':   37,
        'vortzrms':40,
        'brms':     5,
        'tdisp':   44,
        'mdisp':   46,
    },
}

tavg_filename = 'tavg_bflux.dat'

# ---------------------------------------------------------------------------
# Output format strings and headers
# Tavg file: 31 columns (one row per simulation)
# ---------------------------------------------------------------------------
# tavg: Re Pe B lb ub
#       ux_rms ux_err uy_rms uy_err uz_rms uz_err
#       vortz_rms vortz_err brms brms_err
#       tdisp tdisp_err mdisp mdisp_err
#       vturb vturb_err
#       turb_uzrms turb_uzrms_err lam_uzrms lam_uzrms_err
#       turb_brms turb_brms_err lam_brms lam_brms_err
#       uh_rms uh_rms_err
tavg_ncols = 31
tavg_fmt_str = "{:.4e}    " * tavg_ncols
tavg_header_string = ("#" + "{:<s}    " * tavg_ncols).format(
    'Re', 'Pe', 'B', 'lb', 'ub',
    'ux_rms', 'ux_err', 'uy_rms', 'uy_err',
    'uz_rms', 'uz_err',
    'vortz_rms', 'vortz_err',
    'brms', 'brms_err',
    'tdisp', 'tdisp_err',
    'mdisp', 'mdisp_err',
    'vturb', 'vturb_err',
    'turb_uzrms', 'turb_uzrms_err',
    'lam_uzrms', 'lam_uzrms_err',
    'turb_brms', 'turb_brms_err',
    'lam_brms', 'lam_brms_err',
    'uh_rms', 'uh_rms_err')

tavg_file = open(tavg_filename, 'w')
tavg_file.write(tavg_header_string)
tavg_file.write('\n')

# write column key file
key_filename = tavg_filename.replace('.dat', '_column_key.txt')
with open(key_filename, 'w') as key_file:
    tavg_cols = [
        'Re', 'Pe', 'B', 'lb', 'ub',
        'ux_rms', 'ux_err', 'uy_rms', 'uy_err',
        'uz_rms', 'uz_err',
        'vortz_rms', 'vortz_err',
        'brms', 'brms_err',
        'tdisp', 'tdisp_err',
        'mdisp', 'mdisp_err',
        'vturb', 'vturb_err',
        'turb_uzrms', 'turb_uzrms_err',
        'lam_uzrms', 'lam_uzrms_err',
        'turb_brms', 'turb_brms_err',
        'lam_brms', 'lam_brms_err',
        'uh_rms', 'uh_rms_err']

    key_file.write('Column key for: {}\n'.format(tavg_filename))
    key_file.write('(time-averaged output — one row per simulation)\n\n')
    for i, name in enumerate(tavg_cols):
        key_file.write('  Col {:3d}: {}\n'.format(i + 1, name))

index_counter = 0


for m, sim_set in enumerate(simulations):
    # loops over each parameter set (i.e. Re = 600, Pe = 60 is one iteration
    # of the loop)
    tavg_file.write('# Index {:03d}\n'.format(index_counter))
    index_counter += 1
    for n, dir_path in enumerate(sim_set):
        Re_p, Pe_p, B_p = sim_params[dir_path]
        print("Starting sim: Re {}, Pe {}, B {}".format(Re_p, Pe_p, B_p))
        # loops over each B value in that parameter suite
        directory = os.listdir(dir_path)
        simdat_files = fnmatch.filter(directory, 'simdat*.cdf')
        simdat_files = [dir_path + x for x in simdat_files]
        lb, ub = bounds[m][n]

        # ----------------------------------------------------------------
        # OUT* file pre-processing
        # Load here so that uhrms_out is available for the lam/turb
        # threshold in the simdat extraction below.
        # ----------------------------------------------------------------
        out_files = sorted(fnmatch.filter(directory, 'OUT*'))
        out_files = [dir_path + f for f in out_files]

        out_data_list = []
        if len(out_files) > 0:
            # peek at first data line (skip comments/blanks) to count columns
            with open(out_files[0]) as _f:
                for _line in _f:
                    _line = _line.strip()
                    if _line and not _line.startswith('#'):
                        ncols_out = len(_line.split())
                        break
            cols = OUT_COLS['steady'] if ncols_out == 36 else OUT_COLS['stoch']
            # only load the columns we actually use
            needed_cols = sorted(set(cols.values()))
            col_map = {orig: new for new, orig in enumerate(needed_cols)}
            cols = {k: col_map[v] for k, v in cols.items()}

            for of in out_files:
                try:
                    out_data_list.append(np.loadtxt(of, usecols=needed_cols))
                except Exception as e:
                    print('Warning: could not read {}: {}'.format(of, e))

        if len(out_data_list) == 0:
            print('Warning: no OUT* files found in {} — OUT quantities set to NaN'
                  .format(dir_path))
            nan = float('nan')
            uhrms_out    = uhrms_err    = nan
            uxrms_out    = uxrms_err    = nan
            uyrms_out    = uyrms_err    = nan
            uzrms_out    = uzrms_err    = nan
            vortzrms_out = vortzrms_err = nan
            brms_out     = brms_err     = nan
            tdisp_out    = tdisp_err    = nan
            mdisp_out    = mdisp_err    = nan
        else:
            out_data = np.vstack(out_data_list)

            # deduplicate and sort by time
            t_out = out_data[:, cols['t']]
            _, out_idx = np.unique(t_out, return_index=True)
            out_data = out_data[out_idx, :]
            t_out = out_data[:, cols['t']]

            tidx_out = np.where((lb <= t_out) & (t_out <= ub))

            def _tavg_col(col_idx):
                return db.tavg(out_data[:, col_idx], tidx_out)

            uxrms_out,    uxrms_err    = _tavg_col(cols['uxrms'])
            uyrms_out,    uyrms_err    = _tavg_col(cols['uyrms'])
            uzrms_out,    uzrms_err    = _tavg_col(cols['uzrms'])
            vortzrms_out, vortzrms_err = _tavg_col(cols['vortzrms'])
            brms_out,     brms_err     = _tavg_col(cols['brms'])
            tdisp_out,    tdisp_err    = _tavg_col(cols['tdisp'])
            mdisp_out,    mdisp_err    = _tavg_col(cols['mdisp'])
            uhrms_ts  = np.sqrt(out_data[:, cols['uxrms']]**2 + out_data[:, cols['uyrms']]**2)
            uhrms_out,    uhrms_err    = db.tavg(uhrms_ts, tidx_out)

        B = 0
        Re = 0
        Pe = 0
        t_tot = np.array([])

        avg_lam_brms = np.array([])
        avg_turb_brms = np.array([])

        avg_uzrms_lam = np.array([])
        avg_uzrms_turb = np.array([])

        vturb = np.array([])
        vlam = np.array([])   # needed internally for discounted_tavg
        np_tot = 0

        for j, fn in enumerate(simdat_files):
            cdf_file = Dataset(fn, exclude='Chem')

            x = np.array(cdf_file.variables["x"])
            y = np.array(cdf_file.variables["y"])
            z = np.array(cdf_file.variables["z"])
            t = np.array(cdf_file.variables["t"])
            Nx = len(x)
            Ny = len(y)
            Nz = len(z)
            np_tot = Nx * Ny * Nz
            Nt = len(t)
            #print('Number of Timesteps: ', Nt)
            del x
            del y
            del z
            gx = cdf_file.variables["Gammax"][0]
            gy = cdf_file.variables["Gammay"][0]
            gz = cdf_file.variables["Gammaz"][0]
            dx = gx / Nx
            dy = gy / Ny
            dz = gz / Nz
            Re = 1. / cdf_file.variables["D_visc"][0]
            Pe = 1. / cdf_file.variables["D_therm"][0]
            B = cdf_file.variables["B_therm"][0]
            Fr = 1. / np.sqrt(B)
            #print("Parameters: Re - ", Re, ", Pe - ", Pe, ', B - ', B)

            ux = np.array(cdf_file.variables["ux"])
            uy = np.array(cdf_file.variables["uy"])
            uz = np.array(cdf_file.variables["uz"])
            temp = np.array(cdf_file.variables["Temp"])

            cdf_file.close()

            brms = temp
            del temp

            vortz = db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
            del ux
            del uy

            uzrms = uz
            del uz

            for i in range(Nt):
                t_tot = np.append(t_tot, t[i])

                # lam/turb threshold: vortz**2 <= uhrms_out / Fr
                idx_lam = np.where(vortz[i, :, :, :]**2 <= uhrms_out / Fr)
                idx_turb = np.where(vortz[i, :, :, :]**2 > uhrms_out / Fr)

                # brms (temperature/buoyancy rms) in lam/turb regions
                avg_lam_brms = np.append(avg_lam_brms,
                    db.rms(brms[i, idx_lam[0], idx_lam[1], idx_lam[2]]))
                avg_turb_brms = np.append(avg_turb_brms,
                    db.rms(brms[i, idx_turb[0], idx_turb[1], idx_turb[2]]))

                # uz rms in lam/turb regions
                avg_uzrms_lam = np.append(avg_uzrms_lam,
                    db.rms(uzrms[i, idx_lam[0], idx_lam[1], idx_lam[2]]))
                avg_uzrms_turb = np.append(avg_uzrms_turb,
                    db.rms(uzrms[i, idx_turb[0], idx_turb[1], idx_turb[2]]))

                vlam = np.append(vlam, len(idx_lam[0]) / np_tot)
                vturb = np.append(vturb, len(idx_turb[0]) / np_tot)

            del uzrms
            del brms
            del vortz

        # remove repeated timesteps and sort chronologically
        t, indices = np.unique(t_tot, return_index=True)
        del t_tot

        avg_lam_brms = avg_lam_brms[indices]
        avg_turb_brms = avg_turb_brms[indices]

        avg_uzrms_lam = avg_uzrms_lam[indices]
        avg_uzrms_turb = avg_uzrms_turb[indices]

        vlam = vlam[indices]
        vturb = vturb[indices]

        Nt = len(t)
        tavg_file.write('# Re = {:6.2f}, Pe = {:6.2f}, B = {:6.2f}\n'
            .format(Re, Pe, B))

        # temporal averages
        tidx = np.where((lb <= t) & (t <= ub))

        avg_lam_brms, err_lam_brms = db.discounted_tavg(avg_lam_brms, vlam, tidx)
        avg_turb_brms, err_turb_brms = db.discounted_tavg(avg_turb_brms, vturb, tidx)

        avg_uzrms_lam, err_uzrms_lam = db.discounted_tavg(avg_uzrms_lam, vlam, tidx)
        avg_uzrms_turb, err_uzrms_turb = db.discounted_tavg(avg_uzrms_turb, vturb, tidx)
        avg_uzrms_lam  /= uhrms_out
        err_uzrms_lam  /= uhrms_out
        avg_uzrms_turb /= uhrms_out
        err_uzrms_turb /= uhrms_out

        vturb_avg, vturb_err = db.tavg(vturb, tidx)

        # shared scalar quantities come from the OUT file time averages;
        # lam/turb decompositions and vturb come from the simdat extraction
        tavg_file.write(tavg_fmt_str.format(
            Re, Pe, B, lb, ub,
            uxrms_out, uxrms_err, uyrms_out, uyrms_err,
            uzrms_out, uzrms_err,
            vortzrms_out, vortzrms_err,
            brms_out, brms_err,
            tdisp_out, tdisp_err,
            mdisp_out, mdisp_err,
            vturb_avg, vturb_err,
            avg_uzrms_turb, err_uzrms_turb,
            avg_uzrms_lam, err_uzrms_lam,
            avg_turb_brms, err_turb_brms,
            avg_lam_brms, err_lam_brms,
            uhrms_out, uhrms_err))
        tavg_file.write('\n')

    tavg_file.write('\n\n\n')

tavg_file.close()
