import os
import sys
sys.path.append(os.environ['AURORA_REPOS'])

import numpy as np
import pandas as pd
import math
from flagships.gs_solver.fs_flagships import RunFlagships, LambdaAndBetaPol1Params, FluxConstrainedLambdaAndBetaPol1Params
from flagships.gs_solver.fs_curves import PiecewiseCurve
from flagships.gs_solver.fs_profile import LambdaCurveProfile
from flagships.gs_solver.fs_gs_solver import FSFluxNormalizedSolver
from flagships.gs_solver.fs_init_guess import PrevRunInitialGuess

INPUT_DATA_PATH = 'input_data'
EXT_PSI_DIR = os.path.join(INPUT_DATA_PATH, 'ext_psi_files')
GEOM_DIR = os.path.join(INPUT_DATA_PATH, 'geom_files')
SCALARS_TABLE_PATH = os.path.join(INPUT_DATA_PATH, 'M3_1005D_0.10Wb-n2e19-lampeak-eps0.6-t300-r0.csv')

EXT_PSI_PREFEX = 'ext_psi_'
GEOM_PREFIX = 'geom_'

PSI_AXIS_SETPOINT = 100e-3 #Wb

'''
lambda = lambda0*(1-psibar^2), lambda0 = 150mWb
get Ishaft from scalars
Te = 50 + 250*psibar
ne = 4e19*( 1 - 0.67*psibar - 0.33*psibar​^4 )

'''

num_psibar = 1000
psibar = np.linspace(0, 1, num_psibar)

lambda_0_guess = 8
lambda_values = lambda_0_guess*(1 - psibar**2)
lambda_curve = PiecewiseCurve(psibar, lambda_values) 

ext_psi_scale_loc= [0.9999, 0.0] # Shouldn't matter as just setting scale to 1
out_folder = os.path.join('out')
soakscale = 0 
psi_lim_in = 0

scalars_df = pd.read_csv(SCALARS_TABLE_PATH)

prev_file = None #'csim_027c_equil/lm26_test_0.002296.hdf5'
all_run_params = []

for i_row in range(1, len(scalars_df)):
    row = scalars_df.iloc[i_row]
    time = row["time(s)"]
    timestep_str = f'{time:.6f}'
    ext_psi_file = os.path.join(EXT_PSI_DIR, EXT_PSI_PREFEX + timestep_str + '.csv')
    geom_file = os.path.join(GEOM_DIR, GEOM_PREFIX + timestep_str + '.csv')
    out_file = f'm3v3_{time:.6f}.hdf5'

    res_factor = (5 / (np.e - 1)) * (np.exp((time / 0.002910)**6) - 1)
    res_factor = max(res_factor, 1) 

    if i_row == 0:
        res_factor *= 2

    expected_opoint_r = 0.5 - time * 120
    res = res_factor * int(30 * math.sqrt(0.35 / (expected_opoint_r - 0.05)))
    res = min(200, res)

    print(time, res)

    run_params = FluxConstrainedLambdaAndBetaPol1Params(out_folder=out_folder, out_file=out_file, geom_file=geom_file,
                                        ext_psi_scale_loc=ext_psi_scale_loc, gun_femm_file=ext_psi_file, flux_setpoint=PSI_AXIS_SETPOINT,
                                        soak_file=ext_psi_file, soakscale=soakscale, lambda_curve=lambda_curve,
                                        Ishaft_in=1e6*row['Ishaft(MA)'], Ipl_in=1e6*row['Ipl'], beta_pol1_in=row['betap1'],
                                        psi_lim_in=psi_lim_in, expected_opoint=[expected_opoint_r, 0], use_csharp_solver=True,
                                        mesh_resolution=res)
    
    if prev_file is not None and os.path.exists(prev_file):
        run_params.init_guess = PrevRunInitialGuess(prev_file)

    prev_file = os.path.join(out_folder, out_file)

    try:
        RunFlagships(run_params, write_history=True)
    except Exception as e:
        print(e)