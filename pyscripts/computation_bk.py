import os
from classy import Class
import numpy as np
import argparse
from powerbispectrum import ComputePowerBiSpectrum

parser = argparse.ArgumentParser(description='Computation of bispectrum model')
parser.add_argument('-ells', type=str, default='000', required=False)
cmdline = parser.parse_args()

h_fid = .6736
omega_b_fid = .02237
omega_cdm_fid = .12
n_s_fid = .9649
A_s_fid = 2.083e-9

params_cosmo = {
            'output': 'mPk',
            'h': h_fid,
            'omega_b': omega_b_fid,
            'omega_cdm': omega_cdm_fid,
            'n_s': n_s_fid,
            'A_s': A_s_fid,
            'tau_reio': 0.0544,
            'N_ncdm': 1.,
            'm_ncdm': 0.06,
            'N_ur': 2.0328,
            'z_max_pk': 3.,
            'P_k_max_h/Mpc': 50.,
            }
z = .8

kbin = np.logspace(-4, np.log10(0.5), num=256, base=10)

cl = ComputePowerBiSpectrum(params_cosmo, z, diag=False)
cl.initial_power_spectrum()

ELL = int(cmdline.ells[2])
ell1 = int(cmdline.ells[0])
ell2 = int(cmdline.ells[1])

max_a_post = {'b1': 1.9943274e+00,
 'b2': -2.8493848e-01,
 'bG2': -2.9738286e-01,
 'c1': 0,
 'c2': 0,
 'ch': 3.9950281e+02,
 'Pshot': 5.6773871e-01,
 'Bshot': 2.9556961e+00,
             }

cl.calc_B(
        kbin, ell1, ell2, ELL,
        alpha_perp=1, alpha_parallel=1, b1=max_a_post['b1'], b2=max_a_post['b2'], bG2=max_a_post['bG2'], c1=max_a_post['c1'], c2=max_a_post['c2'], 
        Pshot=max_a_post['Pshot']*1000, Bshot=max_a_post['Bshot']*1000,
        integrand='tree',
        ks=.05,
    )
tree = cl.BK
cl.calc_B(
        kbin, ell1, ell2, ELL,
        alpha_perp=1, alpha_parallel=1, b1=max_a_post['b1'], b2=max_a_post['b2'], bG2=max_a_post['bG2'],
        c1=max_a_post['c1'], c2=max_a_post['c2'], knl=.3,
        Pshot=max_a_post['Pshot']*1000, Bshot=max_a_post['Bshot']*1000,
        integrand='SN',
        ks=.05,
    )
SN = cl.BK
BK = tree['K']+SN['K']

np.savetxt(f'/home/rneveux/to_mike/model_example/klog/B_{ell1}{ell2}{ELL}.txt',BK)