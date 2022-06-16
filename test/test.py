import numpy as np
from matplotlib import pyplot as plt

from classy import Class

import initial
import bispec

# Cosmological parameters
h = .676
Omega_m = .31
omega_b = .022
sigma8_0 = .8
n_s = .96
m_ncdm = .06

# CLASS cosmo
params_cosmo = {
    'output': 'tCl mPk',
    'h': h,
    'omega_b': omega_b,
    'omega_cdm': Omega_m*h**2 - omega_b,
    'n_s': n_s,
    'm_ncdm': m_ncdm,
    'N_ncdm': 1.,
    'ln10^{10}A_s': 3.094,
    'tau_reio': 0.0826026,
    'N_ur': 2.,
    'z_max_pk': 3.,
    'P_k_max_h/Mpc': 50.,
    }

cosmo = Class()
cosmo.set(params_cosmo)
cosmo.compute()

#redshift
z = 1.48

# Initial power spectra
initial_cosmo = initial.InputPowerSpectrum(z, cosmo, params_fid=params_cosmo)
initial_cosmo.calcMatterPowerSpectrum()
initial_cosmo.calcPrimordialPowerSpectrum(params_cosmo['ln10^{10}A_s'], n_s)
initial_cosmo.calcFiducialHubbleAndDiameterDistance()

#hardcode: k_in = k_pri = np.logspace(np.log(2e-5), np.log(50), 500, base=np.e)
k_in, pk_in = initial_cosmo.getMatterPowerSpectrum()
k_pri, pk_pri = initial_cosmo.getPrimordialPowerSpectrum()

alpha_perp = initial_cosmo.getAlphaPerp()
alpha_parallel = initial_cosmo.getAlphaParallel()
sigma8 = initial_cosmo.getSigma8z(sigma8_0)
sigma8_norm = initial_cosmo.getSigma8ForNormalization()
D_of_z = initial_cosmo.getGrowthFactor()
f_of_z = initial_cosmo.getGrowthRate()

params = {
    'alpha_perp': alpha_perp,
    'alpha_parallel': alpha_parallel,
    'sigma8': sigma8,
    'fz': f_of_z,
    'b1': 2,
    'b2': 0,
    'b3': 0,
    'bK2': 0,
    'bK3': 0,
    'bDK': 0,
    'bO':  0,
    'c1': 1,
    }

kbin = np.linspace(0.01, 0.3, 5)

#Set bispectrum
bispectrum = bispec.ClassBiSpectrum()
bispectrum.set_params(params)
bispectrum.set_input_pk(k_in, pk_in)
bispectrum.set_normalization(sigma8_norm)

for name in ['Tree_NoWiggle']:
    if name is 'Tree_DampIvanov': flag=True
    else: flag=False
    B = bispectrum.calc_B(name='Tree_BAO_b1_b1_f',flag_BAO=True, sigma8_fid=params['sigma8'], fz_fid=params['fz'],
            kbin=kbin,ell1=0,ell2=0,ELL=0)
    print(B)

