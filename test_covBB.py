import numpy as np
from classy import Class
import initial
import cov_BB_diag
import cov_BB
import cov_PP

h = .6777
Omega_m = .307115
omega_b = 0.048206*h**2
sigma8_0 = .8288

# CLASS cosmo
params_cosmo = {
    'output': 'tCl mPk',
    'h': h,
    'omega_b': omega_b,
    'omega_cdm': Omega_m*h**2 - omega_b,
    'sigma8':sigma8_0,
    'tau_reio': 0.0826026,
    'N_ur': 2.,
    'z_max_pk': 3.,
    'P_k_max_h/Mpc': 50.,
    }

cosmo = Class()
cosmo.set(params_cosmo)
cosmo.compute()

initial_cosmo = initial.InputPowerSpectrum(.51, cosmo, params_fid=params_cosmo)
initial_cosmo.calcMatterPowerSpectrum()
k_in, pk_in = initial_cosmo.getMatterPowerSpectrum()
sigma8 = initial_cosmo.getSigma8z(sigma8_0)
f_of_z = initial_cosmo.getGrowthRate()
print(sigma8)
print(f_of_z)

cov = cov_BB_diag.ClassCovarianceBBDiag()
#cov = cov_BB.ClassCovarianceBB()
#cov = cov_PP.ClassCovariancePP()
cov.set_input_pk(k_in,pk_in)

cov.set_params({'sigma8':sigma8,'fz':f_of_z,'b1':1.98,'alpha_perp':1,'alpha_parallel':1,'b2':0,'b3':0,'bK2':0,'bK3':0,'bDK':0,'bO':0,})
c = cov.calc_cov_BB( 'cov_BB_G', kbin=np.arange(0.01, 0.4, .01), volume = 1.76e9, nmean = 3.26e-4, deltaK = 0.01)

np.save('test/cov_Patchy_NGC_anal_1_40_diag',c)
