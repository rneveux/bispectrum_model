import os
import sys
from classy import Class
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
import argparse
import classy


parser = argparse.ArgumentParser(description='Computation of power kernels using CLASS-PT')
parser.add_argument('-redshift', type=float, required=True)
parser.add_argument('-cosmo', type=int, help='number of cosmology in Jamie file', required=True)
parser.add_argument('-directory', type=str, help='save directory', required=True)
parser.add_argument('-kmin', type=float, default=.005, required=False)
parser.add_argument('-kmax', type=float, default=.3025, required=False)
parser.add_argument('-Omfid', type=float, default=0, required=False)
cmdline = parser.parse_args()

k=np.arange(cmdline.kmin, cmdline.kmax, .0025)

z=cmdline.redshift
all_kernels =   [
                        15, 21, 16, 22, 17, 23, 1, 30, 31, 32, 33,
                         4, 5, 11, 7, 8, 13,
                        18, 24, 19, 25, 26, 34, 35, 36, 37, 12, 9,
                        20, 27, 28, 29, 38, 39,
                        40, 41, 42
                        ]

'''cosmo = np.load('/home/rneveux/bispectrum/theory/cosmologies/lnAs/eft_cosmologies_noDQ1.npy')[cmdline.cosmo]

h = cosmo[2]
omega_cdm = cosmo[0]
omega_b = cosmo[1]
lnAs = cosmo[3]
ns = cosmo[4]

params_cosmo = {'k_output_values':2.0,'output': 'tCl mPk','z_max_pk': 3.,'P_k_max_h/Mpc': 50., 'z_pk':z,
            'omega_cdm':omega_cdm,'omega_b':omega_b,'h':h,
            'ln10^{10}A_s':lnAs,
            'n_s':ns,
               'N_ur':2.0328, 'N_ncdm':1, 'omega_ncdm':0.0006442}'''

cosmo = np.load('/home/rneveux/bispectrum/theory/cosmologies/forFFcomp/eft_cosmologies_h_lnAs_Ocdm.npy')[cmdline.cosmo]

h = cosmo[0]
omega_cdm = cosmo[2]
lnAs = cosmo[1]

params_cosmo = {
    'output': 'mPk',
    'P_k_max_h/Mpc': 50.,
    'h': h,
    'omega_b': 0.022,
    'omega_cdm': omega_cdm,
    'ln10^{10}A_s': lnAs,
    'n_s': 0.965,
    'N_ur': 2.038,
    'N_ncdm': 1,
    'm_ncdm': 0.1,
    'tau_reio': 0.0544,
    'z_max_pk': 3,
}

c = Class()
c.set(params_cosmo)
if cmdline.Omfid:
    c.set({'output':'mPk',
        'non linear':'PT',
        'IR resummation':'Yes',
        'Bias tracers':'Yes',
        'cb':'Yes',
        'RSD':'Yes',
        'AP':'Yes',
        'Omfid':cmdline.Omfid
       })
else:
    c.set({'output':'mPk',
        'non linear':'PT',
        'IR resummation':'Yes',
        'Bias tracers':'Yes',
        'cb':'Yes',
        'RSD':'Yes',
        'AP':'No',
       })

c.compute()

kh = k*h

fz = c.scale_independent_growth_factor_f(z)
c.initialize_output(kh, z, len(kh))
c_mult = c.get_pk_mult(kh, z, len(kh))

for i in all_kernels[:-3]:
        if i in [11,12,13]:
                c_mult[i] *= h
        else: c_mult[i] *= h**3
        if i==13:
                c_mult[40] = c_mult[i] * fz**2 * k**2
                c_mult[41] = c_mult[i] * fz**3 * k**2
                c_mult[42] = c_mult[i] * fz**4 * k**2

c_mult = c_mult[all_kernels]

np.savetxt(os.path.join(cmdline.directory,f'pk_kernels_{cmdline.cosmo}.txt'),c_mult)