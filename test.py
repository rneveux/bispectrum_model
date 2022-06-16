import os
import sys

from classy import Class
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
import argparse

import initial
from powerbispectrum import ComputePowerBiSpectrum

cosmo = np.load('/home/rneveux/bispectrum/theory/eft_cosmologies_noDQ1.npy')[5000]

params_cosmo = {'output': 'tCl mPk','z_max_pk': 3.,'P_k_max_h/Mpc': 50.,
                'omega_cdm':cosmo[0],'omega_b':cosmo[1],'h':cosmo[2],
                'ln10^{10}A_s':cosmo[3],'n_s':cosmo[4]}
z = .8
k = np.arange(0.005,0.2025,.0025)

kernels = ['fnlloc_b1_b1_b1', 'fnlloc_b1_b1_f', 'fnlloc_b1_f_f', 'fnlloc_f_f_f',
            'fnlequi_b1_b1_b1', 'fnlequi_b1_b1_f', 'fnlequi_b1_f_f', 'fnlequi_f_f_f',
            'fnlortho_b1_b1_b1', 'fnlortho_b1_b1_f', 'fnlortho_b1_f_f', 'fnlortho_f_f_f'
           ]

c=ComputePowerBiSpectrum(params_cosmo,z)
c.initial_power_spectrum()

c.kernel_computation_Bk(kernels[0],k,0,0,0,integrand='PNG')

