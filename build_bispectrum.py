import sys
import numpy as np

from compute_kernels import ComputeKernels


k1 = np.logspace(np.log(3.0e-4), np.log(.005), 4, base=np.e)
k2 = np.arange(.00501,.301,.0025)
k3 = np.logspace(np.log(.301), np.log(10.0), 25, base=np.e)
k = np.hstack([k1, k2, k3])

z=.56

#cosmo_nb = int(sys.argv[1])
cosmo_nb = 'test'

#cosmo = np.load('/home/rneveux/bispectrum/theory/eft_cosmologies_noDQ1.npy')[cosmo_nb]

params_cosmo = {'output': 'tCl mPk','z_max_pk': 3.,'P_k_max_h/Mpc': 50.,
				'omega_cdm':0.119,'omega_b':0.022,'h':.678,
                'sigma8':.829,'n_s':.963}
#                'omega_cdm':cosmo[0],'omega_b':cosmo[1],'h':cosmo[2],
#                'ln10^{10}A_s':cosmo[3],'n_s':cosmo[4]}

c=ComputeKernels(k,params_cosmo,z)
c.intial_power_spectrum()

all_kernels = ['b1_b1_b1', 'b1_b1_b2','b1_b1_bK2','b1_b1_f','b1_b1_b1_f','b1_b1_f_b1',
'b1_b2_f','b1_bK2_f','b1_f_f','b1_f_f_f','b2_f_f','bK2_f_f','f_f_f','f_f_f_f',
'c1_b1_b1','c1_b1_b2','c1_b1_bK2','c1_b1_f','c1_b1_b1_f','c1_b1_f_f','c1_b2_f',
'c1_bK2_f','c1_f_f','c1_f_f_f','c1_c1_b1','c1_c1_b2','c1_c1_bK2','c1_c1_f',
'c1_c1_b1_f','c1_c1_f_f','c2_b1_b1','c2_b1_b2','c2_b1_bK2','c2_b1_f','c2_b1_b1_f',
'c2_b1_f_f','c2_b2_f','c2_bK2_f','c2_f_f','c2_f_f_f','c2_c1_b1','c2_c1_b2',
'c2_c1_bK2','c2_c1_f','c2_c1_b1_f','c2_c1_f_f','c2_c2_b1','c2_c2_b2','c2_c2_bK2',
'c2_c2_f','c2_c2_b1_f','c2_c2_f_f']

bias = {'b1':2,'b2':0,'bK2':0,'c1':0,'c2':0,}

c.load_kernels(all_kernels,cosmo_nb)
c.bispectrum_model(**bias)
c.plot_bispectrum('test.pdf')