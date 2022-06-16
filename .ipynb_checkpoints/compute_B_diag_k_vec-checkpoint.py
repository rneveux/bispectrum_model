from classy import Class
import os
import pickle as pkl
import numpy as np
from triumvirate.fieldmesh import record_binned_vectors
from triumvirate.dataobjs import Binning
from powerbispectrum import ComputePowerBiSpectrum

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

cl = ComputePowerBiSpectrum(params_cosmo, z, diag=True)
cl.initial_power_spectrum()

max_a_post = {'b1': 1.9943274e+00,
 'b2': -2.8493848e-01,
 'bG2': -2.9738286e-01,
 'bGamma3': 6.8452853e-01,
 'c0': 1.1104226e+01,
 'c2pp': 3.0399105e+01,
 'c1': 5.9360800e+00,
 'c2': -1.6427413e-01,
 'ch': 3.9950281e+02,
 'Pshot': 5.6773871e-01,
 'a0': -4.3189287e+00,
 'Bshot': 2.9556961e+00}

nb_bin = 8
binning = Binning('fourier', 'lin', bin_min=0.005, bin_max=0.085, num_bins=nb_bin)

kmodes = record_binned_vectors(binning, boxsize=2000, ngrid=768)

ks = []
for i in range(nb_bin):
    mask = kmodes['index']==i
    ks.append(np.array((kmodes[mask]['vecx'],kmodes[mask]['vecy'],kmodes[mask]['vecz']), dtype='double').T)
    
np.save('/home/rneveux/test/test_kbin/bk/k_vector', ks)
tot_bk_disc = []
for i in range(nb_bin):
    for part in ['tree','SN']:
        path = f'/home/rneveux/test/test_kbin/bk/k_vector_kbin{i}_{part}.npy'
        if os.path.isfile(path): continue
        cl.calc_B_diag_k_vector(
                ks[i].copy(order='C'), los=np.array([0,0,1], dtype='double'),
                alpha_perp=1, alpha_parallel=1, b1=max_a_post['b1'], b2=max_a_post['b2'], bG2=max_a_post['bG2'], c1=max_a_post['c1'],
                c2=max_a_post['c2'],
                knl=.3,
                Pshot=max_a_post['Pshot'], Bshot=max_a_post['Bshot'],
                ks=.05,
                part = part,
            )
        with open(path,'wb') as f:
            pkl.dump(cl.BK['K'], f, protocol=4)
