import sys

sys.path.append('../')

import numpy as np
from classy import Class

import initial
import bispec


NAME = "Tree"


if __name__ == '__main__':
    # Measurement parameters
    ell1 = 0
    ell2 = 0
    ELL  = 0

    z = 0.51

    # Cosmological parameters (BOSS fiducial)
    h = 0.676
    Omega_m = 0.31
    Omega_b = 0.022 / h**2
    sigma8_0 = 0.824
    n_s = 0.96
    ln_10to10_A_s = 3.094

    # Bias parameters
    b1 = 1.0
    b2 = 0.0
    b3 = 0.0

    bK2 = 0.0
    bK3 = 0.0
    bDK = 0.0
    bO = 0.0

    # Cosmological model
    params_cosmo = {
        'output': 'tCl mPk',
        'h': h,
        'omega_b': Omega_b * h**2,
        'omega_cdm': (Omega_m - Omega_b) * h**2,
        'n_s': n_s,
        'ln10^{10}A_s': ln_10to10_A_s,
        'tau_reio': 0.0826026,
        'N_ncdm': 1.,
        'm_ncdm': 0.06,
        'N_ur': 2.,
        'z_max_pk': 3.,
        'P_k_max_h/Mpc': 50.,
    }

    cosmo = Class()
    cosmo.set(params_cosmo)
    cosmo.compute()

    # Initial cosmological model
    initial_cosmo = initial.InputPowerSpectrum(
        z, cosmo, params_fid=params_cosmo
    )
    initial_cosmo.calcFiducialHubbleAndDiameterDistance()
    initial_cosmo.calcMatterPowerSpectrum()
    # initial_cosmo.calcPrimordialPowerSpectrum(ln_10to10_A_s, n_s)

    k_in, pk_in = initial_cosmo.getMatterPowerSpectrum()
    k_pri, pk_pri = initial_cosmo.getPrimordialPowerSpectrum()
    D_of_z = initial_cosmo.getGrowthFactor()
    f_of_z = initial_cosmo.getGrowthRate()
    alpha_perp = initial_cosmo.getAlphaPerp()
    alpha_parallel = initial_cosmo.getAlphaParallel()
    sigma8 = initial_cosmo.getSigma8z(sigma8_0)
    sigma8_norm = initial_cosmo.getSigma8ForNormalization()

    params = {
        'alpha_perp': alpha_perp,
        'alpha_parallel': alpha_parallel,
        'sigma8': sigma8,
        'fz': f_of_z,
        'b1': b1,
        'b2': b2,
        'b3': b3,
        'bK2': bK2,
        'bK3': bK3,
        'bDK': bDK,
        'bO':  bO,
    }

    # Bispectrum prediction
    bispectrum = bispec.ClassBiSpectrum()
    bispectrum.set_params(params)
    bispectrum.set_input_pk(k_in, pk_in)
    bispectrum.set_normalization(sigma8_norm)

    # Binning
    kbin = np.linspace(0.01, 0.2, 20)
    rbin = np.linspace(0.5, 200, 400)

    # Transformations
    bk_dict = bispectrum.calc_B(
        name=NAME, kbin=kbin,
        ell1=ell1, ell2=ell2, ELL=ELL,
        flag_3pcf=True
    )
    zeta_dict = bispectrum.calc_B_to_3PCF(bk_dict, rbin=rbin)
    bk_reverse_dict = bispectrum.calc_3PCF_to_B(zeta_dict, kbin=kbin)

    # Output
    N_fft = len(bk_dict["kbin1"])
    X = np.array([
        bk_dict["kbin1"].reshape(N_fft**2),
        bk_dict["kbin2"].reshape(N_fft**2),
        bk_dict["B"].reshape(N_fft**2)
    ]).T
    np.savetxt(
        "results/bk{:d}{:d}{:d}_{}.dat".format(ell1, ell2, ELL, NAME),
        X, fmt="%.7e"
    )

    N_rbin = len(rbin)
    X = np.array([
        zeta_dict["rbin1"].reshape(N_rbin**2),
        zeta_dict["rbin2"].reshape(N_rbin**2),
        zeta_dict["3pcf"].reshape(N_rbin**2)
    ]).T
    np.savetxt(
        "results/zeta{:d}{:d}{:d}_{}.dat".format(ell1, ell2, ELL, NAME),
        X, fmt="%.7e"
    )

    X = np.array([
        bk_reverse_dict["kbin1"].reshape(N_fft**2),
        bk_reverse_dict["kbin2"].reshape(N_fft**2),
        bk_reverse_dict["B"].reshape(N_fft**2)
    ]).T
    np.savetxt(
        "results/bk%d%d%d_%s_from3PCF.dat" % (ell1, ell2, ELL, NAME),
        X, fmt="%.7e"
    )



