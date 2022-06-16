#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os

class ClassCovariancePP():

    def __init__(self):
        self.initialize()

    def initialize(self, patchy=True):
        self.k_temp = np.zeros(1)
        self.P_temp = np.zeros(1)

        self.sigma8_norm = 1.0

        self.alpha_perp = 1.0
        self.alpha_parallel = 1.0
        self.sigma8 = 0.0
        self.fz = 0.0

        self.b1 = 0.0
        self.b2 = 0.0
        self.b3 = 0.0

        self.bK2 = 0.0
        self.bK3 = 0.0
        self.bDK = 0.0
        self.bO = 0.0
        
        self.params_cosmo = {}
        
        if patchy:
            self.params_cosmo['h'] = .678
            self.params_cosmo['Omega_b'] = .048
            self.params_cosmo['Omega_m'] = .307
            self.params_cosmo['Tcmb'] = 2.7255
            self.params_cosmo['n_s'] = .961

    def set_input_pk(self, k_in, P_in):
        self.k_temp = k_in
        self.P_temp = P_in

    def set_normalization(self, sigma8_norm):
        self.sigma8_norm = sigma8_norm

    def set_params(self, params):
        self.alpha_perp = params["alpha_perp"]
        self.alpha_parallel = params["alpha_parallel"]
        self.sigma8 = params["sigma8"]
        self.fz = params["fz"]

        self.b1 = params["b1"]
        self.b2 = params["b2"]
        self.b3 = params["b3"]

        self.bK2 = params["bK2"]
        self.bK3 = params["bK3"]
        self.bDK = params["bDK"]
        self.bO = params["bO"]

    def select_C(self, name):

        n_kbin = len(self.kbin)

        if name == "cov_PP_G":
            return hitomipy.integrand_cov_PP_G_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)
        
        elif name == "cov_PP_G_NL":
            return hitomipy.integrand_cov_PP_G_NL_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.sigma2_perp, self.sigma2_para,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG":
            return hitomipy.integrand_cov_PP_NG_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.b2, self.b3, self.bK2, self.bK3, self.bDK, self.bO,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_b1":
            return hitomipy.integrand_cov_PP_NG_b1_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_b2":
            return hitomipy.integrand_cov_PP_NG_b2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_bK2":
            return hitomipy.integrand_cov_PP_NG_bK2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_b2_b2":
            return hitomipy.integrand_cov_PP_NG_b2_b2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_b2_bK2":
            return hitomipy.integrand_cov_PP_NG_b2_bK2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_bK2_bK2":
            return hitomipy.integrand_cov_PP_NG_bK2_bK2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_b3":
            return hitomipy.integrand_cov_PP_NG_b3_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_bK3":
            return hitomipy.integrand_cov_PP_NG_bK3_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_bDK":
            return hitomipy.integrand_cov_PP_NG_bDK_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_bO":
            return hitomipy.integrand_cov_PP_NG_bO_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        
        elif name == "cov_PP_NG_BeatCoupling":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.b2, self.b3, self.bK2, self.bK3, self.bDK, self.bO,
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_BeatCoupling_b1":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_b1_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

      
        elif name == "cov_PP_NG_BeatCoupling_b2":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_b2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_bK2":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_bK2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_b2_b2":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_b2_b2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_b2_bK2":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_b2_bK2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_bK2_bK2":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_bK2_bK2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_b3":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_b3_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_bK3":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_bK3_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_bDK":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_bDK_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_BeatCoupling_bO":
            return hitomipy.integrand_cov_PP_NG_BeatCoupling_bO_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_LocalMean":
            return hitomipy.integrand_cov_PP_NG_LocalMean_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.b2, self.b3, self.bK2, self.bK3, self.bDK, self.bO,
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_LocalMean_NL":
            return hitomipy.integrand_cov_PP_NG_LocalMean_NL_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.b2, self.b3, self.bK2, self.bK3, self.bDK, self.bO,
                    self.sigma2_perp, self.sigma2_para, 
                    self.deltaK, self.nmean, self.volume)

        elif name == "cov_PP_NG_LocalMean_NL_Sigma2B":
            return hitomipy.integrand_cov_PP_NG_LocalMean_NL_Sigma2B_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.b2, self.b3, self.bK2, self.bK3, self.bDK, self.bO,
                    self.sigma2_perp, self.sigma2_para, 
                    self.deltaK, self.nmean, self.volume, self.sigma2_b)

        elif name == "cov_PP_NG_LocalMean_NL_b1":
            return hitomipy.integrand_cov_PP_NG_LocalMean_NL_b1_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.sigma2_perp, self.sigma2_para, 
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_LocalMean_NL_b2":
            return hitomipy.integrand_cov_PP_NG_LocalMean_NL_b2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.sigma2_perp, self.sigma2_para, 
                    self.deltaK, self.nmean, self.volume)


        elif name == "cov_PP_NG_LocalMean_NL_bK2":
            return hitomipy.integrand_cov_PP_NG_LocalMean_NL_bK2_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.sigma2_perp, self.sigma2_para, 
                    self.deltaK, self.nmean, self.volume)

        else:
            print("select_C: ERROR")

        return 0.0

    def select_ndim(self, name):

        name_dim3 = []
        name_dim3.append("cov_PP_NG")
        name_dim3.append("cov_PP_NG_b1")
        name_dim3.append("cov_PP_NG_b2")
        name_dim3.append("cov_PP_NG_bK2")
        name_dim3.append("cov_PP_NG_b2_b2")
        name_dim3.append("cov_PP_NG_b2_bK2")
        name_dim3.append("cov_PP_NG_bK2_bK2")
        name_dim3.append("cov_PP_NG_b3")
        name_dim3.append("cov_PP_NG_bK3")
        name_dim3.append("cov_PP_NG_bDK")
        name_dim3.append("cov_PP_NG_bO")


        name_dim6 = []
        name_dim6.append("cov_PP_NG_BeatCoupling")
        name_dim6.append("cov_PP_NG_BeatCoupling_b1")
        name_dim6.append("cov_PP_NG_BeatCoupling_b2")
        name_dim6.append("cov_PP_NG_BeatCoupling_bK2")
        name_dim6.append("cov_PP_NG_BeatCoupling_b2_b2")
        name_dim6.append("cov_PP_NG_BeatCoupling_b2_bK2")
        name_dim6.append("cov_PP_NG_BeatCoupling_bK2_bK2")
        name_dim6.append("cov_PP_NG_BeatCoupling_b3")
        name_dim6.append("cov_PP_NG_BeatCoupling_bK3")
        name_dim6.append("cov_PP_NG_BeatCoupling_bDK")
        name_dim6.append("cov_PP_NG_BeatCoupling_bO")

        name_dim6.append("cov_PP_NG_LocalMean")
        name_dim6.append("cov_PP_NG_LocalMean_NL")
        name_dim6.append("cov_PP_NG_LocalMean_NL_Sigma2B")
        name_dim6.append("cov_PP_NG_LocalMean_NL_b1")
        name_dim6.append("cov_PP_NG_LocalMean_NL_b2")
        name_dim6.append("cov_PP_NG_LocalMean_NL_bK2")

        if name in name_dim3:
            return 3
        elif name in name_dim6:
            return 6
        else:
            return 2

    def check_flag_window(self, name, flag_window, volume):

        flag = 0
        name_window = []

        if name in name_window:
            if flag_window and volume > 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_BAO(self, name, flag_BAO, sigma8_fid, fz_fid):

        flag = 0
        name_BAO = []
        name_BAO.append("cov_PP_G_NL")
        
        if name in name_BAO:
            if flag_BAO and sigma8_fid > 0.0 and fz_fid > 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_nl(self, name, flag_BAO, sigma8_fid, fz_fid, flag_nl, kbin_nl, pk_nl_0, pk_nl_2):

        flag = 0
        name_nl = []
        name_nl.append("cov_PP_G_NL")
        name_nl.append("cov_PP_NG_LocalMean_NL")
        name_nl.append("cov_PP_NG_LocalMean_NL_Sigma2B")
        name_nl.append("cov_PP_NG_LocalMean_NL_b1")
        name_nl.append("cov_PP_NG_LocalMean_NL_b2")
        name_nl.append("cov_PP_NG_LocalMean_NL_bK2")

        if name in name_nl:
            if flag_BAO and sigma8_fid > 0.0 and fz_fid > 0.0 and flag_nl and np.sum(kbin_nl) > 0.0 and np.sum(pk_nl_0) > 0.0 and np.sum(pk_nl_2) > 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_cov(self, flag_cov, volume, nmean, deltaK):

        flag = 0
        if flag_cov and volume > 0.0 and nmean > 0.0 and deltaK > 0.0:
            flag = 0
        else:
            flag = -1
        return flag


    def Integrand_cov_PP(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])

        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]
            self.ndim = ndim[0]
            self.ncomp = ncomp[0]

        self.select_C(self.name)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]

        return 0

    def calc_cov_PP( self,  name,
                     kbin=np.linspace(0.01, 0.2, 20), ELL = 0, ELL_dash = 0,
                     flag_BAO=False, sigma8_fid=-1.0, fz_fid=- 1.0,
                     flag_nl=False, kbin_nl=np.zeros(20), pk_nl_0=np.zeros(20), pk_nl_2=np.zeros(20), 
                     flag_window=False, volume=500.0**3,
                     flag_cov=True, nmean = 3.0e-4, deltaK = 0.01, sigma2_b = 0.0
                     ):

        self.sigma2_b = sigma2_b

        ## flags ##

        (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
        
        output_dict_ini = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "cov_PP": np.zeros((len(kbin), len(kbin))),
                "ELL": ELL,
                "ELL_dash": ELL_dash,
                "flag_BAO": flag_BAO,
                "flag_nl": flag_nl,
                "flag_window": flag_window}

        check_bao = self.check_flag_BAO(name, flag_BAO, sigma8_fid, fz_fid)
        check_nl = self.check_flag_nl(name, flag_BAO, sigma8_fid, fz_fid, flag_nl, kbin_nl, pk_nl_0, pk_nl_2)
        check_window = self.check_flag_window(name, flag_window, volume)
        check_cov = self.check_flag_cov(flag_cov, volume, nmean, deltaK)

        if check_bao < 0:
            print("FLAG_BAO: ERROR")
            return output_dict_ini

        if check_nl < 0:
            print("FLAG_NL: ERROR")
            return output_dict_ini

        if check_window < 0:
            print("FLAG_WINDOW ERROR")
            return output_dict_ini

        if check_cov < 0:
            print("FLAG_COV ERROR")
            return output_dict_ini

        if flag_window:
            self.volume = volume

        if flag_cov:
            self.volume = volume
            self.nmean  = nmean
            self.deltaK = deltaK

        ## type of powerspectra, e.g,, "Tree", "NoWiggle" and etc. ##
        self.name = name

        ## set kbin ##
        self.kbin = kbin

        ## set multipole indices ##
        self.ELL = ELL
        self.ELL_dash = ELL_dash

        ## initialization ##
        hitomipy.initializeInputPowerSpectrum_py()

        ## read linear power spectrum ##
        hitomipy.readInputPowerSpectrum_py(
            self.k_temp, self.P_temp, len(self.k_temp))

        ## read non-linear power spectra ##
        if flag_nl:

            kbin_nl = np.array(kbin_nl)
            pk_nl_0 = np.array(pk_nl_0)
            pk_nl_2 = np.array(pk_nl_2)

            hitomipy.readNonLinearPowerSpectrum_py(kbin_nl, pk_nl_0, pk_nl_2, len(kbin_nl))


        ## normalization ##
        hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
        hitomipy.calcNormalizationNoWiggle_py(1.0, self.params_cosmo['h'], self.params_cosmo['Omega_b'], 
                                              self.params_cosmo['Omega_m'], self.params_cosmo['Tcmb'], 
                                              self.params_cosmo['n_s'])

        ## sigma2_perp and sigma2_para ##
        if flag_BAO:
            self.sigma2_perp = hitomipy.calcSigma_dd_py(sigma8_fid)
            self.sigma2_para = (1.0 + fz_fid) * (1.0 + fz_fid) * self.sigma2_perp

        ## compute power spectra ##
        NDIM = self.select_ndim(self.name)
        NCOMP = len(self.kbin)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

##########################

        if self.name == "cov_PP_G" or self.name == "cov_PP_G_NL":
            AA = pycuba.Cuhre(
                    self.Integrand_cov_PP,
                    NDIM,
                    ncomp=NCOMP,
                    key=0,
                    verbose=0 | 4)

            result = np.zeros(NCOMP)
           
            for i in range(NCOMP):
                result[i] = AA["results"][i]['integral']
            
            result = np.diag(result)

        else:

            AA = []
            for i in range(NCOMP):
                print("k1 = ", self.kbin[i], "h/Mpc")
                self.kmag1 = self.kbin[i]
                if NDIM <= 3:

                    AA.append(
                        pycuba.Cuhre(
                            self.Integrand_cov_PP, 
                            NDIM,
                            ncomp=NCOMP,
                            key=0, 
                            verbose=0 | 4)["results"]
                    )

                else:

                    NNEW = 50000
                    NMIN = 2
                    FLATNESS = 50
                    MAXEVAL = 50000
                    
                    AA.append(
                        pycuba.Suave(
                            self.Integrand_cov_PP, 
                            NDIM, 
                            NNEW, NMIN, FLATNESS, 
                            ncomp=NCOMP,
                            maxeval = MAXEVAL,
                            verbose=0 | 4)["results"]
                    )
                    
            result = np.zeros((NCOMP,NCOMP))
            
            for i in range(NCOMP):
                for j in range(NCOMP):
                    result[i,j] = AA[i][j]["integral"]
            
        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        ## output dict ##
        
        output_dict = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "cov_PP": result,
                "ELL": self.ELL,
                "ELL_dash": self.ELL_dash,
                "flag_BAO": flag_BAO,
                "flag_nl": flag_nl,
                "flag_window": flag_window}

        return output_dict


