#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os

class ClassCovarianceBB():

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

        if name == "cov_BB_G":
            return hitomipy.integrand_cov_BB_G_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ell1, self.ell2, self.ELL, self.ell1_dash, self.ell2_dash, self.ELL_dash, self.kmag1, self.kmag1_dash,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        else:
            print("select_C: ERROR")

        return 0.0

    def select_ndim(self, name):

        name_dim3 = []
        name_dim3.append("cov_BB_G")
        
        name_dim6 = []

        if name in name_dim3:
            return 3
        elif name in name_dim6:
            return 6
        else:
            return 3

    def Integrand_cov_BB(self, ndim, xx, ncomp, ff, userdata):

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

    def calc_cov_BB( self,  name,
                     kbin=np.linspace(0.01, 0.2, 20), ell1 = 0, ell2 = 0, ELL = 0, ell1_dash = 0, ell2_dash = 0, ELL_dash = 0,
                     volume = 500.0**3, nmean = 3.0e-4, deltaK = 0.01
                     ):

        ## flags ##

        (kbin2_out, kbin1_out, kbin1_dash_out, kbin2_dash_out) = np.meshgrid(kbin, kbin, kbin, kbin)
        
        output_dict_ini = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "kbin1_dash": kbin1_dash_out,
                "kbin2_dash": kbin2_dash_out,
                "cov_BB": np.zeros((len(kbin), len(kbin), len(kbin), len(kbin))),
                "ell1": ell1,
                "ell2": ell2,
                "ELL": ELL,
                "ell1_dash": ell1_dash,
                "ell2_dash": ell2_dash,
                "ELL_dash": ELL_dash
                }

        self.volume = volume
        self.nmean  = nmean
        self.deltaK = deltaK

        ## type of powerspectra, e.g,, "Tree", "NoWiggle" and etc. ##
        self.name = name

        ## set kbin ##
        self.kbin = kbin

        ## set multipole indices ##
        self.ell1 = ell1 
        self.ell2 = ell2 
        self.ELL = ELL
        self.ell1_dash = ell1_dash
        self.ell2_dash = ell2_dash
        self.ELL_dash  = ELL_dash

        ## initialization ##
        hitomipy.initializeInputPowerSpectrum_py()

        ## calc. wigner 3j symbols ##
        hitomipy.setWigner3j_py()
 
        ## read linear power spectrum ##
        hitomipy.readInputPowerSpectrum_py(self.k_temp, self.P_temp, len(self.k_temp))

        ## normalization ##
        hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
        hitomipy.calcNormalizationNoWiggle_py(1.0, self.params_cosmo['h'], self.params_cosmo['Omega_b'], 
                                              self.params_cosmo['Omega_m'], self.params_cosmo['Tcmb'], 
                                              self.params_cosmo['n_s'])
            
        ## number of k-bins
        num_kbin = len(self.kbin)
        
        ## compute power spectra ##
        NDIM = self.select_ndim(self.name)
        NCOMP = num_kbin**2
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

##########################

        AA = []
        for i in range(num_kbin):
            for k in range(num_kbin):
                print("k1 = ", self.kbin[i], "h/Mpc")
                print("k1_dash = ", self.kbin[k], "h/Mpc")
                self.kmag1 = self.kbin[i]
                self.kmag1_dash = self.kbin[k]
                
                if NDIM <= 3:
        
                    AA.append(pycuba.Cuhre(
                            self.Integrand_cov_BB, 
                            NDIM,
                            ncomp=NCOMP,
                            key=0, 
                            verbose=0 | 4)["results"])
        
                else:
                
                    NNEW = 50000
                    NMIN = 2
                    FLATNESS = 50
                    MAXEVAL = 50000
                    
                    AA.append(pycuba.Suave(
                            self.Integrand_cov_BB,
                            NDIM, 
                            NNEW, NMIN, FLATNESS, 
                            ncomp=NCOMP,
                            maxeval = MAXEVAL,
                            verbose=0 | 4)["results"])
                
        result = np.zeros((num_kbin, num_kbin, num_kbin, num_kbin))
        
        for i in range(num_kbin):
            for j in range(num_kbin):
                for k in range(num_kbin):
                    for l in range(num_kbin):
                        II = i * num_kbin + k
                        JJ = j * num_kbin + l
                        result[i,j,k,l] = AA[II][JJ]["integral"]
        
        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        ## output dict ##
        
        output_dict = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "kbin1_dash": kbin1_dash_out,
                "kbin2_dash": kbin2_dash_out,
                "cov_BB": result,
                "ell1": self.ell1,
                "ell2": self.ell2,
                "ELL": self.ELL,
                "ell1_dash": self.ell1_dash,
                "ell2_dash": self.ell2_dash,
                "ELL_dash": self.ELL_dash}

        return output_dict


