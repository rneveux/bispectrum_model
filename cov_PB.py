#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os
from cov_BB_diag import ClassCovarianceBBDiag

class ClassCovariancePB(ClassCovarianceBBDiag):

    def select_C(self, name):

        n_kbin = len(self.kbin)

        if name == "cov_PB_NG_PB_diag":
            return hitomipy.integrand_cov_PB_NG_PB_diag_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ell1_dash, self.ell2_dash, self.ELL_dash,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1, self.b2, self.bK2,
                    self.deltaK, self.nmean, self.volume)

        if name == "cov_PB_NG_P5_diag":
            return hitomipy.integrand_cov_PB_NG_P5_diag_py(
                    self.xx_in, self.ndim, self.ff_out, self.ncomp,
                    self.kbin, n_kbin, self.ELL, self.ell1_dash, self.ell2_dash, self.ELL_dash, self.kmag1,
                    self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                    self.deltaK, self.nmean, self.volume)

        else:
            print("select_C: ERROR")

        return 0.0

    def select_ndim(self, name):

        name_dim3 = []
        name_dim3.append("cov_PB_NG_PB_diag")

        name_dim5 = []
        name_dim5.append("cov_PB_NG_P5_diag")

        if name in name_dim3:
            return 3
        elif name in name_dim5:
            return 5

    def Integrand_cov_PB(self, ndim, xx, ncomp, ff, userdata):

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

    def calc_cov_PB( self,  name,
                     kbin=np.linspace(0.01, 0.2, 20), ELL = 0, ell1_dash = 0, ell2_dash = 0, ELL_dash = 0,
                     volume = 500.0**3, nmean = 3.0e-4, deltaK = 0.01
                     ):

        ## flags ##

        (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
        
        output_dict_ini = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "cov_PB": np.zeros((len(kbin), len(kbin))),
                "ELL": ELL,
                "ell1_dash": ell1_dash,
                "ell2_dash": ell2_dash,
                "ELL_dash": ELL_dash
                }

        self.volume = volume
        self.nmean  = nmean
        self.deltaK = deltaK

        ## type of powerspectra ##
        self.name = name

        ## set kbin ##
        self.kbin = kbin

        ## set multipole indices ##
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
        NCOMP = num_kbin
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

##########################

        AA = []
        for i in range(num_kbin):
            print("k1 = ", self.kbin[i], "h/Mpc")
            self.kmag1 = self.kbin[i]
            
            if NDIM <= 3:
    
                AA.append(pycuba.Cuhre(
                        self.Integrand_cov_PB, 
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
                        self.Integrand_cov_PB,
                        NDIM, 
                        NNEW, NMIN, FLATNESS, 
                        ncomp=NCOMP,
                        maxeval = MAXEVAL,
                        verbose=0 | 4)["results"])
                
        result = np.zeros((num_kbin, num_kbin))
        
        for i in range(num_kbin):
            for j in range(num_kbin):
                result[i,j] = AA[i][j]["integral"]
        
        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        ## output dict ##
        
        output_dict = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "cov_PB": result,
                "ELL": self.ELL,
                "ell1_dash": self.ell1_dash,
                "ell2_dash": self.ell2_dash,
                "ELL_dash": self.ELL_dash}
        
        print(output_dict)

        return output_dict

