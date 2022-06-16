#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os

class ClassPowerSpectrum():

    def __init__(self):
        self.initialize()

    def set_fitting_parameters(self, A20, A11, A02, A30, A21, A12, A03, A40, A31, A22, A13, A04):

        self.A20 = A20
        self.A11 = A11
        self.A02 = A02

        self.A30 = A30
        self.A21 = A21
        self.A12 = A12
        self.A03 = A03

        self.A40 = A40
        self.A31 = A31
        self.A22 = A22
        self.A13 = A13
        self.A04 = A04

    def initialize(self):
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

    def select_P(self, name):

        n_kbin = len(self.kbin)

        if name == "Tree":
            return hitomipy.integrand_P_Tree_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1
                )

        elif name == "Tree_NoWiggle":
            return hitomipy.integrand_P_Tree_NoWiggle_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1
                )

        elif name == "Tree_BAO":
            return hitomipy.integrand_P_Tree_BAO_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                self.sigma2_perp, self.sigma2_para
                )

        elif name == "Tree_BAO_Fitting":
            return hitomipy.integrand_P_Tree_BAO_Fitting_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                self.sigma2_perp, self.sigma2_para,
                self.A20, self.A11, self.A02, 
                self.A30, self.A21, self.A12, self.A03,
                self.A40, self.A31, self.A22, self.A13, self.A04)

        elif name == "NonLinearFitting":
            return hitomipy.integrand_P_NonLinearFitting_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                self.sigma2_perp, self.sigma2_para
                )

        elif name == "NonLinearFitting_Window":
            return hitomipy.integrand_P_NonLinearFitting_Window_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
                self.sigma2_perp, self.sigma2_para,
                self.volume
                )

        elif name == "LocalMean":
            return hitomipy.integrand_P_LocalMean_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1, self.b2, self.bK2,
                self.sigma2_perp, self.sigma2_para,
                self.nmean, self.volume
                )
        
        elif name == "Tree_BAO_b1_b1":
            return hitomipy.integrand_P_Tree_BAO_b1_b1_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel,
                self.sigma2_perp, self.sigma2_para
                )

        elif name == "Tree_BAO_b1_f":
            return hitomipy.integrand_P_Tree_BAO_b1_f_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel,
                self.sigma2_perp, self.sigma2_para
                )

        elif name == "Tree_BAO_f_f":
            return hitomipy.integrand_P_Tree_BAO_f_f_py(
                self.xx_in, self.ndim, self.ff_out, self.ncomp,
                self.kbin, n_kbin, self.ELL,
                self.alpha_perp, self.alpha_parallel,
                self.sigma2_perp, self.sigma2_para
                )

        else:
            print("select_P: ERROR")

        return 0.0

    def select_ndim(self, name):

        name_dim = []
        name_dim.append("Tree_Window")
        name_dim.append("Tree_Window_IC")
        name_dim.append("Tree_Window_IC_Approx")
        name_dim.append("NonLinearFitting_Window")
        name_dim.append("LocalMean")

        if name in name_dim:
            return 4
        elif name == "SPT_2loop":
            return 7
        else:
            return 2

    def check_flag_window(self, name, flag_window, volume):

        flag = 0
        name_window = []
        name_window.append("Tree_Window")
        name_window.append("Tree_Window_IC")
        name_window.append("Tree_Window_IC_Approx")
        name_window.append("NonLinearFitting_Window")
        name_window.append("LocalMean")

        if name in name_window:
            if flag_window and volume > 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_BAO(self, name, flag_BAO, sigma8_fid, fz_fid):

        flag = 0
        name_BAO = []
        name_BAO.append("Tree_BAO")
        name_BAO.append("Tree_BAO_b1_b1")
        name_BAO.append("Tree_BAO_b1_f")
        name_BAO.append("Tree_BAO_f_f")
        name_BAO.append("Tree_BAO_Fitting")

        if name in name_BAO:
            if flag_BAO and sigma8_fid >= 0.0 and fz_fid >= 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_nl(self, name, flag_BAO, sigma8_fid, fz_fid, flag_nl, kbin_nl, pk_nl_0, pk_nl_2):

        flag = 0
        name_nl = []
        name_nl.append("NonLinearFitting")
        name_nl.append("NonLinearFitting_Window")
        name_nl.append("LocalMean")
        if name in name_nl:
            if flag_BAO and sigma8_fid >= 0.0 and fz_fid >= 0.0 and flag_nl and np.sum(kbin_nl) > 0.0 and np.sum(pk_nl_0) > 0.0 and np.sum(pk_nl_2) > 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_Recon(self, name, flag_Recon, b1_fid, b1_input, R):

        flag = 0
        name_Recon = []
        
        if name in name_Recon:
            if flag_Recon and b1_fid > 0.0 and R >= 0.0 and b1_input > 0.0:
                flag = 0
            else:
                flag = -1
                    
        return flag


    def Integrand_P(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])

        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]
        
        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        self.select_P(self.name)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]

        return 0

    def calc_P( self,  name,
                kbin=np.linspace(0.01, 0.2, 20), ELL=0,
                flag_2pcf=False,
                flag_BAO=False, sigma8_fid=-1.0, fz_fid=- 1.0,
                flag_nl=False, kbin_nl=np.zeros(20), pk_nl_0=np.zeros(20), pk_nl_2=np.zeros(20), 
                flag_window=False, volume=500.0**3,
                flag_LM=False, nmean = 3.0e-4,
                flag_sigmaB=False,
                flag_Recon = False, b1_fid = - 1.0, R = - 1.0, b1_input = -1.0,
                flag_Damping = False, sigma2_perp = 0.0, sigma2_para = 0.0
                ):

        ## flags ##
        output_dict_ini = {
                "kbin": kbin,
                "P": np.zeros(len(kbin)),
                "kbin_fft": kbin,
                "P_fft": np.zeros(len(kbin)),
                "ELL": ELL,
                "flag_2pcf": flag_2pcf,
                "flag_BAO": flag_BAO,
                "flag_nl": flag_nl,
                "flag_window": flag_window}

        check_bao = self.check_flag_BAO(name, flag_BAO, sigma8_fid, fz_fid)
        check_nl = self.check_flag_nl(name, flag_BAO, sigma8_fid, fz_fid, flag_nl, kbin_nl, pk_nl_0, pk_nl_2)
        check_window = self.check_flag_window(name, flag_window, volume)
        check_recon = self.check_flag_Recon(name, flag_Recon, b1_fid, b1_input, R)

        if check_bao < 0:
            print("FLAG_BAO: ERROR")
            return output_dict_ini

        if check_nl < 0:
            print("FLAG_NL: ERROR")
            return output_dict_ini

        if check_window < 0:
            print("FLAG_WINDOW ERROR")
            return output_dict_ini

        if check_recon < 0:
            print("FLAG_RECON: ERROR")
            return output_dict_ini

        if flag_window:
            self.volume = volume

        if flag_LM:
            self.nmean = nmean

        if flag_BAO:
            self.fz_fid = fz_fid
            self.sigma8_fid = sigma8_fid

        if flag_Recon:
            self.b1_fid = b1_fid
            self.R = R
            self.b1_input = b1_input

        ## type of powerspectra, e.g,, "Tree", "NoWiggle" and etc. ##
        self.name = name

        ## set kbin ##
        if not flag_2pcf:
            self.kbin = kbin
        elif flag_2pcf:
            kbin0 = np.logspace(np.log(3.0e-4), np.log(0.2), 100, base=np.e)
            kbin1 = np.logspace(np.log(0.201), np.log(10.0), 50, base=np.e)
            self.kbin = np.hstack([kbin0, kbin1])

        ## set multipole indices ##
        self.ELL = ELL

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
        hitomipy.calcNormalizationNoWiggle_py(1.0)

        ## sigma2_perp and sigma2_para ##
        if flag_BAO:

            if flag_Damping:

                self.sigma2_perp = sigma2_perp
                self.sigma2_para = sigma2_para
            
            else:

                self.sigma2_perp = hitomipy.calcSigma_dd_py(self.sigma8_fid)
                self.sigma2_para = (1.0 + self.fz_fid) * (1.0 + self.fz_fid) * self.sigma2_perp

#                print("sigma2_perp = ", self.sigma2_perp)
#                print("sigma2_para = ", self.sigma2_para)

#            print("Normal=", self.sigma2_perp, self.sigma2_para)
            
            if flag_Recon:

                self.sigma2_perp = pycuba.Cuhre(
                        self.Integrand_P_sigma2_perp_Reconstructed, 
                        2, 
                        ncomp=1, 
                        key=0, verbose=0 | 4)["results"][0]['integral']


                self.sigma2_para = pycuba.Cuhre(
                        self.Integrand_P_sigma2_para_Reconstructed, 
                        2, 
                        ncomp=1, 
                        key=0, verbose=0 | 4)["results"][0]['integral']

#                print("Recon=", self.sigma2_perp, self.sigma2_para, flag_Recon)

        ## compute power spectra ##
        NDIM = self.select_ndim(self.name)
        NCOMP = len(self.kbin)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

        AA = pycuba.Cuhre(
                self.Integrand_P,
                NDIM,
                ncomp=NCOMP,
                key=0,
                verbose=0 | 4)

        if flag_sigmaB:
            BB = pycuba.Cuhre(
                    self.Integrand_P_SigmaB, 
                    2, 
                    ncomp=1, 
                    key=0, verbose=0 | 4, 
                    maxeval=100000)
            
            sigmaB = BB["results"][0]['integral']

#        NNEW = 2000
#        NMIN = 2
#        FLATNESS = 50
#        MAXEVAL = 2000
#        print("Suave")
#        AA = pycuba.Suave(self.Integrand_P, NDIM, NNEW, NMIN, FLATNESS,\
#                          ncomp=NCOMP, maxeval = MAXEVAL, verbose=0 | 4)

        pk_temp = np.zeros(NCOMP)
        for i in range(NCOMP):
            pk_temp[i] = AA["results"][i]['integral']

        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        if flag_2pcf:
            f_pk = interpolate.interp1d(
                self.kbin, pk_temp, fill_value="extrapolate", kind="cubic")
            kbin_out = kbin
            pk_out = f_pk(kbin)
            kbin_fft = self.kbin
            pk_fft = pk_temp
        else:
            kbin_out = self.kbin
            pk_out = pk_temp
            kbin_fft = self.kbin
            pk_fft = pk_temp

        ## output dict ##
        output_dict = {
                "kbin": kbin_out,
                "P": pk_out,
                "kbin_fft": kbin_fft,
                "P_fft": pk_fft,
                "ELL": self.ELL,
                "flag_2pcf": flag_2pcf,
                "flag_BAO": flag_BAO,
                "flag_nl": flag_nl,
                "flag_window": flag_window}

        if flag_sigmaB:
            return sigmaB
        else:
            return output_dict

    def calc_P_to_2PCF(self, pk_in, rbin = np.linspace(0.0, 200, 41), N_fftlog = 1000):

        if not pk_in["flag_2pcf"]:

            output_dict_ini = {
                "rbin": rbin, 
                "2pcf": np.zeros(len(rbin)), 
                "rbin_fft": rbin, 
                "2pcf_fft": np.zeros(len(rbin)), 
                "ELL": pk_in["ELL"],
                "flag_2pcf": pk_in["flag_2pcf"], 
                "flag_BAO": pk_in["flag_BAO"], 
                "flag_nl": pk_in["flag_nl"], 
                "flag_window": pk_in["flag_window"], 
                "N_fftlog": N_fftlog}
            
            return output_dict_ini

        ## set bispec. ##
        BB = pk_in["P_fft"]

        ## set kbin ##
        self.kbin = pk_in["kbin_fft"]
        
        ## set multipole indices ##
        self.ELL  = pk_in["ELL"]
        
        ## set rbin ##
        self.rbin = rbin
        
        ## input parameter for fftlog ##
        NNN = N_fftlog
        
        ## compute 2PCF ##
        kbin_for_xi = np.logspace(np.log(self.kbin[0]), np.log(self.kbin[-1]), NNN, base=np.e)
        
        f_pk = interpolate.interp1d(self.kbin, BB[:], fill_value = "extrapolate", kind="cubic")
        pk_for_xi = f_pk(kbin_for_xi)
        
        r_temp = np.zeros(NNN)
        xi_temp = np.zeros(NNN)
        hitomipy.hankel_py(self.ELL, 2, NNN, kbin_for_xi, pk_for_xi, r_temp, xi_temp)
        
        f_xi = interpolate.interp1d(r_temp, xi_temp, fill_value = "extrapolate", kind="cubic")
        
        CC = np.zeros(len(self.rbin))
        CC[:] = f_xi(self.rbin[:])
        
        CC_for_pk = np.zeros(len(r_temp))
        CC_for_pk[:] = f_xi(r_temp[:])
        
        sign = np.real(1.0j**(self.ELL))
        xi_out = sign * CC
        xi_fft = sign * CC_for_pk
        
        ## rbin_out ##
        rbin_out = self.rbin

        ## output dict ##
        output_dict = {"rbin": rbin_out, 
                       "2pcf": xi_out, 
                       "rbin_fft": r_temp, 
                       "2pcf_fft": xi_fft, 
                       "ELL": self.ELL, 
                       "flag_2pcf": pk_in["flag_2pcf"], 
                       "flag_BAO": pk_in["flag_BAO"], 
                       "flag_nl": pk_in["flag_nl"], 
                       "flag_window": pk_in["flag_window"], 
                       "N_fftlog": N_fftlog}

        return output_dict

    def calc_2PCF_to_P(self, xi_in, kbin = np.linspace(0.01, 0.2, 20)):

        if not xi_in["flag_2pcf"]:
            
            output_dict_ini = {
                "kbin": kbin,
                "P": np.zeros(len(kbin)),
                "kbin_fft": kbin,
                "P_fft": np.zeros(len(kbin)),
                "ELL": xi_in["ELL"],
                "flag_2pcf": xi_in["flag_2pcf"],
                "flag_BAO": xi_in["flag_BAO"],
                "flag_nl": xi_in["flag_nl"],
                "flag_window": xi_in["flag_window"],
                "N_fftlog": xi_in["N_fftlog"]}
    
            return output_dict_ini

        ## set powerspec. ##
        BB = xi_in["2pcf_fft"]
        
        ## set rbin ##
        self.rbin = xi_in["rbin_fft"]
        
        ## set multipole indices ##
        self.ELL  = xi_in["ELL"]
        
        ## set kbin ##
        self.kbin = kbin
        
        ## input parameter for fftlog ##
        NNN = xi_in["N_fftlog"]
        
        ## compute 3PCF ##
        rbin_for_pk = np.logspace(np.log(self.rbin[0]), np.log(self.rbin[-1]), NNN, base=np.e)
        CC = np.zeros(len(self.kbin))
        
        f_xi = interpolate.interp1d(self.rbin, BB[:], fill_value = "extrapolate", kind="cubic")
        xi_for_pk = f_xi(rbin_for_pk)
        
        k_temp = np.zeros(NNN)
        pk_temp = np.zeros(NNN)
        hitomipy.hankel_py(self.ELL, 2, NNN, rbin_for_pk, xi_for_pk, k_temp, pk_temp)
        
        f_pk = interpolate.interp1d(k_temp, pk_temp, fill_value = "extrapolate", kind="cubic")
        CC[:] = f_pk(self.kbin[:])
        
        sign = np.real((-1.0j)**(self.ELL)) * (2.0*np.pi)**3
        pk_out = sign * CC
        pk_temp = sign * pk_temp
        
        ## rbin_out ##
        kbin_out = self.kbin
        
        ## output dict ##
        output_dict = {
            "kbin": kbin_out, 
            "P": pk_out, 
            "kbin_fft": k_temp, 
            "P_fft": pk_temp, 
            "ELL": self.ELL,
            "flag_2pcf": xi_in["flag_2pcf"],
            "flag_BAO": xi_in["flag_BAO"],
            "flag_nl": xi_in["flag_nl"],
            "flag_window": xi_in["flag_window"],
            "N_fftlog": xi_in["N_fftlog"]}
        
        return output_dict


    def Integrand_P_SigmaB(self, ndim, xx, ncomp, ff, userdata):
        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])
        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]
           
        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        hitomipy.integrand_P_SigmaB_py(
            self.xx_in, self.ndim, self.ff_out, self.ncomp, 
            self.kbin, len(self.kbin), 
            self.alpha_perp, self.alpha_parallel, self.sigma8, self.fz, self.b1,
            self.sigma2_perp, self.sigma2_para, 
            self.nmean, self.volume)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]
        
        return 0

    def Integrand_P_sigma2_perp_Reconstructed(self, ndim, xx, ncomp, ff, userdata):
        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])
        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        hitomipy.integrand_P_sigma2_perp_Reconstructed_py(
            self.xx_in, self.ndim, self.ff_out, self.ncomp, 
            self.kbin, len(self.kbin), 
            self.sigma8_fid, self.fz_fid, self.b1_input,
            self.b1_fid, self.R)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]
        
        return 0

    def Integrand_P_sigma2_para_Reconstructed(self, ndim, xx, ncomp, ff, userdata):
        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])
        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        hitomipy.integrand_P_sigma2_para_Reconstructed_py(
            self.xx_in, self.ndim, self.ff_out, self.ncomp, 
            self.kbin, len(self.kbin), 
            self.sigma8_fid, self.fz_fid, self.b1_input,
            self.b1_fid, self.R)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]
        
        return 0




#    def calc_P_SigmaB(self, nmean, volume, flag_BAO = "False", sigma8_fid = -1.0, fz_fid = - 1.0, flag_nl = "False", kbin_nl = np.zeros(20), pk_nl_0=np.zeros(20), pk_nl_2=np.zeros(20)):
#	self.nmean = nmean
#	self.volume = volume
#	self.kbin = np.zeros(20)
# check = self.check_flag_BAO(name, flag_BAO, sigma8_fid, fz_fid)
# if check < 0:
# print "FLAG ERROR"
# output_dict = {"kbin": kbin, "P": np.nan * np.ones(len(kbin)), "ELL": ELL, "flag_2pcf": flag_2pcf, "flag_BAO": flag_BAO, "flag_nl": flag_nl}
# return output_dict
#
# check = self.check_flag_nl(name, flag_nl, kbin_nl, pk_nl_0, pk_nl_2)
# if check < 0:
# print "FLAG ERROR"
# output_dict = {"kbin": kbin, "P": np.nan * np.ones(len(kbin)), "ELL": ELL, "flag_2pcf": flag_2pcf, "flag_BAO": flag_BAO, "flag_nl": flag_nl}
# return output_dict
#
#	## type of bispectra, e.g,, "Tree", "NoWiggle" and etc. ##
# self.name = name
#
#	## initialization ##
#	hitomipy.initializeInputPowerSpectrum_py()
#
#	## read linear power spectrum ##
#	hitomipy.readInputPowerSpectrum_py(self.k_temp, self.P_temp, len(self.k_temp))
#
#	kbin_nl_ver2 = np.zeros(len(kbin_nl))
#	pk_nl_0_ver2 = np.zeros(len(kbin_nl))
#	pk_nl_2_ver2 = np.zeros(len(kbin_nl))
#	for i in range(len(kbin_nl)):
#	    kbin_nl_ver2[i] = kbin_nl[i]
#	    pk_nl_0_ver2[i] = pk_nl_0[i]
#	    pk_nl_2_ver2[i] = pk_nl_2[i]
#
# hitomipy.readNonLinearPowerSpectrum_py(kbin_nl, pk_nl_0, pk_nl_2, len(kbin_nl))
#	hitomipy.readNonLinearPowerSpectrum_py(kbin_nl_ver2, pk_nl_0_ver2, pk_nl_2_ver2, len(kbin_nl))
#
#
#	hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
#	hitomipy.calcNormalizationNoWiggle_py(1.0)
#
#	## sigma2_perp and sigma2_para ##
#	if flag_BAO == "True":
#	    self.sigma2_perp = hitomipy.calcSigma_dd_py(sigma8_fid)
#	    self.sigma2_para = ( 1.0 + fz_fid ) * ( 1.0 + fz_fid ) * self.sigma2_perp
#
#	## compute power spectra ##
#	NDIM = 2
#	NCOMP = 1
#	if NCOMP>1024:
#	    print "# of NCOMP should be <= 1024, otherwise results become zero."
#
#	AA = pycuba.Cuhre(self.Integrand_P_SigmaB, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4, maxeval=100000)
#
#	sigmaB = AA["results"][0]['integral']
#
#	## finalize parameters ##
#	hitomipy.finalizeInputPowerSpectrum_py()
#
##	## kbin_out ##
# kbin_out = self.kbin
##
##	## output dict ##
# output_dict = {"kbin": kbin_out, "P": pk_out, "ELL": self.ELL, "flag_2pcf": flag_2pcf, "flag_BAO": flag_BAO}
#
#	return sigmaB
#
#    def calc_cov_PP(self, name, kbin = np.linspace(0.01, 0.2, 20), ELL = 0, ELL_dash = 0,\
#	            flag_2pcf="False",\
#	            flag_BAO = "False", sigma8_fid = -1.0, fz_fid = - 1.0,\
#	            flag_nl = "False", kbin_nl = np.zeros(20), pk_nl_0=np.zeros(20), pk_nl_2=np.zeros(20),\
#	            flag_window = "False", volume = None,\
#	            flag_cov = False, nmean = None, DeltaK = None):
#
#	check = self.check_flag_BAO(name, flag_BAO, sigma8_fid, fz_fid)
#	if check < 0:
#	    print "FLAG ERROR", "flag_BAO"
#	    output_dict = {"kbin": kbin, "P": np.nan * np.ones(len(kbin)), "ELL": ELL, "flag_2pcf": flag_2pcf, "flag_BAO": flag_BAO, "flag_nl": flag_nl}
#	    return output_dict
#
#	check = self.check_flag_nl(name, flag_nl, kbin_nl, pk_nl_0, pk_nl_2)
#	if check < 0:
#	    print "FLAG ERROR", "flag_nl"
#	    output_dict = {"kbin": kbin, "P": np.nan * np.ones(len(kbin)), "ELL": ELL, "flag_2pcf": flag_2pcf, "flag_BAO": flag_BAO, "flag_nl": flag_nl}
#	    return output_dict
#
#	check = self.check_flag_window(name, flag_window, volume)
#	if check < 0:
#	    print "FLAG ERROR", "flag_window"
#	    output_dict = {"kbin": kbin, "P": np.nan * np.ones(len(kbin)), "ELL": ELL, "flag_2pcf": flag_2pcf, "flag_BAO": flag_BAO, "flag_nl": flag_nl}
#	    return output_dict
#
#	check = self.check_flag_cov(name, flag_cov, DeltaK, nmean, volume)
#	if check < 0:
#	    return -1
#
#	if flag_window == "True":
#	    self.volume = volume
#
#	if flag_cov == True:
#	    self.volume = volume
#	    self.DeltaK = DeltaK
#	    self.nmean = nmean
#
#	## type of bispectra, e.g,, "Tree", "NoWiggle" and etc. ##
#	self.name = name
#
#	## set kbin ##
#	if flag_2pcf == "False":
#	    self.kbin = kbin
#	elif flag_2pcf == "True":
#	    kbin0 = np.logspace(np.log(3.0e-4), np.log(0.2), 100, base=np.e)
#	    kbin1 = np.logspace(np.log(0.201), np.log(10.0), 50, base=np.e)
#	    self.kbin= np.hstack([kbin0, kbin1])
#
#	## set multipole indices ##
#	self.ELL  = ELL
#	self.ELL_dash  = ELL_dash
#
#	## initialization ##
#	hitomipy.initializeInputPowerSpectrum_py()
#
#	## read linear power spectrum ##
#	hitomipy.readInputPowerSpectrum_py(self.k_temp, self.P_temp, len(self.k_temp))
#
#	if flag_nl == "True":
#
#	    kbin_nl_ver2 = np.zeros(len(kbin_nl))
#	    pk_nl_0_ver2 = np.zeros(len(kbin_nl))
#	    pk_nl_2_ver2 = np.zeros(len(kbin_nl))
#	    for i in range(len(kbin_nl)):
#		kbin_nl_ver2[i] = kbin_nl[i]
#		pk_nl_0_ver2[i] = pk_nl_0[i]
#		pk_nl_2_ver2[i] = pk_nl_2[i]
#
# hitomipy.readNonLinearPowerSpectrum_py(kbin_nl, pk_nl_0, pk_nl_2, len(kbin_nl))
#	    hitomipy.readNonLinearPowerSpectrum_py(kbin_nl_ver2, pk_nl_0_ver2, pk_nl_2_ver2, len(kbin_nl))
#
#	hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
#	hitomipy.calcNormalizationNoWiggle_py(1.0)
#
#	## sigma2_perp and sigma2_para ##
#	if flag_BAO == "True":
#	    self.sigma2_perp = hitomipy.calcSigma_dd_py(sigma8_fid)
#	    self.sigma2_para = ( 1.0 + fz_fid ) * ( 1.0 + fz_fid ) * self.sigma2_perp
#
#	## compute power spectra ##
#	NDIM = self.select_ndim(self.name)
#	NCOMP = len(self.kbin)
#	if NCOMP>1024:
#	    print "# of NCOMP should be <= 1024, otherwise results become zero."
#
#	if self.name == "cov_PP_G" or self.name == "cov_PP_G_NL":
#	    AA = pycuba.Cuhre(self.Integrand_P, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4)
#	    result = np.zeros(NCOMP)
#	    for i in range(NCOMP):
#	        result[i] = AA["results"][i]['integral']
#	    result = np.diag(result)
#	else:
#	    AA = []
#	    for i in range(NCOMP):
#	        print "k1 = ", self.kbin[i], "h/Mpc"
#	        self.kmag1 = self.kbin[i]
#		if NDIM <= 3:
# if NDIM <= 2:
#		    AA.append(pycuba.Cuhre(self.Integrand_P, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4)["results"])
#		else:
#		    NNEW = 10000
#		    NMIN = 2
#		    FLATNESS = 50
#		    MAXEVAL = 10000
#		    print "Suave", NDIM
#		    AA.append(pycuba.Suave(self.Integrand_P, NDIM, NNEW, NMIN, FLATNESS, ncomp=NCOMP, maxeval = MAXEVAL, verbose=0 | 4)["results"])
#
#	    result = np.zeros((NCOMP,NCOMP))
#	    for i in range(NCOMP):
#	        for j in range(NCOMP):
#		   result[i,j] = AA[i][j]["integral"]
#
#	## finalize parameters ##
#	hitomipy.finalizeInputPowerSpectrum_py()
#
#	(kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
#
#	## output dict ##
#	output_dict = {"kbin1": kbin1_out, "kbin2": kbin2_out, "cov_PP": result, "ELL": self.ELL, "ELL_dash": self.ELL_dash, "flag_BAO": flag_BAO}
#
#	return output_dict
#
#
