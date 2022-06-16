#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from classy import Class
import hitomipy

class InputPowerSpectrum():

    def __init__(self, redshift, cosmo, params_fid=None):
        self.cosmo = cosmo
        self.redshift = redshift
        if redshift < 1.0e-10:
            self.redshift = 1.0e-10
        ## matter power spectrum ##
        self.kmax = 50.0
        self.kmin = 2.0e-5
        self.num_kbin = 500
        self.k_pivot = 0.05
        ln_kmin = np.log(self.kmin)
        ln_kmax = np.log(self.kmax)
        self.k = np.logspace(ln_kmin, ln_kmax, self.num_kbin, base=np.e)
        self.pk = np.zeros(len(self.k))
        self.pk_no_wiggle = np.zeros(len(self.k))
        self.pk_pri = np.zeros(len(self.k))
        self.Dz = 1.0
        self.fz = 0.0
        self.sigma8_norm = 1.0
        self.sigma8_0 = 0.0
        self.H = 0.0
        self.Da = 0.0
        self.H_fid = 0.0
        self.Da_fid = 0.0
        self._c_ = 2.99792458e5  # [ km/s ] ##

        if params_fid is None:
            h_fid = 0.6727
            omega_b_fid = 0.0492 * h_fid**2
            omega_cdm_fid = (0.3156 - 0.0492) * h_fid**2
            n_s_fid = 0.9645
            ln10A_s10_fid = 3.094
            tau_reio_fid = 0.0826026

            self.params_fid = {
                    'h': h_fid,
                    'omega_b': omega_b_fid,
                    'omega_cdm': omega_cdm_fid,
                    'n_s': n_s_fid,
                    'ln10^{10}A_s': ln10A_s10_fid,
                    'tau_reio': tau_reio_fid,
            }
        else:
            self.params_fid = params_fid

    def calcFiducialHubbleAndDiameterDistance(self):

        cosmo_temp = Class()

        ##########################################
        ########## fiducial parameter ############
        ##########################################

        h = self.params_fid["h"]
        ### fiducial Dz and H ###
        cosmo_temp.set(self.params_fid)
        cosmo_temp.compute()
        self.Da_fid = cosmo_temp.angular_distance(self.redshift) * h
        self.H_fid = cosmo_temp.Hubble(self.redshift) * self._c_
        cosmo_temp.struct_cleanup()
        cosmo_temp.empty()

        return self.Da_fid,self.H_fid

    def calcMatterPowerSpectrum(self, no_wiggle=True):

        h = self.cosmo.h()
        Omega_b = self.params_fid['omega_b']/h**2
        Omega_m = (self.params_fid['omega_b']+self.params_fid['omega_cdm'])/h**2
        Tcmb = 2.7255
        n_s = self.params_fid['n_s']

        ## compute linear matter power spectrum ##
        for i in range(self.num_kbin):
            self.pk[i] = self.cosmo.pk_lin(self.k[i] * h, self.redshift) * h**3
            if no_wiggle:
                self.pk_no_wiggle[i] = hitomipy.f_pk_no_wiggle_integrand_py(self.k[i] , h, Omega_b, Omega_m, Tcmb, n_s) #* h**3

        self.Da = self.cosmo.angular_distance(self.redshift) * h
        self.H = self.cosmo.Hubble(self.redshift) * self._c_
        self.Dz = self.cosmo.scale_independent_growth_factor(self.redshift)
        self.fz = self.cosmo.scale_independent_growth_factor_f(self.redshift)
        self.sigma8_norm = self.cosmo.sigma(8.0 / h, self.redshift)

    def calcPrimordialPowerSpectrum(self, ln10A_s10, n_s):

        h = self.cosmo.h()

        def f_pk_pri(k):
            A_s = 1.0e-10 * np.exp(ln10A_s10)
            P_pri = (3.0 / 5.0)**2 * (2.0 * np.pi**2 / k**3) * \
                    A_s * (k / (self.k_pivot * h))**(n_s - 1.0)
            return P_pri

        for i in range(self.num_kbin):
            self.pk_pri[i] = f_pk_pri(self.k[i] * h) * h**3

    def getMatterPowerSpectrum(self):
        return self.k, self.pk

    def getNoWigglePowerSpectrum(self):
        return self.pk_no_wiggle

    def getPrimordialPowerSpectrum(self):
        return self.pk_pri

    def getGrowthRate(self):
        return self.fz

    def getGrowthFactor(self):
        return self.Dz

    def getSigma8z(self, sigma8_0=-1.0):
        if sigma8_0 < 0.0:
            return self.sigma8_norm
        else:
            return sigma8_0 * self.Dz

    def getSigma8ForNormalization(self):
        return self.sigma8_norm

    def getAlphaPerp(self):
        if self.Da_fid < 1.0e-10:
            return 0.0
        else:
            return self.Da / self.Da_fid

    def getAlphaParallel(self):
        if self.H_fid < 1.0e-10:
            return 0.0
        else:
            return self.H_fid / self.H

    def getAngularDiameterDistance(self):
        return self.Da

    def getAngularDiameterDistanceFiducial(self):
        return self.Da_fid

    def getHubbleParameter(self):
        return self.H

    def getHubbleParameterFiducial(self):
        return self.H_fid
