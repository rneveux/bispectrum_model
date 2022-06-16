import os
import numpy as np
from classy import Class
import pycuba
import hitomipy
import initial

class ComputeKernels():

    def __init__(self, k, params_cosmo, z, ell1=0, ell2=0, ELL=0, diag=True):

        self.params_cosmo = params_cosmo
        self.k = k
        self.z = z
        self.ell1 = ell1
        self.ell2 = ell2
        self.ELL = ELL
        self.diag = diag

    def intial_power_spectrum(self):

        cosmo = Class()
        cosmo.set(self.params_cosmo)
        cosmo.compute()

        self.initial_cosmo = initial.InputPowerSpectrum(self.z, cosmo, params_fid=self.params_cosmo)
        self.initial_cosmo.calcMatterPowerSpectrum()
        self.k_in, self.pk_in = self.initial_cosmo.getMatterPowerSpectrum()
        self.pk_in_no_wiggle = self.initial_cosmo.getNoWigglePowerSpectrum()

        self.f = self.initial_cosmo.getGrowthRate()
        self.rs_drag = cosmo.rs_drag()
        self.sigma8_norm = self.initial_cosmo.getSigma8ForNormalization()

    def kernel_computation(self,kernel_name,to_save=False):

        self.kernel_name = kernel_name

        if self.diag: self.calc_K_diag()
        else: self.calc_K_full()

        if to_save: 
            self.K['params_cosmo'] = self.params_cosmo
            self.K['z'] = self.z
            self.save(to_save)

    def calc_K_full(
        self,
        integrand='tree',
        ks=.05, #flexible - no impact for physical scales
    ):

        self.integrand = integrand
        (kbin2_out, kbin1_out) = np.meshgrid(self.k, self.k)
        ## flags ##
        output_dict_ini = {
            "kbin1": kbin1_out,
            "kbin2": kbin2_out,
            "K": np.zeros((len(self.k), len(self.k))),
            "ell1": self.ell1,
            "ell2": self.ell2,
            "ELL": self.ELL,
        }

        ## initialization ##
        hitomipy.initializeInputPowerSpectrum_py()

        ## calc. wigner 3j symbols ##
        hitomipy.setWigner3j_py()

        ## read linear power spectrum ##
        hitomipy.readInputPowerSpectrum_py(self.k_in, self.pk_in, len(self.k_in))
        hitomipy.readInputNoWigglePowerSpectrum_py(self.k_in, self.pk_in_no_wiggle, len(self.k_in))

        ## normalization ##
        hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
        h = self.params_cosmo['h']
        Omega_b = self.params_cosmo['omega_b']/h**2
        Omega_m = (self.params_cosmo['omega_b']+self.params_cosmo['omega_cdm'])/h**2
        Tcmb = 2.7255
        n_s = self.params_cosmo['n_s']
        hitomipy.calcNormalizationNoWiggle_py(1.0, h, Omega_b, Omega_m, Tcmb, n_s)

        ## IR resummation damping term ##
        self.Sigma2 = hitomipy.Sig2_py(self.rs_drag,ks)
        self.dSigma2 = hitomipy.dSig2_py(self.rs_drag,ks)
        
        ## compute bispectra ##
        NDIM = 3
        NCOMP = len(self.k)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

        AA = []
        for i in range(NCOMP):
            print("k1 = ", self.k[i], "h/Mpc")
            self.kmag1 = self.k[i]
            AA.append(
                pycuba.Cuhre(
                    self.Integrand_K, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4
                )["results"]
                )

        bk_temp = np.zeros((NCOMP, NCOMP))
        for i in range(NCOMP):
            for j in range(NCOMP):
                bk_temp[i, j] = AA[i][j]["integral"]

        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        #########################
        (kbin2_out, kbin1_out) = np.meshgrid(self.k, self.k)

        ## sigma8 normalisation ##
        if self.integrand == 'tree': bk_out = bk_temp*self.sigma8_norm**4
        elif self.integrand == 'SN': bk_out = bk_temp*self.sigma8_norm**2
        elif self.integrand == 'PNG': bk_out = bk_temp*self.sigma8_norm**3

        ## output dict ##
        self.K = {
            "kbin1": kbin1_out,
            "kbin2": kbin2_out,
            "K": bk_out,
            "ell1": self.ell1,
            "ell2": self.ell2,
            "ELL": self.ELL,
            "kernel": self.kernel_name
        }

    def calc_K_diag(
        self,
        integrand='tree',
        ks=.05, #flexible - no impact for physical scale
    ):

        self.integrand = integrand

        ## flags ##
        output_dict_ini = {
            "kbin": self.k,
            "K": np.zeros((len(self.k))),
            "ell1": self.ell1,
            "ell2": self.ell2,
            "ELL": self.ELL,
        }

        ## initialization ##
        hitomipy.initializeInputPowerSpectrum_py()

        ## calc. wigner 3j symbols ##
        hitomipy.setWigner3j_py()

        ## read linear power spectrum ##
        hitomipy.readInputPowerSpectrum_py(self.k_in, self.pk_in, len(self.k_in))
        hitomipy.readInputNoWigglePowerSpectrum_py(self.k_in, self.pk_in_no_wiggle, len(self.k_in))

        ## normalization ##
        hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
        h = self.params_cosmo['h']
        Omega_b = self.params_cosmo['omega_b']/h**2
        Omega_m = (self.params_cosmo['omega_b']+self.params_cosmo['omega_cdm'])/h**2
        Tcmb = 2.7255
        n_s = self.params_cosmo['n_s']
        hitomipy.calcNormalizationNoWiggle_py(1.0, h, Omega_b, Omega_m, Tcmb, n_s)

        # IR resummation damping term
        self.Sigma2 = hitomipy.Sig2_py(self.rs_drag,ks)
        self.dSigma2 = hitomipy.dSig2_py(self.rs_drag,ks)

        ## compute bispectra ##
        NDIM = 3
        NCOMP = len(self.k)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

        AA = pycuba.Cuhre(
                    self.Integrand_K, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4
                )["results"]

        bk_temp = np.zeros((NCOMP))
        for i in range(NCOMP):
            bk_temp[i] = AA[i]["integral"]

        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        #########################

        ## sigma8 normalisation ##
        if self.integrand == 'tree': bk_out = bk_temp*self.sigma8_norm**4
        elif self.integrand == 'SN': bk_out = bk_temp*self.sigma8_norm**2
        elif self.integrand == 'PNG': bk_out = bk_temp*self.sigma8_norm**3

        ## output dict ##
        self.K = {
            "kbin": self.k,
            "K": bk_out,
            "ell1": self.ell1,
            "ell2": self.ell2,
            "ELL": self.ELL,
            "kernel": self.kernel_name
        }

    def Integrand_K(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])


        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        if self.integrand=='tree':
            if self.diag:
                hitomipy.integrand_B_Kernel_diag_py(
                    self.xx_in,
                    self.ndim,
                    self.ff_out,
                    self.ncomp,
                    self.k,
                    len(self.k),
                    self.ell1,
                    self.ell2,
                    self.ELL,
                    self.f,
                    self.Sigma2,
                    self.dSigma2,
                    self.kernel_name
                )
            else:
                hitomipy.integrand_B_Kernel_py(
                        self.xx_in,
                        self.ndim,
                        self.ff_out,
                        self.ncomp,
                        self.k,
                        len(self.k),
                        self.ell1,
                        self.ell2,
                        self.ELL,
                        self.kmag1,
                        self.f,
                        self.Sigma2, 
                        self.dSigma2, 
                        self.kernel_name
                    )

        elif self.integrand=='SN':
            if self.diag:
                hitomipy.integrand_B_Kernel_SN_diag_py(
                    self.xx_in,
                    self.ndim,
                    self.ff_out,
                    self.ncomp,
                    self.kbin,
                    len(self.kbin),
                    self.ell1,
                    self.ell2,
                    self.ELL,
                    self.f,
                    self.Sigma2,
                    self.dSigma2,
                    self.kernel_name
                )
            else:
                hitomipy.integrand_B_Kernel_SN_py(
                    self.xx_in,
                    self.ndim,
                    self.ff_out,
                    self.ncomp,
                    self.k,
                    len(self.k),
                    self.ell1,
                    self.ell2,
                    self.ELL,
                    self.kmag1,
                    self.f,
                    self.Sigma2, 
                    self.dSigma2, 
                    self.kernel_name
                )

        elif self.integrand=='PNG':
            hitomipy.integrand_B_Kernel_PNG_diag_py(
                    self.xx_in,
                    self.ndim,
                    self.ff_out,
                    self.ncomp,
                    self.kbin,
                    len(self.kbin),
                    self.ell1,
                    self.ell2,
                    self.ELL,
                    self.f,
                    self.Sigma2,
                    self.dSigma2,
                    self.kernel_name
                )

        for i in range(ncomp[0]):
                ff[i] = self.ff_out[i]

        return 0

    def save(self,to_save):

        os.makedirs(os.path.dirname(to_save), exist_ok=True)
        np.save(to_save,self.K)

    def load(self,to_load):

        return np.load(to_load,allow_pickle=True).item()




