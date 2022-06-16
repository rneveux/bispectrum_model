import os
import numpy as np
from classy import Class
import pycuba
import hitomipy
import initial

class ComputeBispectrum():
    r"""Computation of the tree level, the shot noise and the primordial non 
    gaussianity contribution of the bispectrum, either split into kernels or 
    not.
    All computation are done using the tripolar spherical harmonic 
    decomposition formalised in ref. [1]_.
    .. [1] Sugiyama N. et al., 2019. **MNRAS** 484(1), 364-384.
    [arXiv: `1803.02132 <https://arxiv.org/abs/1803.02132>`_]
    The bispectrum is here caracterised by two wavenumbers :math:`k_1` and 
    :math:`k_1` and three multipoles :math:`l_1`, :math:`l_2`, :math:`L`, the 
    latest states for the LOS anisotropies.

    Parameters
    ----------
    params_cosmo : dict
        Cosmological parameters to compute linear power spectrum with Class.
    z : float
        Redshift.
    diag : bool, optional
        If `True` (default), model computed on the diagonal part of the bispectrum.
    Attributes
    ----------
    Da : float
        Angular distance at the given cosmology and redshift.
    H : float
        Hubble parameter at the given cosmology and redshift.
    K : dict
        Result of the computation stored with parameters used (wavenumbers, 
        multipoles, kernels).
    """

    def __init__(self, params_cosmo, z, diag=True):

        self.params_cosmo = params_cosmo
        self.z = z
        self.diag = diag

    def intial_power_spectrum(self,cosmo=False):
        """Computation of the linear power spectrum, the no wiggle part of it 
        and derived cosmological parameters.
        ----------
        cosmo : :class:`classy.Class`, optional
            If `False` (default), Class is instanciate with params_cosmo.
        """

        if not cosmo:
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
        self.Da,self.H = self.initial_cosmo.calcFiducialHubbleAndDiameterDistance()

    def kernel_computation(self, kernel_name, k, ell1, ell2, ELL, integrand='tree', to_save=False):
        """Call the diagonal or full part of bispectrum model kernel.
        ----------
        kernel_name : str
            bias parameters associated to the kernel to compute.
        k : list
            wavevenumbers at which compute the kernel.
        ell1 : int
            first multipole of the tripoSH decomposition.
        ell2 : int
            second multipole of the tripoSH decomposition.
        ELL : int
            third multipole of the tripoSH decomposition, LOS anisotropies.
        integrand : str, optional
            part of bispectrum to compute, default is `tree`.
        to_save : str, optional
            If `False` (default), kernel not saved.
        """

        self.kernel_name = kernel_name
        self.k = k
        self.ell1 = ell1
        self.ell2 = ell2
        self.ELL = ELL

        if self.diag: self.calc_K_diag(integrand=integrand)
        else: self.calc_K_full(integrand=integrand)

        if to_save: 
            self.K['params_cosmo'] = self.params_cosmo
            self.K['z'] = self.z
            self.save(to_save)

    def calc_K_full(
        self,
        integrand='tree',
        ks=.05,
    ):
        """Computation of the full bispectrum kernel.
        ----------
        integrand : str, optional
            part of the bispectrum model to compute. `tree` (default), 
            `SN` (shot noise) or `PNG` (primordial non gaussianity).
        ks : float, optional
            Separation scale for IR resummation (default is 0.05), no impact on 
            the calculation for physical scales.
        """

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
        ks=.05,
    ):
        """Computation of the diagonal part of the bispectrum kernel.
        ----------
        integrand : str, optional
            part of the bispectrum model to compute. `tree` (default), 
            `SN` (shot noise) or `PNG` (primordial non gaussianity).
        ks : float, optional
            Separation scale for IR resummation (default is 0.05), no impact on 
            the calculation for physical scales.
        """

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
        """Call of the function to integrate for kernel part.
        ----------
        No direct call of this function needed.
        """

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

        for i in range(ncomp[0]):
                ff[i] = self.ff_out[i]

        return 0

    def save(self,to_save):
        """Save the result dictionary of the computation.
        ----------
        to_save : str
            name of the file to save the result.
        """

        os.makedirs(os.path.dirname(to_save), exist_ok=True)
        np.save(to_save,self.K)

    def load(self,to_load):
        """Loading of dictionary saved in .npy file.
        ----------
        to_load : str
            name of the file to load.
        """

        return np.load(to_load,allow_pickle=True).item()

    def calc_B_diag(
        self,
        k, ell1, ell2, ELL,
        alpha_perp=1, alpha_parallel=1, b1=2, b2=0, bG2=0, c1=0, c2=0, knl=.3,
        integrand='tree',
        ks=.05,
    ):
        """Computation of the diagonal part of the bispectrum.
        ----------
        k : list
            wavevenumbers at which compute the kernel.
        ell1 : int
            first multipole of the tripoSH decomposition.
        ell2 : int
            second multipole of the tripoSH decomposition.
        ELL : int
            third multipole of the tripoSH decomposition, LOS anisotropies.
        integrand : str, optional
            part of the bispectrum model to compute. `tree` (default), 
            `SN` (shot noise) or `PNG` (primordial non gaussianity).
        ks : float, optional
            Separation scale for IR resummation (default is 0.05), no impact on 
            the calculation for physical scales.
        """

        self.k = k
        self.ell1 = ell1
        self.ell2 = ell2
        self.ELL = ELL
        self.alpha_perp = alpha_perp
        self.alpha_parallel = alpha_parallel
        self.b1 = b1
        self.b2 = b2
        self.bG2 = bG2
        self.c1 = c1
        self.c2 = c2
        self.knl = knl
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
                    self.Integrand_B, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4
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
        }

    def Integrand_B(self, ndim, xx, ncomp, ff, userdata):
        """Call of the function to integrate for bispectrum model.
        ----------
        No direct call of this function needed.
        """

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])


        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        if self.integrand=='tree':
            if self.diag:
                hitomipy.integrand_B_FoG_Damping_diag_py(
                    self.xx_in,
                    self.ndim,
                    self.ff_out,
                    self.ncomp,
                    self.k,
                    len(self.k),
                    self.ell1,
                    self.ell2,
                    self.ELL,
                    self.alpha_perp,
                    self.alpha_parallel,
                    self.f,
                    self.b1,
                    self.b2,
                    self.bG2,
                    self.c1,
                    self.c2,
                    self.knl,
                    self.Sigma2,
                    self.dSigma2,
                )
            else:
                pass
                '''hitomipy.integrand_B_FoG_py(
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
                    )'''

        elif self.integrand=='tree+SN':
            if self.diag:
                pass
                '''hitomipy.integrand_B_Kernel_SN_diag_py(
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
                )'''
            else:
                pass
                '''hitomipy.integrand_B_Kernel_SN_py(
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
                )'''

        elif self.integrand=='PNG':
            pass
            '''hitomipy.integrand_B_Kernel_PNG_diag_py(
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
                )'''

        for i in range(ncomp[0]):
                ff[i] = self.ff_out[i]

        return 0




