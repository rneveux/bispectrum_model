# -*- coding: utf-8 -*-
import pyximport

import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os


class ClassBiSpectrum:
    def __init__(self):
        self.initialize()

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

        self.c1 = 0.0

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

        self.c1 = params["c1"]

    def select_B(self, name):

        n_kbin = len(self.kbin)

        if name == "Tree":
            return hitomipy.integrand_B_Tree_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
                self.b2,
                self.bK2,
            )

        elif name == "Tree_FoG":
            return hitomipy.integrand_B_Tree_FoG_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
                self.b2,
                self.bK2,
                self.c1,
                self.knl,
            )

        elif name == "Tree_DampIvanov":
            return hitomipy.integrand_B_Tree_DampIvanov_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
                self.b2,
                self.bK2,
                self.rbao,
                self.ks,
            )

        elif name == "Tree_NoWiggle":
            return hitomipy.integrand_B_Tree_NoWiggle_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
                self.b2,
                self.bK2,
            )

        elif name == "Tree_NoWiggle_RealSpace_DarkMatter_Growth":
            return hitomipy.integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Growth_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.sigma8,
            )

        elif name == "Tree_NoWiggle_RealSpace_DarkMatter_Shift":
            return hitomipy.integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Shift_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.sigma8,
            )

        elif name == "Tree_NoWiggle_RealSpace_DarkMatter_Tidal":
            return hitomipy.integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Tidal_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.sigma8,
            )

        elif name == "Tree_BAO_RealSpace_DarkMatter_Growth":
            return hitomipy.integrand_B_Tree_BAO_RealSpace_DarkMatter_Growth_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.sigma8,
                self.sigma2_perp,
            )

        elif name == "Tree_BAO_RealSpace_DarkMatter_Shift":
            return hitomipy.integrand_B_Tree_BAO_RealSpace_DarkMatter_Shift_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.sigma8,
                self.sigma2_perp,
            )

        elif name == "Tree_BAO_RealSpace_DarkMatter_Tidal":
            return hitomipy.integrand_B_Tree_BAO_RealSpace_DarkMatter_Tidal_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.sigma8,
                self.sigma2_perp,
            )

        elif name == "Tree_BAO":
            return hitomipy.integrand_B_Tree_BAO_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
                self.b2,
                self.bK2,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_Reconstructed":
            return hitomipy.integrand_B_Tree_Reconstructed_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
                self.b2,
                self.bK2,
                self.b1_fid,
                self.R,
            )

        elif name == "Tree_NoWiggle_Reconstructed":
            return hitomipy.integrand_B_Tree_NoWiggle_Reconstructed_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
                self.b2,
                self.bK2,
                self.b1_fid,
                self.R,
            )

        ###############

        elif name == "Tree_BAO_Template":
            return hitomipy.integrand_B_Tree_BAO_Template_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
                self.parameters,
            )

        elif name == "Tree_BAO_b1_b1_b1":
            return hitomipy.integrand_B_Tree_BAO_b1_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b1_b1_f":
            return hitomipy.integrand_B_Tree_BAO_b1_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b1_f_f":
            return hitomipy.integrand_B_Tree_BAO_b1_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b2_b1_b1":
            return hitomipy.integrand_B_Tree_BAO_b2_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b2_b1_f":
            return hitomipy.integrand_B_Tree_BAO_b2_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b2_f_f":
            return hitomipy.integrand_B_Tree_BAO_b2_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_bK2_b1_b1":
            return hitomipy.integrand_B_Tree_BAO_bK2_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_bK2_b1_f":
            return hitomipy.integrand_B_Tree_BAO_bK2_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_bK2_f_f":
            return hitomipy.integrand_B_Tree_BAO_bK2_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b1f_b1_b1":
            return hitomipy.integrand_B_Tree_BAO_b1f_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b1f_b1_f":
            return hitomipy.integrand_B_Tree_BAO_b1f_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_b1f_f_f":
            return hitomipy.integrand_B_Tree_BAO_b1f_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_ff_f_f":
            return hitomipy.integrand_B_Tree_BAO_ff_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        elif name == "Tree_BAO_f_f_f":
            return hitomipy.integrand_B_Tree_BAO_f_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma2_perp,
                self.sigma2_para,
            )

        ###############

        elif name == "Tree_b1_b1_b1":
            return hitomipy.integrand_B_Tree_b1_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b1_b1_f":
            return hitomipy.integrand_B_Tree_b1_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b1_f_f":
            return hitomipy.integrand_B_Tree_b1_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b2_b1_b1":
            return hitomipy.integrand_B_Tree_b2_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b2_b1_f":
            return hitomipy.integrand_B_Tree_b2_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b2_f_f":
            return hitomipy.integrand_B_Tree_b2_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_bK2_b1_b1":
            return hitomipy.integrand_B_Tree_bK2_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_bK2_b1_f":
            return hitomipy.integrand_B_Tree_bK2_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_bK2_f_f":
            return hitomipy.integrand_B_Tree_bK2_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b1f_b1_b1":
            return hitomipy.integrand_B_Tree_b1f_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b1f_b1_f":
            return hitomipy.integrand_B_Tree_b1f_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_b1f_f_f":
            return hitomipy.integrand_B_Tree_b1f_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_ff_f_f":
            return hitomipy.integrand_B_Tree_ff_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_f_f_f":
            return hitomipy.integrand_B_Tree_f_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        ###############

        elif name == "Tree_NoWiggle_b1_b1_b1":
            return hitomipy.integrand_B_Tree_NoWiggle_b1_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b1_b1_f":
            return hitomipy.integrand_B_Tree_NoWiggle_b1_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b1_f_f":
            return hitomipy.integrand_B_Tree_NoWiggle_b1_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b2_b1_b1":
            return hitomipy.integrand_B_Tree_NoWiggle_b2_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b2_b1_f":
            return hitomipy.integrand_B_Tree_NoWiggle_b2_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b2_f_f":
            return hitomipy.integrand_B_Tree_NoWiggle_b2_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_bK2_b1_b1":
            return hitomipy.integrand_B_Tree_NoWiggle_bK2_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_bK2_b1_f":
            return hitomipy.integrand_B_Tree_NoWiggle_bK2_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_bK2_f_f":
            return hitomipy.integrand_B_Tree_NoWiggle_bK2_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b1f_b1_b1":
            return hitomipy.integrand_B_Tree_NoWiggle_b1f_b1_b1_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b1f_b1_f":
            return hitomipy.integrand_B_Tree_NoWiggle_b1f_b1_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_b1f_f_f":
            return hitomipy.integrand_B_Tree_NoWiggle_b1f_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_ff_f_f":
            return hitomipy.integrand_B_Tree_NoWiggle_ff_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "Tree_NoWiggle_f_f_f":
            return hitomipy.integrand_B_Tree_NoWiggle_f_f_f_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
            )

        elif name == "fNL_Local":
            return hitomipy.integrand_B_NonGaussian_Local_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
            )

        elif name == "fNL_Equilateral":
            return hitomipy.integrand_B_NonGaussian_Equilateral_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
            )

        elif name == "fNL_Orthogonal":
            return hitomipy.integrand_B_NonGaussian_Orthogonal_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.alpha_perp,
                self.alpha_parallel,
                self.sigma8,
                self.fz,
                self.b1,
            )

        elif name == "KernelTemplate":
            return hitomipy.integrand_KernelTemplate_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.kbin,
                n_kbin,
                self.ell1,
                self.ell2,
                self.ELL,
                self.kmag1,
                self.parameters,
            )

        else:
            print("select_B: ERROR")

        return 0.0

    def select_ndim(self, name):

        name_dim = []
        name_dim_single = []
        name_dim_single.append("Tree_NoWiggle_RealSpace_DarkMatter_Growth")
        name_dim_single.append("Tree_NoWiggle_RealSpace_DarkMatter_Shift")
        name_dim_single.append("Tree_NoWiggle_RealSpace_DarkMatter_Tidal")
        name_dim_single.append("Tree_BAO_RealSpace_DarkMatter_Growth")
        name_dim_single.append("Tree_BAO_RealSpace_DarkMatter_Shift")
        name_dim_single.append("Tree_BAO_RealSpace_DarkMatter_Tidal")

        if name in name_dim:
            return 6
        if name in name_dim_single:
            return 2
        else:
            return 3

    def check_flag_BAO(self, name, flag_BAO, sigma8_fid, fz_fid):

        flag = 0
        name_BAO = []
        name_BAO.append("Tree_BAO_Template")

        name_BAO.append("Tree_BAO")

        name_BAO.append("Tree_BAO_b1_b1_b1")
        name_BAO.append("Tree_BAO_b1_b1_f")
        name_BAO.append("Tree_BAO_b1_f_f")
        name_BAO.append("Tree_BAO_b2_b1_b1")
        name_BAO.append("Tree_BAO_b2_b1_f")
        name_BAO.append("Tree_BAO_b2_f_f")
        name_BAO.append("Tree_BAO_bK2_b1_b1")
        name_BAO.append("Tree_BAO_bK2_b1_f")
        name_BAO.append("Tree_BAO_bK2_f_f")
        name_BAO.append("Tree_BAO_b1f_b1_b1")
        name_BAO.append("Tree_BAO_b1f_b1_f")
        name_BAO.append("Tree_BAO_b1f_f_f")
        name_BAO.append("Tree_BAO_ff_f_f")
        name_BAO.append("Tree_BAO_f_f_f")

        if name in name_BAO:
            if flag_BAO and sigma8_fid >= 0.0 and fz_fid >= 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_Recon(self, name, flag_Recon, b1_fid, R):

        flag = 0
        name_Recon = []
        name_Recon.append("Tree_Reconstructed")
        name_Recon.append("Tree_NoWiggle_Reconstructed")

        if name in name_Recon:
            if flag_Recon and b1_fid > 0.0 and R >= 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def check_flag_PNG(self, name, flag_PNG, k_pri, pk_pri):

        flag = 0
        name_PNG = []

        if name in name_PNG:
            if flag_PNG and len(k_pri) > 1 and np.sum(pk_pri) > 0.0:
                flag = 0
            else:
                flag = -1

        return flag

    def Integrand_B(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])

        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]
        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        self.select_B(self.name)

        for i in range(ncomp[0]):
            ff[i] = self.ff_out[i]

        return 0

    def calc_B(
        self,
        name,
        kbin=np.linspace(0.01, 0.2, 20),
        ell1=0,
        ell2=0,
        ELL=0,
        flag_3pcf=False,
        flag_BAO=False,
        sigma8_fid=-1.0,
        fz_fid=-1.0,
        flag_Recon=False,
        b1_fid=-1.0,
        R=-1.0,
        flag_PNG=False,
        k_pri=np.zeros(1),
        pk_pri=np.zeros(1),
        flag_FoG=False,        
        knl=.3,
        flag_DampIvanov=False,
        ks=.05,
        rbao=147,
        parameters=None,
    ):

        (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
        ## flags ##
        output_dict_ini = {
            "kbin1": kbin1_out,
            "kbin2": kbin2_out,
            "B": np.zeros((len(kbin), len(kbin))),
            "kbin1_fft": kbin1_out,
            "kbin2_fft": kbin2_out,
            "B_fft": np.zeros((len(kbin), len(kbin))),
            "ell1": ell1,
            "ell2": ell2,
            "ELL": ELL,
            "flag_3pcf": flag_3pcf,
            "flag_BAO": flag_BAO,
            "flag_Recon": flag_Recon,
            "flag_FoG": flag_FoG,
            "flag_DampIvanov": flag_DampIvanov,
            "flag_PNG": flag_PNG,
        }

        check_bao = self.check_flag_BAO(name, flag_BAO, sigma8_fid, fz_fid)
        check_recon = self.check_flag_Recon(name, flag_Recon, b1_fid, R)
        check_png = self.check_flag_PNG(name, flag_PNG, k_pri, pk_pri)

        if check_bao < 0:
            print("FLAG_BAO: ERROR")
            return output_dict_ini

        if check_recon < 0:
            print("FLAG_RECON: ERROR")
            return output_dict_ini

        if check_png < 0:
            print("FLAG_PNG ERROR")
            return output_dict_ini

        if flag_BAO:
            self.fz_fid = fz_fid
            self.sigma8_fid = sigma8_fid

        ## type of powerspectra, e.g,, "Tree", "NoWiggle" and etc. ##
        self.name = name
        self.parameters = parameters

        ## set kbin ##
        if not flag_3pcf:
            self.kbin = kbin
        elif flag_3pcf:
            kbin0 = np.logspace(np.log(3.0e-4), np.log(0.2), 100, base=np.e)
            kbin1 = np.logspace(np.log(0.201), np.log(10.0), 50, base=np.e)
            self.kbin = np.hstack([kbin0, kbin1])

        ## set multipole indices ##
        self.ell1 = ell1
        self.ell2 = ell2
        self.ELL = ELL

        ## initialization ##
        hitomipy.initializeInputPowerSpectrum_py()

        ## calc. wigner 3j symbols ##
        hitomipy.setWigner3j_py()

        ## read linear power spectrum ##
        hitomipy.readInputPowerSpectrum_py(self.k_temp, self.P_temp, len(self.k_temp))
        if flag_PNG:
            hitomipy.readInputPrimordialPowerSpectrum_py(k_pri, pk_pri, len(k_pri))

        ## normalization ##
        hitomipy.calcNormalizationUsingSigma8_py(self.sigma8_norm)
        hitomipy.calcNormalizationNoWiggle_py(1.0)

        ## sigma2_perp and sigma2_para ##
        if flag_BAO:
            self.sigma2_perp = hitomipy.calcSigma_dd_py(self.sigma8_fid)
            self.sigma2_para = (
                (1.0 + self.fz_fid) * (1.0 + self.fz_fid) * self.sigma2_perp
            )

        if flag_Recon:
            self.b1_fid = b1_fid
            self.R = R

        if flag_FoG:
            self.knl = knl

        if flag_DampIvanov:
            self.rbao = rbao
            self.ks = ks

        ## compute bispectra ##
        NDIM = self.select_ndim(self.name)
        NCOMP = len(self.kbin)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

        AA = []
        for i in range(NCOMP):
            print("k1 = ", self.kbin[i], "h/Mpc")
            self.kmag1 = self.kbin[i]
            if NDIM > 3:
                NNEW = 5000
                NMIN = 2
                FLATNESS = 50
                MAXEVAL = 5000
                AA.append(
                    pycuba.Suave(
                        self.Integrand_B,
                        NDIM,
                        NNEW,
                        NMIN,
                        FLATNESS,
                        ncomp=NCOMP,
                        maxeval=MAXEVAL,
                        verbose=0 | 4,
                    )["results"]
                )
            else:
                AA.append(
                    pycuba.Cuhre(
                        self.Integrand_B, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4
                    )["results"]
                )

        bk_temp = np.zeros((NCOMP, NCOMP))
        for i in range(NCOMP):
            for j in range(NCOMP):
                bk_temp[i, j] = AA[i][j]["integral"]

        ## finalize parameters ##
        hitomipy.finalizeInputPowerSpectrum_py()

        #########################
        if flag_3pcf:
            f_bk = interpolate.interp2d(self.kbin, self.kbin, bk_temp, kind="cubic")
            (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
            bk_out = f_bk(kbin, kbin)
            (kbin2_fft, kbin1_fft) = np.meshgrid(self.kbin, self.kbin)
            bk_fft = bk_temp
        else:
            (kbin2_out, kbin1_out) = np.meshgrid(self.kbin, self.kbin)
            bk_out = bk_temp
            (kbin2_fft, kbin1_fft) = np.meshgrid(self.kbin, self.kbin)
            bk_fft = bk_temp

        ## output dict ##
        output_dict = {
            "kbin1": kbin1_out,
            "kbin2": kbin2_out,
            "B": bk_out,
            "kbin1_fft": kbin1_fft,
            "kbin2_fft": kbin2_fft,
            "B_fft": bk_fft,
            "ell1": ell1,
            "ell2": ell2,
            "ELL": ELL,
            "flag_3pcf": flag_3pcf,
            "flag_BAO": flag_BAO,
            "flag_Recon": flag_Recon,
            "flag_PNG": flag_PNG,
            "flag_FoG": flag_FoG,
            "flag_DampIvanov": flag_DampIvanov,
        }

        return output_dict

    def calc_B_to_3PCF(self, bk_in, rbin=np.linspace(0.0, 200, 41), N_fftlog=1000):

        if not bk_in["flag_3pcf"]:
            (rbin2_out, rbin1_out) = np.meshgrid(rbin, rbin)
            output_dict_ini = {
                "rbin1": rbin1_out,
                "rbin2": rbin2_out,
                "3pcf": np.zeros((len(rbin), len(rbin))),
                "rbin1_fft": rbin1_out,
                "rbin2_fft": rbin2_out,
                "3pcf_fft": np.zeros((len(rbin), len(rbin))),
                "ell1": bk_in["ell1"],
                "ell1": bk_in["ell2"],
                "ELL": bk_in["ELL"],
                "flag_3pcf": bk_in["flag_3pcf"],
                "flag_BAO": bk_in["flag_BAO"],
                "flag_Recon": bk_in["flag_Recon"],
                "flag_PNG": bk_in["flag_PNG"],
                "flag_FoG": bk_in["flag_FoG"],
                "flag_DampIvanov": bk_in["flag_DampIvanov"],
                "N_fftlog": N_fftlog,
            }

            return output_dict_ini

        ## set bispec. ##
        BB = bk_in["B_fft"]

        ## set kbin ##
        self.kbin = bk_in["kbin1_fft"][:, 0]

        ## set multipole indices ##
        self.ell1 = bk_in["ell1"]
        self.ell2 = bk_in["ell2"]
        self.ELL = bk_in["ELL"]

        ## set rbin ##
        self.rbin = rbin

        ## input parameter for fftlog ##
        NNN = N_fftlog

        ## compute 3PCF ##
        kbin_for_zeta = np.logspace(
            np.log(self.kbin[0]), np.log(self.kbin[-1]), NNN, base=np.e
        )
        CC = np.zeros((len(self.kbin), NNN))
        for i in range(len(self.kbin)):

            f_bk = interpolate.interp1d(
                self.kbin, BB[i, :], fill_value="extrapolate", kind="cubic"
            )
            bk_for_zeta = f_bk(kbin_for_zeta)

            r_temp = np.zeros(NNN)
            zeta_temp = np.zeros(NNN)
            hitomipy.hankel_py(
                self.ell2, 2, NNN, kbin_for_zeta, bk_for_zeta, r_temp, zeta_temp
            )

            f_zeta = interpolate.interp1d(
                r_temp, zeta_temp, fill_value="extrapolate", kind="cubic"
            )
            CC[i, :] = f_zeta(r_temp[:])

        DD = np.zeros((NNN, NNN))
        for j in range(NNN):

            f_bk = interpolate.interp1d(
                self.kbin, CC[:, j], fill_value="extrapolate", kind="cubic"
            )
            bk_for_zeta = f_bk(kbin_for_zeta)

            r_temp = np.zeros(NNN)
            zeta_temp = np.zeros(NNN)
            hitomipy.hankel_py(
                self.ell1, 2, NNN, kbin_for_zeta, bk_for_zeta, r_temp, zeta_temp
            )

            f_zeta = interpolate.interp1d(
                r_temp, zeta_temp, fill_value="extrapolate", kind="cubic"
            )
            DD[:, j] = f_zeta(r_temp[:])

        sign = np.real(1.0j ** (self.ell1 + self.ell2))
        zeta_fft = sign * DD

        f_zeta = interpolate.interp2d(r_temp, r_temp, zeta_fft, kind="cubic")
        zeta_out = np.zeros((len(self.rbin), len(self.rbin)))
        zeta_out[:, :] = f_zeta(self.rbin[:], self.rbin[:])

        ## rbin_out ##
        (rbin2_out, rbin1_out) = np.meshgrid(self.rbin, self.rbin)
        (rbin2_fft, rbin1_fft) = np.meshgrid(r_temp, r_temp)

        ## output dict ##
        output_dict = {
            "rbin1": rbin1_out,
            "rbin2": rbin2_out,
            "3pcf": zeta_out,
            "rbin1_fft": rbin1_fft,
            "rbin2_fft": rbin2_fft,
            "3pcf_fft": zeta_fft,
            "ell1": bk_in["ell1"],
            "ell2": bk_in["ell2"],
            "ELL": bk_in["ELL"],
            "flag_3pcf": bk_in["flag_3pcf"],
            "flag_BAO": bk_in["flag_BAO"],
            "flag_Recon": bk_in["flag_Recon"],
            "flag_PNG": bk_in["flag_PNG"],
            "flag_FoG": bk_in["flag_FoG"],
            "flag_DampIvanov": bk_in["flag_DampIvanov"],
            "N_fftlog": N_fftlog,
        }

        return output_dict

    def calc_3PCF_to_B(self, zeta_in, kbin=np.linspace(0.01, 0.2, 20)):

        if not zeta_in["flag_3pcf"]:
            (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
            output_dict_ini = {
                "kbin1": kbin1_out,
                "kbin2": kbin2_out,
                "B": np.zeros((len(kbin), len(kbin))),
                "ell1": zeta_in["ell1"],
                "ell2": zeta_in["ell2"],
                "ELL": zeta_in["ELL"],
                "flag_3pcf": zeta_in["flag_3pcf"],
                "flag_BAO": zeta_in["flag_BAO"],
                "flag_Recon": zeta_in["flag_Recon"],
                "flag_PNG": zeta_in["flag_PNG"],
                "flag_FoG": zeta_in["flag_FoG"],
                "flag_DampIvanov": zeta_in["flag_DampIvanov"],
                "N_fftlog": zeta_in["N_fftlog"],
            }

            return output_dict_ini

        ## set 3pcf ##
        BB = zeta_in["3pcf_fft"]

        ## set rbin ##
        self.rbin = zeta_in["rbin1_fft"][:, 0]

        ## set multipole indices ##
        self.ell1 = zeta_in["ell1"]
        self.ell2 = zeta_in["ell2"]
        self.ELL = zeta_in["ELL"]

        ## set kbin ##
        self.kbin = kbin

        ## input parameter for fftlog ##
        NNN = zeta_in["N_fftlog"]

        ## compute bispec ##
        rbin_for_bk = np.logspace(
            np.log(self.rbin[0]), np.log(self.rbin[-1]), NNN, base=np.e
        )
        CC = np.zeros((len(self.rbin), len(self.kbin)))
        for i in range(len(self.rbin)):

            f_zeta = interpolate.interp1d(
                self.rbin, BB[i, :], fill_value="extrapolate", kind="cubic"
            )
            zeta_for_bk = f_zeta(rbin_for_bk)

            k_temp = np.zeros(NNN)
            bk_temp = np.zeros(NNN)
            hitomipy.hankel_py(
                self.ell2, 2, NNN, rbin_for_bk, zeta_for_bk, k_temp, bk_temp
            )

            f_bk = interpolate.interp1d(
                k_temp, bk_temp, fill_value="extrapolate", kind="cubic"
            )
            CC[i, :] = f_bk(self.kbin[:])

        DD = np.zeros((len(self.kbin), len(self.kbin)))
        for j in range(len(self.kbin)):

            f_zeta = interpolate.interp1d(
                self.rbin, CC[:, j], fill_value="extrapolate", kind="cubic"
            )
            zeta_for_bk = f_zeta(rbin_for_bk)

            k_temp = np.zeros(NNN)
            bk_temp = np.zeros(NNN)
            hitomipy.hankel_py(
                self.ell1, 2, NNN, rbin_for_bk, zeta_for_bk, k_temp, bk_temp
            )

            f_bk = interpolate.interp1d(
                k_temp, bk_temp, fill_value="extrapolate", kind="cubic"
            )
            DD[:, j] = f_bk(self.kbin[:])

        sign = np.real((-1.0j) ** (self.ell1 + self.ell2)) * (2.0 * np.pi) ** 6
        bk_out = sign * DD

        ## rbin_out ##
        (kbin2_out, kbin1_out) = np.meshgrid(self.kbin, self.kbin)

        ## output dict ##
        output_dict = {
            "kbin1": kbin1_out,
            "kbin2": kbin2_out,
            "B": bk_out,
            "ell1": self.ell1,
            "ell2": self.ell2,
            "ELL": self.ELL,
            "flag_3pcf": zeta_in["flag_3pcf"],
            "flag_BAO": zeta_in["flag_BAO"],
            "flag_Recon": zeta_in["flag_Recon"],
            "flag_PNG": zeta_in["flag_PNG"],
            "flag_FoG": zeta_in["flag_FoG"],
            "flag_DampIvanov": zeta_in["flag_DampIvanov"],
            "N_fftlog": zeta_in["N_fftlog"],
        }

        return output_dict
        return output_dict
        return output_dict
