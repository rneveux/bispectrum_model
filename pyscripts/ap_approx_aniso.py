# -*- coding: utf-8 -*-
import pyximport

import os
import argparse
import numpy as np
from scipy import interpolate
import hitomipy
import pycuba
import os


class ClassApproxAniso:

    def __init__(self, epsilon, estimator, ell1=0, ell2=0, ELL=0, ell1_dash=0, ell2_dash=0, ELL_dash=0, n=0, m=0):

        self.epsilon = epsilon
        self.estimator = estimator
        self.ell1 = ell1
        self.ell2 = ell2
        self.ELL = ELL
        self.ell1_dash = ell1_dash
        self.ell2_dash = ell2_dash
        self.ELL_dash = ELL_dash
        self.n = n
        self.m = m

    def Integrand_T(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])


        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        hitomipy.integrand_SS_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.ell1,
                self.ell2,
                self.ELL,
                self.ell1_dash,
                self.ell2_dash,
                self.ELL_dash,
                self.n,
                self.m, 
                self.epsilon, 
                len(self.epsilon)
            )

        for i in range(ncomp[0]):
                ff[i] = self.ff_out[i]

        return 0

    def Integrand_Tpow(self, ndim, xx, ncomp, ff, userdata):

        self.xx_in = np.zeros(ndim[0])
        self.ff_out = np.zeros(ncomp[0])


        for i in range(ndim[0]):
            self.xx_in[i] = xx[i]

        self.ndim = ndim[0]
        self.ncomp = ncomp[0]

        hitomipy.integrand_SSpow_py(
                self.xx_in,
                self.ndim,
                self.ff_out,
                self.ncomp,
                self.ELL,
                self.ELL_dash,
                self.n,
                self.epsilon, 
                len(self.epsilon)
            )

        for i in range(ncomp[0]):
                ff[i] = self.ff_out[i]

        return 0

    def calc_T(self):

        ## calc. wigner 3j symbols ##
        hitomipy.setWigner3j_py()

        ## compute ##
        NCOMP = len(self.epsilon)
        if NCOMP > 1024:
            print("# of NCOMP should be <= 1024, otherwise results become zero.")
            return output_dict_ini

        if self.estimator=='bk':
            NDIM = 3
            t_temp = pycuba.Cuhre(
                        self.Integrand_T, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4
                    )["results"]
        elif self.estimator=='pk':
            NDIM = 2
            t_temp = pycuba.Cuhre(
                        self.Integrand_Tpow, NDIM, ncomp=NCOMP, key=0, verbose=0 | 4
                    )["results"]

        TT = np.zeros(NCOMP)
        for i in range(NCOMP):
            TT[i] = t_temp[i]['integral']

        #########################

        ## output dict ##
        if self.estimator=='bk':
            output_dict = {
                "epsilon": self.epsilon,
                "T": TT,
                "ell1": self.ell1,
                "ell2": self.ell2,
                "ELL": self.ELL,
                "ell1_dash": self.ell1_dash,
                "ell2_dash": self.ell2_dash,
                "ELL_dash": self.ELL_dash,
                "n": self.n,
                "m": self.m,
            }
        elif self.estimator=='pk':
            output_dict = {
                "epsilon": self.epsilon,
                "T": TT,
                "ELL": self.ELL,
                "ELL_dash": self.ELL_dash,
                "n": self.n,
            }

        return output_dict

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='approximation of AP effect')
    parser.add_argument('-estimator', type=str, help='pk or bk', required=True)
    parser.add_argument('-ell1', type=int, default=0, help='multipoles of target bispectrum - ell1', required=False)
    parser.add_argument('-ell2', type=int, default=0, help='multipoles of target bispectrum - ell2', required=False)
    parser.add_argument('-ELL', type=int, default=0, help='multipoles of target bispectrum or power spectrum - ELL', required=False)
    parser.add_argument('-ell1_dash', type=int, default=0, help='multipoles of contributing bispectrum - ell1', required=False)
    parser.add_argument('-ell2_dash', type=int, default=0, help='multipoles of contributing bispectrum - ell2', required=False)
    parser.add_argument('-ELL_dash', type=int, default=0, help='multipoles of contributing bispectrum - ELL', required=False)
    parser.add_argument('-n', type=int, default=0, help='expansion power of k1 or k', required=False)
    parser.add_argument('-m', type=int, default=0, help='expansion power of k2', required=False)
    cmdline = parser.parse_args()

    epsilon = np.arange(-.04,.04,.00025)
    c = ClassApproxAniso(epsilon,cmdline.estimator,cmdline.ell1,cmdline.ell2,cmdline.ELL,
                            cmdline.ell1_dash,cmdline.ell2_dash,cmdline.ELL_dash,
                            cmdline.n,cmdline.m)
    out = c.calc_T()
    if cmdline.estimator=='bk':
        to_save='bispectrum/theory/approx_epsilon/{}{}{}/{}{}{}/{}{}'.format(cmdline.ell1,cmdline.ell2,cmdline.ELL,
                                                                    cmdline.ell1_dash,cmdline.ell2_dash,cmdline.ELL_dash,
                                                                    cmdline.n,cmdline.m)
    else:
        to_save='bispectrum/theory/approx_epsilon/{}/{}/{}'.format(cmdline.ELL,cmdline.ELL_dash,cmdline.n)
    os.makedirs(os.path.dirname(to_save), exist_ok=True)
    print(out)
    np.save(to_save,out)