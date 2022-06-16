# cython: c_string_type=unicode, c_string_encoding=utf8

import numpy as np
cimport numpy as np

# distutils: language=["c++"]
# distutils: sources=["cpp/fftlog.hpp"]
cdef extern from "cpp/fftlog.hpp":
    int pk2xi(int N, double * k, double * pk, double * r, double * xi)
    int xi2pk(int N, double * r, double * xi, double * k, double * pk)
    int fftlog_ComputeXiLM(int l, int m, int N, double * k, double * pk, double * r, double * xi)

def pk2xi_py(
        int N,
        np.ndarray[double, ndim=1, mode="c"] k not None,
        np.ndarray[double, ndim=1, mode="c"] pk not None,
        np.ndarray[double, ndim=1, mode="c"] r not None,
        np.ndarray[double, ndim=1, mode="c"] xi not None
        ):
    return pk2xi(N, & k[0], & pk[0], & r[0], & xi[0])


def xi2pk_py(
        int N,
        np.ndarray[double, ndim=1, mode="c"] r not None,
        np.ndarray[double, ndim=1, mode="c"] xi not None,
        np.ndarray[double, ndim=1, mode="c"] k not None,
        np.ndarray[double, ndim=1, mode="c"] pk not None
        ):
    return xi2pk(N, & r[0], & xi[0], & k[0], & pk[0])


def hankel_py(
        int l, int m, int N,
        np.ndarray[double, ndim=1, mode="c"] k not None,
        np.ndarray[double, ndim=1, mode="c"] pk not None,
        np.ndarray[double, ndim=1, mode="c"] r not None,
        np.ndarray[double, ndim=1, mode="c"] xi not None
        ):
    return fftlog_ComputeXiLM(l, m, N, & k[0], & pk[0], & r[0], & xi[0])


# distutils: sources=["cpp/wigner.hpp"]
cdef extern from "cpp/wigner.hpp":
    int setWigner3j()

def setWigner3j_py():
    return setWigner3j()

# distutils: sources=["cpp/pk_lin.hpp"]
cdef extern from "cpp/pk_lin.hpp":
    int readInputPowerSpectrum(double * kbin_in, double * pk_in, int pk_num_in)
    int readInputNoWigglePowerSpectrum(double * kbin_in, double * pk_in, int pk_num_in)
    int readInputPrimordialPowerSpectrum(double * kbin_in, double * pk_in, int pk_num_in)

    int readNonLinearPowerSpectrum(double * kmag_in, double * pk_0_in, double * pk_2_in, int pk_num_in)

    void initializeInputPowerSpectrum()
    void finalizeInputPowerSpectrum()
    void calcNormalizationUsingSigma8(double sigma8)
    void calcNormalizationNoWiggle(double sigma8, double h, double omegab, double omega0, double Tcmb, double n_s)
    double set_kmin(double kmin)
    double set_kmax(double kmax)
    double calcSigma_dd(double sigma8)
    double f_pk_no_wiggle_integrand(double kmag, double h, double omegab, double omega0, double Tcmb, double n_s)


def f_pk_no_wiggle_integrand_py(double kmag, double h, double omegab, double omega0, double Tcmb, double n_s):
    return f_pk_no_wiggle_integrand(kmag, h, omegab, omega0, Tcmb, n_s)

def calcSigma_dd_py(double sigma8):
    return calcSigma_dd(sigma8)


def readInputPowerSpectrum_py(
        np.ndarray[double, ndim=1, mode="c"] kbin_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_in not None,
        int pk_num_in
        ):
    return readInputPowerSpectrum( & kbin_in[0], & pk_in[0], pk_num_in)

def readInputNoWigglePowerSpectrum_py(
        np.ndarray[double, ndim=1, mode="c"] kbin_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_in not None,
        int pk_num_in
        ):
    return readInputNoWigglePowerSpectrum( & kbin_in[0], & pk_in[0], pk_num_in)


def readInputPrimordialPowerSpectrum_py(
        np.ndarray[double, ndim=1, mode="c"] kbin_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_in not None,
        int pk_num_in
        ):
    return readInputPrimordialPowerSpectrum( & kbin_in[0], & pk_in[0], pk_num_in)


def readNonLinearPowerSpectrum_py(
        np.ndarray[double, ndim=1, mode="c"] kbin_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_0_in not None,
        np.ndarray[double, ndim=1, mode="c"] pk_2_in not None,
        int pk_num_in
        ):
    return readNonLinearPowerSpectrum(& kbin_in[0], & pk_0_in[0], & pk_2_in[0], pk_num_in)


def initializeInputPowerSpectrum_py():
    return initializeInputPowerSpectrum()


def finalizeInputPowerSpectrum_py():
    return finalizeInputPowerSpectrum()


def calcNormalizationUsingSigma8_py(double sigma8):
    return calcNormalizationUsingSigma8(sigma8)


def calcNormalizationNoWiggle_py(double sigma8, double h, double omegab, double omega0, double Tcmb, double n_s):
    return calcNormalizationNoWiggle(sigma8, h, omegab, omega0, Tcmb, n_s)


def set_kmin_py(double kmin):
    return set_kmin(kmin)


def set_kmax_py(double kmax):
    return set_kmax(kmax)


# distutils: sources=["cpp/calc_P.hpp"]
cdef extern from "cpp/calc_P.hpp":

    int integrand_P_Tree(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1)

    int integrand_P_Damping_Tree(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, double dSigma2)
    
    int integrand_P_Damping_Tree_for_1loop(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, double dSigma2)

    int integrand_P_Damping_Tree_gm(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, double dSigma2)

    int integrand_P_Counterterm(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double c0, double c1, double c2, double ch, double knl)

    int integrand_P_Damping_Counterterm(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double c0, double c1, double c2, double ch, double knl,
            double Sigma2, double dSigma2)

    int integrand_P_1loop(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2, double bGamma3)

    int integrand_P_Damping_1loop(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2, double bGamma3,
            double Sigma2, double dSigma2)
    
    int integrand_P_Damping_1loop_k_vector(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kx, double * ky, double * kz, int num_k_bin,
            double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2, double bGamma3,
            double Sigma2, double dSigma2)
    
    int integrand_P_Damping_1loop_nw(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2, double bGamma3,
            double Sigma2, double dSigma2)

    int integrand_P_Damping_1loop_gm(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2, double bGamma3,
            double Sigma2, double dSigma2)

    int integrand_P_Kernel_Tree(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double f, double Sigma2, double dSigma2,char * parameters)

    int integrand_P_Kernel_Counterterm(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double f, double Sigma2, double dSigma2,char * parameters)

    int integrand_P_Kernel_1loop_22(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double f, double Sigma2, double dSigma2,char * parameters)

    int integrand_P_Kernel_1loop_22_norm(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double f, double Sigma2, double dSigma2)

    int integrand_P_Kernel_1loop_13(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double f, double Sigma2, double dSigma2,char * parameters)

    int integrand_P_Tree_NoWiggle(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1)

    int integrand_P_Tree_BAO(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double sigma2_perp, double sigma2_para)

    int integrand_P_Tree_BAO_Fitting(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double sigma2_perp, double sigma2_para,
            double A20, double A11, double A02, 
            double A30, double A21, double A12, double A03,
            double A40, double A31, double A22, double A13, double A04)

    int integrand_P_Tree_Window(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double volume)

    int integrand_P_Tree_Window_IC(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double volume)

    int integrand_P_Tree_Window_IC_Approx(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double volume)


    int integrand_P_NonLinearFitting(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double sigma2_perp, double sigma2_para)

    int integrand_P_NonLinearFitting_Window(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double sigma2_perp, double sigma2_para,
            double volume)

    int integrand_P_LocalMean(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2,
            double sigma2_perp, double sigma2_para, 
            double nmean, double volume)

    int integrand_P_SigmaB(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double sigma2_perp, double sigma2_para, 
            double nmean, double volume)

    int integrand_P_sigma2_perp_Reconstructed(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, 
            double sigma8, double fz, double b1,
            double b1_fid, double R)

    int integrand_P_sigma2_para_Reconstructed(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, 
            double sigma8, double fz, double b1,
            double b1_fid, double R)

    int integrand_P_Tree_BAO_b1_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_P_Tree_BAO_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)
 
    int integrand_P_Tree_BAO_f_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ELL,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)
 
# distutils: sources=["cpp/calc_cov_PP.hpp"]
cdef extern from "cpp/calc_cov_PP.hpp":

    int integrand_cov_PP_G(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_G_NL(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double sigma2_perp, double sigma2_para, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double b2, double b3, double bK2, double bK3, double bDK, double bO,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_b2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_bK2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_b2_b2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_b2_bK2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_bK2_bK2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_b3(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_bK3(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_bDK(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_bO(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

########################################
    
    int integrand_cov_PP_NG_BeatCoupling(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double b3, double bK2, double bK3, double bDK, double bO,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_b2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_bK2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_b2_b2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_b2_bK2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_bK2_bK2(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_b3(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_bK3(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_bDK(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_BeatCoupling_bO(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double DeltaK, double nmean, double volume)

##################################################

    int integrand_cov_PP_NG_LocalMean(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double b3, double bK2, double bK3, double bDK, double bO, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_LocalMean_NL(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double b3, double bK2, double bK3, double bDK, double bO, 
            double sigma2_perp, double sigma2_para, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_LocalMean_NL_Sigma2B(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double b3, double bK2, double bK3, double bDK, double bO, 
            double sigma2_perp, double sigma2_para,
            double DeltaK, double nmean, double volume, double sigma2_b)

    int integrand_cov_PP_NG_LocalMean_NL_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double sigma2_perp, double sigma2_para,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_LocalMean_NL_b2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double sigma2_perp, double sigma2_para,
            double DeltaK, double nmean, double volume)

    int integrand_cov_PP_NG_LocalMean_NL_bK2(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double sigma2_perp, double sigma2_para,
            double DeltaK, double nmean, double volume)

# distutils: sources=["cpp/calc_cov_BB.hpp"]
cdef extern from "cpp/calc_cov_BB.hpp":

    int integrand_cov_BB_G(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, double k1, double k1_dash,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume)

    int integrand_cov_BB_G_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, double k1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume)

# distutils: sources=["cpp/calc_cov_PB.hpp"]
cdef extern from "cpp/calc_cov_PB.hpp":

    int integrand_cov_PB_NG_PB_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, 
            double DeltaK, double nmean, double volume)

    int integrand_cov_PB_NG_P5_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume)


############
############
############

# distutils: sources=["cpp/calc_B.hpp"]
cdef extern from "cpp/calc_B.hpp":

    int integrand_B_Kernel(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, double f, 
            double Sigma2, double dSigma2, char * parameters)

    int integrand_B_Kernel_SN(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, double f, 
            double Sigma2, double dSigma2, char * parameters)

    int integrand_B_Kernel_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double f,
            double Sigma2, double dSigma2, double alpha_perp, double alpha_parallel, char * parameters)

    int integrand_B_Kernel_SN_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double f,
            double Sigma2, double dSigma2, double alpha_perp, double alpha_parallel, char * parameters)

    int integrand_B_Kernel_PNG_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double f,
            double Sigma2, double dSigma2, char * parameters)

    int integrand_B_NonGaussian_Local(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1)

    int integrand_B_NonGaussian_Equilateral(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1)

    int integrand_B_NonGaussian_Orthogonal(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1)


    int integrand_B_Tree(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double bK2)

    int integrand_B_Tree_FoG(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double bK2,
            double c1, double c2, double knl)

    int integrand_B_Tree_DampIvanov(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double bK2,
            double rbao, double ks)

    int integrand_B_FoG_Damping_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1, 
            double b2, double bG2,
            double c1, double c2, double knl,
            double Sigma2, double dSigma2)

    int integrand_B_FoG_Damping(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, double f, double b1, 
            double b2, double bG2,
            double c1, double c2, double knl,
            double Sigma2, double dSigma2)

    int integrand_B_SN_FoG_Damping_diag(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double alpha_perp, double alpha_parallel, double f, double b1,
            double c1, double c2, double knl, double Pshot, double Bshot,
            double Sigma2, double dSigma2)
    
    int integrand_B_SN_FoG_Damping(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, double f, double b1, 
            double c1, double c2, double knl, double Pshot, double Bshot,
            double Sigma2, double dSigma2)
 
    int integrand_B_Tree_NoWiggle(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double b2, double bK2)


    int integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Growth(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double sigma8)

    int integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Shift(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double sigma8)

    int integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Tidal(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double sigma8)


    int integrand_B_Tree_BAO_RealSpace_DarkMatter_Growth(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double sigma8, double sigma2_perp)

    int integrand_B_Tree_BAO_RealSpace_DarkMatter_Shift(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double sigma8, double sigma2_perp)

    int integrand_B_Tree_BAO_RealSpace_DarkMatter_Tidal(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double sigma8, double sigma2_perp)




    int integrand_B_Tree_Reconstructed(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
            double b2, double bK2,
            double b1_fid, double R)

    int integrand_B_Tree_BAO(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double bK2, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_NoWiggle_Reconstructed(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
            double b2, double bK2, 
            double b1_fid, double R)

    # Tree BAO for fitting #
    int integrand_B_Tree_BAO_Template( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para,
            char * parameters)

    int integrand_B_Tree_BAO_b1_b1_b1( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)
    
    int integrand_B_Tree_BAO_b1_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel,
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_b1_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_b2_b1_b1( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_b2_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel,
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_b2_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_bK2_b1_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_bK2_b1_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_bK2_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_b1f_b1_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_b1f_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_b1f_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double kmag1, double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_ff_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double kmag1, double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

    int integrand_B_Tree_BAO_f_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel, 
            double sigma2_perp, double sigma2_para)

############
############
############

    # Tree for fitting #
    int integrand_B_Tree_b1_b1_b1( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)
    
    int integrand_B_Tree_b1_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_b1_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_b2_b1_b1( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_b2_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_b2_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_bK2_b1_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_bK2_b1_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_bK2_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_b1f_b1_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_b1f_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_b1f_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double kmag1, double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_ff_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double kmag1, double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_f_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

############
############
############

    # Tree NoWiggle for fitting #
    int integrand_B_Tree_NoWiggle_b1_b1_b1( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)
    
    int integrand_B_Tree_NoWiggle_b1_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_b1_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_b2_b1_b1( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_b2_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_b2_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_bK2_b1_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_bK2_b1_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_bK2_f_f( 
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_b1f_b1_b1(
            double * xx_in, int ndim, double * ff_out, int ncomp,
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_b1f_b1_f(
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
            double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_b1f_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double kmag1, double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_ff_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
            double kmag1, double alpha_perp, double alpha_parallel)

    int integrand_B_Tree_NoWiggle_f_f_f(  
            double * xx_in, int ndim, double * ff_out, int ncomp, 
            double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
            double alpha_perp, double alpha_parallel)

    int integrand_SS(
            double * xx_in, int ndim, double * ff_out, int ncomp,  
            int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
            int n, int m,
            double * epsilon, int num_epsilon) 

    int integrand_SSpow(
            double * xx_in, int ndim, double * ff_out, int ncomp,  
            int ELL, int ELL_dash, 
            int n,
            double * epsilon, int num_epsilon) 

# distutils: sources=["cpp/kernel.hpp"]
cdef extern from "cpp/kernel.hpp":

    double Sig2(double rbao, double ks)

    double dSig2(double rbao, double ks)

    double Bispectrum_Tree_FoG_Damping(
            double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2,
            double c1, double c2, double knl, double Sigma2, double dSigma2)
    
    double Bispectrum_SN_FoG_Damping(
            double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1,
            double c1, double c2, double knl, double Pshot, double Bshot, double Sigma2, double dSigma2)
    
    double Powerspectrum_Tree_Damping_for_1loop(
            double * kvec1_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1,
            double Sigma2, double dSigma2)   
    
    double Powerspectrum_Counterterm_Damping(
            double * kvec1_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1,
            double c0, double c1, double c2, double ch, double knl, 
            double Sigma2, double dSigma2)


# distutils: sources=["cpp/sph.hpp"]
cdef extern from "cpp/sph.hpp":

    double calcYYY(
            double * kvec1, double * kvec2, double * los, int ell1, int ell2, int ELL)

def Sig2_py(
        double rbao, double ks):
    return Sig2(rbao, ks)

def dSig2_py(
        double rbao, double ks):
    return dSig2(rbao, ks)

def integrand_P_LocalMean_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2,
        double sigma2_perp, double sigma2_para, 
        double nmean, double volume
        ):
    return integrand_P_LocalMean(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2,
            sigma2_perp, sigma2_para, 
            nmean, volume)

def integrand_P_Tree_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1
        ):
    return integrand_P_Tree(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1)

def integrand_P_Damping_Tree_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, dSigma2
        ):
    return integrand_P_Damping_Tree(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, Sigma2, dSigma2)

def integrand_P_Damping_Tree_for_1loop_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, dSigma2
        ):
    return integrand_P_Damping_Tree_for_1loop(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, Sigma2, dSigma2)

'''def integrand_P_Damping_Tree_gm_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, dSigma2
        ):
    return integrand_P_Damping_Tree_gm(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, Sigma2, dSigma2)

def integrand_P_Counterterm_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double c0, double c1, double c2, double knl,
        ):
    return integrand_P_Counterterm(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, 
            c0, c1, c2, knl)
'''
def integrand_P_Damping_Counterterm_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double c0, double c1, double c2, double ch, double knl,
        double Sigma2, dSigma2,
        ):
    return integrand_P_Damping_Counterterm(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, 
            c0, c1, c2, ch, knl, 
            Sigma2, dSigma2)
'''
def integrand_P_Damping_Counterterm_gm_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double c0, double c1, double c2, double knl,
        double Sigma2, dSigma2,
        ):
    return integrand_P_Damping_Counterterm_gm(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, 
            c0, c1, c2, knl, 
            Sigma2, dSigma2)'''

def integrand_P_1loop_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double b2, double bG2, double bGamma3,
        ):
    return integrand_P_1loop(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, 
            b2, bG2, bGamma3)

def integrand_P_Damping_1loop_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double b2, double bG2, double bGamma3,
        double Sigma2, dSigma2,
        ):
    return integrand_P_Damping_1loop(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, 
            b2, bG2, bGamma3,
            Sigma2, dSigma2)

def integrand_P_Damping_1loop_k_vector_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kx not None, np.ndarray[double, ndim=1, mode="c"] ky not None, 
        np.ndarray[double, ndim=1, mode="c"] kz not None, int num_k,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double b2, double bG2, double bGamma3,
        double Sigma2, dSigma2,
        ):
    return integrand_P_Damping_1loop_k_vector(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kx[0], &ky[0], &kz[0], num_k,
            alpha_perp, alpha_parallel, f, b1, 
            b2, bG2, bGamma3,
            Sigma2, dSigma2)

def integrand_P_Damping_1loop_nw_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double b2, double bG2, double bGamma3,
        double Sigma2, dSigma2,
        ):
    return integrand_P_Damping_1loop_nw(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, 
            b2, bG2, bGamma3,
            Sigma2, dSigma2)

'''def integrand_P_Damping_1loop_gm_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double f, double b1, 
        double b2, double bG2, double bGamma3,
        double Sigma2, dSigma2,
        ):
    return integrand_P_Damping_1loop_gm(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, f, b1, 
            b2, bG2, bGamma3,
            Sigma2, dSigma2)'''

def integrand_P_Kernel_Tree_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double f, double Sigma2, dSigma2,  char * parameters
        ):
    return integrand_P_Kernel_Tree(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            f, Sigma2, dSigma2,
            &parameters[0])

def integrand_P_Kernel_Counterterm_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double f, double Sigma2, dSigma2,  char * parameters
        ):
    return integrand_P_Kernel_Counterterm(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            f, Sigma2, dSigma2,
            &parameters[0])

def integrand_P_Kernel_1loop_22_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double f, double Sigma2, dSigma2,  char * parameters
        ):
    return integrand_P_Kernel_1loop_22(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            f, Sigma2, dSigma2,
            &parameters[0])

def integrand_P_Kernel_1loop_22_norm_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double f, double Sigma2, dSigma2
        ):
    return integrand_P_Kernel_1loop_22_norm(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            f, Sigma2, dSigma2)

def integrand_P_Kernel_1loop_13_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double f, double Sigma2, dSigma2,  char * parameters
        ):
    return integrand_P_Kernel_1loop_13(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            f, Sigma2, dSigma2,
            &parameters[0])

def integrand_P_Tree_NoWiggle_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1
        ):
    return integrand_P_Tree_NoWiggle(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, sigma8, fz, b1)

def integrand_P_Tree_BAO_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double sigma2_perp, double sigma2_para):
    return integrand_P_Tree_BAO( 
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            sigma2_perp, sigma2_para)

def integrand_P_Tree_BAO_Fitting_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double sigma2_perp, double sigma2_para,
        double A20, double A11, double A02, 
        double A30, double A21, double A12, double A03,
        double A40, double A31, double A22, double A13, double A04):
    return integrand_P_Tree_BAO_Fitting( 
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            sigma2_perp, sigma2_para,
            A20, A11, A02, 
            A30, A21, A12, A03,
            A40, A31, A22, A13, A04)


def integrand_P_NonLinearFitting_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double sigma2_perp, double sigma2_para
        ):
    return integrand_P_NonLinearFitting(
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            sigma2_perp, sigma2_para)


def integrand_P_NonLinearFitting_Window_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double sigma2_perp, double sigma2_para, 
        double volume
        ):
    return integrand_P_NonLinearFitting_Window(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            sigma2_perp, sigma2_para, 
            volume)

def integrand_P_SigmaB_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double sigma2_perp, double sigma2_para, 
        double nmean, double volume
        ):
    return integrand_P_SigmaB(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            sigma2_perp, sigma2_para, 
            nmean, volume)


def integrand_P_sigma2_perp_Reconstructed_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, 
        double sigma8, double fz, double b1,
        double b1_fid, double R
        ):
    return integrand_P_sigma2_perp_Reconstructed(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, 
            sigma8, fz, b1,
            b1_fid, R)


def integrand_P_sigma2_para_Reconstructed_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, 
        double sigma8, double fz, double b1,
        double b1_fid, double R
        ):
    return integrand_P_sigma2_para_Reconstructed(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, 
            sigma8, fz, b1,
            b1_fid, R)



# def integrand_P_Tree_b1_b1_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                              np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel):
#    return integrand_P_Tree_b1_b1( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel)
#
#
# def integrand_P_Tree_b1_f_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                             np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel):
#    return integrand_P_Tree_b1_f( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel)
#
#
# def integrand_P_Tree_f_f_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                            np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel):
#    return integrand_P_Tree_f_f( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel)
#
#
def integrand_P_Tree_BAO_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL,
        double alpha_perp, double alpha_parallel,
        double sigma2_perp, double sigma2_para
        ):
    return integrand_P_Tree_BAO_b1_b1( 
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)


def integrand_P_Tree_BAO_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_P_Tree_BAO_b1_f(
            &xx_in[0], ndim, & ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)


def integrand_P_Tree_BAO_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, 
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_P_Tree_BAO_f_f( 
            &xx_in[0], ndim, & ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)


# def integrand_P_Tree_NoWiggle_b1_b1_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                       np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel):
#    return integrand_P_Tree_NoWiggle_b1_b1( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel)
#
#
# def integrand_P_Tree_NoWiggle_b1_f_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                      np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel):
#    return integrand_P_Tree_NoWiggle_b1_f( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel)
#
#
# def integrand_P_Tree_NoWiggle_f_f_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                     np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel):
#    return integrand_P_Tree_NoWiggle_f_f( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel)
#
#

# def integrand_P_Tree_Aniso_real_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                   np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int _L_, int _M_, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double g20, double g21_real, double g21_imag, double g22_real, double g22_imag):
#    return integrand_P_Tree_Aniso_real( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ell1, ell2, _L_, _M_, alpha_perp, alpha_parallel, sigma8, fz, b1, g20, g21_real, g21_imag, g22_real, g22_imag)
#
#
# def integrand_P_Tree_Aniso_imag_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                   np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int _L_, int _M_, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double g20, double g21_real, double g21_imag, double g22_real, double g22_imag):
#    return integrand_P_Tree_Aniso_imag( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ell1, ell2, _L_, _M_, alpha_perp, alpha_parallel, sigma8, fz, b1, g20, g21_real, g21_imag, g22_real, g22_imag)
#
#

# def integrand_P_Tree_Window_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double volume):
#    return integrand_P_Tree_Window( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, volume)
#
#
# def integrand_P_Tree_Window_IC_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                  np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double volume):
#    return integrand_P_Tree_Window_IC( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, volume)
#
#
# def integrand_P_Tree_Window_IC_Approx_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                         np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double volume):
#    return integrand_P_Tree_Window_IC_Approx( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, volume)
#
#
# def integrand_P_Tree_KSZ_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                            np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                            aH_tau_T0_over_c):
#    return integrand_P_Tree_KSZ( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, aH_tau_T0_over_c)
#
#
# def integrand_P_Tree_A_perp_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_A_perp( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_A_para_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_A_para( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_B_perp_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_B_perp( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_B_para_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_B_para( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_B_perp_para_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                                    np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                                    double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_B_perp_para( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_Poly_1_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_Poly_1( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_Poly_2_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_Poly_2( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_Poly_3_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_Poly_3( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_Tree_Poly_4_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                               np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                               double sigma2_perp, double sigma2_para):
#    return integrand_P_Tree_Poly_4( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para)
#
#
# def integrand_P_SPT1loop_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                            np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                            double b2, double b3, double bK2, double bK3, double bDK, double bO,
#                            double sigma2_perp, double sigma2_para):
#    return integrand_P_SPT1loop( & xx_in[0], ndim, & ff_out[0], ncomp, & kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, sigma2_perp, sigma2_para)
#
#
# def integrand_P_SPT2loop_py(np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
#                            np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
#                            double b2, double b3, double bK2, double bK3, double bDK, double bO,
#                            double sigma2_perp, double sigma2_para):
# return integrand_P_SPT2loop( & xx_in[0], ndim, & ff_out[0], ncomp, &
# kbin[0], num_k_bin, ELL, alpha_perp, alpha_parallel, sigma8, fz, b1, b2,
# b3, bK2, bK3, bDK, bO, sigma2_perp, sigma2_para)


#    int integrand_P_Tree_KSZ(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,\
#                         double aH_tau_T0_over_c)
#    int integrand_P_Tree_Aniso_real(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int _L_, int _M_, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double g20, double g21_real, double g21_imag, double g22_real, double g22_imag)
#    int integrand_P_Tree_Aniso_imag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int _L_, int _M_, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double g20, double g21_real, double g21_imag, double g22_real, double g22_imag)
#
#    int integrand_P_Tree_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel)
#    int integrand_P_Tree_b1_f( double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel)
#    int integrand_P_Tree_f_f(  double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel)
#
#    int integrand_P_Tree_BAO_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL,\
#                               double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para)
#    int integrand_P_Tree_BAO_b1_f( double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL,\
#                               double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para)
#    int integrand_P_Tree_BAO_f_f(  double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL,\
#                               double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para)
#
#    int integrand_P_Tree_NoWiggle_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel)
#    int integrand_P_Tree_NoWiggle_b1_f( double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel)
#    int integrand_P_Tree_NoWiggle_f_f(  double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel)
#
#
#    int integrand_P_Tree_A_perp(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,\
#                         double sigma2_perp, double sigma2_para)
#
#    int integrand_P_Tree_A_para(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,\
#                         double sigma2_perp, double sigma2_para)
#
#    int integrand_P_Tree_B_perp(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,\
#                         double sigma2_perp, double sigma2_para)
#
#    int integrand_P_Tree_B_para(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,\
#                         double sigma2_perp, double sigma2_para)
#
#    int integrand_P_Tree_B_perp_para(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,\
#                         double sigma2_perp, double sigma2_para)
#
#    int integrand_P_SPT1loop(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin,
#                      int ELL, double alpha_perp, double alpha_parallel,
#                  double sigma8, double fz, double b1, double b2, double b3, double bK2,
#             double bK3, double bDK, double bO, double sigma2_perp, double sigma2_para)
#    int integrand_P_SPT2loop(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin,
#                      int ELL, double alpha_perp, double alpha_parallel,
#                  double sigma8, double fz, double b1, double b2, double b3, double bK2,
#             double bK3, double bDK, double bO, double sigma2_perp, double sigma2_para)
#


def integrand_cov_PP_G_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_G(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash,
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG( 
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            b2, b3, bK2, bK3, bDK, bO,
            DeltaK, nmean, volume)

def integrand_cov_PP_G_NL_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double sigma2_perp, double sigma2_para, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_G_NL(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            sigma2_perp, sigma2_para, 
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_b2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_b2(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_bK2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_bK2(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_b2_b2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_b2_b2( 
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_b2_bK2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_b2_bK2(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_bK2_bK2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_bK2_bK2( 
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_b3_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_b3( 
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_bK3_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_bK3( 
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_bDK_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_bDK( 
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_bO_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_bO( 
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)

#####################################

def integrand_cov_PP_NG_BeatCoupling_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            b2, b3, bK2, bK3, bDK, bO,
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_BeatCoupling_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_BeatCoupling_b2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_b2(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_BeatCoupling_bK2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_bK2(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_BeatCoupling_b2_b2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_b2_b2(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_BeatCoupling_b2_bK2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_b2_bK2(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_BeatCoupling_bK2_bK2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_bK2_bK2(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_BeatCoupling_b3_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_b3(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_BeatCoupling_bK3_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_bK3(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_BeatCoupling_bDK_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_bDK(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_BeatCoupling_bO_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_BeatCoupling_bO(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)

#################################

def integrand_cov_PP_NG_LocalMean_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_LocalMean(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            b2, b3, bK2, bK3, bDK, bO,
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_LocalMean_NL_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double sigma2_perp, double sigma2_para,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_LocalMean_NL(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            b2, b3, bK2, bK3, bDK, bO,
            sigma2_perp, sigma2_para, 
            DeltaK, nmean, volume)

def integrand_cov_PP_NG_LocalMean_NL_Sigma2B_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double sigma2_perp, double sigma2_para,
        double DeltaK, double nmean, double volume, double sigma2_b
        ):
    return integrand_cov_PP_NG_LocalMean_NL_Sigma2B(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            b2, b3, bK2, bK3, bDK, bO,
            sigma2_perp, sigma2_para, 
            DeltaK, nmean, volume, sigma2_b)

def integrand_cov_PP_NG_LocalMean_NL_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double sigma2_perp, double sigma2_para,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_LocalMean_NL_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            sigma2_perp, sigma2_para, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_LocalMean_NL_b2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1, 
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double sigma2_perp, double sigma2_para,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_LocalMean_NL_b2(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            sigma2_perp, sigma2_para, 
            DeltaK, nmean, volume)


def integrand_cov_PP_NG_LocalMean_NL_bK2_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double sigma2_perp, double sigma2_para,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PP_NG_LocalMean_NL_bK2(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ELL_dash, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            sigma2_perp, sigma2_para, 
            DeltaK, nmean, volume)

################################
################################
################################

def integrand_B_NonGaussian_Local_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1):
    return integrand_B_NonGaussian_Local(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1)

def integrand_B_NonGaussian_Equilateral_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1):
    return integrand_B_NonGaussian_Equilateral(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1)


def integrand_B_NonGaussian_Orthogonal_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1):
    return integrand_B_NonGaussian_Orthogonal(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1)


def integrand_B_Tree_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2
        ):
    return integrand_B_Tree(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            b2, bK2)

def integrand_B_Tree_FoG_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2,
        double c1, double c2, double knl,
        ):
    return integrand_B_Tree_FoG(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            b2, bK2,
            c1, c2, knl)

def integrand_B_Tree_DampIvanov_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2,
        double rbao, double ks,
        ):
    return integrand_B_Tree_DampIvanov(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            b2, bK2,
            rbao, ks)

def integrand_B_FoG_Damping_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, 
        double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2,
        double c1, double c2, double knl,
        double Sigma2, double dSigma2
        ):
    return integrand_B_FoG_Damping_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, alpha_perp, alpha_parallel, f, b1,
            b2, bG2,
            c1, c2, knl,
            Sigma2, dSigma2)

def integrand_B_FoG_Damping_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2,
        double c1, double c2, double knl,
        double Sigma2, double dSigma2
        ):
    return integrand_B_FoG_Damping(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, alpha_perp, alpha_parallel, f, b1,
            b2, bG2,
            c1, c2, knl,
            Sigma2, dSigma2)

def integrand_B_SN_FoG_Damping_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, 
        double alpha_perp, double alpha_parallel, double f, double b1,
        double c1, double c2, double knl, double Pshot, double Bshot,
        double Sigma2, double dSigma2
        ):
    return integrand_B_SN_FoG_Damping_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, alpha_perp, alpha_parallel, f, b1,
            c1, c2, knl, Pshot, Bshot,
            Sigma2, dSigma2)

def integrand_B_SN_FoG_Damping_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double f, double b1,
        double c1, double c2, double knl, double Pshot, double Bshot,
        double Sigma2, double dSigma2
        ):
    return integrand_B_SN_FoG_Damping(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, alpha_perp, alpha_parallel, f, b1,
            c1, c2, knl,  Pshot, Bshot,
            Sigma2, dSigma2)

def integrand_B_Tree_NoWiggle_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double b2, double bK2
        ):
    return integrand_B_Tree_NoWiggle(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            b2, bK2)

def integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Growth_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double sigma8):
    return integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Growth(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, sigma8)

def integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Shift_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double sigma8):
    return integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Shift(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, sigma8)

def integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Tidal_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double sigma8):
    return integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Tidal(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, sigma8)

def integrand_B_Tree_BAO_RealSpace_DarkMatter_Growth_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double sigma8, double sigma2_perp):
    return integrand_B_Tree_BAO_RealSpace_DarkMatter_Growth(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, sigma8, sigma2_perp)

def integrand_B_Tree_BAO_RealSpace_DarkMatter_Shift_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double sigma8, double sigma2_perp):
    return integrand_B_Tree_BAO_RealSpace_DarkMatter_Shift(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, sigma8, sigma2_perp)


def integrand_B_Tree_BAO_RealSpace_DarkMatter_Tidal_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double sigma8, double sigma2_perp):
    return integrand_B_Tree_BAO_RealSpace_DarkMatter_Tidal(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, sigma8, sigma2_perp)


def integrand_B_Tree_BAO_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double b2, double bK2, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            b2, bK2,
            sigma2_perp, sigma2_para)

def integrand_B_Tree_Reconstructed_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double b2, double bK2, 
        double b1_fid, double R
        ):
    return integrand_B_Tree_Reconstructed(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            b2, bK2,
            b1_fid, R)

def integrand_B_Tree_NoWiggle_Reconstructed_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
        double b2, double bK2, 
        double b1_fid, double R
        ):
    return integrand_B_Tree_NoWiggle_Reconstructed(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            b2, bK2,
            b1_fid, R)

def integrand_B_Tree_BAO_Template_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para,
        char * parameters):
    return integrand_B_Tree_BAO_Template(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para,
            &parameters[0])


def integrand_B_Tree_BAO_b1_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b1_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b1_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b1_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel,
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b1_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b1_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b2_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b2_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b2_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b2_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b2_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b2_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_bK2_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_bK2_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_bK2_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_bK2_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_bK2_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_bK2_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b1f_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b1f_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b1f_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b1f_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_b1f_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_b1f_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_ff_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel,
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_ff_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel, 
            sigma2_perp, sigma2_para)

def integrand_B_Tree_BAO_f_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel, 
        double sigma2_perp, double sigma2_para
        ):
    return integrand_B_Tree_BAO_f_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel,
            sigma2_perp, sigma2_para)

###################
###################
###################

def integrand_B_Tree_b1_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b1_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b1_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b1_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b1_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b1_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b2_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b2_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b2_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b2_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b2_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b2_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_bK2_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel):
    return integrand_B_Tree_bK2_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_bK2_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_bK2_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_bK2_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_bK2_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b1f_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b1f_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b1f_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b1f_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_b1f_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_b1f_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_ff_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_ff_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_f_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_f_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)


###################
###################
###################

def integrand_B_Tree_NoWiggle_b1_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b1_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b1_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b1_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b1_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b1_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b2_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b2_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b2_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b2_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b2_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b2_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_bK2_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel):
    return integrand_B_Tree_NoWiggle_bK2_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_bK2_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_bK2_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_bK2_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_bK2_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b1f_b1_b1_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b1f_b1_b1(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b1f_b1_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b1f_b1_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1,
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_b1f_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_b1f_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_ff_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_ff_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)

def integrand_B_Tree_NoWiggle_f_f_f_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
        double alpha_perp, double alpha_parallel
        ):
    return integrand_B_Tree_NoWiggle_f_f_f(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, 
            alpha_perp, alpha_parallel)


def integrand_SS_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
        int n, int m,
        np.ndarray[double, ndim=1, mode="c"] epsilon not None, int num_epsilon
        ):
    return integrand_SS(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            ell1, ell2, ELL, ell1_dash, ell2_dash, ELL_dash,
            n, m,
            &epsilon[0], num_epsilon)

def integrand_SSpow_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        int ELL, int ELL_dash, 
        int n,
        np.ndarray[double, ndim=1, mode="c"] epsilon not None, int num_epsilon
        ):
    return integrand_SSpow(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            ELL, ELL_dash,
            n,
            &epsilon[0], num_epsilon)

def integrand_cov_BB_G_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, double k1, double k1_dash,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_BB_G(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, ell1_dash, ell2_dash, ELL_dash, k1, k1_dash,
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)

def integrand_cov_BB_G_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, double k1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_BB_G_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ell1, ell2, ELL, ell1_dash, ell2_dash, ELL_dash, k1,
            alpha_perp, alpha_parallel, sigma8, fz, b1, 
            DeltaK, nmean, volume)

def integrand_cov_PB_NG_PB_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PB_NG_PB_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ell1_dash, ell2_dash, ELL_dash,
            alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2,
            DeltaK, nmean, volume)

def integrand_cov_PB_NG_P5_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double DeltaK, double nmean, double volume
        ):
    return integrand_cov_PB_NG_P5_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp, 
            &kbin[0], num_k_bin, ELL, ell1_dash, ell2_dash, ELL_dash,
            alpha_perp, alpha_parallel, sigma8, fz, b1,
            DeltaK, nmean, volume)


################################
################################

def integrand_B_Kernel_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, double f,
        double Sigma2, double dSigma2, char * parameters
        ):
    return integrand_B_Kernel(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, f, Sigma2, dSigma2,
            &parameters[0])

################################

def integrand_B_Kernel_SN_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, double f,
        double Sigma2, double dSigma2, char * parameters
        ):
    return integrand_B_Kernel_SN(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, kmag1, f, Sigma2, dSigma2,
            &parameters[0])

################################

def integrand_B_Kernel_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double f,
        double Sigma2, double dSigma2, double alpha_perp, double alpha_parallel, char * parameters
        ):
    return integrand_B_Kernel_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, f, Sigma2, dSigma2, alpha_perp, alpha_parallel,
            &parameters[0])

################################

def integrand_B_Kernel_SN_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double f,
        double Sigma2, double dSigma2, double alpha_perp, double alpha_parallel, char * parameters
        ):
    return integrand_B_Kernel_SN_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, f, Sigma2, dSigma2, alpha_perp, alpha_parallel,
            &parameters[0])

def integrand_B_Kernel_PNG_diag_py(
        np.ndarray[double, ndim=1, mode="c"] xx_in not None, int ndim, np.ndarray[double, ndim=1, mode="c"] ff_out not None, int ncomp,
        np.ndarray[double, ndim=1, mode="c"] kbin not None, int num_k_bin, int ell1, int ell2, int ELL, double f,
        double Sigma2, double dSigma2, char * parameters
        ):
    return integrand_B_Kernel_PNG_diag(
            &xx_in[0], ndim, &ff_out[0], ncomp,
            &kbin[0], num_k_bin, ell1, ell2, ELL, f, Sigma2, dSigma2,
            &parameters[0])

def B_Tree_FoG_Damping_diag_k_vector_py(
        np.ndarray[double, ndim=1, mode="c"] kvec1_in not None, np.ndarray[double, ndim=1, mode="c"] kvec2_in not None, np.ndarray[double, ndim=1, mode="c"] kvec3_in not None, 
        np.ndarray[double, ndim=1, mode="c"] los not None, 
        double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2, double c1, double c2, double knl,
        double Sigma2, double dSigma2
        ):
    return Bispectrum_Tree_FoG_Damping(
            &kvec1_in[0], &kvec2_in[0], &kvec3_in[0], &los[0], alpha_perp, alpha_parallel, f, b1, b2, bG2,
            c1, c2, knl, Sigma2, dSigma2)

def B_SN_FoG_Damping_diag_k_vector_py(
        np.ndarray[double, ndim=1, mode="c"] kvec1_in not None, np.ndarray[double, ndim=1, mode="c"] kvec2_in not None, np.ndarray[double, ndim=1, mode="c"] kvec3_in not None, 
        np.ndarray[double, ndim=1, mode="c"] los not None, 
        double alpha_perp, double alpha_parallel, double f, double b1, double c1, double c2, double knl, double Pshot, double Bshot,
        double Sigma2, double dSigma2
        ):
    return Bispectrum_SN_FoG_Damping(
            &kvec1_in[0], &kvec2_in[0], &kvec3_in[0], &los[0], alpha_perp, alpha_parallel, f, b1,
            c1, c2, knl, Pshot, Bshot, Sigma2, dSigma2)

def P_Tree_Damping_k_vector_py(
        np.ndarray[double, ndim=1, mode="c"] kvec1_in not None,
        np.ndarray[double, ndim=1, mode="c"] los not None, 
        double alpha_perp, double alpha_parallel, double f, double b1,
        double Sigma2, double dSigma2
        ):
    return Powerspectrum_Tree_Damping_for_1loop(
            &kvec1_in[0],&los[0], alpha_perp, alpha_parallel, f, b1, Sigma2, dSigma2)

def P_Counterterm_Damping_k_vector_py(
        np.ndarray[double, ndim=1, mode="c"] kvec1_in not None,
        np.ndarray[double, ndim=1, mode="c"] los not None, 
        double alpha_perp, double alpha_parallel, double f, double b1,
        double c0, double c1, double c2, double ch, double knl, 
        double Sigma2, double dSigma2
        ):
    return Powerspectrum_Counterterm_Damping(
            &kvec1_in[0],&los[0], alpha_perp, alpha_parallel, f, b1, c0, c1, c2, ch, knl, Sigma2, dSigma2)

def calcYYY_py(
        np.ndarray[double, ndim=1, mode="c"] kvec1 not None, np.ndarray[double, ndim=1, mode="c"] kvec2 not None, np.ndarray[double, ndim=1, mode="c"] los not None, 
        int ell1, int ell2, int ELL
        ):
    return calcYYY(
            &kvec1[0], &kvec2[0], &los[0], ell1, ell2, ELL)



