#ifndef __calc_cov_PP__
#define __calc_cov_PP__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __pk_lin__
#include "pk_lin.hpp"
#endif

#ifndef __kernel__
#include "kernel.hpp"
#endif

#ifndef __sph__
#include "sph.hpp"
#endif

/* covariance of P */
int integrand_cov_PP_G(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double alpha_perp, double alpha_parallel,
                       double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu  = -1.0 + 2.0 * xx_in[0];
    double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
    phi = phi;
    
    for(int i = 0; i < num_kbin; i++) {
    
	double kvec[3] = {0.0, 0.0, kbin[i]};
	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
	    
	double Delta = calcDeltaFunction(kbin[i], kbin[i], DeltaK);
	
	double L1 = LegendrePolynomials(ELL, mu);
	double L2 = LegendrePolynomials(ELL_dash, mu);
	double Nlll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0) * 2.0;
	
	double result = Nlll * L1 * L2 * Delta * Cov_PP_G(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	ff_out[i] = result * jacobian;
    
    }
    return 0;

}

int integrand_cov_PP_G_NL(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double alpha_perp, double alpha_parallel,
                          double sigma8, double fz, double b1, double sigma2_perp, double sigma2_para, double DeltaK, double nmean, double volume) {

    double mu  = -1.0 + 2.0 * xx_in[0];
    double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
    phi = phi;
    
    for(int i = 0; i < num_kbin; i++) {
    
	double kvec[3] = {0.0, 0.0, kbin[i]};
	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
	    
	double Delta = calcDeltaFunction(kbin[i], kbin[i], DeltaK);
	
	double L1 = LegendrePolynomials(ELL, mu);
	double L2 = LegendrePolynomials(ELL_dash, mu);
	double Nlll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0) * 2.0;
	
	double result = Nlll * L1 * L2 * Delta * Cov_PP_G_NL(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para, nmean, volume);
	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	ff_out[i] = result * jacobian;
    
    }
    return 0;
}


/* covariance of P */
int integrand_cov_PP_NG(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                        double sigma8, double fz, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	    double mu1_los = MU(kvec1, los);
	    double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	    double result = Nll * L1 * L2 * Cov_PP_NG(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nmean, volume);
    
    	double jacobian = 1.0;
    	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_b2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_b2(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_bK2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_bK2(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_b2_b2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_b2_b2(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_b2_bK2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_b2_bK2(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_bK2_bK2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_bK2_bK2(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_b3(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_b3(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_bK3(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_bK3(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_bDK(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_bDK(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_bO(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int j = 0; j < num_kbin; j++) {
    
	/********/
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	/********/
    
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_bO(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}



int integrand_cov_PP_NG_BeatCoupling(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_BeatCoupling_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                        double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	    double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
    	double mu1_los = MU(kvec1, los);
    	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
    	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, nmean, volume);
    
    	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
    	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_BeatCoupling_b2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_b2(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_BeatCoupling_bK2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_bK2(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_BeatCoupling_b2_b2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_b2_b2(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_BeatCoupling_b2_bK2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_b2_bK2(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}



int integrand_cov_PP_NG_BeatCoupling_bK2_bK2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_bK2_bK2(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}



int integrand_cov_PP_NG_BeatCoupling_b3(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_b3(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_BeatCoupling_bK3(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_bK3(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_BeatCoupling_bDK(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_bDK(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_BeatCoupling_bO(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                     double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};

	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
   
	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_BeatCoupling_bO(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_LocalMean(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                        double sigma8, double fz, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
	double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};

    	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_LocalMean(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_LocalMean_NL(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double sigma2_perp, double sigma2_para,
        double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
        double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
        double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
        
        double mu1_los = MU(kvec1, los);
        double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
        double L2 = LegendrePolynomials(ELL_dash, mu2_los);
        double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
        
        double result = Nll * L1 * L2 * Cov_PP_NG_LocalMean_NL(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, sigma2_perp, sigma2_para, nmean, volume);
        
        double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
        ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_LocalMean_NL_Sigma2B(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1,
        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double sigma2_perp, double sigma2_para,
        double DeltaK, double nmean, double volume, double sigma2_b) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
        double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
        double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
        double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
        
        double mu1_los = MU(kvec1, los);
        double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
        double L2 = LegendrePolynomials(ELL_dash, mu2_los);
        double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
        
        double result = Nll * L1 * L2 * Cov_PP_NG_LocalMean_NL_Sigma2B(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, sigma2_perp, sigma2_para, nmean, volume, sigma2_b);
        
        double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
        ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_LocalMean_NL_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                        double sigma8, double fz, double b1, double sigma2_perp, double sigma2_para, double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
	    double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
    	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
    	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};

    	double mu1_los = MU(kvec1, los);
    	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
    	double result = Nll * L1 * L2 * Cov_PP_NG_LocalMean_NL(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sigma2_perp, sigma2_para, nmean, volume);
    
    	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
    	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_NG_LocalMean_NL_b2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                        double sigma8, double fz, double b1, double sigma2_perp, double sigma2_para,
				        double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
	double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};

    	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_LocalMean_NL_b2(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}

int integrand_cov_PP_NG_LocalMean_NL_bK2(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double kmag1, double alpha_perp, double alpha_parallel,
                                        double sigma8, double fz, double b1, double sigma2_perp, double sigma2_para,
				        double DeltaK, double nmean, double volume) {

    double mu2 = - 1.0 + 2.0 * xx_in[0];
    double mu  = - 1.0 + 2.0 * xx_in[1];
    double phi = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int j = 0; j < num_kbin; j++) {
    
	double kvec1[3] = {0.0, 0.0, kmag1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};

    	double mu1_los = MU(kvec1, los);
	double mu2_los = MU(kvec2, los);
        double L1 = LegendrePolynomials(ELL, mu1_los);
      	double L2 = LegendrePolynomials(ELL_dash, mu2_los);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_LocalMean_NL_bK2(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[j] = result * jacobian;
    
    }
    
    return 0;
}


int integrand_cov_PP_SN(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double alpha_perp, double alpha_parallel,
                        double sigma8, double fz, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume) {

    double mu1  = - 1.0 + 2.0 * xx_in[0];
    double mu2  = - 1.0 + 2.0 * xx_in[1];
    double phi2 = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
    
        int index = i * num_kbin + j;
        double kvec1[3] = {kbin[i] * sqrt(1.0 - mu1 * mu1), 0.0, kbin[i] * mu1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2) * cos(phi2), kbin[j] * sqrt(1.0 - mu2 * mu2) * sin(phi2), kbin[j] * mu2};
        double los[3] = {0.0, 0.0, 1.0};
	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
    
        double L1 = LegendrePolynomials(ELL, mu1);
      	double L2 = LegendrePolynomials(ELL_dash, mu2);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_SN(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[index] = result * jacobian;
    
    }}
    
    return 0;
}

int integrand_cov_PP_SN_LocalMean(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double alpha_perp, double alpha_parallel,
                        double sigma8, double fz, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume) {

    double mu1  = - 1.0 + 2.0 * xx_in[0];
    double mu2  = - 1.0 + 2.0 * xx_in[1];
    double phi2 = 0.0 + 2.0 * M_PI * xx_in[2];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[3] );
    double mu_e = -1.0 + 2.0 * xx_in[4];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[5];   

    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
    
        int index = i * num_kbin + j;
        double kvec1[3] = {kbin[i] * sqrt(1.0 - mu1 * mu1), 0.0, kbin[i] * mu1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2) * cos(phi2), kbin[j] * sqrt(1.0 - mu2 * mu2) * sin(phi2), kbin[j] * mu2};
        double los[3] = {0.0, 0.0, 1.0};
	double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
    
        double L1 = LegendrePolynomials(ELL, mu1);
      	double L2 = LegendrePolynomials(ELL_dash, mu2);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_SN_LocalMean(kvec1, kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
    
	double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
	ff_out[index] = result * jacobian;
    
    }}
    
    return 0;
}

/* covariance of P */
int integrand_cov_PP_NG_Reconstructed(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double alpha_perp, double alpha_parallel,
                                      double sigma8, double fz, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume,
				      double b1_fid, double R) {

    double mu1  = - 1.0 + 2.0 * xx_in[0];
    double mu2  = - 1.0 + 2.0 * xx_in[1];
    double phi2 = 0.0 + 2.0 * M_PI * xx_in[2];
    
    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
    
        int index = i * num_kbin + j;
        double kvec1[3] = {kbin[i] * sqrt(1.0 - mu1 * mu1), 0.0, kbin[i] * mu1};
        double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2) * cos(phi2), kbin[j] * sqrt(1.0 - mu2 * mu2) * sin(phi2), kbin[j] * mu2};
        double los[3] = {0.0, 0.0, 1.0};
    
        double L1 = LegendrePolynomials(ELL, mu1);
      	double L2 = LegendrePolynomials(ELL_dash, mu2);
      	double Nll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0);
    
	double result = Nll * L1 * L2 * Cov_PP_NG_Reconstructed(kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume, b1_fid, R);
    
	double jacobian = 1.0;
	ff_out[index] = result * jacobian;
    
    }}
    
    return 0;
}

/* covariance of P */
int integrand_cov_PP_KSZ_G(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double alpha_perp, double alpha_parallel,
                           double sigma8, double fz, double b1, double DeltaK, double nmean, double volume, double aH_tau_T0_over_c1, double aH_tau_T0_over_c2, double sigma2_vv, double R2_N) {

    double mu  = -1.0 + 2.0 * xx_in[0];
    double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
    phi = phi;
    
    for(int i = 0; i < num_kbin; i++) {
    
	double kvec[3] = {kbin[i] * sqrt(1.0 - mu * mu), 0.0, kbin[i] * mu};
	double los[3] = {0.0, 0.0, 1.0};
	    
	double Delta = calcDeltaFunction(kbin[i], kbin[i], DeltaK);
	
	double L1 = LegendrePolynomials(ELL, mu);
	double L2 = LegendrePolynomials(ELL_dash, mu);
	double Nlll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0) * 2.0;
	
	double result = Nlll * L1 * L2 * Delta * Cov_PP_KSZ_G(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume,
              	                                              aH_tau_T0_over_c1, aH_tau_T0_over_c2, sigma2_vv, R2_N);
	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	ff_out[i] = result * jacobian;
    
    }
    return 0;

}

/* covariance of P */
int integrand_cov_PP_GALAXY_KSZ_G(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ELL_dash, double alpha_perp, double alpha_parallel,
                                  double sigma8, double fz, double b1, double DeltaK, double nmean, double volume, double aH_tau_T0_over_c) {

    double mu  = -1.0 + 2.0 * xx_in[0];
    double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
    phi = phi;
    
    for(int i = 0; i < num_kbin; i++) {
    
	double kvec[3] = {kbin[i] * sqrt(1.0 - mu * mu), 0.0, kbin[i] * mu};
	double los[3] = {0.0, 0.0, 1.0};
	    
	double Delta = calcDeltaFunction(kbin[i], kbin[i], DeltaK);
	
	double L1 = LegendrePolynomials(ELL, mu);
	double L2 = LegendrePolynomials(ELL_dash, mu);
	double Nlll = (2.0*double(ELL)+1.0) * (2.0*double(ELL_dash)+1.0) * 2.0;
	
	double result = Nlll * L1 * L2 * Delta * Cov_PP_GALAXY_KSZ_G(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume, aH_tau_T0_over_c);
	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	ff_out[i] = result * jacobian;
    
    }
    return 0;

}





#endif
