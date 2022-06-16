#ifndef __calc_cov_BB__
#define __calc_cov_BB__

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

int integrand_cov_BB_G(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash,
			           double k1, double k1_dash,
                 	   double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
//	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[1];
	double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
//	double mu = 1.0;
//	double phi = 0.0;

    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {

		double kmag1 = k1;
		double kmag2 = kbin[i];
		double kmag1_dash = k1_dash;
		double kmag2_dash = kbin[j];

	
		int index = i * num_kbin + j;
 
		/********/
		double kvec1[3] = {k1 * sqrt(1.0 - mu1 * mu1), 0.0, k1 * mu1};
		double kvec2[3] = {kbin[i] * sqrt(1.0 - mu2 * mu2) * cos(phi2), kbin[i] * sqrt(1.0 - mu2 * mu2) * sin(phi2), kbin[i] * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {0.0, 0.0, 1.0};
		/********/
	

		double kmag3 = NORM(kvec3);
	
		double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
		double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
	
		double SllL = calcYYY(kvec1, kvec2, los, ell1, ell2, ELL);

		double cov = Cov_BB_G_PPP(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
	
		double fac = calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			       * calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash)
			   
			       + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			       * calcYYY(kvec2, kvec1, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			       * calcYYY(kvec1, kvec3, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			       * calcYYY(kvec3, kvec1, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			       * calcYYY(kvec2, kvec3, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			       * calcYYY(kvec3, kvec2, los, ell1_dash, ell2_dash, ELL_dash);
	
		double result = Nlll * Nlll_dash * SllL * fac * cov;
    
        double jacobian = 1.0;
	
		ff_out[index] = result * jacobian;

	}}
	
	return 0;
}

int integrand_cov_BB_G_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, double k1,
                 	   double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
//	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[1];
	double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
//	double mu = 1.0;
//	double phi = 0.0;

    for(int j = 0; j < num_kbin; j++) {

		double kmag1 = k1;
		double kmag2 = k1;
		double kmag1_dash = kbin[j];
		double kmag2_dash = kbin[j];

	
		/********/
		double kvec1[3] = {k1 * sqrt(1.0 - mu1 * mu1), 0.0, k1 * mu1};
		double kvec2[3] = {k1 * sqrt(1.0 - mu2 * mu2) * cos(phi2), k1 * sqrt(1.0 - mu2 * mu2) * sin(phi2), k1 * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {0.0, 0.0, 1.0};
		/********/
	

		double kmag3 = NORM(kvec3);
	
		double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
		double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
	
		double SllL = calcYYY(kvec1, kvec2, los, ell1, ell2, ELL);

		double cov = Cov_BB_G_PPP(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
	
		double fac = calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			       * calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash)
			   
			       + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			       * calcYYY(kvec2, kvec1, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			       * calcYYY(kvec1, kvec3, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			       * calcYYY(kvec3, kvec1, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			       * calcYYY(kvec2, kvec3, los, ell1_dash, ell2_dash, ELL_dash)
	
			       + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			       * calcYYY(kvec3, kvec2, los, ell1_dash, ell2_dash, ELL_dash);
	
		double result = Nlll * Nlll_dash * SllL * fac * cov;
    
        double jacobian = 1.0;
	
		ff_out[j] = result * jacobian;

	}
	
	return 0;
}

int integrand_cov_BB_NG_BB_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, double DeltaK, double nmean, double volume) {

	double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
	double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
	
	for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {
		
		int index = i * num_kbin + j;
 
		double kmag1 = kbin[i];
		double kmag2 = kbin[i];
		double kmag1_dash = kbin[j];
		double kmag2_dash = kbin[j];
	
		double result = 0.0;
		/**********/
		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu2  = - 1.0 + 2.0 * xx_in[1];
			double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
			double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
			double kvec3_dash_A[3] = { - kvec1[0] - kvec2_dash[0], - kvec1[1] - kvec2_dash[1], - kvec1[2] - kvec2_dash[2] };
			double kvec3_dash_B[3] = { - kvec2[0] - kvec2_dash[0], - kvec2[1] - kvec2_dash[1], - kvec2[2] - kvec2_dash[2] };
			double kvec3_dash_C[3] = { - kvec3[0] - kvec2_dash[0], - kvec3[1] - kvec2_dash[1], - kvec3[2] - kvec2_dash[2] };

			double kmag3 = NORM(kvec3);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec1, kvec2, kvec3, kvec1, kvec2_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec2, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec2, kvec1, kvec3, kvec2, kvec2_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec3, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3, kvec1, kvec2, kvec3, kvec2_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu2  = - 1.0 + 2.0 * xx_in[1];
			double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu1_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
			double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
		    
			double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
			double kvec3_dash_A[3] = { - kvec1_dash[0] - kvec1[0], - kvec1_dash[1] - kvec1[1], - kvec1_dash[2] - kvec1[2] };
			double kvec3_dash_B[3] = { - kvec1_dash[0] - kvec2[0], - kvec1_dash[1] - kvec2[1], - kvec1_dash[2] - kvec2[2] };
			double kvec3_dash_C[3] = { - kvec1_dash[0] - kvec3[0], - kvec1_dash[1] - kvec3[1], - kvec1_dash[2] - kvec3[2] };

			double kmag3 = NORM(kvec3);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec1, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
				* Cov_BB_NG_BB(kvec1, kvec2, kvec3, kvec1, kvec1_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec2, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
				* Cov_BB_NG_BB(kvec2, kvec1, kvec3, kvec2, kvec1_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec3, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3, kvec1, kvec2, kvec3, kvec1_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu2  = - 1.0 + 2.0 * xx_in[0];
			double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
			double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2), 0.0, kmag2 * mu2 };
			double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
			double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
			double kvec3_A[3] = { - kvec3_dash[0] - kvec2[0], - kvec3_dash[1] - kvec2[1], - kvec3_dash[2] - kvec2[2] };

			double kmag3_dash = NORM(kvec3_dash);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec3_dash, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag1, kmag3_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3_dash, kvec2, kvec3_A, kvec3_dash, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
			double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
			double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
			double kvec3_B[3] = { - kvec1[0] - kvec3_dash[0], - kvec1[1] - kvec3_dash[1], - kvec1[2] - kvec3_dash[2] };

			double kmag3_dash = NORM(kvec3_dash);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec3_dash, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag2, kmag3_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3_dash, kvec1, kvec3_B, kvec3_dash, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu2  = - 1.0 + 2.0 * xx_in[1];
			double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };

			double kvec_alpha[3] = { kvec1[0] + kvec2[0] - kvec2_dash[0], kvec1[1] + kvec2[1] - kvec2_dash[1], kvec1[2] + kvec2[2] - kvec2_dash[2]};
			double kmag_alpha = NORM(kvec_alpha);
		    
			double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };

			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec_alpha, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag_alpha, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3, kvec1, kvec2, kvec3, kvec_alpha, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

	/**********/

	    double jacobian = 1.0;

	    ff_out[index] = result * jacobian;

	}}

	return 0;
}

int integrand_cov_BB_NG_PT_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, 
				double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume) {

    double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
    double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
    
    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
		
	int index = i * num_kbin + j;
	
	double kmag1 = kbin[i];
	double kmag2 = kbin[i];
	double kmag1_dash = kbin[j];
	double kmag2_dash = kbin[j];
	
	double result = 0.0;
	/**********/
	if(1) {
	    
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu2  = - 1.0 + 2.0 * xx_in[1];
	    double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];
	    
	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
	    
	    double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	    double kvec3_dash_A[3] = { - kvec2[0] - kvec3[0] - kvec2_dash[0],  - kvec2[1] - kvec3[1] - kvec2_dash[1],  - kvec2[2] - kvec3[2] - kvec2_dash[2]};
	    double kvec3_dash_B[3] = { - kvec1[0] - kvec3[0] - kvec2_dash[0],  - kvec1[1] - kvec3[1] - kvec2_dash[1],  - kvec1[2] - kvec3[2] - kvec2_dash[2]};
	    double kvec3_dash_C[3] = { - kvec1[0] - kvec2[0] - kvec2_dash[0],  - kvec1[1] - kvec2[1] - kvec2_dash[1],  - kvec1[2] - kvec2[2] - kvec2_dash[2]};
	    
	    double kmag3 = NORM(kvec3);
	    
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec1[3] = MINUS_XVEC(kvec1);
	    double M_kvec2[3] = MINUS_XVEC(kvec2);
	    double M_kvec3[3] = MINUS_XVEC(kvec3);
	    
	    result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec1, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec1, kvec2, kvec3, kvec2_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	    
	            + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec2, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec2, kvec1, kvec3, kvec2_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	    
	    	    + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec3, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3, kvec1, kvec2, kvec2_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	    
	}

	if(1) {
	
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu2  = - 1.0 + 2.0 * xx_in[1];
	    double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu1_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[4];
	    
	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
	    double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
	    
	    double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	    double kvec3_dash_A[3] = { - kvec2[0] - kvec3[0] - kvec1_dash[0],  - kvec2[1] - kvec3[1] - kvec1_dash[1],  - kvec2[2] - kvec3[2] - kvec1_dash[2]};
	    double kvec3_dash_B[3] = { - kvec1[0] - kvec3[0] - kvec1_dash[0],  - kvec1[1] - kvec3[1] - kvec1_dash[1],  - kvec1[2] - kvec3[2] - kvec1_dash[2]};
	    double kvec3_dash_C[3] = { - kvec1[0] - kvec2[0] - kvec1_dash[0],  - kvec1[1] - kvec2[1] - kvec1_dash[1],  - kvec1[2] - kvec2[2] - kvec1_dash[2]};
	    
	    double kmag3 = NORM(kvec3);
	    
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec1[3] = MINUS_XVEC(kvec1);
	    double M_kvec2[3] = MINUS_XVEC(kvec2);
	    double M_kvec3[3] = MINUS_XVEC(kvec3);
	    
	    result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, M_kvec1, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec1, kvec2, kvec3, kvec1_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	
		    + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, M_kvec2, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec2, kvec1, kvec3, kvec1_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	
		    + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, M_kvec3, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3, kvec1, kvec2, kvec1_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}

	if(1) {
	
	    double mu2  = - 1.0 + 2.0 * xx_in[0];
	    double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
	    double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2), 0.0, kmag2 * mu2 };
	    double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
	    double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
	    double kvec3_A[3] = { - kvec2[0] - kvec1_dash[0] - kvec2_dash[0], - kvec2[1] - kvec1_dash[1] - kvec2_dash[1], - kvec2[2] - kvec1_dash[2] - kvec2_dash[2] };

	    double kmag3_dash = NORM(kvec3_dash);
	
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec3_dash[3] = MINUS_XVEC(kvec3_dash);
	    
	    result += Nlll * Nlll_dash * calcYYY(M_kvec3_dash, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag1, kmag3_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3_dash, kvec2, kvec3_A, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}


	if(1) {
	
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
	    double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
	    double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
	    double kvec3_B[3] = { - kvec1[0] - kvec1_dash[0] - kvec2_dash[0], - kvec1[1] - kvec1_dash[1] - kvec2_dash[1], - kvec1[2] - kvec1_dash[2] - kvec2_dash[2] };

	    double kmag3_dash = NORM(kvec3_dash);
	
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec3_dash[3] = MINUS_XVEC(kvec3_dash);
	    
	    result += Nlll * Nlll_dash * calcYYY(kvec1, M_kvec3_dash, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag2, kmag3_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3_dash, kvec1, kvec3_B, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}


	if(1) {
	
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu2  = - 1.0 + 2.0 * xx_in[1];
	    double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];
	    
	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
	    
	    double kvec_beta[3] = { kvec1[0] + kvec2[0] + kvec2_dash[0], kvec1[1] + kvec2[1] + kvec2_dash[1], kvec1[2] + kvec2[2] + kvec2_dash[2]};
	    double kmag_beta = NORM(kvec_beta);
	    
	    double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	    
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec_beta[3] = MINUS_XVEC(kvec_beta);

	    result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec_beta, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag_beta, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3, kvec1, kvec2, M_kvec_beta, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}

	/**********/

	double jacobian = 1.0;
	
	ff_out[index] = result * jacobian;

    }}

    return 0;
}

int integrand_cov_BB_NG_P6_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu1  = - 1.0 + 2.0 * xx_in[0];
    double mu2  = - 1.0 + 2.0 * xx_in[1];
    double phi2 = 0.0 + 2.0 * M_PI * xx_in[2];
    double mu1_dash  = - 1.0 + 2.0 * xx_in[3];
    double phi1_dash = 0.0 + 2.0 * M_PI * xx_in[4];
    double mu2_dash  = - 1.0 + 2.0 * xx_in[5];
    double phi2_dash = 0.0 + 2.0 * M_PI * xx_in[6];
     
    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
    
	double kmag1 = kbin[i];
	double kmag2 = kbin[i];
	double kmag1_dash = kbin[j];
	double kmag2_dash = kbin[j];
	
        int index = i * num_kbin + j;
        double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1};
        double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};
        double kvec1_dash[3] = {kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash};
        double kvec2_dash[3] = {kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash};

	double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };

        double los[3] = {0.0, 0.0, 1.0};

	double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
	double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
	
	double result = Nlll * Nlll_dash 
	              * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	              * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
                      * Cov_BB_NG_P6(kvec1, kvec2, kvec3, kvec1_dash, kvec2_dash, kvec3_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[index] = result * jacobian;
    
    }}
    
    return 0;
}

int integrand_cov_BB_TripoSH_G_full_real(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ell12, int ELL, int _J_, int _MJ_, int ell1_dash, int ell2_dash, int ell12_dash, int ELL_dash, int _J_dash_, int _MJ_dash_,
			    double k1, double k1_dash,
                 	    double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
	double phi1 = 0.0 + 2.0 * M_PI * xx_in[1];
	double mu2  = - 1.0 + 2.0 * xx_in[2];
	double phi2 = 0.0 + 2.0 * M_PI * xx_in[3];
	double mu  =  - 1.0 + 2.0 * xx_in[4]; 
	double phi =  0.0 + 2.0 * M_PI * xx_in[5];

    	for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {
		
		int index = i * num_kbin + j;
 
		/********/
		double kvec1[3] = {k1 * sqrt(1.0 - mu1 * mu1) * cos(phi1), k1 * sqrt(1.0-mu1*mu1) * sin(phi1), k1 * mu1};
		double kvec2[3] = {kbin[i] * sqrt(1.0 - mu2 * mu2) * cos(phi2), kbin[i] * sqrt(1.0 - mu2 * mu2) * sin(phi2), kbin[i] * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0-mu*mu)*sin(phi), mu};
		/********/
	
		double kmag1 = k1;
		double kmag2 = kbin[i];
		double kmag3 = NORM(kvec3);
		double kmag1_dash = k1_dash;
		double kmag2_dash = kbin[j];
	
		double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
		double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
		std::complex<double> SllL = calcTripoSH(kvec1, kvec2, los, ell1, ell2, ell12, ELL, _J_, _MJ_);
		double cov = Cov_BB_G_PPP(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
	
			//result += Nlll * Nlll_dash * SllL * SllL_dash * Delta1 * Delta2 * cov;
		std::complex<double> fac = calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec1, kvec2, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
			   
			   + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec2, kvec1, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec1, kvec3, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec3, kvec1, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec2, kvec3, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec3, kvec2, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_);
	
		std::complex<double> result = Nlll * Nlll_dash * SllL * fac * cov;
	
		double jacobian = 1.0;
	
		ff_out[index] = result.real() * jacobian;

	}}
	
	return 0;
}


int integrand_cov_BB_TripoSH_G_full_imag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ell12, int ELL, int _J_, int _MJ_, int ell1_dash, int ell2_dash, int ell12_dash, int ELL_dash, int _J_dash_, int _MJ_dash_,
			    double k1, double k1_dash,
                 	    double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
	double phi1 = 0.0 + 2.0 * M_PI * xx_in[1];
	double mu2  = - 1.0 + 2.0 * xx_in[2];
	double phi2 = 0.0 + 2.0 * M_PI * xx_in[3];
	double mu  =  - 1.0 + 2.0 * xx_in[4]; 
	double phi =  0.0 + 2.0 * M_PI * xx_in[5];

    	for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {
		
		int index = i * num_kbin + j;
 
		/********/
		double kvec1[3] = {k1 * sqrt(1.0 - mu1 * mu1) * cos(phi1), k1 * sqrt(1.0-mu1*mu1) * sin(phi1), k1 * mu1};
		double kvec2[3] = {kbin[i] * sqrt(1.0 - mu2 * mu2) * cos(phi2), kbin[i] * sqrt(1.0 - mu2 * mu2) * sin(phi2), kbin[i] * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0-mu*mu)*sin(phi), mu};
		/********/
	
		double kmag1 = k1;
		double kmag2 = kbin[i];
		double kmag3 = NORM(kvec3);
		double kmag1_dash = k1_dash;
		double kmag2_dash = kbin[j];
	
		double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
		double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
		std::complex<double> SllL = calcTripoSH(kvec1, kvec2, los, ell1, ell2, ell12, ELL, _J_, _MJ_);
		double cov = Cov_BB_G_PPP(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
	
			//result += Nlll * Nlll_dash * SllL * SllL_dash * Delta1 * Delta2 * cov;
		std::complex<double> fac = calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec1, kvec2, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
			   
			   + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec2, kvec1, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec1, kvec3, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec3, kvec1, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec2, kvec3, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_)
	
			   + calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
		           * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
			   * calcTripoSH(kvec3, kvec2, los, ell1_dash, ell2_dash, ell12_dash, ELL_dash, _J_dash_, _MJ_dash_);
	
		std::complex<double> result = Nlll * Nlll_dash * SllL * fac * cov;
	
		double jacobian = 1.0;
	
		ff_out[index] = result.imag() * jacobian;

	}}
	
	return 0;
}

int integrand_cov_BB_NG_BB_full(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
			        double k1, double k1_dash,
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, double DeltaK, double nmean, double volume) {

	double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
	double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
	
    	for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {

		double kmag1 = k1;
		double kmag2 = kbin[i];
		double kmag1_dash = k1_dash;
		double kmag2_dash = kbin[j];
	
		int index = i * num_kbin + j;

		double result = 0.0;
		/**********/
		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu2  = - 1.0 + 2.0 * xx_in[1];
			double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
			double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
			double kvec3_dash_A[3] = { - kvec1[0] - kvec2_dash[0], - kvec1[1] - kvec2_dash[1], - kvec1[2] - kvec2_dash[2] };
			double kvec3_dash_B[3] = { - kvec2[0] - kvec2_dash[0], - kvec2[1] - kvec2_dash[1], - kvec2[2] - kvec2_dash[2] };
			double kvec3_dash_C[3] = { - kvec3[0] - kvec2_dash[0], - kvec3[1] - kvec2_dash[1], - kvec3[2] - kvec2_dash[2] };

			double kmag3 = NORM(kvec3);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec1, kvec2, kvec3, kvec1, kvec2_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec2, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec2, kvec1, kvec3, kvec2, kvec2_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec3, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3, kvec1, kvec2, kvec3, kvec2_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu2  = - 1.0 + 2.0 * xx_in[1];
			double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu1_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
			double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
		    
			double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
			double kvec3_dash_A[3] = { - kvec1_dash[0] - kvec1[0], - kvec1_dash[1] - kvec1[1], - kvec1_dash[2] - kvec1[2] };
			double kvec3_dash_B[3] = { - kvec1_dash[0] - kvec2[0], - kvec1_dash[1] - kvec2[1], - kvec1_dash[2] - kvec2[2] };
			double kvec3_dash_C[3] = { - kvec1_dash[0] - kvec3[0], - kvec1_dash[1] - kvec3[1], - kvec1_dash[2] - kvec3[2] };

			double kmag3 = NORM(kvec3);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec1, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
				* Cov_BB_NG_BB(kvec1, kvec2, kvec3, kvec1, kvec1_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec2, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
				* Cov_BB_NG_BB(kvec2, kvec1, kvec3, kvec2, kvec1_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

				+ Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec3, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3, kvec1, kvec2, kvec3, kvec1_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu2  = - 1.0 + 2.0 * xx_in[0];
			double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
			double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2), 0.0, kmag2 * mu2 };
			double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
			double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
			double kvec3_A[3] = { - kvec3_dash[0] - kvec2[0], - kvec3_dash[1] - kvec2[1], - kvec3_dash[2] - kvec2[2] };

			double kmag3_dash = NORM(kvec3_dash);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec3_dash, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag1, kmag3_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3_dash, kvec2, kvec3_A, kvec3_dash, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
			double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
			double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
			double kvec3_B[3] = { - kvec1[0] - kvec3_dash[0], - kvec1[1] - kvec3_dash[1], - kvec1[2] - kvec3_dash[2] };

			double kmag3_dash = NORM(kvec3_dash);
		    
			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec3_dash, los, ell1, ell2, ELL)
			        * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag2, kmag3_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3_dash, kvec1, kvec3_B, kvec3_dash, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

		if(1) {
	
			double mu1  = - 1.0 + 2.0 * xx_in[0];
			double mu2  = - 1.0 + 2.0 * xx_in[1];
			double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
			double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
			double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

			double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
			double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
			double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };

			double kvec_alpha[3] = { kvec1[0] + kvec2[0] - kvec2_dash[0], kvec1[1] + kvec2[1] - kvec2_dash[1], kvec1[2] + kvec2[2] - kvec2_dash[2]};
			double kmag_alpha = NORM(kvec_alpha);
		    
			double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };

			double los[3] = { 0.0, 0.0, 1.0 };
	
			result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
			        * calcYYY(kvec_alpha, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
				* calcDeltaFunction(kmag_alpha, kmag1_dash, DeltaK)
				* Cov_BB_NG_BB(kvec3, kvec1, kvec2, kvec3, kvec_alpha, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		}

	/**********/

	    double jacobian = 1.0;

	    ff_out[index] = result * jacobian;

	}}

	return 0;
}

int integrand_cov_BB_NG_PT_full(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
			        double k1, double k1_dash,
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, 
				double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double DeltaK, double nmean, double volume) {

    double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
    double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
    
    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
    
    	double kmag1 = k1;
    	double kmag2 = kbin[i];
    	double kmag1_dash = k1_dash;
    	double kmag2_dash = kbin[j];
	
	int index = i * num_kbin + j;

	double result = 0.0;
	/**********/
	if(1) {
	    
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu2  = - 1.0 + 2.0 * xx_in[1];
	    double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];
	    
	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
	    
	    double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	    double kvec3_dash_A[3] = { - kvec2[0] - kvec3[0] - kvec2_dash[0],  - kvec2[1] - kvec3[1] - kvec2_dash[1],  - kvec2[2] - kvec3[2] - kvec2_dash[2]};
	    double kvec3_dash_B[3] = { - kvec1[0] - kvec3[0] - kvec2_dash[0],  - kvec1[1] - kvec3[1] - kvec2_dash[1],  - kvec1[2] - kvec3[2] - kvec2_dash[2]};
	    double kvec3_dash_C[3] = { - kvec1[0] - kvec2[0] - kvec2_dash[0],  - kvec1[1] - kvec2[1] - kvec2_dash[1],  - kvec1[2] - kvec2[2] - kvec2_dash[2]};
	    
	    double kmag3 = NORM(kvec3);
	    
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec1[3] = MINUS_XVEC(kvec1);
	    double M_kvec2[3] = MINUS_XVEC(kvec2);
	    double M_kvec3[3] = MINUS_XVEC(kvec3);
	    
	    result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec1, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag1, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec1, kvec2, kvec3, kvec2_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	    
	            + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec2, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag2, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec2, kvec1, kvec3, kvec2_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	    
	    	    + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec3, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag3, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3, kvec1, kvec2, kvec2_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	    
	}

	if(1) {
	
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu2  = - 1.0 + 2.0 * xx_in[1];
	    double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu1_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[4];
	    
	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
	    double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
	    
	    double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	    double kvec3_dash_A[3] = { - kvec2[0] - kvec3[0] - kvec1_dash[0],  - kvec2[1] - kvec3[1] - kvec1_dash[1],  - kvec2[2] - kvec3[2] - kvec1_dash[2]};
	    double kvec3_dash_B[3] = { - kvec1[0] - kvec3[0] - kvec1_dash[0],  - kvec1[1] - kvec3[1] - kvec1_dash[1],  - kvec1[2] - kvec3[2] - kvec1_dash[2]};
	    double kvec3_dash_C[3] = { - kvec1[0] - kvec2[0] - kvec1_dash[0],  - kvec1[1] - kvec2[1] - kvec1_dash[1],  - kvec1[2] - kvec2[2] - kvec1_dash[2]};
	    
	    double kmag3 = NORM(kvec3);
	    
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec1[3] = MINUS_XVEC(kvec1);
	    double M_kvec2[3] = MINUS_XVEC(kvec2);
	    double M_kvec3[3] = MINUS_XVEC(kvec3);
	    
	    result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, M_kvec1, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag1, kmag2_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec1, kvec2, kvec3, kvec1_dash, kvec3_dash_A, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	
		    + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, M_kvec2, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag2, kmag2_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec2, kvec1, kvec3, kvec1_dash, kvec3_dash_B, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume)
	
		    + Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, M_kvec3, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag3, kmag2_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3, kvec1, kvec2, kvec1_dash, kvec3_dash_C, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}

	if(1) {
	
	    double mu2  = - 1.0 + 2.0 * xx_in[0];
	    double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
	    double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2), 0.0, kmag2 * mu2 };
	    double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
	    double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
	    double kvec3_A[3] = { - kvec2[0] - kvec1_dash[0] - kvec2_dash[0], - kvec2[1] - kvec1_dash[1] - kvec2_dash[1], - kvec2[2] - kvec1_dash[2] - kvec2_dash[2] };

	    double kmag3_dash = NORM(kvec3_dash);
	
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec3_dash[3] = MINUS_XVEC(kvec3_dash);
	    
	    result += Nlll * Nlll_dash * calcYYY(M_kvec3_dash, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag1, kmag3_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3_dash, kvec2, kvec3_A, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}


	if(1) {
	
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu1_dash  = - 1.0 + 2.0 * xx_in[1];
	    double phi1_dash =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];

	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec1_dash[3] = { kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
		    
	    double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };
	    double kvec3_B[3] = { - kvec1[0] - kvec1_dash[0] - kvec2_dash[0], - kvec1[1] - kvec1_dash[1] - kvec2_dash[1], - kvec1[2] - kvec1_dash[2] - kvec2_dash[2] };

	    double kmag3_dash = NORM(kvec3_dash);
	
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec3_dash[3] = MINUS_XVEC(kvec3_dash);
	    
	    result += Nlll * Nlll_dash * calcYYY(kvec1, M_kvec3_dash, los, ell1, ell2, ELL)
	            * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag2, kmag3_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3_dash, kvec1, kvec3_B, kvec1_dash, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}


	if(1) {
	
	    double mu1  = - 1.0 + 2.0 * xx_in[0];
	    double mu2  = - 1.0 + 2.0 * xx_in[1];
	    double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
	    double mu2_dash  = - 1.0 + 2.0 * xx_in[3];
	    double phi2_dash =   0.0 + 2.0 * M_PI * xx_in[4];
	    
	    double kvec1[3] = { kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1 };
	    double kvec2[3] = { kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2 };
	    double kvec2_dash[3] = { kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash };
	    
	    double kvec_beta[3] = { kvec1[0] + kvec2[0] + kvec2_dash[0], kvec1[1] + kvec2[1] + kvec2_dash[1], kvec1[2] + kvec2[2] + kvec2_dash[2]};
	    double kmag_beta = NORM(kvec_beta);
	    
	    double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	    
	    double los[3] = { 0.0, 0.0, 1.0 };
	    
	    double M_kvec_beta[3] = MINUS_XVEC(kvec_beta);

	    result += Nlll * Nlll_dash * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	            * calcYYY(M_kvec_beta, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
	    	    * calcDeltaFunction(kmag_beta, kmag1_dash, DeltaK)
	    	    * Cov_BB_NG_PT(kvec3, kvec1, kvec2, M_kvec_beta, kvec2_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, nmean, volume);
	
	}

	/**********/

	double jacobian = 1.0;
	
	ff_out[index] = result * jacobian;

    }}

    return 0;
}

int integrand_cov_BB_NG_P6_full(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
				double k1, double k1_dash,
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu1  = - 1.0 + 2.0 * xx_in[0];
    double mu2  = - 1.0 + 2.0 * xx_in[1];
    double phi2 = 0.0 + 2.0 * M_PI * xx_in[2];
    double mu1_dash  = - 1.0 + 2.0 * xx_in[3];
    double phi1_dash = 0.0 + 2.0 * M_PI * xx_in[4];
    double mu2_dash  = - 1.0 + 2.0 * xx_in[5];
    double phi2_dash = 0.0 + 2.0 * M_PI * xx_in[6];
     
    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
    
    	double kmag1 = k1;
    	double kmag2 = kbin[i];
    	double kmag1_dash = k1_dash;
    	double kmag2_dash = kbin[j];
	
        int index = i * num_kbin + j;
        double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1};
        double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};
        double kvec1_dash[3] = {kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * cos(phi1_dash), kmag1_dash * sqrt(1.0 - mu1_dash * mu1_dash) * sin(phi1_dash), kmag1_dash * mu1_dash};
        double kvec2_dash[3] = {kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * cos(phi2_dash), kmag2_dash * sqrt(1.0 - mu2_dash * mu2_dash) * sin(phi2_dash), kmag2_dash * mu2_dash};

	double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	double kvec3_dash[3] = { - kvec1_dash[0] - kvec2_dash[0], - kvec1_dash[1] - kvec2_dash[1], - kvec1_dash[2] - kvec2_dash[2] };

        double los[3] = {0.0, 0.0, 1.0};

	double Nlll      = (2.0*double(ELL)+1.0) * (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
	double Nlll_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
	
	double result = Nlll * Nlll_dash 
	              * calcYYY(kvec1, kvec2, los, ell1, ell2, ELL)
	              * calcYYY(kvec1_dash, kvec2_dash, los, ell1_dash, ell2_dash, ELL_dash)
                      * Cov_BB_NG_P6(kvec1, kvec2, kvec3, kvec1_dash, kvec2_dash, kvec3_dash, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[index] = result * jacobian;
    
    }}
    
    return 0;
}


#endif
