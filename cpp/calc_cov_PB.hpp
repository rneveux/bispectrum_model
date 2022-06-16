#ifndef __calc_cov_PB__
#define __calc_cov_PB__

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


int integrand_cov_PB_NG_PB_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash,
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, double DeltaK, double nmean, double volume) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
	double mu2  = - 1.0 + 2.0 * xx_in[1];
	double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];
    
        for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {
 
		double kmag = kbin[i];
		double kmag1 = kbin[j];
		double kmag2 = kbin[j];
            
            int index = i * num_kbin + j;
	
		double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1};
		double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {0.0, 0.0, 1.0};

		double kmag3 = NORM(kvec3);

		double NL = 2.0 * double(ELL) + 1.0;
		double NllL_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
		double SllL_dash = calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash);

		double result = calcDeltaFunction(kmag, kmag1, DeltaK)
			      * LegendrePolynomials(ELL, mu1) 
			      * Cov_PB_NG_PB(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

			      + calcDeltaFunction(kmag, kmag2, DeltaK)
			      * LegendrePolynomials(ELL, mu2) 
			      * Cov_PB_NG_PB(kvec2, kvec1, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

			      + calcDeltaFunction(kmag, kmag3, DeltaK)
			      * LegendrePolynomials(ELL, MU(kvec3, los)) 
			      * Cov_PB_NG_PB(kvec3, kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		result *= (NL * NllL_dash * SllL_dash);
	
		double jacobian = 1.0;
	
		ff_out[index] = result * jacobian;

	}}
	
	return 0;
}

int integrand_cov_PB_NG_P5_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash,
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu         = - 1.0 + 2.0 * xx_in[0];
    double mu1  = - 1.0 + 2.0 * xx_in[1];
    double phi1 = 0.0 + 2.0 * M_PI * xx_in[2];
    double mu2  = - 1.0 + 2.0 * xx_in[3];
    double phi2 = 0.0 + 2.0 * M_PI * xx_in[4];
     
    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
    
	double kmag = kbin[i];
	double kmag1 = kbin[j];
	double kmag2 = kbin[j];
	
        int index = i * num_kbin + j;

        double kvec[3] = {kmag * sqrt(1.0 - mu * mu), 0.0, kmag * mu};
        double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1) * cos(phi1), kmag1 * sqrt(1.0 - mu1 * mu1) * sin(phi1), kmag1 * mu1};
        double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};

	double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };

        double los[3] = {0.0, 0.0, 1.0};

	double NL = 2.0 * double(ELL) + 1.0;
	double NllL_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);

	double result = NL * NllL_dash 
	              * LegendrePolynomials(ELL, mu) 
	              * calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash)
                      * Cov_PB_NG_P5(kvec, kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[index] = result * jacobian;
    
    }}
    
    return 0;
}

int integrand_cov_PB_NG_PB_full(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
			        double k1_in, 
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, double DeltaK, double nmean, double volume) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
	double mu2  = - 1.0 + 2.0 * xx_in[1];
	double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];

    	for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {
		
		int index = i * num_kbin + j;
 
		double kmag = kbin[i];
		double kmag1 = k1_in;
		double kmag2 = kbin[j];

		double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1};
		double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {0.0, 0.0, 1.0};

		double kmag3 = NORM(kvec3);

		double NL = 2.0 * double(ELL) + 1.0;
		double NllL_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
		double SllL_dash = calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash);

		double result = calcDeltaFunction(kmag, kmag1, DeltaK)
			      * LegendrePolynomials(ELL, mu1) 
			      * Cov_PB_NG_PB(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

			      + calcDeltaFunction(kmag, kmag2, DeltaK)
			      * LegendrePolynomials(ELL, mu2) 
			      * Cov_PB_NG_PB(kvec2, kvec1, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume)

			      + calcDeltaFunction(kmag, kmag3, DeltaK)
			      * LegendrePolynomials(ELL, MU(kvec3, los)) 
			      * Cov_PB_NG_PB(kvec3, kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, nmean, volume);

		result *= (NL * NllL_dash * SllL_dash);
	
		double jacobian = 1.0;
	
		ff_out[index] = result * jacobian;

	}}
	
	return 0;
}

int integrand_cov_PB_NG_P5_full(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
			        double k1_in, 
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume) {

    double mu         = - 1.0 + 2.0 * xx_in[0];
    double mu1  = - 1.0 + 2.0 * xx_in[1];
    double phi1 = 0.0 + 2.0 * M_PI * xx_in[2];
    double mu2  = - 1.0 + 2.0 * xx_in[3];
    double phi2 = 0.0 + 2.0 * M_PI * xx_in[4];
     
    for(int i = 0; i < num_kbin; i++) {
    for(int j = 0; j < num_kbin; j++) {
	
	int index = i * num_kbin + j;

	double kmag = kbin[i];
	double kmag1 = k1_in;
	double kmag2 = kbin[j];

        double kvec[3] = {kmag * sqrt(1.0 - mu * mu), 0.0, kmag * mu};
        double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1) * cos(phi1), kmag1 * sqrt(1.0 - mu1 * mu1) * sin(phi1), kmag1 * mu1};
        double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};

	double kvec3[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };

        double los[3] = {0.0, 0.0, 1.0};

	double NL = 2.0 * double(ELL) + 1.0;
	double NllL_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);

	double result = NL * NllL_dash 
	              * LegendrePolynomials(ELL, mu) 
	              * calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash)
                      * Cov_PB_NG_P5(kvec, kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume);
    
	double jacobian = 1.0;
	ff_out[index] = result * jacobian;
    
    }}
    
    return 0;
}


int integrand_cov_PB_NG_PB_diag_Reconstructed(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume,
				double b1_fid, double R) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
	double mu2  = - 1.0 + 2.0 * xx_in[1];
	double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];

    	for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {
		
		int index = i * num_kbin + j;
 
		double kmag = kbin[i];
		double kmag1 = kbin[j];
		double kmag2 = kbin[j];
	
		double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1};
		double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {0.0, 0.0, 1.0};

		double kmag3 = NORM(kvec3);

		double NL = 2.0 * double(ELL) + 1.0;
		double NllL_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
		double SllL_dash = calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash);

		double result = calcDeltaFunction(kmag, kmag1, DeltaK)
			      * LegendrePolynomials(ELL, mu1) 
			      * Cov_PB_NG_PB_Reconstructed(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume, b1_fid, R)

			      + calcDeltaFunction(kmag, kmag2, DeltaK)
			      * LegendrePolynomials(ELL, mu2) 
			      * Cov_PB_NG_PB_Reconstructed(kvec2, kvec1, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume, b1_fid, R)

			      + calcDeltaFunction(kmag, kmag3, DeltaK)
			      * LegendrePolynomials(ELL, MU(kvec3, los)) 
			      * Cov_PB_NG_PB_Reconstructed(kvec3, kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume, b1_fid, R);

		result *= (NL * NllL_dash * SllL_dash);
	
		double jacobian = 1.0;
	
		ff_out[index] = result * jacobian;

	}}
	
	return 0;

}

int integrand_cov_PB_NG_PB_full_Reconstructed(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
					      double k1_in, 
	                                      double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double DeltaK, double nmean, double volume,
				              double b1_fid, double R) {

	double mu1  = - 1.0 + 2.0 * xx_in[0];
	double mu2  = - 1.0 + 2.0 * xx_in[1];
	double phi2 =   0.0 + 2.0 * M_PI * xx_in[2];

    	for(int i = 0; i < num_kbin; i++) {
    	for(int j = 0; j < num_kbin; j++) {
		
		int index = i * num_kbin + j;
 
		double kmag = kbin[i];
		double kmag1 = k1_in;
		double kmag2 = kbin[j];
	
		double kvec1[3] = {kmag1 * sqrt(1.0 - mu1 * mu1), 0.0, kmag1 * mu1};
		double kvec2[3] = {kmag2 * sqrt(1.0 - mu2 * mu2) * cos(phi2), kmag2 * sqrt(1.0 - mu2 * mu2) * sin(phi2), kmag2 * mu2};
		double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
		double los[3]   = {0.0, 0.0, 1.0};

		double kmag3 = NORM(kvec3);

		double NL = 2.0 * double(ELL) + 1.0;
		double NllL_dash = (2.0*double(ELL_dash)+1.0) * (2.0*double(ell1_dash)+1.0) * (2.0*double(ell2_dash)+1.0);
		double SllL_dash = calcYYY(kvec1, kvec2, los, ell1_dash, ell2_dash, ELL_dash);

		double result = calcDeltaFunction(kmag, kmag1, DeltaK)
			      * LegendrePolynomials(ELL, mu1) 
			      * Cov_PB_NG_PB_Reconstructed(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume, b1_fid, R)

			      + calcDeltaFunction(kmag, kmag2, DeltaK)
			      * LegendrePolynomials(ELL, mu2) 
			      * Cov_PB_NG_PB_Reconstructed(kvec2, kvec1, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume, b1_fid, R)

			      + calcDeltaFunction(kmag, kmag3, DeltaK)
			      * LegendrePolynomials(ELL, MU(kvec3, los)) 
			      * Cov_PB_NG_PB_Reconstructed(kvec3, kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, fz, b1, nmean, volume, b1_fid, R);

		result *= (NL * NllL_dash * SllL_dash);
	
		double jacobian = 1.0;
	
		ff_out[index] = result * jacobian;

	}}
	
	return 0;

}

#endif

