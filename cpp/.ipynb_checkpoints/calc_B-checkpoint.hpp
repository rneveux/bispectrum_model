#ifndef __calc_B__
#define __calc_B__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __pk_lin__
#include "pk_lin.hpp"
#endif

#ifndef __kernel__
#include "kernel.hpp"
#endif

#ifndef __kernel2__
#include "kernel2.hpp"
#endif

#ifndef __sph__
#include "sph.hpp"
#endif

int integrand_B_Tree(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		     double kmag1,                  
		     double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];

	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_FoG(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
                     double kmag1,
                     double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2,
		     double c1, double c2, double knl) {

        /********/
        double mu1  = 1.0;
        double phi1 = 0.0;
        double mu2  = - 1.0 + 2.0 * xx_in[0];
        double phi2 =   0.0;
        double mu   = - 1.0 + 2.0 * xx_in[1];
        double phi  =   0.0 + 2.0 * M_PI * xx_in[2];

        /********/

        for(int j = 0; j < num_k_bin; j++) {

            /********/
            double kvec1[3] = {0.0, 0.0, kmag1};
            double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
            double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
            double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
            /********/

            /********/
            double bispec = Bispectrum_Tree_FoG(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, c1, c2, knl);
            /********/

            /********/
            double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
            double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
            double result = Nlll * Slll * bispec;
            /********/

            double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
            ff_out[j] = result * jacobian;

        }
        return 0;

}

int integrand_B_Tree_DampIvanov(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL,
                     double kmag1,
                     double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2,
		     double rbao, double ks) {

        /********/
        double mu1  = 1.0;
        double phi1 = 0.0;
        double mu2  = - 1.0 + 2.0 * xx_in[0];
        double phi2 =   0.0;
        double mu   = - 1.0 + 2.0 * xx_in[1];
        double phi  =   0.0 + 2.0 * M_PI * xx_in[2];

        /********/

	double S2 = Sig2(rbao, ks);
	double dS2 = dSig2(rbao, ks);

        for(int j = 0; j < num_k_bin; j++) {

            /********/
            double kvec1[3] = {0.0, 0.0, kmag1};
            double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
            double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
            double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
            /********/

            /********/
            double bispec = Bispectrum_Tree_Ivanov_Damping(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, S2, dS2);
            /********/

            /********/
            double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
            double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
            double result = Nlll * Slll * bispec;
            /********/

            double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
            ff_out[j] = result * jacobian;

        }
        return 0;

}

int integrand_B_FoG_Damping_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double alpha_perp, double alpha_parallel, double f,
		            double b1, double b2, double bG2,
		     		double c1, double c2, double knl,
		            double Sigma2, double dSigma2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kbin[j]};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
        double bispec = Bispectrum_Tree_FoG_Damping(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, c1, c2, knl, Sigma2, dSigma2);
        /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_FoG_Damping(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
		            double alpha_perp, double alpha_parallel, double f,
		            double b1, double b2, double bG2,
		     		double c1, double c2, double knl,
		            double Sigma2, double dSigma2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
        double bispec = Bispectrum_Tree_FoG_Damping(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, c1, c2, knl, Sigma2, dSigma2);
        /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_SN_FoG_Damping_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double alpha_perp, double alpha_parallel, double f,
		            double b1,
		     		double c1, double c2, double knl, double Pshot, double Bshot,
		            double Sigma2, double dSigma2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kbin[j]};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
        double bispec = Bispectrum_SN_FoG_Damping(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, f, b1, c1, c2, knl, Pshot, Bshot, Sigma2, dSigma2);
        /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_SN_FoG_Damping(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
		            double alpha_perp, double alpha_parallel, double f,
		            double b1,
		     		double c1, double c2, double knl, double Pshot, double Bshot,
		            double Sigma2, double dSigma2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
        double bispec = Bispectrum_SN_FoG_Damping(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, f, b1, c1, c2, knl, Pshot, Bshot, Sigma2, dSigma2);
        /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

double Bispectrum_Kernel(double * kvec1, double * kvec2, double * los, double f, char * parameters_in) {

    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "b1_b1_b1") {
        return Bispectrum_Kernel_b1_b1_b1(kvec1, kvec2);
    } else if (parameters == "b1_b1_b2") {
        return Bispectrum_Kernel_b1_b1_b2(kvec1, kvec2);
    } else if (parameters == "b1_b1_bG2") {
        return Bispectrum_Kernel_b1_b1_bG2(kvec1, kvec2);
    } else if (parameters == "b1_b1_f") {
        return Bispectrum_Kernel_b1_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_b1_b1_f") {
        return Bispectrum_Kernel_b1_b1_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_b1_f_f") {
        return Bispectrum_Kernel_b1_b1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_b2_f") {
        return Bispectrum_Kernel_b1_b2_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_bG2_f") {
        return Bispectrum_Kernel_b1_bG2_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_f_f") {
        return Bispectrum_Kernel_b1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_f_f_f") {
        return Bispectrum_Kernel_b1_f_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "b2_f_f") {
        return Bispectrum_Kernel_b2_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "bG2_f_f") {
        return Bispectrum_Kernel_bG2_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "f_f_f") {
        return Bispectrum_Kernel_f_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "f_f_f_f") {
        return Bispectrum_Kernel_f_f_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_b1_b1") {
        return Bispectrum_Kernel_c1_b1_b1(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_b2") {
        return Bispectrum_Kernel_c1_b1_b2(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_bG2") {
        return Bispectrum_Kernel_c1_b1_bG2(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_f") {
        return Bispectrum_Kernel_c1_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_b1_b1_f") {
        return Bispectrum_Kernel_c1_b1_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_b1_f_f") {
        return Bispectrum_Kernel_c1_b1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_b2_f") {
        return Bispectrum_Kernel_c1_b2_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_bG2_f") {
        return Bispectrum_Kernel_c1_bG2_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_f_f") {
        return Bispectrum_Kernel_c1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_f_f_f") {
        return Bispectrum_Kernel_c1_f_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_c1_b1") {
        return Bispectrum_Kernel_c1_c1_b1(kvec1, kvec2, los);
    } else if (parameters == "c1_c1_b2") {
        return Bispectrum_Kernel_c1_c1_b2(kvec1, kvec2, los);
    } else if (parameters == "c1_c1_bG2") {
        return Bispectrum_Kernel_c1_c1_bG2(kvec1, kvec2, los);
    } else if (parameters == "c1_c1_f") {
        return Bispectrum_Kernel_c1_c1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_c1_b1_f") {
        return Bispectrum_Kernel_c1_c1_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c1_c1_f_f") {
        return Bispectrum_Kernel_c1_c1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_b1_b1") {
        return Bispectrum_Kernel_c2_b1_b1(kvec1, kvec2, los);
    } else if (parameters == "c2_b1_b2") {
        return Bispectrum_Kernel_c2_b1_b2(kvec1, kvec2, los);
    } else if (parameters == "c2_b1_bG2") {
        return Bispectrum_Kernel_c2_b1_bG2(kvec1, kvec2, los);
    } else if (parameters == "c2_b1_f") {
        return Bispectrum_Kernel_c2_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_b1_b1_f") {
        return Bispectrum_Kernel_c2_b1_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_b1_f_f") {
        return Bispectrum_Kernel_c2_b1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_b2_f") {
        return Bispectrum_Kernel_c2_b2_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_bG2_f") {
        return Bispectrum_Kernel_c2_bG2_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_f_f") {
        return Bispectrum_Kernel_c2_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_f_f_f") {
        return Bispectrum_Kernel_c2_f_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_c1_b1") {
        return Bispectrum_Kernel_c1_c2_b1(kvec1, kvec2, los);
    } else if (parameters == "c2_c1_b2") {
        return Bispectrum_Kernel_c1_c2_b2(kvec1, kvec2, los);
    } else if (parameters == "c2_c1_bG2") {
        return Bispectrum_Kernel_c1_c2_bG2(kvec1, kvec2, los);
    } else if (parameters == "c2_c1_f") {
        return Bispectrum_Kernel_c1_c2_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_c1_b1_f") {
        return Bispectrum_Kernel_c1_c2_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_c1_f_f") {
        return Bispectrum_Kernel_c1_c2_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_c2_b1") {
        return Bispectrum_Kernel_c2_c2_b1(kvec1, kvec2, los);
    } else if (parameters == "c2_c2_b2") {
        return Bispectrum_Kernel_c2_c2_b2(kvec1, kvec2, los);
    } else if (parameters == "c2_c2_bG2") {
        return Bispectrum_Kernel_c2_c2_bG2(kvec1, kvec2, los);
    } else if (parameters == "c2_c2_f") {
        return Bispectrum_Kernel_c2_c2_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_c2_b1_f") {
        return Bispectrum_Kernel_c2_c2_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "c2_c2_f_f") {
        return Bispectrum_Kernel_c2_c2_f_f(kvec1, kvec2, los, f);
    }  else {
        return 0.0;
    }
}

int integrand_B_Kernel(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, double f,
		            double Sigma2, double dSigma2, char * parameters) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double K12 = Bispectrum_Kernel(kvec1, kvec2, los, f, parameters);
	    double K13 = Bispectrum_Kernel(kvec1, kvec3, los, f, parameters);
	    double K23 = Bispectrum_Kernel(kvec2, kvec3, los, f, parameters);
	    /********/

	    double k1 = NORM(kvec1);
	    double k2 = NORM(kvec2);
	    double k3 = NORM(kvec3);

	    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
            double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
            double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);

	    double bispec = K12 * (BAO1 * BAO2 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2) * ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2)
					+ BAO1 * f_pk_no_wiggle(k2) 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k1) * BAO2 
					* ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2))
				+ K13 * (BAO1 * BAO3 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2) * ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ BAO1 * f_pk_no_wiggle(k3) 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2)
					+ f_pk(k1) * BAO3 
					* ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3))
				+ K23 * (BAO2 * BAO3 
					* ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2) * ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ BAO2 * f_pk_no_wiggle(k3) 
					* ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2)
					+ f_pk(k2) * BAO3 
					* ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3));
        
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Kernel_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double f,
		            double Sigma2, double dSigma2, double alpha_perp, double alpha_parallel, char * parameters) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1_in[3] = {0.0, 0.0, kbin[j]};
	    double kvec2_in[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3_in[3] = {- kvec1_in[0] - kvec2_in[0], - kvec1_in[1] - kvec2_in[1], - kvec1_in[2] - kvec2_in[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
        
        double kvec1[3] = { 0.0, 0.0, 0.0 };
        double kvec2[3] = { 0.0, 0.0, 0.0 };
        double kvec3[3] = { 0.0, 0.0, 0.0 };

	    if (alpha_parallel == 0){
	    	kvec1[0] = kvec1_in[0];
			kvec1[1] = kvec1_in[1];
			kvec1[2] = kvec1_in[2];
	    	kvec2[0] = kvec2_in[0];
			kvec2[1] = kvec2_in[1];
			kvec2[2] = kvec2_in[2];
	    	kvec3[0] = kvec3_in[0];
			kvec3[1] = kvec3_in[1];
			kvec3[2] = kvec3_in[2];
	    }
	    else {
	    calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
		calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
		calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
		}
	    
	    /********/
	    double K12 = Bispectrum_Kernel(kvec1, kvec2, los, f, parameters);
	    double K13 = Bispectrum_Kernel(kvec1, kvec3, los, f, parameters);
	    double K23 = Bispectrum_Kernel(kvec2, kvec3, los, f, parameters);
	    /********/

	    double k1 = NORM(kvec1);
	    double k2 = NORM(kvec2);
	    double k3 = NORM(kvec3);

	    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
        double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
        double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
        

	    double bispec = K12 * (BAO1 * BAO2 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2) * ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2)
					+ BAO1 * f_pk_no_wiggle(k2) 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k1) * BAO2 
					* ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2))
				+ K13 * (BAO1 * BAO3 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2) * ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ BAO1 * f_pk_no_wiggle(k3) 
					* ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2)
					+ f_pk(k1) * BAO3 
					* ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3))
				+ K23 * (BAO2 * BAO3 
					* ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2) * ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ BAO2 * f_pk_no_wiggle(k3) 
					* ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2)
					+ f_pk(k2) * BAO3 
					* ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2)
					+ f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3));
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

double Bispectrum_Kernel_SN(double * kvec1, double * los, double f, char * parameters_in) {

    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "Bshot_b1_b1") {
        return Bispectrum_Kernel_Bshot_b1_b1();
    } else if (parameters == "Bshot_b1_f") {
        return Bispectrum_Kernel_Bshot_b1_f(kvec1, los, f);
    } else if (parameters == "Bshot_b1_c1") {
        return Bispectrum_Kernel_Bshot_b1_c1(kvec1, los);
    } else if (parameters == "Bshot_b1_c2") {
        return Bispectrum_Kernel_Bshot_b1_c2(kvec1, los);
    } else if (parameters == "Pshot_f_b1") {
        return Bispectrum_Kernel_Pshot_f_b1(kvec1, los, f);
    } else if (parameters == "Pshot_f_f") {
        return Bispectrum_Kernel_Pshot_f_f(kvec1, los, f);
    } else if (parameters == "Pshot_f_c1") {
        return Bispectrum_Kernel_Pshot_f_c1(kvec1, los, f);
    } else if (parameters == "Pshot_f_c2") {
        return Bispectrum_Kernel_Pshot_f_c2(kvec1, los, f);
    } else {
        return 0.0;
    }
}

int integrand_B_Kernel_SN(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, double f,
		            double Sigma2, double dSigma2, char * parameters) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double K1 = Bispectrum_Kernel_SN(kvec1, los, f, parameters);
	    double K2 = Bispectrum_Kernel_SN(kvec2, los, f, parameters);
	    double K3 = Bispectrum_Kernel_SN(kvec3, los, f, parameters);
	    /********/

	    double k1 = NORM(kvec1);
	    double k2 = NORM(kvec2);
	    double k3 = NORM(kvec3);

	    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
            double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
            double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);

	    double bispec = K1 * (BAO1 * ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2) + f_pk_no_wiggle(k1))
				+ K2 * (BAO2 * ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2) + f_pk_no_wiggle(k2))
				+ K3 * (BAO3 * ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2) + f_pk_no_wiggle(k3));
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Kernel_SN_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double f,
		            double Sigma2, double dSigma2, double alpha_perp, double alpha_parallel, char * parameters) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1_in[3] = {0.0, 0.0, kbin[j]};
	    double kvec2_in[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3_in[3] = {- kvec1_in[0] - kvec2_in[0], - kvec1_in[1] - kvec2_in[1], - kvec1_in[2] - kvec2_in
            [2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	   
        
        double kvec1[3] = { 0.0, 0.0, 0.0 };
        double kvec2[3] = { 0.0, 0.0, 0.0 };
        double kvec3[3] = { 0.0, 0.0, 0.0 };
	    /********/

	    if (alpha_parallel == 0){
	    	kvec1[0] = kvec1_in[0];
			kvec1[1] = kvec1_in[1];
			kvec1[2] = kvec1_in[2];
	    	kvec2[0] = kvec2_in[0];
			kvec2[1] = kvec2_in[1];
			kvec2[2] = kvec2_in[2];
	    	kvec3[0] = kvec3_in[0];
			kvec3[1] = kvec3_in[1];
			kvec3[2] = kvec3_in[2];
	    }
	    else {
	    calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
		calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
		calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
		}
	    
	    /********/
	    double K1 = Bispectrum_Kernel_SN(kvec1, los, f, parameters);
	    double K2 = Bispectrum_Kernel_SN(kvec2, los, f, parameters);
	    double K3 = Bispectrum_Kernel_SN(kvec3, los, f, parameters);
	    /********/

	    double k1 = NORM(kvec1);
	    double k2 = NORM(kvec2);
	    double k3 = NORM(kvec3);

	    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
            double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
            double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);

	    double bispec = K1 * (BAO1 * ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2) + f_pk_no_wiggle(k1))
				+ K2 * (BAO2 * ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2) + f_pk_no_wiggle(k2))
				+ K3 * (BAO3 * ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2) + f_pk_no_wiggle(k3));
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

double Bispectrum_Kernel_PNG(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2, char * parameters_in) {

    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "fnlloc_b1_b1_b1") {
        return Bispectrum_Kernel_PNG_fnlloc_b1_b1_b1(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlloc_b1_b1_f") {
        return Bispectrum_Kernel_PNG_fnlloc_b1_b1_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlloc_b1_f_f") {
        return Bispectrum_Kernel_PNG_fnlloc_b1_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlloc_f_f_f") {
        return Bispectrum_Kernel_PNG_fnlloc_f_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlequi_b1_b1_b1") {
        return Bispectrum_Kernel_PNG_fnlequi_b1_b1_b1(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlequi_b1_b1_f") {
        return Bispectrum_Kernel_PNG_fnlequi_b1_b1_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlequi_b1_f_f") {
        return Bispectrum_Kernel_PNG_fnlequi_b1_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlequi_f_f_f") {
        return Bispectrum_Kernel_PNG_fnlequi_f_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_b1_b1_b1") {
        return Bispectrum_Kernel_PNG_fnlortho_b1_b1_b1(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_b1_b1_f") {
        return Bispectrum_Kernel_PNG_fnlortho_b1_b1_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_b1_f_f") {
        return Bispectrum_Kernel_PNG_fnlortho_b1_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_f_f_f") {
        return Bispectrum_Kernel_PNG_fnlortho_f_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_LSS_b1_b1_b1") {
        return Bispectrum_Kernel_PNG_fnlortho_LSS_b1_b1_b1(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_LSS_b1_b1_f") {
        return Bispectrum_Kernel_PNG_fnlortho_LSS_b1_b1_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_LSS_b1_f_f") {
        return Bispectrum_Kernel_PNG_fnlortho_LSS_b1_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else if (parameters == "fnlortho_LSS_f_f_f") {
        return Bispectrum_Kernel_PNG_fnlortho_LSS_f_f_f(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
    } else {
        return 0.0;
    }
}

int integrand_B_Kernel_PNG_diag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double f,
		            double Sigma2, double dSigma2, char * parameters) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kbin[j]};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};  
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Kernel_PNG(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2, parameters);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}



int integrand_B_Tree_NoWiggle(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1,                  
		              double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Growth(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1,                  
		              double sigma8) {

	/********/
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_RealSpace_DarkMatter_Growth(kvec1, kvec2, kvec3, sigma8);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0);
	    double Slll = LegendrePolynomials(ell1, mu2);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Shift(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1,                  
		              double sigma8) {

	/********/
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_RealSpace_DarkMatter_Shift(kvec1, kvec2, kvec3, sigma8);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0);
	    double Slll = LegendrePolynomials(ell1, mu2);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_NoWiggle_RealSpace_DarkMatter_Tidal(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1,                  
		              double sigma8) {

	/********/
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_RealSpace_DarkMatter_Tidal(kvec1, kvec2, kvec3, sigma8);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0);
	    double Slll = LegendrePolynomials(ell1, mu2);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_Tree_BAO_RealSpace_DarkMatter_Growth(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1,                  
		              double sigma8, double sigma2_perp) {

	/********/
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_RealSpace_DarkMatter_Growth(kvec1, kvec2, kvec3, sigma8, sigma2_perp);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0);
        double Slll = LegendrePolynomials(ell1, mu2);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_Tree_BAO_RealSpace_DarkMatter_Shift(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1,                  
		              double sigma8, double sigma2_perp) {

	/********/
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_RealSpace_DarkMatter_Shift(kvec1, kvec2, kvec3, sigma8, sigma2_perp);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0);
        double Slll = LegendrePolynomials(ell1, mu2);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_BAO_RealSpace_DarkMatter_Tidal(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1,                  
		              double sigma8, double sigma2_perp) {

	/********/
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_RealSpace_DarkMatter_Tidal(kvec1, kvec2, kvec3, sigma8, sigma2_perp);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0);
        double Slll = LegendrePolynomials(ell1, mu2);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_Tree_BAO(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		         double kmag1,                  
		         double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_Reconstructed(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                   double kmag1,                  
		                   double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_Reconstructed(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_Tree_NoWiggle_Reconstructed(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                            double kmag1,                  
		                            double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double b2, double bK2, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_Reconstructed(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_Local(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                  double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Local(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_Equilateral(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                        double kmag1, 
             	                        double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Equilateral(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_Orthogonal(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                       double kmag1, 
             	                       double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Orthogonal(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}



int integrand_B_Tree_b1_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    double bispec = Bispectrum_Tree_b1_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}

	return 0;
}

int integrand_B_Tree_b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b1_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b1_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_b2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b2_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_b2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b2_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_b2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b2_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

/********/

int integrand_B_Tree_bK2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_bK2_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_bK2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_bK2_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_bK2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_bK2_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


/********/

int integrand_B_Tree_b1f_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b1f_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_b1f_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b1f_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_b1f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_b1f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


/********/

int integrand_B_Tree_ff_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_ff_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}



int integrand_B_Tree_Reconstructed_b1b1_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                              double kmag1, 
			                      double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_Reconstructed_b1b1_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_Reconstructed_b1b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_Reconstructed_b1b1_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_Reconstructed_b1b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_Reconstructed_b1b1_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_Reconstructed_b1f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_Reconstructed_b1f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_Reconstructed_ff_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_Reconstructed_ff_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

/***************/
/**  NoWiggle **/
/***************/

int integrand_B_Tree_NoWiggle_b1_b1_b1(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
        double alpha_perp, double alpha_parallel
        ) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b1_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}

	return 0;
}

int integrand_B_Tree_NoWiggle_b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b1_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_NoWiggle_b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b1_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_NoWiggle_b2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b2_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_NoWiggle_b2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b2_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_NoWiggle_b2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b2_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

/********/

int integrand_B_Tree_NoWiggle_bK2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_bK2_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_NoWiggle_bK2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_bK2_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_NoWiggle_bK2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_bK2_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


/********/

int integrand_B_Tree_NoWiggle_b1f_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b1f_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_NoWiggle_b1f_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b1f_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_NoWiggle_b1f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_b1f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


/********/

int integrand_B_Tree_NoWiggle_ff_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_ff_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_NoWiggle_f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}



int integrand_B_Tree_NoWiggle_Reconstructed_b1b1_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                              double kmag1, 
			                      double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_Reconstructed_b1b1_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_NoWiggle_Reconstructed_b1b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_Reconstructed_b1b1_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_NoWiggle_Reconstructed_b1b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_Reconstructed_b1b1_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_NoWiggle_Reconstructed_b1f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_Reconstructed_b1f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_NoWiggle_Reconstructed_ff_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		                             double kmag1, 
			                     double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_NoWiggle_Reconstructed_ff_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, b1_fid, R);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}








































//int integrand_B_DDD_SS(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, 
//		       int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
//		       double kmag1, 
//		       double alpha_perp, double alpha_parallel) {
//
//	/********/
//	double mu1  = 1.0;
//	double phi1 = 0.0;
//	double mu2  = - 1.0 + 2.0 * xx_in[0];
//	double phi2 =   0.0;
//	double mu   = - 1.0 + 2.0 * xx_in[1];      
//	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
//	/********/
//
//	for(int j = 0; j < num_k_bin; j++) {
//	
//	    /********/
//	    double kvec1[3] = {0.0, 0.0, kmag1};
//	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
//	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
//
//	    double kvec1_dash[3] = { 0.0, 0.0, 0.0, };
//	    double kvec2_dash[3] = { 0.0, 0.0, 0.0, };
//	    /********/
//
//	    /********/
//	    calcTrueWavevector(kvec1, los, alpha_perp, alpha_parallel, kvec1_dash);
//	    calcTrueWavevector(kvec2, los, alpha_perp, alpha_parallel, kvec2_dash);
//	    /********/
//
//	    /********/
//	    double YYY_YYY = calcYYY_YYY(kvec1, kvec2, kvec1_dash, kvec2_dash, los, ell1, ell2, ELL, ell1_dash, ell2_dash, ELL_dash);
//	    /********/
//	    
//	    /********/
//	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
//	    double result = Nlll * YYY_YYY;
//	    /********/
//	    
//	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
//	    ff_out[j] = result * jacobian;
//	
//	}
//
//	return 0;
//}


double get_nm(int n) {

    if(0) {
    } else if (n == 0) {
	return 1.0;
    } else if (n == 1) {
	return 1.0;
    } else if (n == 2) {
	return 2.0;
    } else if (n == 3) {
	return 6.0;
    } else if (n == 4) {
	return 24.0;
    } else {
	return 1.0e-10;
    }

}

int integrand_SS(
        double * xx_in, int ndim, double * ff_out, int ncomp,  
        int ell1, int ell2, int ELL, int ell1_dash, int ell2_dash, int ELL_dash, 
        int n, int m,
        double * epsilon, int num_epsilon) {

	/********/
//	double mu1  = 1.0;
//	double phi1 = 0.0;
	double mu1  = - 1.0 + 2.0 * xx_in[0];      
	double phi1 =   0.0 + 2.0 * M_PI * xx_in[1];
	double mu2  = - 1.0 + 2.0 * xx_in[2];
//	double phi2 =   0.0 + 2.0 * M_PI * xx_in[3];
	double phi2 =   0.0;
	/********/

    for(int i = 0; i < num_epsilon; i++) {

        /********/
//    	double kvec1_hat[3] = {0.0, 0.0, 1.0};
//    	double kvec2_hat[3] = {sqrt(1.0 - mu2 * mu2), 0.0, mu2};
//    	double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
    	double kvec1_hat[3] = {sqrt(1.0 - mu1 * mu1) * cos(phi1), sqrt(1.0 - mu1 * mu1) * sin(phi1), mu1};
    	double kvec2_hat[3] = {sqrt(1.0 - mu2 * mu2) * cos(phi2), sqrt(1.0 - mu2 * mu2) * sin(phi2), mu2};
//    	double kvec2_hat[3] = {sqrt(1.0 - mu2 * mu2), 0.0, mu2};
    	double los[3]   = {0.0, 0.0, 1.0};
     
    	double kvec1_hat_dash[3] = { 0.0, 0.0, 1.0 };
    	double kvec2_hat_dash[3] = { 0.0, 0.0, 1.0 };
    	calcTrueWavevectorHat(kvec1_hat, los, epsilon[i], kvec1_hat_dash);
    	calcTrueWavevectorHat(kvec2_hat, los, epsilon[i], kvec2_hat_dash);
    	/********/
    
    	/********/
    	double mu1_n = MU(kvec1_hat, los); 
    	double mu2_m = MU(kvec2_hat, los); 
    	double L1 = (3.0 * mu1_n * mu1_n - 1.0) / 2.0;
    	double L2 = (3.0 * mu2_m * mu2_m - 1.0) / 2.0;
    
    	double fac1_a =  (2.0/3.0) * L1 * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
    	double fac1_b = 1.0 + (1.0/3.0) * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
    	double Delta1 = pow( 1.0 + fac1_a / fac1_b, 0.5 ) - 1.0;
    	double Delta1_n = pow( Delta1, double(n));
    
    	double fac2_a =  (2.0/3.0) * L2 * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
    	double fac2_b = 1.0 + (1.0/3.0) * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
    	double Delta2 = pow( 1.0 + fac2_a / fac2_b, 0.5 ) - 1.0;
    	double Delta2_m = pow( Delta2, double(m));
    	/********/
    
    	/********/
    	double YYY_YYY = calcYYY_YYY(kvec1_hat, kvec2_hat, kvec1_hat_dash, kvec2_hat_dash, los, ell1, ell2, ELL, ell1_dash, ell2_dash, ELL_dash);
    	/********/
    	
    	/********/
    	double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
    	double NN = get_nm(n);
    	double MM = get_nm(m);
    	double result = Nlll * YYY_YYY * Delta1_n * Delta2_m / (NN * MM);
    	/********/
    	
    	double jacobian = 1.0;
    	ff_out[i] = result * jacobian;

    }

	return 0;

}

int integrand_SSpow(
        double * xx_in, int ndim, double * ff_out, int ncomp,  
        int ELL, int ELL_dash, 
        int n,
        double * epsilon, int num_epsilon) {

	/********/
	double mu  = - 1.0 + 2.0 * xx_in[0];      
	double phi =   0.0 + 2.0 * M_PI * xx_in[1];
	/********/

    for(int i = 0; i < num_epsilon; i++) {

        /********/  
    	double kvec_hat[3] = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
    	double los[3]   = {0.0, 0.0, 1.0};
     
    	double kvec_hat_dash[3] = { 0.0, 0.0, 1.0 };
    	calcTrueWavevectorHat(kvec_hat, los, epsilon[i], kvec_hat_dash);
    	/********/
    
    	/********/
    	double mu_n = MU(kvec_hat, los); 
    	double mu_dash = MU(kvec_hat_dash, los); 
    	double L = (3.0 * mu_n * mu_n - 1.0) / 2.0;
    
    	double fac_a =  (2.0/3.0) * L * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
    	double fac_b = 1.0 + (1.0/3.0) * ( pow(1.0 + epsilon[i], - 6.0) - 1.0 );
    	double Delta = pow( 1.0 + fac_a / fac_b, 0.5 ) - 1.0;
    	double Delta_n = pow( Delta, double(n));
    	/********/
    
    	/********/
    	double YYY_YYY = LegendrePolynomials(ELL, mu_n) * LegendrePolynomials(ELL_dash, mu_dash);
    	/********/
    	
    	/********/
    	double Nlll = 2.0*double(ELL)+1.0;
    	double NN = get_nm(n);
    	double result = Nlll * YYY_YYY * Delta_n / NN;
    	/********/
    	
    	double jacobian = 1.0;
    	ff_out[i] = result * jacobian;

    }

	return 0;

}

/* PNG from PB */ 

double kmax_integ = 1.0;
//double kmax_integ = 10.0;
double kmin_integ = 1.0e-4;

int integrand_B_NonGaussian_From_PB_Local(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				          double kmag1, 
             	                          double alpha_perp, double alpha_parallel, double sigma8, double fz, 
    					  double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO
					    ) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local(kvec1, kvec2, kvec3, pvec, los, 
                                                                 alpha_perp, alpha_parallel, sigma8, fz, 
								 b1, b2, b3, bK2, bK3, bDK, bO);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Equilateral(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				          double kmag1, 
             	                          double alpha_perp, double alpha_parallel, double sigma8, double fz, 
    					  double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO
					    ) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Equilateral(kvec1, kvec2, kvec3, pvec, los, 
                                                                 alpha_perp, alpha_parallel, sigma8, fz, 
								 b1, b2, b3, bK2, bK3, bDK, bO);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_Reconstructed(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				          double kmag1, 
             	                          double alpha_perp, double alpha_parallel, double sigma8, double fz, 
    					  double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double b1_fid, double R
					    ) {


	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_Reconstructed(kvec1, kvec2, kvec3, pvec, los, 
                                                                 alpha_perp, alpha_parallel, sigma8, fz, 
								 b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Equilateral_Reconstructed(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				          double kmag1, 
             	                          double alpha_perp, double alpha_parallel, double sigma8, double fz, 
    					  double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double b1_fid, double R
					    ) {




	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Equilateral_Reconstructed(kvec1, kvec2, kvec3, pvec, los, 
                                                                 alpha_perp, alpha_parallel, sigma8, fz, 
								 b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R);


	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_Tree_BAO_b1_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b1_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}

	return 0;
}

int integrand_B_Tree_BAO_b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b1_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_BAO_b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b1_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_BAO_b2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b2_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_BAO_b2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b2_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_BAO_b2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b2_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

/********/

int integrand_B_Tree_BAO_bK2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_bK2_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_BAO_bK2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_bK2_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_BAO_bK2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_bK2_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


/********/

int integrand_B_Tree_BAO_b1f_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		              double kmag1, 
			      double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b1f_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_BAO_b1f_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		             double kmag1, 
			     double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b1f_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_Tree_BAO_b1f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_b1f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


/********/

int integrand_B_Tree_BAO_ff_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_ff_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}


int integrand_B_Tree_BAO_f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
		            double kmag1, 
			    double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};	    
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_Tree_BAO_f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_B_NonGaussian_Local_b1_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				           double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Local_b1_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_Local_b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				           double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Local_b1_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_Local_b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				         double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Local_b1_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_Local_f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				         double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Local_f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_Equilateral_b1_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				           double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Equilateral_b1_b1_b1(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_Equilateral_b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				           double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Equilateral_b1_b1_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_Equilateral_b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				         double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Equilateral_b1_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_Equilateral_f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				         double kmag1, 
             	                           double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_Equilateral_f_f_f(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel);;
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

/** PNG from PB decomposed */
int integrand_B_NonGaussian_From_PB_Local_b1_b1_b1_R(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	double R = 10.0;


	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_b1_b1_R(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel, R);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_b1_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_b1_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1_b2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_b2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_b2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_b2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1_bK2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_bK2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_bK2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_bK2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1_b1f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_b1f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_b1f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_b1f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1_ff_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_ff_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_ff_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_ff_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                  double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                         alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                 double kmag1, 
             	                                 double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1_f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                        alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}



int integrand_B_NonGaussian_From_PB_Local_b2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_b1_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_b1_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b2_b2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_b2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b2_b2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_b2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b2_bK2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_bK2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b2_bK2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_bK2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b2_b1f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_b1f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b2_b1f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_b1f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b2_ff_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_ff_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b2_ff_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_ff_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b2_f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                  double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                         alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                 double kmag1, 
             	                                 double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b2_f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                        alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}



int integrand_B_NonGaussian_From_PB_Local_bK2_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_b1_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_bK2_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_b1_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_bK2_b2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_b2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_bK2_b2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_b2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_bK2_bK2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_bK2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_bK2_bK2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_bK2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_bK2_b1f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_b1f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_bK2_b1f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_b1f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_bK2_ff_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_ff_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_bK2_ff_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_ff_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_bK2_f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                  double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                         alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_bK2_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                 double kmag1, 
             	                                 double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_bK2_f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                        alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_b1_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_b1_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1f_b2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_b2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_b2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_b2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1f_bK2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_bK2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_bK2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_bK2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1f_b1f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_b1f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_b1f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_b1f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_b1f_ff_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_ff_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_ff_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_ff_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                  double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                         alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_b1f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                 double kmag1, 
             	                                 double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_b1f_f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                        alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}



int integrand_B_NonGaussian_From_PB_Local_ff_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_b1_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_ff_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_b1_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_ff_b2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_b2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_ff_b2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_b2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_ff_bK2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_bK2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_ff_bK2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_bK2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_ff_b1f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_b1f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_ff_b1f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_b1f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_ff_ff_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_ff_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_ff_ff_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_ff_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_ff_f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                  double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                         alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_ff_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                 double kmag1, 
             	                                 double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_ff_f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                        alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_b1_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_b1_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_f_b2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_b2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_b2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_b2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_f_bK2_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_bK2_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_bK2_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_bK2_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_f_b1f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_b1f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                           alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_b1f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_b1f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}

int integrand_B_NonGaussian_From_PB_Local_f_ff_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                   double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_ff_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_ff_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                   double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_ff_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                          alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_f_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                  double kmag1, 
             	                                  double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_f_b1(kvec1, kvec2, kvec3, pvec, los, 
                                                                         alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}


int integrand_B_NonGaussian_From_PB_Local_f_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, 
				                 double kmag1, 
             	                                 double alpha_perp, double alpha_parallel) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

    	double kmax = kmax_integ;
        double kmin = kmin_integ;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[3] );
        double mu_p = -1.0 + 2.0 * xx_in[4];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[5];   

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/

	    /********/
	    double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
	    /********/
	    
	    /********/
	    double bispec = Bispectrum_NonGaussian_From_PB_Local_f_f_f(kvec1, kvec2, kvec3, pvec, los, 
                                                                        alpha_perp, alpha_parallel);
	    /********/
 
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);;
	    ff_out[j] = result * jacobian;
	
	}
	return 0;

}



int integrand_B_Tree_BAO_Template(
        double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1, 
		double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para,
        char * parameters) {

	/********/
	double mu1  = 1.0;
	double phi1 = 0.0;
	double mu2  = - 1.0 + 2.0 * xx_in[0];
	double phi2 =   0.0;
	double mu   = - 1.0 + 2.0 * xx_in[1];      
	double phi  =   0.0 + 2.0 * M_PI * xx_in[2];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    /********/
	    double kvec1[3] = {0.0, 0.0, kmag1};
	    double kvec2[3] = {kbin[j] * sqrt(1.0 - mu2 * mu2), 0.0, kbin[j] * mu2};
	    double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    double bispec = Bispectrum_Tree_BAO_Template(kvec1, kvec2, kvec3, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para, parameters);
	    /********/
	    
	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    double result = Nlll * Slll * bispec;
	    /********/
	    
	    double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = result * jacobian;
	
	}

	return 0;
}






#endif


