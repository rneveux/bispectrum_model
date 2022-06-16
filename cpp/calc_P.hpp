#ifndef __calc_P__
#define __calc_P__

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

int integrand_P_Tree(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double f, double b1) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree(kvec, los, alpha_perp, alpha_parallel, f, b1);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Damping_Tree(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
				double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, double dSigma2) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_Damping(kvec, los, alpha_perp, alpha_parallel, f, b1, Sigma2, dSigma2);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Damping_Tree_for_1loop(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
				double alpha_perp, double alpha_parallel, double f, double b1, double Sigma2, double dSigma2) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_Damping_for_1loop(kvec, los, alpha_perp, alpha_parallel, f, b1, Sigma2, dSigma2);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Counterterm(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
					double alpha_perp, double alpha_parallel, double f, double b1, double c0, double c1, double c2, double ch, double knl) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Counterterm(kvec, los, alpha_perp, alpha_parallel, f, b1, c0, c1, c2, ch, knl);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}


int integrand_P_Damping_Counterterm(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
					double alpha_perp, double alpha_parallel, double f, double b1, double c0, double c1, double c2, double ch, double knl,
					double Sigma2, double dSigma2) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Counterterm_Damping(kvec, los, alpha_perp, alpha_parallel, f, b1, c0, c1, c2, ch, knl, Sigma2, dSigma2);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_1loop(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
					double alpha_perp, double alpha_parallel, double f, double b1,
					double b2, double bG2, double bGamma3) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 3.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_1loop(kvec, pvec, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, bGamma3);
        	double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Damping_1loop(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
					double alpha_perp, double alpha_parallel, double f, double b1,
					double b2, double bG2, double bGamma3,
					double Sigma2, double dSigma2) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);

        	//double pnw = Powerspectrum_1loop_Damping_NW(kvec, pvec, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, b3, bG3, bG2d, bGamma3);
        	//double pw = Powerspectrum_1loop_Damping_tot(kvec, pvec, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, b3, bG3, bG2d, bGamma3) - pnw;

        	//double D = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);


        	//double Ptot = (pnw + pw * D);
        	
        	double result = Nlll * L * Powerspectrum_1loop_Damping(kvec, pvec, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, bGamma3, Sigma2, dSigma2); //Ptot;
        	double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Damping_1loop_k_vector(double * xx_in, int ndim, double * ff_out, int ncomp, double * kx, double * ky, double * kz, int num_kbin, 
                   double alpha_perp, double alpha_parallel, double f, double b1,
					double b2, double bG2, double bGamma3,
					double Sigma2, double dSigma2) {

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[0] );
        double mu_p = -1.0 + 2.0 * xx_in[1];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[2];

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {kx[i], ky[i], kz[i]};
        	double los[3]  = {0.0, 0.0, 1.0};

            double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};

        	double result = Powerspectrum_1loop_Damping(kvec, pvec, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, bGamma3, Sigma2, dSigma2); 
        
        	double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Damping_1loop_nw(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
					double alpha_perp, double alpha_parallel, double f, double b1,
					double b2, double bG2, double bGamma3,
					double Sigma2, double dSigma2) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);

        	double pnw = Powerspectrum_1loop_Damping_NW(kvec, pvec, los, alpha_perp, alpha_parallel, f, b1, b2, bG2, bGamma3);
        
            double result = Nlll * L * pnw;
        
        	double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

double Power_Spectrum_Kernel_Tree(double * kvec, double * los, double f, char * parameters_in) {

    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "b1_b1") {
        return Power_Spectrum_Kernel_Tree_b1_b1();
    } else if (parameters == "b1_f") {
        return Power_Spectrum_Kernel_Tree_b1_f(kvec, los, f);
    } else if (parameters == "f_f") {
        return Power_Spectrum_Kernel_Tree_f_f(kvec, los, f);
    } else {
        return 0.0;
    }
}

int integrand_P_Kernel_Tree(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, 
		            double f, double Sigma2, double dSigma2, char * parameters) {

	/********/
	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    	/********/
	    	double kvec[3] = {0.0, 0.0, kbin[j]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

        	double k = NORM(kvec);
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);

		double BAO = f_pk(k) - f_pk_no_wiggle(k);

		double D = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);

		double K = Power_Spectrum_Kernel_Tree(kvec, los, f, parameters);

		double powerspec = K * (f_pk_no_wiggle(k) + BAO * D);
	    
		/********/
		double result = Nlll * L * powerspec;
		/********/
		    
		double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

		ff_out[j] = result * jacobian;
	
	}
	return 0;
}

double Power_Spectrum_Kernel_Counterterm(double * kvec, double * los, double f, char * parameters_in) {

    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "c0_b1") {
        return Power_Spectrum_Kernel_Counterterm_c0_b1(kvec);
    } else if (parameters == "c0_f") {
        return Power_Spectrum_Kernel_Counterterm_c0_f(kvec, los, f);
    } else if (parameters == "c1_b1") {
        return Power_Spectrum_Kernel_Counterterm_c1_b1(kvec, los);
    } else if (parameters == "c1_f") {
        return Power_Spectrum_Kernel_Counterterm_c1_f(kvec, los, f);
    } else if (parameters == "c2_b1") {
        return Power_Spectrum_Kernel_Counterterm_c2_b1(kvec, los);
    } else if (parameters == "c2_f") {
        return Power_Spectrum_Kernel_Counterterm_c2_f(kvec, los, f);
    } else {
        return 0.0;
    }
}

int integrand_P_Kernel_Counterterm(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, 
		            double f, double Sigma2, double dSigma2, char * parameters) {

	/********/
	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    	/********/
	    	double kvec[3] = {0.0, 0.0, kbin[j]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

        	double k = NORM(kvec);
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);

		double BAO = f_pk(k) - f_pk_no_wiggle(k);

		double D = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);

		double K = Power_Spectrum_Kernel_Counterterm(kvec, los, f, parameters);

		double powerspec = 2.0 * K * (f_pk_no_wiggle(k) + BAO * D);
	    
		/********/
		double result = Nlll * L * powerspec;
		/********/
		    
		double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));

		ff_out[j] = result * jacobian;
	
	}
	return 0;
}

double Power_Spectrum_Kernel_1loop_22(double * kvec1, double * kvec2, double * los, double f, char * parameters_in) {

    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "b2_b2") {
        return Power_Spectrum_Kernel_1loop_22_b2_b2();
    } else if (parameters == "b2_bG2") {
        return Power_Spectrum_Kernel_1loop_22_b2_bG2(kvec1, kvec2);
    } else if (parameters == "b1_b2") {
        return Power_Spectrum_Kernel_1loop_22_b1_b2(kvec1, kvec2);
    } else if (parameters == "b2_f") {
        return Power_Spectrum_Kernel_1loop_22_b2_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_b2_f") {
        return Power_Spectrum_Kernel_1loop_22_b1_b2_f(kvec1, kvec2, los, f);
    } else if (parameters == "b2_f_f") {
        return Power_Spectrum_Kernel_1loop_22_b2_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "bG2_bG2") {
        return Power_Spectrum_Kernel_1loop_22_bG2_bG2(kvec1, kvec2);
    } else if (parameters == "b1_bG2") {
        return Power_Spectrum_Kernel_1loop_22_b1_bG2(kvec1, kvec2);
    } else if (parameters == "bG2_f") {
        return Power_Spectrum_Kernel_1loop_22_bG2_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_bG2_f") {
        return Power_Spectrum_Kernel_1loop_22_b1_bG2_f(kvec1, kvec2, los, f);
    } else if (parameters == "bG2_f_f") {
        return Power_Spectrum_Kernel_1loop_22_bG2_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_b1") {
        return Power_Spectrum_Kernel_1loop_22_b1_b1(kvec1, kvec2);
    } else if (parameters == "b1_f") {
        return Power_Spectrum_Kernel_1loop_22_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_b1_f") {
        return Power_Spectrum_Kernel_1loop_22_b1_b1_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_f_f") {
        return Power_Spectrum_Kernel_1loop_22_b1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "f_f") {
        return Power_Spectrum_Kernel_1loop_22_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "f_f_f") {
        return Power_Spectrum_Kernel_1loop_22_f_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_b1_f_f") {
        return Power_Spectrum_Kernel_1loop_22_b1_b1_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "b1_f_f_f") {
        return Power_Spectrum_Kernel_1loop_22_b1_f_f_f(kvec1, kvec2, los, f);
    } else if (parameters == "f_f_f_f") {
        return Power_Spectrum_Kernel_1loop_22_f_f_f_f(kvec1, kvec2, los, f);
    } else {
        return 0.0;
    }
}

int integrand_P_Kernel_1loop_22(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, 
		            double f, double Sigma2, double dSigma2, char * parameters) {

	/********/
	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    	/********/
	    	double kvec[3] = {0.0, 0.0, kbin[j]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

        	double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
		double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
		double k_M_p = NORM(kvec_M_pvec);
		double p = NORM(pvec);

		double BAOp = f_pk(p) - f_pk_no_wiggle(p);
		double BAOkMp = f_pk(k_M_p) - f_pk_no_wiggle(k_M_p);

		double Dp = ExpDamping_Ivanov(pvec, los, f, Sigma2, dSigma2);
		double DkMp = ExpDamping_Ivanov(kvec_M_pvec, los, f, Sigma2, dSigma2);
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);

		double K = Power_Spectrum_Kernel_1loop_22(pvec, kvec_M_pvec, los, f, parameters);

		double powerspec = 2.0 * K * (f_pk_no_wiggle(k_M_p) + BAOkMp * DkMp) * (f_pk_no_wiggle(p) + BAOp * Dp);
	    
		/********/
		double result = Nlll * L * powerspec;
		/********/
		    
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);

		ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_P_Kernel_1loop_22_norm(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, 
		            double f, double Sigma2, double dSigma2) {

	/********/
	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
    	/********/
    	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

    	double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
		double p = NORM(pvec);

		double BAOp = f_pk(p) - f_pk_no_wiggle(p);

		double Dp = ExpDamping_Ivanov(pvec, los, f, Sigma2, dSigma2);
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);

		double powerspec = - 2.0 / 4.0 * (f_pk_no_wiggle(p) + BAOp * Dp) * (f_pk_no_wiggle(p) + BAOp * Dp);
	    
		/********/
		double result = Nlll * L * powerspec;
		/********/
		    
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);

		ff_out[j] = result * jacobian;
	
	}
	return 0;
}

double Power_Spectrum_Kernel_1loop_13(double * kvec1, double * kvec2, double * kvec3, double * los, double f, char * parameters_in) {

    std::string parameters = parameters_in;

    if(0) {
    } else if (parameters == "b1_b3") {
        return Power_Spectrum_Kernel_1loop_13_b1_b3();
    } else if (parameters == "b1_bG3") {
        return Power_Spectrum_Kernel_1loop_13_b1_bG3(kvec1, kvec2, kvec3);
    } else if (parameters == "b1_bG2d") {
        return Power_Spectrum_Kernel_1loop_13_b1_bG2d(kvec1, kvec2, kvec3);
    } else if (parameters == "b1_bGamma3") {
        return Power_Spectrum_Kernel_1loop_13_b1_bGamma3(kvec1, kvec2, kvec3);
    } else if (parameters == "b1_b1") {
        return Power_Spectrum_Kernel_1loop_13_b1_b1(kvec1, kvec2, kvec3);
    } else if (parameters == "b1_f") {
        return Power_Spectrum_Kernel_1loop_13_b1_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b1_b1_f_f") {
        return Power_Spectrum_Kernel_1loop_13_b1_b1_f_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b1_f_f_f") {
        return Power_Spectrum_Kernel_1loop_13_b1_f_f_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b1_b1_f") {
        return Power_Spectrum_Kernel_1loop_13_b1_b1_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b1_f_f") {
        return Power_Spectrum_Kernel_1loop_13_b1_f_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b1_b2") {
        return Power_Spectrum_Kernel_1loop_13_b1_b2(kvec1, kvec2, kvec3);
    } else if (parameters == "b1_bG2") {
        return Power_Spectrum_Kernel_1loop_13_b1_bG2(kvec1, kvec2, kvec3);
    } else if (parameters == "b1_b2_f") {
        return Power_Spectrum_Kernel_1loop_13_b1_b2_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b1_bG2_f") {
        return Power_Spectrum_Kernel_1loop_13_b1_bG2_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b3_f") {
        return Power_Spectrum_Kernel_1loop_13_b3_f(kvec3, los, f);
    } else if (parameters == "bG3_f") {
        return Power_Spectrum_Kernel_1loop_13_bG3_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "bG2d_f") {
        return Power_Spectrum_Kernel_1loop_13_bG2d_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "bGamma3_f") {
        return Power_Spectrum_Kernel_1loop_13_bGamma3_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "f_f") {
        return Power_Spectrum_Kernel_1loop_13_f_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "f_f_f_f") {
        return Power_Spectrum_Kernel_1loop_13_f_f_f_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "f_f_f") {
        return Power_Spectrum_Kernel_1loop_13_f_f_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b2_f") {
        return Power_Spectrum_Kernel_1loop_13_b2_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "bG2_f") {
        return Power_Spectrum_Kernel_1loop_13_bG2_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "b2_f_f") {
        return Power_Spectrum_Kernel_1loop_13_b2_f_f(kvec1, kvec2, kvec3, los, f);
    } else if (parameters == "bG2_f_f") {
        return Power_Spectrum_Kernel_1loop_13_bG2_f_f(kvec1, kvec2, kvec3, los, f);
    } else {
        return 0.0;
    }
}

int integrand_P_Kernel_1loop_13(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ELL, 
		            double f, double Sigma2, double dSigma2, char * parameters) {

	/********/
	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];
	/********/

	for(int j = 0; j < num_k_bin; j++) {
	
	    	/********/
	    	double kvec[3] = {0.0, 0.0, kbin[j]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

        	double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
		double k = NORM(kvec);
		double p = NORM(pvec);

		double BAOp = f_pk(p) - f_pk_no_wiggle(p);
		double BAOk = f_pk(k) - f_pk_no_wiggle(k);

		double Dp = ExpDamping_Ivanov(pvec, los, f, Sigma2, dSigma2);
		double Dk = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);

		double K = Power_Spectrum_Kernel_1loop_13(pvec, M_pvec, kvec, los, f, parameters);

		double powerspec = 6.0 * K * (f_pk_no_wiggle(k) + BAOk * Dk) * (f_pk_no_wiggle(p) + BAOp * Dp);
	    
		/********/
		double result = Nlll * L * powerspec;
		/********/
		    
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);

		ff_out[j] = result * jacobian;
	
	}
	return 0;
}

int integrand_P_Tree_NoWiggle(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
         
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_NoWiggle(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}
int integrand_P_Tree_BAO(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_BAO(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_BAO_Fitting(
        double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	    double sigma2_perp, double sigma2_para, 
        double A20, double A11, double A02, 
        double A30, double A21, double A12, double A03,
        double A40, double A31, double A22, double A13, double A04) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        double kvec[3] = {0.0, 0.0, kbin[i]};
        double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        double L = LegendrePolynomials(ELL, mu);
        double Nlll = (2.0*double(ELL)+1.0);
        
        double result = Nlll * L * Powerspectrum_Tree_BAO_Fitting(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para, A20, A11, A02, A30, A21, A12, A03, A40, A31, A22, A13, A04);
        double jacobian = 1.0;
        
        ff_out[i] = result * jacobian;
	}

	return 0;
}


int integrand_P_Tree_Window(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, 
	                    double sigma8, double fz, double b1, double volume) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];   

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_Window(kvec, pvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, volume);
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_Window_IC(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, 
	                    double sigma8, double fz, double b1, double volume) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];   

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_Window_IC(kvec, pvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, volume);
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_Window_IC_Approx(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, 
	                    double sigma8, double fz, double b1, double volume) {



	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];   

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
       		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
  
        	double L = 1.0;
		if(ELL != 0 ) {
		    L = 0.0;
		}
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_Window_IC_Approx(kvec, pvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, volume);
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        	ff_out[i] = result * jacobian;
	}
	return 0;

}

int integrand_P_Tree_KSZ(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, 
	                 double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, double aH_tau_T0_over_c) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
               	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
 
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_KSZ(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, aH_tau_T0_over_c);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_Aniso_real(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, 
	                   int ell1, int ell2, int _L_, int _M_, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1, 
	                   double g20, double g21_real, double g21_imag, double g22_real, double g22_imag) {

	double mu1  = -1.0 + 2.0 * xx_in[0];
	double phi1  = 0.0 + 2.0 * M_PI * xx_in[1];
	double mu2  = -1.0 + 2.0 * xx_in[2];
	double phi2  = 0.0 + 2.0 * M_PI * xx_in[3];

	for(int i = 0; i < num_kbin; i++) {
		double kvec[3] = {kbin[i] * sqrt(1.0 - mu1 * mu1) * cos(phi1), kbin[i] * sqrt(1.0 - mu1 * mu1) * sin(phi1), kbin[i] * mu1};
		double los[3] = {sqrt(1.0 - mu2 * mu2) * cos(phi2), sqrt(1.0 - mu2 * mu2) * sin(phi2), mu2};
        
		std::complex<double> Sll = calcBipoSH(mu1, phi1, mu2, phi2, ell1, ell2, _L_, _M_);
        	double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
        	
		std::complex<double> result = Nlll * Sll * Powerspectrum_Tree_Aniso(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, g20, g21_real, g21_imag, g22_real, g22_imag);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result.real() * jacobian;
	}

	return 0;
}


int integrand_P_Tree_Aniso_imag(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, 
	                        int ell1, int ell2, int _L_, int _M_, double alpha_perp, double alpha_parallel, 
				double sigma8, double fz, double b1, 
	                        double g20, double g21_real, double g21_imag, double g22_real, double g22_imag) {

	double mu1  = -1.0 + 2.0 * xx_in[0];
	double phi1  = 0.0 + 2.0 * M_PI * xx_in[1];
	double mu2  = -1.0 + 2.0 * xx_in[2];
	double phi2  = 0.0 + 2.0 * M_PI * xx_in[3];

	for(int i = 0; i < num_kbin; i++) {
		double kvec[3] = {kbin[i] * sqrt(1.0 - mu1 * mu1) * cos(phi1), kbin[i] * sqrt(1.0 - mu1 * mu1) * sin(phi1), kbin[i] * mu1};
		double los[3] = {sqrt(1.0 - mu2 * mu2) * cos(phi2), sqrt(1.0 - mu2 * mu2) * sin(phi2), mu2};
        
		std::complex<double> Sll = calcBipoSH(mu1, phi1, mu2, phi2, ell1, ell2, _L_, _M_);
        	double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0);
        	
		std::complex<double> result = Nlll * Sll * Powerspectrum_Tree_Aniso(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, g20, g21_real, g21_imag, g22_real, g22_imag);

        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result.imag() * jacobian;
	}
	return 0;
}

//int integrand_P_Tree_BipoSH(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ell1, int ell2, int ELL, int M, 
//	                    double sigma8, double fz, double alpha_perp, double alpha_parallel, double b1) {
//
//	double mu  = -1.0 + 2.0 * xx_in[0];
//	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
//	phi = phi;
//
//	for(int i = 0; i < num_kbin; i++) {
//        	double kvec[3] = {kbin[i] * sqrt(1.0 - mu * mu), 0.0, kbin[i] * mu};
//        	double los[3] = {0.0, 0.0, 1.0};
//        
//        	double L = calcBipoSH(mu1, phi1, mu2, phi2, ell1, ell2, ELL, M);
//
//		   
//		    LegendrePolynomials(ELL, mu);
//        	double Nlll = (2.0*double(ELL)+1.0);
//        	
//        	double result = Nlll * L * Powerspectrum_Tree(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1);
//        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
//        
//        	ff_out[i] = result * jacobian;
//	}
//
//	return 0;
//}
//

int integrand_P_Tree_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_b1_b1(kvec, los, alpha_perp, alpha_parallel);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}


int integrand_P_Tree_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_b1_f(kvec, los, alpha_perp, alpha_parallel);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}

int integrand_P_Tree_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_f_f(kvec, los, alpha_perp, alpha_parallel);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}


int integrand_P_Tree_BAO_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_BAO_b1_b1(kvec, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}

int integrand_P_Tree_BAO_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_BAO_b1_f(kvec, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}

int integrand_P_Tree_BAO_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_BAO_f_f(kvec, los, alpha_perp, alpha_parallel, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}


int integrand_P_Tree_NoWiggle_b1_b1(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_NoWiggle_b1_b1(kvec, los, alpha_perp, alpha_parallel);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}


int integrand_P_Tree_NoWiggle_b1_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_NoWiggle_b1_f(kvec, los, alpha_perp, alpha_parallel);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}

int integrand_P_Tree_NoWiggle_f_f(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_NoWiggle_f_f(kvec, los, alpha_perp, alpha_parallel);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}
	return 0;
}

int integrand_P_NonLinearFitting(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_NonLinearFitting(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_NonLinearFitting_Window(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para, double volume) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double pmag = exp( log(kmin) + dlnk * xx_in[1] );
    double mu_p = -1.0 + 2.0 * xx_in[2];
    double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];   

	for(int i = 0; i < num_kbin; i++) {
        double kvec[3] = {0.0, 0.0, kbin[i]};
        double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        
        double L = LegendrePolynomials(ELL, mu);
        double Nlll = (2.0*double(ELL)+1.0);
        
        double result = Nlll * L * Powerspectrum_NonLinearFitting_Window(kvec, pvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para, volume);
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_SigmaB(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	               double sigma2_perp, double sigma2_para, double nmean, double volume) {

    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double kmag = exp( log(kmin) + dlnk * xx_in[0] );
    double mu = -1.0 + 2.0 * xx_in[1];
    double phi = 0.0;   
    phi = phi;
    
    double kvec[3] = {0.0, 0.0, kmag};
    double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double result = W2(kmag, R) * ( Powerspectrum_NonLinearFitting(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para) + 1.0/nmean);
    double jacobian = dlnk * pow(kmag,3) / (2.0 * M_PI * M_PI);
        
    ff_out[0] = result * jacobian;
    
    return 0;
}

int integrand_P_sigma2_perp_Reconstructed(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, 
        double sigma8, double fz, double b1,
        double b1_fid, double R) {

    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double kmag = exp( log(kmin) + dlnk * xx_in[0] );
    double mu = -1.0 + 2.0 * xx_in[1];
    mu = mu;
    double kmag2 = kmag * kmag;
    
    double W = exp( - pow(kmag * R, 2) / 2.0 );
    double result = (1.0/3.0) * (pow(sigma8, 2) * f_pk(kmag) / kmag2) * 
                   ( ( 1.0 - (W*b1 / b1_fid) ) * ( 1.0 - (W*b1 / b1_fid) )
                   + ( 1.0 - (W*b1 / b1_fid) ) * (-W*fz/b1_fid) * (2.0/5.0)
                   + (-W*fz/b1_fid) * (-W*fz/b1_fid) * (3.0/35.0) );
    
    double jacobian = dlnk * pow(kmag,3) / (2.0 * M_PI * M_PI);
    ff_out[0] = result * jacobian;
    
    return 0;
}

int integrand_P_sigma2_para_Reconstructed(
        double * xx_in, int ndim, double * ff_out, int ncomp, 
        double * kbin, int num_kbin, 
        double sigma8, double fz, double b1,
        double b1_fid, double R) {

    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double kmag = exp( log(kmin) + dlnk * xx_in[0] );
    double mu = -1.0 + 2.0 * xx_in[1];
    mu = mu;
    double kmag2 = kmag * kmag;
    
    double W = exp( - pow(kmag * R, 2) / 2.0 );
    double result = (1.0/3.0) * (pow(sigma8, 2) * f_pk(kmag) / kmag2) * 
                   ( ( ( 1.0 - (W*b1 / b1_fid) ) + fz ) * ( ( 1.0 - (W*b1 / b1_fid) ) + fz )
                   + ( 1.0 - (W*b1 / b1_fid) ) * (-W*fz/b1_fid) * (6.0/5.0)
                   + (-W*fz/b1_fid) * (6.0 * fz / 5.0)
                   + (-W*fz/b1_fid) * (-W*fz/b1_fid) * (3.0/7.0)
                   );

    double jacobian = dlnk * pow(kmag,3) / (2.0 * M_PI * M_PI);
    ff_out[0] = result * jacobian;
    
    return 0;

}

int integrand_P_Tree_A_perp(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_A_perp(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_A_para(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_A_para(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_B_perp(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_B_perp(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_B_para(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_B_para(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_Tree_B_perp_para(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, double sigma8, double fz, double b1,
	                 double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];
	double phi  = 0.0 + 2.0 * M_PI * xx_in[1];
	phi = phi;

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_Tree_B_perp_para(kvec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, sigma2_perp, sigma2_para);
        	double jacobian = 1.0; //(2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
        
        	ff_out[i] = result * jacobian;
	}

	return 0;
}


int integrand_P_SPT1loop(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin,
                      int ELL, double alpha_perp, double alpha_parallel, 
	              double sigma8, double fz, double b1, double b2, double b3, double bK2,
		      double bK3, double bDK, double bO, double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p = -1.0 + 2.0 * xx_in[2];
        double phi_p = 0.0 + 2.0 * M_PI * xx_in[3];   

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p) * cos(phi_p), pmag * sqrt(1.0 - mu_p * mu_p) * sin(phi_p), pmag * mu_p};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_SPT1loop(kvec, pvec, los, alpha_perp, alpha_parallel,
		         	sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, sigma2_perp, sigma2_para);
		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
        	ff_out[i] = result * jacobian;
	}

//	double mu_p  = -1.0 + 2.0 * xx_in[0];
//
//    	double kmax = 5.0;
//        double kmin = 1.0e-4;
//        double dlnk = log(kmax) - log(kmin);
//        double pmag = exp( log(kmin) + dlnk * xx_in[1] );
//        double mu = -1.0 + 2.0 * xx_in[2];
//        double phi = 0.0 + 2.0 * M_PI * xx_in[3];   
//
//	for(int i = 0; i < num_kbin; i++) {
//
//		double kvec[3] = {0.0, 0.0, kbin[i]};
//        
//        	double los[3]  = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
//
//		double pvec[3] = {pmag * sqrt(1.0 - mu_p * mu_p), 0.0, pmag * mu_p};
//        
//        	double L = LegendrePolynomials(ELL, mu);
//        	double Nlll = (2.0*double(ELL)+1.0);
//        	
//        	double result = Nlll * L * Powerspectrum_SPT1loop(kvec, pvec, los, alpha_perp, alpha_parallel,
//		         	sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, sigma2_perp, sigma2_para);
//		double jacobian = dlnk * pow(pmag,3) / (2.0 * M_PI * M_PI);
//        	ff_out[i] = result * jacobian;
//	}

	return 0;
}

int integrand_P_SPT2loop(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin,
                      int ELL, double alpha_perp, double alpha_parallel, 
	              double sigma8, double fz, double b1, double b2, double b3, double bK2,
		      double bK3, double bDK, double bO, double sigma2_perp, double sigma2_para) {

	double mu  = -1.0 + 2.0 * xx_in[0];

    	double kmax = 5.0;
        double kmin = 1.0e-4;
        double dlnk = log(kmax) - log(kmin);
        double pmag1 = exp( log(kmin) + dlnk * xx_in[1] );
        double mu_p1 = -1.0 + 2.0 * xx_in[2];
        double phi_p1 = 0.0 + 2.0 * M_PI * xx_in[3];   
        double pmag2 = exp( log(kmin) + dlnk * xx_in[4] );
        double mu_p2 = -1.0 + 2.0 * xx_in[5];
        double phi_p2 = 0.0 + 2.0 * M_PI * xx_in[6];   

	for(int i = 0; i < num_kbin; i++) {
        	double kvec[3] = {0.0, 0.0, kbin[i]};
        	double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};

		double pvec1[3] = {pmag1 * sqrt(1.0 - mu_p1 * mu_p1) * cos(phi_p1), pmag1 * sqrt(1.0 - mu_p1 * mu_p1) * sin(phi_p1), pmag1 * mu_p1};
		double pvec2[3] = {pmag2 * sqrt(1.0 - mu_p2 * mu_p2) * cos(phi_p2), pmag2 * sqrt(1.0 - mu_p2 * mu_p2) * sin(phi_p2), pmag2 * mu_p2};
        
        	double L = LegendrePolynomials(ELL, mu);
        	double Nlll = (2.0*double(ELL)+1.0);
        	
        	double result = Nlll * L * Powerspectrum_SPT2loop(kvec, pvec1, pvec2, los, alpha_perp, alpha_parallel,
		         	sigma8, fz, b1, b2, b3, bK2, bK3, bDK, bO, sigma2_perp, sigma2_para);
		double jacobian = (dlnk * pow(pmag1,3) / (2.0 * M_PI * M_PI)) * (dlnk * pow(pmag2,3) / (2.0 * M_PI * M_PI));
        	ff_out[i] = result * jacobian;
	}

	return 0;
}

int integrand_P_LocalMean(double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_kbin, int ELL, double alpha_perp, double alpha_parallel, 
	                      double sigma8, double fz, double b1, double b2, double bK2, double sigma2_perp, double sigma2_para, double nmean, double volume) {

    double mu  = -1.0 + 2.0 * xx_in[0];
    
    double kmax = 5.0;
    double kmin = 1.0e-4;
    double dlnk = log(kmax) - log(kmin);
    double emag = exp( log(kmin) + dlnk * xx_in[1] );
    double mu_e = -1.0 + 2.0 * xx_in[2];
    double phi_e = 0.0 + 2.0 * M_PI * xx_in[3];   
    
	for(int i = 0; i < num_kbin; i++) {

        double kvec[3] = {0.0, 0.0, kbin[i]};
        double los[3]  = {sqrt(1.0 - mu * mu), 0.0, mu};
        
        double evec[3] = {emag * sqrt(1.0 - mu_e * mu_e) * cos(phi_e), emag * sqrt(1.0 - mu_e * mu_e) * sin(phi_e), emag * mu_e};
                
        double L = LegendrePolynomials(ELL, mu);
        double Nlll = (2.0*double(ELL)+1.0);
        
        double result = Nlll * L * Powerspectrum_LocalMean(kvec, evec, los, alpha_perp, alpha_parallel, sigma8, fz, b1, b2, bK2, sigma2_perp, sigma2_para, nmean, volume);
        double jacobian = dlnk * pow(emag,3) / (2.0 * M_PI * M_PI);
        ff_out[i] = result * jacobian;

	}

	return 0;
}



#endif

