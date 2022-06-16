#ifndef __precompute_kernel__
#define __precompute_kernel__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __pk_lin__
#include "pk_lin.hpp"
#endif

#ifndef __sph__
#include "sph.hpp"
#endif

#ifndef __kernel__
#include "kernel.hpp"
#endif

/********************************/
/*        Basic kernels         */
/********************************/



double F2t(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	return 5.0/7.0 + mu / 2.0 * (k1/k2 + k2/k1) + (2.0/7.0) * mu * mu;
}

double G2t(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	return 3.0/7.0 + mu / 2.0 * (k1/k2 + k2/k1) + (4.0/7.0) * mu * mu;
}

double U1() {
	return 1.0;
}

double V1t(double * kvec1, double * los) {
	double mu = MU(kvec1, los);
	return mu * mu;
}

double KFOGt(double * kvec1, double * los) {
	double mu = MU(kvec1, los);
	double k1 = NORM(kvec1);
	return - mu * mu * k1 * k1 / 0.3 / 0.3; // k_nl^r = O.3 h/Mpc following Ivanov et al.
}

double U2() {
	return 1.0 / 2.0;
}

double K2(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
	return (mu * mu - 1.0/3.0);
}

double V2t(double * kvec1, double * kvec2, double * los) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double mu = MU(kvec12, los);
	return mu * mu * G2(kvec1, kvec2);
}

double DVt(double * kvec1, double * kvec2, double * los) {
	
	double kvec12[3] = PLUS(kvec1, kvec2);
	double mu = MU(kvec12, los);
	double kn = NORM(kvec12);

	double mu1 = MU(kvec1, los);
	double k1 = NORM(kvec1);

	double mu2 = MU(kvec2, los);
	double k2 = NORM(kvec2);
	
	return mu * kn / 2.0 * (mu1/k1 + mu2/k2);
}

double V11t(double * kvec1, double * kvec2, double * los) {

	double kvec12[3] = PLUS(kvec1, kvec2);
	double mu = MU(kvec12, los);
	double kn = NORM(kvec12);

	double mu1 = MU(kvec1, los);
	double k1 = NORM(kvec1);

	double mu2 = MU(kvec2, los);
	double k2 = NORM(kvec2);

	return mu * kn / 2.0 * (mu1/k1 * mu2 * mu2 + mu2/k2 * mu1 * mu1);
}

/********************************/
/*     Pre-computed kernels     */
/********************************/

double b1_3(double * kvec1, double * kvec2) {
    return F2t(kvec1, kvec2);
}

double b1_2_b2() {
    return U2();
}

double b1_2_bk2(double * kvec1, double * kvec2) {
    return K2(kvec1, kvec2);
}

double b1_2_f(double * kvec1, double * kvec2, double * los) {
    return V2t(kvec1, kvec2, los) + (V1t(kvec1, los) + V1t(kvec2, los)) * F2t(kvec1, kvec2);
}

double b1_3_f(double * kvec1, double * kvec2, double * los) {
    return DVt(kvec1, kvec2, los);
}

double b1_2_f_2(double * kvec1, double * kvec2, double * los) {
    return V11t(kvec1, kvec2, los) + (V1t(kvec1, los) + V1t(kvec2, los)) * DVt(kvec1, kvec2, los);
}





double b1_b2_f(double * kvec1, double * kvec2, double * los) {
    return U2() * (V1t(kvec1, los) + V1t(kvec2, los));
}

double b1_bk2_f(double * kvec1, double * kvec2, double * los) {
    return K2(kvec1, kvec2) * (V1t(kvec1, los) + V1t(kvec2, los));
}

double b1_f_2(double * kvec1, double * kvec2, double * los) {
    return (V1t(kvec1, los) + V1t(kvec2, los)) * V2t(kvec1, kvec2, los) + V1t(kvec1, los) * V1t(kvec2, los) * F2t(kvec1, kvec2);
}

double b1_f_3(double * kvec1, double * kvec2, double * los) {
    return (V1t(kvec1, los) + V1t(kvec2, los)) * V11(kvec1, kvec2, los) + V1t(kvec1, los) * V1t(kvec2, los) * DVt(kvec1, kvec2, los);
}





double b2_f_2(double * kvec1, double * kvec2, double * los) {
    return U2() * V1t(kvec1, los) * V1t(kvec2, los);
}

double bk2_f_2(double * kvec1, double * kvec2, double * los) {
    return K2(kvec1, kvec2) * V1t(kvec1, los) * V1t(kvec2, los);
}

double f_3(double * kvec1, double * kvec2, double * los) {
    return V1t(kvec1, los) * V1t(kvec2, los) * V2t(kvec1, kvec2, los);
}

double f_4(double * kvec1, double * kvec2, double * los) {
    return V1t(kvec1, los) * V1t(kvec2, los) * V11t(kvec1, kvec2, los);
}

/********************************/
/* Kernels for FoG counterterm  */
/********************************/

double c1_b1_2(double * kvec1, double * kvec2, double * los) {
    return (KFOGt(kvec1, los) + KFOGt(kvec2, los)) * F2t(kvec1, kvec2);
}

double c1_b1_b2(double * kvec1, double * kvec2, double * los) {
    return (KFOGt(kvec1, los) + KFOGt(kvec2, los)) * U2();
}

double c1_b1_bk2(double * kvec1, double * kvec2, double * los) {
    return (KFOGt(kvec1, los) + KFOGt(kvec2, los)) * K2(kvec1, kvec2);
}

double c1_b1_f(double * kvec1, double * kvec2, double * los) {
    return (KFOGt(kvec1, los) + KFOGt(kvec2, los)) * V2t(kvec1, kvec2, los) + KFOGt(kvec1, los) * V1t(kvec2, los) * F2t(kvec1, kvec2);
}

double c1_b1_2_f(double * kvec1, double * kvec2, double * los) {
    return (KFOGt(kvec1, los) + KFOGt(kvec2, los)) * DVt(kvec1, kvec2, los);
}

double c1_b1_f_2(double * kvec1, double * kvec2, double * los) {
    return (KFOGt(kvec1, los) + KFOGt(kvec2, los)) * V11t(kvec1, kvec2, los) + KFOGt(kvec1, los) * V1t(kvec2, los) * DVt(kvec1, kvec2, los);
}





double c1_b2_f(double * kvec1, double * kvec2, double * los) {
    return U2() * KFOGt(kvec1, los) * V1t(kvec2, los);
}

double c1_bk2_f(double * kvec1, double * kvec2, double * los) {
    return K2(kvec1, kvec2) * KFOGt(kvec1, los) * V1t(kvec2, los);
}

double c1_f_2(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * V1t(kvec2, los) * V2t(kvec1, kvec2, los);
}

double c1_f_3(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * V1t(kvec2, los) * V11t(kvec1, kvec2, los);
}





double c1_2_b1(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * KFOGt(kvec2, los) * F2t(kvec1, kvec2);
}

double c1_2_b2(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * KFOGt(kvec2, los) * U2();
}

double c1_2_bk2(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * KFOGt(kvec2, los) * K2(kvec1, kvec2);
}

double c1_2_f(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * KFOGt(kvec2, los) * V2t(kvec1, kvec2, los);
}

double c1_2_b1_f(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * KFOGt(kvec2, los) * DVt(kvec1, kvec2, los);
}

double c1_2_f_2(double * kvec1, double * kvec2, double * los) {
    return KFOGt(kvec1, los) * KFOGt(kvec2, los) * V11t(kvec1, kvec2, los);
}


double KernelTemplate(double * kvec1, double * kvec2, double * los, char * parameters_in) {
    
    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "b1_3") {
        return b1_3(kvec1, kvec2);
    } else if (parameters == "b1_2_b2") {
        return b1_2_b2();
    } else if (parameters == "b1_2_bk2") {
        return b1_2_bk2(kvec1, kvec2);
    } else if (parameters == "b1_2_f") {
        return b1_2_f(kvec1, kvec2, los);
    } else if (parameters == "b1_3_f") {
        return b1_3_f(kvec1, kvec2, los);
    } else if (parameters == "b1_2_f_2") {
        return b1_2_f_2(kvec1, kvec2, los);
    } else if (parameters == "b1_b2_f") {
        return b1_b2_f(kvec1, kvec2, los);
    } else if (parameters == "b1_bk2_f") {
        return b1_bk2_f(kvec1, kvec2, los);
    } else if (parameters == "b1_f_2") {
        return b1_f_2(kvec1, kvec2, los);
    } else if (parameters == "b1_f_3") {
        return b1_f_3(kvec1, kvec2, los);
    } else if (parameters == "b2_f_2") {
        return b2_f_2(kvec1, kvec2, los);
    } else if (parameters == "bk2_f_2") {
        return bk2_f_2(kvec1, kvec2, los);
    } else if (parameters == "f_3") {
        return f_3(kvec1, kvec2, los);
    } else if (parameters == "f_4") {
        return f_4(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_2") {
        return c1_b1_2(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_b2") {
        return c1_b1_b2(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_bk2") {
        return c1_b1_bk2(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_f") {
        return c1_b1_f(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_2_f") {
        return c1_b1_2_f(kvec1, kvec2, los);
    } else if (parameters == "c1_b1_f_2") {
        return c1_b1_f_2(kvec1, kvec2, los);
    } else if (parameters == "c1_b2_f") {
        return c1_b2_f(kvec1, kvec2, los);
    } else if (parameters == "c1_bk2_f") {
        return c1_bk2_f(kvec1, kvec2, los);
    } else if (parameters == "c1_f_2") {
        return c1_f_2(kvec1, kvec2, los);
    } else if (parameters == "c1_f_3") {
        return c1_f_3(kvec1, kvec2, los);
    } else if (parameters == "c1_2_b1") {
        return c1_2_b1(kvec1, kvec2, los);
    } else if (parameters == "c1_2_b2") {
        return c1_2_b2(kvec1, kvec2, los);
    } else if (parameters == "c1_2_bk2") {
        return c1_2_bk2(kvec1, kvec2, los);
    } else if (parameters == "c1_2_f") {
        return c1_2_f(kvec1, kvec2, los);
    } else if (parameters == "c1_2_b1_f") {
        return c1_2_b1_f(kvec1, kvec2, los);
    } else if (parameters == "c1_2_f_2") {
        return c1_2_f_2(kvec1, kvec2, los);
    } else {
        return 0.0;
    }
}





int integrand_KernelTemplate(
        double * xx_in, int ndim, double * ff_out, int ncomp, double * kbin, int num_k_bin, int ell1, int ell2, int ELL, double kmag1,
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
	    // double kvec3[3] = {- kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2]};
	    double los[3]   = {sqrt(1.0 - mu * mu) * cos(phi), sqrt(1.0 - mu * mu) * sin(phi), mu};
	    /********/
	    double result = KernelTemplate(kvec1, kvec2, los, parameters);
            // double result_1 = KernelTemplate(kvec1, kvec3, los, parameters);
            // double result_2 = KernelTemplate(kvec2, kvec3, los, parameters);
	    /********/

	    /********/
	    double Nlll = (2.0*double(ell1)+1.0) * (2.0*double(ell2)+1.0) * (2.0*double(ELL)+1.0);
	    double Slll = calcYYY(mu1, phi1, mu2, phi2, mu, phi, ell1, ell2, ELL);
	    /********/

	    // double jacobian = (2.0 / 2.0) * (2.0 / 2.0) * (2.0 * M_PI / (2.0 * M_PI));
	    ff_out[j] = Nlll * Slll * result; //* jacobian;
	    // ff_out_1[j] = Slll * result_1;
	    // ff_out_2[j] = Slll * result_2;

	}

	return 0;
}

#endif
