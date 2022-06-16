#ifndef __kernel__
#define __kernel__

#ifndef __common__
#include "common.hpp"
#endif

#ifndef __pk_lin__
#include "pk_lin.hpp"
#endif

#ifndef __sph__
#include "sph.hpp"
#endif

using namespace std;

double MU(double * kvec1, double * kvec2) {
 	double k1 = NORM(kvec1);
 	double k2 = NORM(kvec2);
 	double result = 0.0;
 	if( (k1 > pk_kmin) && (k2 > pk_kmin)) {
 		double m = DOT(kvec1, kvec2) / k1 / k2;
 		if(m > 1.0) {
 			m = 1.0;
 		} else if (m < -1.0) {
 			m = - 1.0;
 		}
 		result += m;
 	}
 	return result;
}

double W1(double kmag, double R) {

	double kr = R * kmag;
	double kr2 = kr * kr;
	double kr3 = kr2 * kr;
  	double w = 3.0 * (sin(kr) / kr3 - cos(kr) / kr2);
	return w;
}


double W2(double kmag, double R) {

	double kr = R * kmag;
	double kr2 = kr * kr;
	double kr3 = kr2 * kr;
  	double w = 3.0 * (sin(kr) / kr3 - cos(kr) / kr2);
	double w2 = w * w;

	return w2;
}

int calcTrueWavevector(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double * kvec_out) {

	double kn = DOT(kvec_in, los);

	double kvec_para[3] = { kn * los[0] / alpha_parallel, kn * los[1] / alpha_parallel, kn * los[2] / alpha_parallel };
	double kvec_perp[3] = { (kvec_in[0] - kn * los[0]) / alpha_perp, (kvec_in[1] - kn * los[1]) / alpha_perp, (kvec_in[2] - kn * los[2]) / alpha_perp };

	double kvec[3] = PLUS(kvec_para, kvec_perp);
	
	kvec_out[0] = kvec[0];
	kvec_out[1] = kvec[1];
	kvec_out[2] = kvec[2];

	return 0;

}

int calcTrueRelativeDistance(double * rvec_in, double * los, double alpha_perp, double alpha_parallel, double * rvec_out) {

	double rn = DOT(rvec_in, los);

	double rvec_para[3] = { rn * los[0] * alpha_parallel, rn * los[1] * alpha_parallel, rn * los[2] * alpha_parallel };
	double rvec_perp[3] = { (rvec_in[0] - rn * los[0]) * alpha_perp, (rvec_in[1] - rn * los[1]) * alpha_perp, (rvec_in[2] - rn * los[2]) * alpha_perp };

	double rvec[3] = PLUS(rvec_para, rvec_perp);
	
	rvec_out[0] = rvec[0];
	rvec_out[1] = rvec[1];
	rvec_out[2] = rvec[2];

	return 0;

}


int calcTrueWavevectorHat(double * kvec_hat_in, double * los, double epsilon, double * kvec_hat_out) {

	double kn = DOT(kvec_hat_in, los);

	double factor2 = 1.0 + (kn*kn) * ( pow(1.0 + epsilon, - 6.0) - 1.0 );
	double factor = pow(factor2, 1.0/2.0);

	kvec_hat_out[0] = ( kvec_hat_in[0] + kn * los[0] * ( pow(1.0 + epsilon, - 3.0) - 1.0 ) ) / factor;
	kvec_hat_out[1] = ( kvec_hat_in[1] + kn * los[1] * ( pow(1.0 + epsilon, - 3.0) - 1.0 ) ) / factor;
	kvec_hat_out[2] = ( kvec_hat_in[2] + kn * los[2] * ( pow(1.0 + epsilon, - 3.0) - 1.0 ) ) / factor;

	return 0;

}

int calcTrueRelativeDistanceUnitVector(double * rvec_hat_in, double * los, double epsilon, double * rvec_hat_out) {

	double rn = DOT(rvec_hat_in, los);

	double factor2 = 1.0 + (rn*rn) * ( pow(1.0 + epsilon, 6.0) - 1.0 );
	double factor = pow(factor2, 1.0/2.0);

	rvec_hat_out[0] = ( rvec_hat_in[0] + rn * los[0] * ( pow(1.0 + epsilon, 3.0) - 1.0 ) ) / factor;
	rvec_hat_out[1] = ( rvec_hat_in[1] + rn * los[1] * ( pow(1.0 + epsilon, 3.0) - 1.0 ) ) / factor;
	rvec_hat_out[2] = ( rvec_hat_in[2] + rn * los[2] * ( pow(1.0 + epsilon, 3.0) - 1.0 ) ) / factor;

	return 0;

}



/***********************************************/
/* Kernel functions for dark matter: Fn and Gn */
/***********************************************/
/***********************************************/
double alpha(double * kvec1, double * kvec2) {
	double k1 = NORM(kvec1);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = 0.0;
	if( (k1 > pk_kmin) ) {
		result += DOT(kvec12, kvec1) / k1 / k1;
	}
	
	return result;
}

double beta(double * kvec1, double * kvec2) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = 0.0;
	if( (k1 > pk_kmin) && (k2 > pk_kmin)) {
		result += DOT(kvec12, kvec12) * DOT(kvec1, kvec2) / 2.0 / pow(k1,2) / pow(k2,2);
	}
	return result;
}

/***********************************************/

double F2(double * kvec1, double * kvec2) {
	return (5.0/14.0) * alpha(kvec1, kvec2) + (5.0/14.0) * alpha(kvec2, kvec1) + (2.0/7.0) * beta(kvec1, kvec2);
} 

double G2(double * kvec1, double * kvec2) {
	return (3.0/14.0) * alpha(kvec1, kvec2) + (3.0/14.0) * alpha(kvec2, kvec1) + (4.0/7.0) * beta(kvec1, kvec2);
}

double F2_Growth(double * kvec1, double * kvec2) {
	return (17.0/21.0);
} 

double F2_Shift(double * kvec1, double * kvec2) {
    double mu = MU(kvec1, kvec2);
    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
	return (1.0/2.0) * mu * (k1/k2 + k2/k1);
} 

double F2_Tidal(double * kvec1, double * kvec2) {
    double mu = MU(kvec1, kvec2);
	return (2.0/7.0) * (mu * mu - (1.0/3.0));
} 


/***********************************************/

double F3_temp(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double F = (7.0/18.0) * alpha(kvec1, kvec23) * F2(kvec2, kvec3) + (1.0/9.0) * beta(kvec1, kvec23) * G2(kvec2, kvec3)
             + (7.0/18.0) * alpha(kvec12, kvec3) * G2(kvec1, kvec2) + (1.0/9.0) * beta(kvec12, kvec3) * G2(kvec1, kvec2); 
	return F;
} 

double F3(double * kvec1, double * kvec2, double * kvec3) {
	return ( F3_temp(kvec1, kvec2, kvec3) + F3_temp(kvec2, kvec3, kvec1) + F3_temp(kvec3, kvec1, kvec2) ) / 3.0;
} 

double G3_temp(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double G = (1.0/6.0) * alpha(kvec1, kvec23) * F2(kvec2, kvec3) + (1.0/3.0) * beta(kvec1, kvec23) * G2(kvec2, kvec3)
                 + (1.0/6.0) * alpha(kvec12, kvec3) * G2(kvec1, kvec2) + (1.0/3.0) * beta(kvec12, kvec3) * G2(kvec1, kvec2); 
	return G;
} 

double G3(double * kvec1, double * kvec2, double * kvec3) {
	return ( G3_temp(kvec1, kvec2, kvec3) + G3_temp(kvec2, kvec3, kvec1) + G3_temp(kvec3, kvec1, kvec2) ) / 3.0;
}


double F4_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec234[3] = PLUS3(kvec2, kvec3, kvec4);

	double F = (9.0) * alpha(kvec1, kvec234) * F3(kvec2, kvec3, kvec4) 
		 + (2.0) *  beta(kvec1, kvec234) * G3(kvec2, kvec3, kvec4)
		 + (9.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec4)
		 + (2.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec4);

	return F / 33.0;
} 

double F4_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec34[3] = PLUS(kvec3, kvec4);

	double F = (9.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec34) * F2(kvec3, kvec4) 
		 + (2.0) * G2(kvec1, kvec2) * beta(kvec12, kvec34) * G2(kvec3, kvec4);

	return F / 33.0;
} 

double F4(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {

	double Ftemp1 = F4_temp_1(kvec1, kvec2, kvec3, kvec4) 
	              + F4_temp_1(kvec2, kvec3, kvec4, kvec1) 
	              + F4_temp_1(kvec3, kvec4, kvec1, kvec2) 
	              + F4_temp_1(kvec4, kvec1, kvec2, kvec3);

	double Ftemp2 = F4_temp_2(kvec1, kvec2, kvec3, kvec4)
		      + F4_temp_2(kvec1, kvec3, kvec2, kvec4)
		      + F4_temp_2(kvec1, kvec4, kvec2, kvec3)
		      + F4_temp_2(kvec2, kvec3, kvec1, kvec4)
		      + F4_temp_2(kvec2, kvec4, kvec1, kvec3)
		      + F4_temp_2(kvec3, kvec4, kvec1, kvec2);

	return Ftemp1 / 4.0 + Ftemp2 / 6.0;
}


double G4_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec234[3] = PLUS3(kvec2, kvec3, kvec4);

	double G = (3.0) * alpha(kvec1, kvec234) * F3(kvec2, kvec3, kvec4) 
                 + (8.0) *  beta(kvec1, kvec234) * G3(kvec2, kvec3, kvec4)
		 + (3.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec4)
		 + (8.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec4);

	return G / 33.0;
} 

double G4_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec34[3] = PLUS(kvec3, kvec4);

	double G = (3.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec34) * F2(kvec3, kvec4) 
		 + (8.0) * G2(kvec1, kvec2) * beta(kvec12, kvec34) * G2(kvec3, kvec4);

	return G / 33.0;
} 


double G4(double * kvec1, double * kvec2, double * kvec3, double * kvec4) {

	double Gtemp1 = G4_temp_1(kvec1, kvec2, kvec3, kvec4) 
		      + G4_temp_1(kvec2, kvec3, kvec4, kvec1) 
		      + G4_temp_1(kvec3, kvec4, kvec1, kvec2) 
		      + G4_temp_1(kvec4, kvec1, kvec2, kvec3);
	double Gtemp2 = G4_temp_2(kvec1, kvec2, kvec3, kvec4)
		      + G4_temp_2(kvec1, kvec3, kvec2, kvec4)
		      + G4_temp_2(kvec1, kvec4, kvec2, kvec3)
		      + G4_temp_2(kvec2, kvec3, kvec1, kvec4)
		      + G4_temp_2(kvec2, kvec4, kvec1, kvec3)
		      + G4_temp_2(kvec3, kvec4, kvec1, kvec2);

	return Gtemp1 / 4.0 + Gtemp2 / 6.0;
}

double F5_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double kvec2345[3] = PLUS4(kvec2, kvec3, kvec4, kvec5);

	double F = (3.0) * alpha(kvec1, kvec2345) * F4(kvec2, kvec3, kvec4, kvec5)
		 + (10.0) * beta(kvec1, kvec2345) * G4(kvec2, kvec3, kvec4, kvec5)
		 + (3.0) * alpha(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4)
		 + (10.0) * beta(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4);

	return F / 52.0;
} 

double F5_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec345[3] = PLUS3(kvec3, kvec4, kvec5);

	double F = (3.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec345) * F3(kvec3, kvec4, kvec5)
		 + (10.0) * G2(kvec1, kvec2) * beta(kvec12, kvec345) * G3(kvec3, kvec4, kvec5);
	return F / 52.0;
} 

double F5_temp_3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec45[3] = PLUS(kvec4, kvec5);

	double F = (3.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec45) * F2(kvec4, kvec5)
		 + (10.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec45) * G2(kvec4, kvec5);
	return F / 52.0;
} 

double F5(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {

	double Ftemp1 = F5_temp_1(kvec1, kvec2, kvec3, kvec4, kvec5)
		      + F5_temp_1(kvec2, kvec3, kvec4, kvec5, kvec1)
		      + F5_temp_1(kvec3, kvec4, kvec5, kvec1, kvec2)
		      + F5_temp_1(kvec4, kvec5, kvec1, kvec2, kvec3)
		      + F5_temp_1(kvec5, kvec1, kvec2, kvec3, kvec4);

	double Ftemp2 = F5_temp_2(kvec1, kvec2, kvec3, kvec4, kvec5)
	              + F5_temp_2(kvec1, kvec3, kvec2, kvec4, kvec5)
	              + F5_temp_2(kvec1, kvec4, kvec2, kvec3, kvec5)
	              + F5_temp_2(kvec1, kvec5, kvec2, kvec3, kvec4)
	              + F5_temp_2(kvec2, kvec3, kvec1, kvec4, kvec5)
	              + F5_temp_2(kvec2, kvec4, kvec1, kvec3, kvec5)
	              + F5_temp_2(kvec2, kvec5, kvec1, kvec3, kvec4)
	              + F5_temp_2(kvec3, kvec4, kvec1, kvec2, kvec5)
	              + F5_temp_2(kvec3, kvec5, kvec1, kvec2, kvec4)
	              + F5_temp_2(kvec4, kvec5, kvec1, kvec2, kvec3);


	double Ftemp3 = F5_temp_3(kvec3, kvec4, kvec5, kvec1, kvec2)
	              + F5_temp_3(kvec2, kvec4, kvec5, kvec1, kvec3)
	              + F5_temp_3(kvec2, kvec3, kvec5, kvec1, kvec4)
	              + F5_temp_3(kvec2, kvec3, kvec4, kvec1, kvec5)
	              + F5_temp_3(kvec1, kvec4, kvec5, kvec2, kvec3)
	              + F5_temp_3(kvec1, kvec3, kvec5, kvec2, kvec4)
	              + F5_temp_3(kvec1, kvec3, kvec4, kvec2, kvec5)
	              + F5_temp_3(kvec1, kvec2, kvec5, kvec3, kvec4)
	              + F5_temp_3(kvec1, kvec2, kvec4, kvec3, kvec5)
	              + F5_temp_3(kvec1, kvec2, kvec3, kvec4, kvec5);

	return Ftemp1 / 5.0 + Ftemp2 / 10.0 + Ftemp3 / 10.0;
}

double G5_temp_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double kvec2345[3] = PLUS4(kvec2, kvec3, kvec4, kvec5);

	double G = (3.0) * alpha(kvec1, kvec2345) * F4(kvec2, kvec3, kvec4, kvec5)
		 + (10.0) * beta(kvec1, kvec2345) * G4(kvec2, kvec3, kvec4, kvec5)
		 + (3.0) * alpha(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4)
		 + (10.0) * beta(kvec1234, kvec5) * G4(kvec1, kvec2, kvec3, kvec4);

	return G / 52.0;
} 

double G5_temp_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec345[3] = PLUS3(kvec3, kvec4, kvec5);

	double G = (3.0) * G2(kvec1, kvec2) * alpha(kvec12, kvec345) * F3(kvec3, kvec4, kvec5)
		 + (10.0) * G2(kvec1, kvec2) * beta(kvec12, kvec345) * G3(kvec3, kvec4, kvec5);
	return G / 52.0;
} 

double G5_temp_3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double kvec45[3] = PLUS(kvec4, kvec5);

	double G = (3.0) * G3(kvec1, kvec2, kvec3) * alpha(kvec123, kvec45) * F2(kvec4, kvec5)
		 + (10.0) * G3(kvec1, kvec2, kvec3) * beta(kvec123, kvec45) * G2(kvec4, kvec5);
	return G / 52.0;
} 

double G5(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5) {

	double Gtemp1 = G5_temp_1(kvec1, kvec2, kvec3, kvec4, kvec5)
		      + G5_temp_1(kvec2, kvec3, kvec4, kvec5, kvec1)
		      + G5_temp_1(kvec3, kvec4, kvec5, kvec1, kvec2)
		      + G5_temp_1(kvec4, kvec5, kvec1, kvec2, kvec3)
		      + G5_temp_1(kvec5, kvec1, kvec2, kvec3, kvec4);

	double Gtemp2 = G5_temp_2(kvec1, kvec2, kvec3, kvec4, kvec5)
	              + G5_temp_2(kvec1, kvec3, kvec2, kvec4, kvec5)
	              + G5_temp_2(kvec1, kvec4, kvec2, kvec3, kvec5)
	              + G5_temp_2(kvec1, kvec5, kvec2, kvec3, kvec4)
	              + G5_temp_2(kvec2, kvec3, kvec1, kvec4, kvec5)
	              + G5_temp_2(kvec2, kvec4, kvec1, kvec3, kvec5)
	              + G5_temp_2(kvec2, kvec5, kvec1, kvec3, kvec4)
	              + G5_temp_2(kvec3, kvec4, kvec1, kvec2, kvec5)
	              + G5_temp_2(kvec3, kvec5, kvec1, kvec2, kvec4)
	              + G5_temp_2(kvec4, kvec5, kvec1, kvec2, kvec3);


	double Gtemp3 = G5_temp_3(kvec3, kvec4, kvec5, kvec1, kvec2)
	              + G5_temp_3(kvec2, kvec4, kvec5, kvec1, kvec3)
	              + G5_temp_3(kvec2, kvec3, kvec5, kvec1, kvec4)
	              + G5_temp_3(kvec2, kvec3, kvec4, kvec1, kvec5)
	              + G5_temp_3(kvec1, kvec4, kvec5, kvec2, kvec3)
	              + G5_temp_3(kvec1, kvec3, kvec5, kvec2, kvec4)
	              + G5_temp_3(kvec1, kvec3, kvec4, kvec2, kvec5)
	              + G5_temp_3(kvec1, kvec2, kvec5, kvec3, kvec4)
	              + G5_temp_3(kvec1, kvec2, kvec4, kvec3, kvec5)
	              + G5_temp_3(kvec1, kvec2, kvec3, kvec4, kvec5);

	return Gtemp1 / 5.0 + Gtemp2 / 10.0 + Gtemp3 / 10.0;
}

/***********************************************/

/*****************************/
/* kernel functions for bias */
/*****************************/

/* second order */
double D1D1(double * kvec1, double * kvec2) {
	return 1.0;
}

double K1K1(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
	return mu * mu - 1.0/3.0;
}

double G1G1(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
	return mu * mu - 1.0;
}

double F2_Bias(double * kvec1, double * kvec2, double b1, double b2, double bK2) {
	return   b1 * F2(kvec1, kvec2) 
	      +  (b2/2.0) * D1D1(kvec1, kvec2) 
	      +  bK2 * K1K1(kvec1, kvec2);
}

double F2_Bias_G2(double * kvec1, double * kvec2, double b1, double b2, double bG2) {
	return   b1 * F2(kvec1, kvec2) 
	      +  (b2/2.0) * D1D1(kvec1, kvec2) 
	      +  bG2 * G1G1(kvec1, kvec2);
}


/* third order */

double D1D2(double * kvec1, double * kvec2, double * kvec3) {
	return ( F2(kvec1, kvec2) + F2(kvec1,kvec3) + F2(kvec2,kvec3) ) / 3.0;
}

double K1K2(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);

	return ( K1K1(kvec1, kvec23) * F2(kvec2, kvec3) 
	       + K1K1(kvec2, kvec13) * F2(kvec1, kvec3)
	       + K1K1(kvec3, kvec12) * F2(kvec1, kvec2) ) / 3.0; 
}

double D1D1D1(double * kvec1, double * kvec2, double * kvec3) {
	return 1.0;
}

double K1K1K1(double * kvec1, double * kvec2, double * kvec3) {
	double mu12 = MU(kvec1, kvec2);
	double mu13 = MU(kvec1, kvec3);
	double mu23 = MU(kvec2, kvec3);

	return mu12 * mu13 * mu23 - (mu12*mu12 + mu13*mu13 + mu23*mu23) / 3.0 + (2.0/9.0);

}

double G1G1G1(double * kvec1, double * kvec2, double * kvec3) {
	double k12 = DOT(kvec1, kvec2);
	double k13 = DOT(kvec1, kvec3);
	double k23 = DOT(kvec2, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double kk1 = k1 * k1;
	double kk2 = k2 * k2;
	double kk3 = k3 * k3;


	return 3.0 / 2.0 * (k12 * k12 / kk1 / kk2 + k13 * k13 / kk1 / kk3 + k23 * k23 / kk2 / kk3)
			- k12 * k13 * k23 / kk1 / kk2 / kk3
			- 1.0 / 2.0;

}

double G_F(double * kvec1, double * kvec2, double * kvec3) {
	double k13[3] = PLUS(kvec1, kvec3);
	double k13s = DOT(k13,k13);
	double k23[3] = PLUS(kvec2, kvec3);
	double k23s = DOT(k23,k23);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double G123 = (DOT(kvec1,k23) * DOT(kvec1,k23) / k1 / k1 / k23s - 1.0) * F2(kvec2,kvec3);
	double G213 = (DOT(kvec2,k13) * DOT(kvec2,k13) / k2 / k2 / k13s - 1.0) * F2(kvec1,kvec3);

	return G123 + G213;

}

double G_G(double * kvec1, double * kvec2, double * kvec3) {
	double k13[3] = PLUS(kvec1, kvec3);
	double k13s = DOT(k13,k13);
	double k23[3] = PLUS(kvec2, kvec3);
	double k23s = DOT(k23,k23);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double G123 = (DOT(kvec1,k23) * DOT(kvec1,k23) / k1 / k1 / k23s - 1.0) * G2(kvec2,kvec3);
	double G213 = (DOT(kvec2,k13) * DOT(kvec2,k13) / k2 / k2 / k13s - 1.0) * G2(kvec1,kvec3);

	return G123 + G213;

}

double DVV(double * kvec1, double * kvec2, double * los) {
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	return mu1 * mu2 / k1 / k2;

}

double VVV(double * kvec1, double * kvec2, double * kvec3, double * los) {
	double mu1 = MU(kvec1, los);

	double dvv = DVV(kvec2, kvec3, los);

	return mu1 * mu1 * dvv;

}

double DV2(double * kvec1, double * kvec2, double * kvec3, double * los) {
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);
	double mu3 = MU(kvec3, los);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);

	double mu13 = MU(kvec13, los);
	double mu23 = MU(kvec23, los);

	double k13 = NORM(kvec13);
	double k23 = NORM(kvec23);

	double DV2_1 = mu1 / k1 * F2(kvec2,kvec3) + mu2 / k2 * F2(kvec1,kvec3) +mu3 / k3 * F2(kvec1,kvec2);
	double DV2_2 = mu13 / k13 * G2(kvec1,kvec3) + mu23 / k23 * G2(kvec2,kvec3);

	return DV2_1+ DV2_2;

}

double VV2(double * kvec1, double * kvec2, double * kvec3, double * los) {
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);
	double mu3 = MU(kvec3, los);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);

	double mu12 = MU(kvec12, los);
	double mu13 = MU(kvec13, los);
	double mu23 = MU(kvec23, los);

	double k13 = NORM(kvec13);
	double k23 = NORM(kvec23);

	double VV2_1 = mu1 / k1 * mu23 * mu23 * G2(kvec2,kvec3) + mu2 / k2 * mu13 * mu13 * G2(kvec1,kvec3) + mu3 / k3 * mu12 * mu12 * G2(kvec1,kvec2);
	double VV2_2 = mu1 * mu1 * mu23 / k23 * G2(kvec2,kvec3) + mu2 * mu2 * mu13 / k13 * G2(kvec1,kvec3);

	return VV2_1+ VV2_2;

}


double D1K1K1(double * kvec1, double * kvec2, double * kvec3) {
	double mu12 = MU(kvec1, kvec2);
	double mu13 = MU(kvec1, kvec3);
	double mu23 = MU(kvec2, kvec3);

	return (mu12 * mu12 + mu13 * mu13 + mu23 * mu23) / 3.0  - (1.0/3.0);
}

double O3(double * kvec1, double * kvec2, double * kvec3) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);

	double mu12 = MU(kvec1, kvec2);
	double mu13 = MU(kvec1, kvec3);
	double mu23 = MU(kvec2, kvec3);

	double K23 = K1K1(kvec1, kvec23);
	double K13 = K1K1(kvec2, kvec13);
	double K12 = K1K1(kvec3, kvec12);

	double KDS2 = (4.0/7.0) * ( K23 * (1.0 - mu23 * mu23) + K13 * (1.0 - mu13 * mu13) + K12 * (1.0 - mu12 * mu12) ) / 3.0; 

	return KDS2;

}

double F3_Bias(double * kvec1, double * kvec2, double * kvec3, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	return b1 * F3(kvec1, kvec2, kvec3)
	    +  b2 * D1D2(kvec1, kvec2, kvec3)
	    +  2.0 * bK2 * K1K2(kvec1, kvec2, kvec3)
	    +  (b3 / 6.0) * D1D1D1(kvec1, kvec2, kvec3)
	    +  bK3 * K1K1K1(kvec1, kvec2, kvec3)
	    +  bDK * D1K1K1(kvec1, kvec2, kvec3)
	    +  bO * O3(kvec1, kvec2, kvec3);

}

/* forth order */

double F4_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double b1) {
	return b1 * F4(kvec1, kvec2, kvec3, kvec4);
}

/* fifth order */

double F5_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double b1) {
	return b1 * F5(kvec1, kvec2, kvec3, kvec4, kvec5);
}

/******************************/
/* kernel functions with RSDs */
/******************************/

double LV1(double * kvec1, double * los) {

	double k1 = NORM(kvec1);
	double result = 0.0;
	if( k1 > pk_kmin ) {
		double k1n = DOT(kvec1, los);
		result += k1n / pow(k1,2);
	}
	return result;
}

double LV2(double * kvec1, double * kvec2, double * los) {

	double kvec12[3] = PLUS(kvec1, kvec2);
	double k12 = NORM(kvec12);

	double result = 0.0;
	if( k12 > pk_kmin ) {
		double k12n = DOT(kvec12, los);
		result += k12n * G2(kvec1, kvec2) / pow(k12,2);
	}
	return result;
}

double LV3(double * kvec1, double * kvec2, double * kvec3, double * los) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123 = NORM(kvec123);

	double result = 0.0;
	if( k123 > pk_kmin ) {
		double k123n = DOT(kvec123, los);
		result += k123n * G3(kvec1, kvec2, kvec3) / pow(k123,2);
	}
	return result;
}

double LV4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234 = NORM(kvec1234);

	double result = 0.0;
	if( k1234 > pk_kmin ) {
		double k1234n = DOT(kvec1234, los);
		result += k1234n * G4(kvec1, kvec2, kvec3, kvec4) / pow(k1234,2);
	}
	return result;
}


/* first order */

double Z1_Bias(double * kvec1, double * los, double f, double b1) {
	double mu = MU(kvec1, los);
	return b1 + f * mu * mu;
}

double Z1_Bias_FoG(double * kvec1, double * los, double f, double b1, double c1, double c2, double knl) {
    double mu = MU(kvec1, los);
	double k2 = NORM(kvec1) * NORM(kvec1);
	double fog_c1 = c1 * mu * mu * k2 / pow(knl,2);
	double fog_c2 = c1 * pow(mu,4) * k2 / pow(knl,2);
        return b1 + f * mu * mu - fog_c1 - fog_c2;
}

/* second order */
double V2(double * kvec1, double * kvec2, double * los, double f) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double mu = MU(kvec12, los);
	
	return f * mu * mu * G2(kvec1, kvec2);
}


double V1V1(double * kvec1, double * kvec2, double * los, double f) {

	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	double result = (f * f) * (kn * kn) * LV1(kvec1,los) * LV1(kvec2,los);
	return result;
}

double D1V1(double * kvec1, double * kvec2, double * los, double f) {
	
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	
	double result = (f) * (kn) * LV1(kvec1,los) / 2.0
		      + (f) * (kn) * LV1(kvec2,los) / 2.0;
	return result;
}

double Z2_Bias(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2) {
	return  V2(kvec1, kvec2, los, f)
	      + V1V1(kvec1, kvec2, los, f) / 2.0
	      + b1 * D1V1(kvec1, kvec2, los, f) 
	      + F2_Bias(kvec1, kvec2, b1, b2, bK2);
}

double Z2_Bias_G(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bG2) {
	return  V2(kvec1, kvec2, los, f)
	      + V1V1(kvec1, kvec2, los, f) / 2.0
	      + b1 * D1V1(kvec1, kvec2, los, f) 
	      + F2_Bias_G2(kvec1, kvec2, b1, b2, bG2);
}

/********************************/
/*  decomposed kernel functions */
/********************************/

double D1() {
	return 1.0;
}

double V1(double * kvec1, double * los) {
	double mu = MU(kvec1, los);
	return mu * mu;
}

double KFOG(double * kvec1, double * los) {
	double mu = MU(kvec1, los);
	double k1 = NORM(kvec1);
	return - mu * mu * k1 * k1 / 0.3 /0.3;
}

double KFOG2(double * kvec1, double * los) {
        double mu = MU(kvec1, los);
        double k1 = NORM(kvec1);
        return - mu * mu * mu * mu * k1 * k1 / 0.3 /0.3;
}
    
double V2(double * kvec1, double * kvec2, double * los) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double mu = MU(kvec12, los);
	return mu * mu * G2(kvec1, kvec2);
}

double V1V1(double * kvec1, double * kvec2, double * los) {
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	double result = (kn * kn) * LV1(kvec1,los) * LV1(kvec2,los);
	return result;
}

double D1V1(double * kvec1, double * kvec2, double * los) {
	
	double kvec12[3] = PLUS(kvec1, kvec2);
	double kn  = DOT(kvec12, los); 
	
	double result = (kn) * LV1(kvec1,los) / 2.0
		      + (kn) * LV1(kvec2,los) / 2.0;
	return result;
}

/****************************/

/* third order */

double V3(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double mu = MU(kvec123, los);
	return f * mu * mu * G3(kvec1, kvec2, kvec3);

}


double V1V2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f * f) * (k123n * k123n) * LV1(kvec1, los) * LV2(kvec2, kvec3, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec2, los) * LV2(kvec1, kvec3, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec3, los) * LV2(kvec1, kvec2, los) / 3.0;

	return result;

}

double V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f * f * f) * (k123n * k123n * k123n) * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los);

	return  result; 

}

double D1V1V1(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f * f) * (k123n * k123n) * LV1(kvec1, los) * LV1(kvec2, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec1, los) * LV1(kvec3, los) / 3.0
		      + (f * f) * (k123n * k123n) * LV1(kvec2, los) * LV1(kvec3, los) / 3.0;

	return result;

}


double D1V2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {
	
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f) * (k123n) * LV2(kvec1, kvec2, los)  / 3.0
		      + (f) * (k123n) * LV2(kvec1, kvec3, los)  / 3.0
		      + (f) * (k123n) * LV2(kvec2, kvec3, los)  / 3.0;

	return result;

}

double D2V1(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double b1, double b2, double bK2) {

	double kvec123[3] = PLUS3(kvec1,kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f) * (k123n) * F2_Bias(kvec2,kvec3, b1, b2, bK2) * LV1(kvec1, los) / 3.0
		      + (f) * (k123n) * F2_Bias(kvec1,kvec3, b1, b2, bK2) * LV1(kvec2, los) / 3.0
		      + (f) * (k123n) * F2_Bias(kvec1,kvec2, b1, b2, bK2) * LV1(kvec3, los) / 3.0;

	return result;

}

double Z3_Bias(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	return V3(kvec1, kvec2, kvec3, los, f)
	     + V1V2(kvec1, kvec2, kvec3, los, f)
	     + V1V1V1(kvec1, kvec2, kvec3, los, f) / 6.0
	     + b1 * D1V1V1(kvec1, kvec2, kvec3, los, f) / 2.0
	     + b1 * D1V2(kvec1, kvec2, kvec3, los, f)
	     + D2V1(kvec1, kvec2, kvec3, los, f, b1, b2, bK2)
	     + F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO);

}

double Z3_Bias_G(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double b1, double b2, double bG2, double bGamma3) {

	double k123_vec[3] = PLUS3(kvec1, kvec2, kvec3);
	double k123 = NORM(k123_vec);
	double mu = MU(k123_vec, los);

	//double mu1 = MU(kvec1, los);
	//double mu2 = MU(kvec2, los);
	//double mu3 = MU(kvec3, los);

	//double k1 = NORM(kvec1);
	//double k2 = NORM(kvec2);
	//double k3 = NORM(kvec3);

	double fmk = f * mu * k123;

	// return b3 / 6.0
	//	+ bG3 * G1G1G1(kvec1, kvec2, kvec3) / 3.0
	//	+ bG2d * (G1G1(kvec1, kvec2) + G1G1(kvec1, kvec3) + G1G1(kvec2, kvec3)) / 3.0
	//return bGamma3 * 2.0 * (G_F(kvec1,kvec2,kvec3) - G_G(kvec1,kvec2,kvec3)) / 3.0
	return b1 * F3(kvec1, kvec2, kvec3)
		+ f * mu * mu * G3(kvec1, kvec2, kvec3)
		+ b1 * fmk * fmk / 2 * (DVV(kvec2,kvec3,los) + DVV(kvec1,kvec3,los) + DVV(kvec1,kvec2,los)) / 3.0
		+ f * fmk * fmk / 2 * (VVV(kvec1,kvec2,kvec3,los) + VVV(kvec2,kvec1,kvec3,los) + VVV(kvec3,kvec1,kvec2,los)) / 3.0
		+ b1 * fmk * DV2(kvec1,kvec2,kvec3,los) / 3.0
		+ f * fmk * VV2(kvec1,kvec2,kvec3,los) / 3.0
	//	+ b2 * (F2(kvec1,kvec2) + F2(kvec1,kvec3) + F2(kvec2,kvec3)) / 3.0
		+ (bG2 * 2.0 + bGamma3 * 0.8) * G_F(kvec1,kvec2,kvec3) / 3.0;
	//	+ b2 * fmk / 2 * (mu1 / k1 + mu2 / k2 + mu3 / k3) / 3.0
	//	+ bG2 * fmk * (mu1 / k1 * G1G1(kvec1, kvec2) + mu2 / k2 * G1G1(kvec1, kvec3) + mu3 / k3 * G1G1(kvec2, kvec3)) / 3.0;
}

double D2V1_b2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1,kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f) * (k123n) * (D1D1(kvec2,kvec3)/2.0) * LV1(kvec1, los) / 3.0
		      + (f) * (k123n) * (D1D1(kvec1,kvec3)/2.0) * LV1(kvec2, los) / 3.0
		      + (f) * (k123n) * (D1D1(kvec1,kvec2)/2.0) * LV1(kvec3, los) / 3.0;

	return result;

}

double D2V1_bK2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec123[3] = PLUS3(kvec1,kvec2, kvec3);
	double k123n = DOT(kvec123, los);

	double result = (f) * (k123n) * K1K1(kvec2,kvec3) * LV1(kvec1, los) / 3.0
		      + (f) * (k123n) * K1K1(kvec1,kvec3) * LV1(kvec2, los) / 3.0
		      + (f) * (k123n) * K1K1(kvec1,kvec2) * LV1(kvec3, los) / 3.0;

	return result;

}

double Z3_Bias_b2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {
	return D1D2(kvec1, kvec2, kvec3) + D2V1_b2(kvec1, kvec2, kvec3, los, f);
}

double Z3_Bias_bK2(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {
	return 2.0 * K1K2(kvec1, kvec2, kvec3) + D2V1_bK2(kvec1, kvec2, kvec3, los, f);
}


/********************************/
/********************************/
/********************************/
/* fourth order */
/********************************/
/********************************/
/********************************/

double V4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {
	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double mu = MU(kvec1234, los);
	return f * mu * mu * G4(kvec1, kvec2, kvec3, kvec4);
}

double V2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV2(kvec1, kvec2, los) * LV2(kvec3, kvec4, los)  / 3.0
		      + (f * f) * (k1234n * k1234n) * LV2(kvec1, kvec3, los) * LV2(kvec2, kvec4, los)  / 3.0
		      + (f * f) * (k1234n * k1234n) * LV2(kvec1, kvec4, los) * LV2(kvec2, kvec3, los)  / 3.0;

	return result;

}


double V1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n= DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV3(kvec2, kvec3, kvec4, los) / 4.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV3(kvec1, kvec3, kvec4, los) / 4.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV3(kvec1, kvec2, kvec4, los) / 4.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV3(kvec1, kvec2, kvec3, los) / 4.0;

	return result;

}

double V1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n= DOT(kvec1234, los);

	double result = (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec3, kvec4, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec2, kvec4, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec2, kvec3, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec1, kvec4, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec1, kvec3, los) / 6.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec2, los) / 6.0;
	return result;

}

double V1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f * f * f) * (k1234n * k1234n * k1234n * k1234n)
		      * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los); 

	return result;
}

double D1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) / 4.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) / 4.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) / 4.0
		      + (f * f * f) * (k1234n * k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) / 4.0;

	return result;

}

double D2V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double bK2) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec2, los)  * F2_Bias(kvec3, kvec4, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec3, los)  * F2_Bias(kvec2, kvec4, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV1(kvec4, los)  * F2_Bias(kvec2, kvec3, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec3, los)  * F2_Bias(kvec1, kvec4, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV1(kvec4, los)  * F2_Bias(kvec1, kvec3, b1, b2, bK2) / 6.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV1(kvec4, los)  * F2_Bias(kvec1, kvec2, b1, b2, bK2) / 6.0;

	return result;

}

double D1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n= DOT(kvec1234, los);

	double result = (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV2(kvec2, kvec3, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV2(kvec2, kvec4, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec1, los) * LV2(kvec3, kvec4, los) / 12.0
	
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV2(kvec1, kvec3, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV2(kvec1, kvec4, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec2, los) * LV2(kvec3, kvec4, los) / 12.0

		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV2(kvec1, kvec2, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV2(kvec1, kvec4, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec3, los) * LV2(kvec2, kvec4, los) / 12.0

		      + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV2(kvec1, kvec2, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV2(kvec1, kvec3, los) / 12.0
		      + (f * f) * (k1234n * k1234n) * LV1(kvec4, los) * LV2(kvec2, kvec3, los) / 12.0;

	return result;

}

double D3V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f) * (k1234n) * LV1(kvec1, los) * F3_Bias(kvec2, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0
		      + (f) * (k1234n) * LV1(kvec2, los) * F3_Bias(kvec1, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0
		      + (f) * (k1234n) * LV1(kvec3, los) * F3_Bias(kvec1, kvec2, kvec4, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0
		      + (f) * (k1234n) * LV1(kvec4, los) * F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO) / 4.0;

	return result;

}

double D2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double bK2) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f) * (k1234n) * LV2(kvec1, kvec2, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec1, kvec3, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec1, kvec4, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec2, kvec3, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec2, kvec4, los) * F2_Bias(kvec1, kvec3, b1, b2, bK2) / 6.0
		      + (f) * (k1234n) * LV2(kvec3, kvec4, los) * F2_Bias(kvec1, kvec2, b1, b2, bK2) / 6.0;

	return result;

}

double D1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f) {

	double kvec1234[3] = PLUS4(kvec1, kvec2, kvec3, kvec4);
	double k1234n = DOT(kvec1234, los);

	double result = (f) * (k1234n) * LV3(kvec1, kvec2, kvec3, los) / 4.0
	              + (f) * (k1234n) * LV3(kvec1, kvec2, kvec4, los) / 4.0
	              + (f) * (k1234n) * LV3(kvec1, kvec3, kvec4, los) / 4.0
	              + (f) * (k1234n) * LV3(kvec2, kvec3, kvec4, los) / 4.0;

	return result;

}

double Z4_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	return V4(kvec1, kvec2, kvec3, kvec4, los, f)
	     + V2V2(kvec1, kvec2, kvec3, kvec4, los, f) / 2.0
	     + V1V3(kvec1, kvec2, kvec3, kvec4, los, f)
	     + V1V1V2(kvec1, kvec2, kvec3, kvec4, los, f) / 2.0
	     + V1V1V1V1(kvec1, kvec2, kvec3, kvec4, los, f) / 24.0

	     + b1 * D1V1V1V1(kvec1, kvec2, kvec3, kvec4, los, f) / 6.0
	     + b1 * D1V1V2(kvec1, kvec2, kvec3, kvec4, los, f) 
	     + b1 * D1V3(kvec1, kvec2, kvec3, kvec4, los, f) 

	     + D2V1V1(kvec1, kvec2, kvec3, kvec4, los, f, b1, b2, bK2) / 2.0
	     + D2V2(kvec1, kvec2, kvec3, kvec4, los, f, b1, b2, bK2) 
	     + D3V1(kvec1, kvec2, kvec3, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 

	     + F4_Bias(kvec1, kvec2, kvec3, kvec4, b1);

}
/********************************/
/********************************/
/********************************/
/* fifth order */
/********************************/
/********************************/
/********************************/

double V5(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec12345[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double mu = MU(kvec12345, los);
	
	return (f) * (mu * mu) * G5(kvec1, kvec2, kvec3, kvec4, kvec5);
}

double V1V4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV4(kvec2, kvec3, kvec4, kvec5, los)
		      + LV1(kvec2, los) * LV4(kvec1, kvec3, kvec4, kvec5, los)
		      + LV1(kvec3, los) * LV4(kvec1, kvec2, kvec4, kvec5, los)
		      + LV1(kvec4, los) * LV4(kvec1, kvec2, kvec3, kvec5, los)
		      + LV1(kvec5, los) * LV4(kvec1, kvec2, kvec3, kvec4, los);

	return pow(f, 2) * pow(kn, 2) * result / 5.0;

}

double V2V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV2(kvec1, kvec2, los) * LV3(kvec3, kvec4, kvec5, los)
	              + LV2(kvec1, kvec3, los) * LV3(kvec2, kvec4, kvec5, los)
	              + LV2(kvec1, kvec4, los) * LV3(kvec2, kvec3, kvec5, los)
	              + LV2(kvec1, kvec5, los) * LV3(kvec2, kvec3, kvec4, los)
	              + LV2(kvec2, kvec3, los) * LV3(kvec1, kvec4, kvec5, los)
	              + LV2(kvec2, kvec4, los) * LV3(kvec1, kvec3, kvec5, los)
	              + LV2(kvec2, kvec5, los) * LV3(kvec1, kvec3, kvec4, los)
	              + LV2(kvec3, kvec4, los) * LV3(kvec1, kvec2, kvec5, los)
	              + LV2(kvec3, kvec5, los) * LV3(kvec1, kvec2, kvec4, los)
	              + LV2(kvec4, kvec5, los) * LV3(kvec1, kvec2, kvec3, los);

	return pow(f, 2) * pow(kn, 2) * result / 10.0;

}

double V1V2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV2(kvec2, kvec3, los)  * LV2(kvec4, kvec5, los)
		      + LV1(kvec1, los) * LV2(kvec2, kvec4, los)  * LV2(kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV2(kvec2, kvec5, los)  * LV2(kvec3, kvec4, los)

		      + LV1(kvec2, los) * LV2(kvec1, kvec3, los)  * LV2(kvec4, kvec5, los)
		      + LV1(kvec2, los) * LV2(kvec1, kvec4, los)  * LV2(kvec3, kvec5, los)
		      + LV1(kvec2, los) * LV2(kvec1, kvec5, los)  * LV2(kvec3, kvec4, los)

		      + LV1(kvec3, los) * LV2(kvec1, kvec2, los)  * LV2(kvec4, kvec5, los)
		      + LV1(kvec3, los) * LV2(kvec1, kvec4, los)  * LV2(kvec2, kvec5, los)
		      + LV1(kvec3, los) * LV2(kvec1, kvec5, los)  * LV2(kvec2, kvec4, los)

		      + LV1(kvec4, los) * LV2(kvec1, kvec2, los)  * LV2(kvec3, kvec5, los)
		      + LV1(kvec4, los) * LV2(kvec1, kvec3, los)  * LV2(kvec2, kvec5, los)
		      + LV1(kvec4, los) * LV2(kvec1, kvec5, los)  * LV2(kvec2, kvec3, los)

		      + LV1(kvec5, los) * LV2(kvec1, kvec2, los)  * LV2(kvec3, kvec4, los)
		      + LV1(kvec5, los) * LV2(kvec1, kvec3, los)  * LV2(kvec2, kvec4, los)
		      + LV1(kvec5, los) * LV2(kvec1, kvec4, los)  * LV2(kvec2, kvec3, los);

	return pow(f, 3) * pow(kn, 3) * result / 15.0;

}

double V1V1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV3(kvec3, kvec4, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV3(kvec2, kvec4, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV3(kvec2, kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec5, los) * LV3(kvec2, kvec3, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV3(kvec1, kvec4, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV3(kvec1, kvec3, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec5, los) * LV3(kvec1, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV3(kvec1, kvec2, kvec5, los)
		      + LV1(kvec3, los) * LV1(kvec5, los) * LV3(kvec1, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV1(kvec5, los) * LV3(kvec1, kvec2, kvec3, los);

	return pow(f, 3) * pow(kn, 3) * result / 10.0;

}

double V1V1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los) * LV2(kvec1, kvec2, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV1(kvec5, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec5, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV1(kvec5, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec5, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec2, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec5, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec4, kvec5, los);

	return pow(f, 4) * pow(kn, 4) * result / 10.0;

}

double V1V1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los);

	return pow(f, 5) * pow(kn, 5) * result;

}

double D1V1V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los)
	              + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec5, los)
	              + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) * LV1(kvec5, los)
	              + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los)
	              + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los);

	return pow(f, 4) * pow(kn, 4) * result / 5.0;
}

double D1V1V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec2, los)

		      + LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec3, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec2, kvec5, los)
		      + LV1(kvec1, los) * LV1(kvec5, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec1, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec5, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec5, los) * LV2(kvec1, kvec2, los)

		      + LV1(kvec1, los) * LV1(kvec2, los) * LV2(kvec5, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec5, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec2, kvec5, los)
		      + LV1(kvec2, los) * LV1(kvec5, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec1, kvec5, los)
		      + LV1(kvec5, los) * LV1(kvec4, los) * LV2(kvec1, kvec2, los)

		      + LV1(kvec1, los) * LV1(kvec5, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV2(kvec5, kvec4, los)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV2(kvec5, kvec3, los)
		      + LV1(kvec5, los) * LV1(kvec3, los) * LV2(kvec1, kvec4, los)
		      + LV1(kvec5, los) * LV1(kvec4, los) * LV2(kvec1, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec1, kvec5, los)

		      + LV1(kvec5, los) * LV1(kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV1(kvec5, los) * LV1(kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV1(kvec5, los) * LV1(kvec4, los) * LV2(kvec2, kvec3, los)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV2(kvec5, kvec4, los)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV2(kvec5, kvec3, los)
		      + LV1(kvec3, los) * LV1(kvec4, los) * LV2(kvec5, kvec2, los);

	return pow(f, 3) * pow(kn, 3) * result / 30.0;

}

double D1V1V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV3(kvec2, kvec3, kvec4, los)
		      + LV1(kvec2, los) * LV3(kvec1, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV3(kvec1, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec1, kvec2, kvec3, los)

		      + LV1(kvec1, los) * LV3(kvec2, kvec3, kvec5, los)
		      + LV1(kvec2, los) * LV3(kvec1, kvec3, kvec5, los)
		      + LV1(kvec3, los) * LV3(kvec1, kvec2, kvec5, los)
		      + LV1(kvec5, los) * LV3(kvec1, kvec2, kvec3, los)

		      + LV1(kvec1, los) * LV3(kvec2, kvec5, kvec4, los)
		      + LV1(kvec2, los) * LV3(kvec1, kvec5, kvec4, los)
		      + LV1(kvec5, los) * LV3(kvec1, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec1, kvec2, kvec5, los)

		      + LV1(kvec1, los) * LV3(kvec5, kvec3, kvec4, los)
		      + LV1(kvec5, los) * LV3(kvec1, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV3(kvec1, kvec5, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec1, kvec5, kvec3, los)

		      + LV1(kvec5, los) * LV3(kvec2, kvec3, kvec4, los)
		      + LV1(kvec2, los) * LV3(kvec5, kvec3, kvec4, los)
		      + LV1(kvec3, los) * LV3(kvec5, kvec2, kvec4, los)
		      + LV1(kvec4, los) * LV3(kvec5, kvec2, kvec3, los);

	return pow(f, 2) * pow(kn, 2) * result / 20.0;
}


double D1V2V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV2(kvec1, kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV2(kvec1, kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV2(kvec1, kvec4, los) * LV2(kvec2, kvec3, los)

		      + LV2(kvec1, kvec2, los) * LV2(kvec3, kvec5, los)
		      + LV2(kvec1, kvec3, los) * LV2(kvec2, kvec5, los)
		      + LV2(kvec1, kvec5, los) * LV2(kvec2, kvec3, los)

		      + LV2(kvec1, kvec2, los) * LV2(kvec5, kvec4, los)
		      + LV2(kvec1, kvec5, los) * LV2(kvec2, kvec4, los)
		      + LV2(kvec1, kvec4, los) * LV2(kvec2, kvec5, los)

		      + LV2(kvec1, kvec5, los) * LV2(kvec3, kvec4, los)
		      + LV2(kvec1, kvec3, los) * LV2(kvec5, kvec4, los)
		      + LV2(kvec1, kvec4, los) * LV2(kvec5, kvec3, los)

		      + LV2(kvec5, kvec2, los) * LV2(kvec3, kvec4, los)
		      + LV2(kvec5, kvec3, los) * LV2(kvec2, kvec4, los)
		      + LV2(kvec5, kvec4, los) * LV2(kvec2, kvec3, los);

	return pow(f, 2) * pow(kn, 2) * result / 15.0;

}


double D1V4(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f) {
	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV4(kvec1, kvec2, kvec3, kvec4, los)
		      + LV4(kvec1, kvec2, kvec3, kvec5, los)
		      + LV4(kvec1, kvec2, kvec5, kvec4, los)
		      + LV4(kvec1, kvec5, kvec3, kvec4, los)
	              + LV4(kvec5, kvec2, kvec3, kvec4, los);

	return pow(f, 1) * pow(kn, 1) * result / 5.0;

}

double D2V1V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double bK2) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec3, los) * LV1(kvec4, los) * LV1(kvec5, los) * F2_Bias(kvec1, kvec2, b1, b2, bK2)
		      + LV1(kvec2, los) * LV1(kvec4, los) * LV1(kvec5, los) * F2_Bias(kvec1, kvec3, b1, b2, bK2)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec5, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2)
		      + LV1(kvec2, los) * LV1(kvec3, los) * LV1(kvec4, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec4, los) * LV1(kvec5, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec5, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec3, los) * LV1(kvec4, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec5, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec4, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV1(kvec2, los) * LV1(kvec3, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2);

	return pow(f, 3) * pow(kn, 3) * result / 10.0;
}

double D2V1V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double bK2) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV2(kvec2, kvec3, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec2, kvec4, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec2, kvec5, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec3, kvec4, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec3, kvec5, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec1, los) * LV2(kvec4, kvec5, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2)

		      + LV1(kvec2, los) * LV2(kvec1, kvec3, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec1, kvec4, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec1, kvec5, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec3, kvec4, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec3, kvec5, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2)
		      + LV1(kvec2, los) * LV2(kvec4, kvec5, los) * F2_Bias(kvec1, kvec3, b1, b2, bK2)

		      + LV1(kvec3, los) * LV2(kvec2, kvec1, los) * F2_Bias(kvec4, kvec5, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec2, kvec4, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec2, kvec5, los) * F2_Bias(kvec1, kvec4, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec1, kvec4, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec1, kvec5, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec3, los) * LV2(kvec4, kvec5, los) * F2_Bias(kvec2, kvec1, b1, b2, bK2)

		      + LV1(kvec4, los) * LV2(kvec2, kvec3, los) * F2_Bias(kvec1, kvec5, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec2, kvec1, los) * F2_Bias(kvec3, kvec5, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec2, kvec5, los) * F2_Bias(kvec3, kvec1, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec3, kvec1, los) * F2_Bias(kvec2, kvec5, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec3, kvec5, los) * F2_Bias(kvec2, kvec1, b1, b2, bK2)
		      + LV1(kvec4, los) * LV2(kvec1, kvec5, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2)

		      + LV1(kvec5, los) * LV2(kvec2, kvec3, los) * F2_Bias(kvec4, kvec1, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec2, kvec4, los) * F2_Bias(kvec3, kvec1, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec2, kvec1, los) * F2_Bias(kvec3, kvec4, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec3, kvec4, los) * F2_Bias(kvec2, kvec1, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec3, kvec1, los) * F2_Bias(kvec2, kvec4, b1, b2, bK2)
		      + LV1(kvec5, los) * LV2(kvec4, kvec1, los) * F2_Bias(kvec2, kvec3, b1, b2, bK2);

	return pow(f, 2) * pow(kn, 2) * result / 30.0;
}


double D2V3(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double bK2) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = F2_Bias(kvec1, kvec2, b1, b2, bK2) * LV3(kvec3, kvec4, kvec5, los)
	              + F2_Bias(kvec1, kvec3, b1, b2, bK2) * LV3(kvec2, kvec4, kvec5, los)
	              + F2_Bias(kvec1, kvec4, b1, b2, bK2) * LV3(kvec2, kvec3, kvec5, los)
	              + F2_Bias(kvec1, kvec5, b1, b2, bK2) * LV3(kvec2, kvec3, kvec4, los)
	              + F2_Bias(kvec2, kvec3, b1, b2, bK2) * LV3(kvec1, kvec4, kvec5, los)
	              + F2_Bias(kvec2, kvec4, b1, b2, bK2) * LV3(kvec1, kvec3, kvec5, los)
	              + F2_Bias(kvec2, kvec5, b1, b2, bK2) * LV3(kvec1, kvec3, kvec4, los)
	              + F2_Bias(kvec3, kvec4, b1, b2, bK2) * LV3(kvec1, kvec2, kvec5, los)
	              + F2_Bias(kvec3, kvec5, b1, b2, bK2) * LV3(kvec1, kvec2, kvec4, los)
	              + F2_Bias(kvec4, kvec5, b1, b2, bK2) * LV3(kvec1, kvec2, kvec3, los);

	return pow(f, 1) * pow(kn, 1) * result / 10.0;

}

double D3V1V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * LV1(kvec2, los) * F3_Bias(kvec3, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec1, los) * LV1(kvec3, los) * F3_Bias(kvec2, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec1, los) * LV1(kvec4, los) * F3_Bias(kvec2, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec1, los) * LV1(kvec5, los) * F3_Bias(kvec2, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec2, los) * LV1(kvec3, los) * F3_Bias(kvec1, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec2, los) * LV1(kvec4, los) * F3_Bias(kvec1, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec2, los) * LV1(kvec5, los) * F3_Bias(kvec1, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec3, los) * LV1(kvec4, los) * F3_Bias(kvec1, kvec2, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec3, los) * LV1(kvec5, los) * F3_Bias(kvec1, kvec2, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
		      + LV1(kvec4, los) * LV1(kvec5, los) * F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO);

	return pow(f, 2) * pow(kn, 2) * result / 10.0;


}

double D3V2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV2(kvec1, kvec2, los) * F3_Bias(kvec3, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec1, kvec3, los) * F3_Bias(kvec2, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec1, kvec4, los) * F3_Bias(kvec2, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec1, kvec5, los) * F3_Bias(kvec2, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec2, kvec3, los) * F3_Bias(kvec1, kvec4, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec2, kvec4, los) * F3_Bias(kvec1, kvec3, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec2, kvec5, los) * F3_Bias(kvec1, kvec3, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec3, kvec4, los) * F3_Bias(kvec1, kvec2, kvec5, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec3, kvec5, los) * F3_Bias(kvec1, kvec2, kvec4, b1, b2, b3, bK2, bK3, bDK, bO)
	              + LV2(kvec4, kvec5, los) * F3_Bias(kvec1, kvec2, kvec3, b1, b2, b3, bK2, bK3, bDK, bO);

	return pow(f, 1) * pow(kn, 1) * result / 10.0;

}

double D4V1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, double b1) {

	double kvec[3] = PLUS5(kvec1, kvec2, kvec3, kvec4, kvec5);
	double kn = DOT(kvec, los);

	double result = LV1(kvec1, los) * F4_Bias(kvec2, kvec3, kvec4, kvec5, b1)
		      + LV1(kvec2, los) * F4_Bias(kvec1, kvec3, kvec4, kvec5, b1)
		      + LV1(kvec3, los) * F4_Bias(kvec1, kvec2, kvec4, kvec5, b1)
		      + LV1(kvec4, los) * F4_Bias(kvec1, kvec2, kvec3, kvec5, b1)
		      + LV1(kvec5, los) * F4_Bias(kvec1, kvec2, kvec3, kvec4, b1);

	return pow(f, 1) * pow(kn, 1) * result / 5.0;

}

double Z5_Bias(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double f, 
	       double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double Z = V5(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + V1V4(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + V2V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + V1V2V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + V1V1V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + V1V1V1V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 6.0
		 + V1V1V1V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 120.0

		 + b1 * D1V1V1V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 24.0
		 + b1 * D1V1V1V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + b1 * D1V1V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 + b1 * D1V2V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f) / 2.0
		 + b1 * D1V4(kvec1, kvec2, kvec3, kvec4, kvec5, los, f)
		 
		 + D2V1V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, bK2) / 6.0
		 + D2V1V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, bK2) / 6.0
		 + D2V3(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, bK2) 

		 + D3V1V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, b3, bK2, bK3, bDK, bO) / 2.0
		 + D3V2(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 

		 + D4V1(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1) 

		 + F5_Bias(kvec1, kvec2, kvec3, kvec4, kvec5, b1);

	return Z;

}

/****************************/
/*     Reconstruction       */
/****************************/

double Z1_Bias_D(double * kvec, double * los, double f, double b1, double b1_fid, double R) {
	double k = NORM(kvec);
    	double W = exp( - pow(k * R, 2) / 2.0);
        return (1.0 - (W / b1_fid)) * Z1_Bias(kvec, los, f, b1);
}

double Z1S1(double * kvec1, double * kvec2, double * los, double f, double b1, double b1_fid, double R) {
    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
    double W1 = exp( - pow(k1 * R, 2) / 2.0);
    double W2 = exp( - pow(k2 * R, 2) / 2.0);
    double kvec12[3] = PLUS(kvec1, kvec2);
    double result = (-1.0/2.0) * (1.0/b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1);
    return result;
}

double Z1S1_b1_b1(double * kvec1, double * kvec2, double b1_fid, double R) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = (-1.0/2.0) * (1.0/b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) * D1() * D1();
	return result;
}

double Z1S1_b1_f(double * kvec1, double * kvec2, double * los, double b1_fid, double R) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = (-1.0/2.0) * (1.0/b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) 
	              * ( D1() * V1(kvec1, los) + V1(kvec2, los) * D1() );
	return result;
}

double Z1S1_f_f(double * kvec1, double * kvec2, double * los, double b1_fid, double R) {
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double kvec12[3] = PLUS(kvec1, kvec2);
	double result = (-1.0/2.0) * (1.0/b1_fid) * ( LV1(kvec1, kvec12) * W1 + LV1(kvec2, kvec12) * W2 ) 
	              * V1(kvec1, los) * V1(kvec2, los);
	return result;
}


double Z2_Bias_Reconstructed(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2, double b1_fid, double R) {
	return Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) 
	     + Z1S1(kvec1, kvec2, los, f, b1, b1_fid, R);
} 

double Z2_Bias_D(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2, double b1_fid, double R) {

    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
    double W1 = exp( - pow(k1 * R, 2) / 2.0);
    double W2 = exp( - pow(k2 * R, 2) / 2.0);
    double kvec12[3] = PLUS(kvec1, kvec2);
    double k12 = NORM(kvec12);
    double W12 = exp( - pow(k12 * R, 2) / 2.0);
    
    double Z1_1 = Z1_Bias(kvec1, los, f, b1);
    double Z1_2 = Z1_Bias(kvec2, los, f, b1);
    
    double Z_rec = Z2_Bias_Reconstructed(kvec1, kvec2, los, f, b1, b2, bK2, b1_fid, R);
    
    double Z_a = (-1.0) * (W12/b1_fid) * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2);
    double Z_b = (1.0/2.0) * ( LV1(kvec1, kvec12) * (W1/b1_fid) * LV1(kvec2, kvec12) * (W2/b1_fid) ) * Z1_1 * Z1_2;
    
    return Z_rec + Z_a + Z_b;
    
} 

double Z3_Bias_Reconstructed(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double b1_fid, double R) {

	double kvec12[3] = PLUS(kvec1, kvec2);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec123[3] = PLUS3(kvec1, kvec2, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double k12 = NORM(kvec12);
	double k13 = NORM(kvec13);
	double k23 = NORM(kvec23);

    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
	double W2 = exp( - pow(k2 * R, 2) / 2.0);
	double W3 = exp( - pow(k3 * R, 2) / 2.0);

    	double W12 = exp( - pow(k12 * R, 2) / 2.0);
	double W23 = exp( - pow(k23 * R, 2) / 2.0);
	double W13 = exp( - pow(k13 * R, 2) / 2.0);

	double Z3 = Z3_Bias(kvec1, kvec2, kvec3, los, f, b1, b2, b3, bK2, bK3, bDK, bO);

	double Z1_1 = Z1_Bias(kvec1, los, f, b1);
	double Z1_2 = Z1_Bias(kvec2, los, f, b1);
	double Z1_3 = Z1_Bias(kvec3, los, f, b1);

	double Z2_12 = Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2);
	double Z2_13 = Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2);
	double Z2_23 = Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2);

	double R_A = ( LV1(kvec23, kvec123) * (W23 / b1_fid) + LV1(kvec1, kvec123) *  (W1 / b1_fid) ) * Z1_1 * Z2_23
		   + ( LV1(kvec13, kvec123) * (W13 / b1_fid) + LV1(kvec2, kvec123) *  (W2 / b1_fid) ) * Z1_2 * Z2_13
		   + ( LV1(kvec12, kvec123) * (W12 / b1_fid) + LV1(kvec3, kvec123) *  (W3 / b1_fid) ) * Z1_3 * Z2_12;

	double R_B = ( LV1(kvec1, kvec123) * (W1/b1_fid) * LV1(kvec2, kvec123) * (W2/b1_fid)
		   +   LV1(kvec1, kvec123) * (W1/b1_fid) * LV1(kvec3, kvec123) * (W3/b1_fid)
		   +   LV1(kvec2, kvec123) * (W2/b1_fid) * LV1(kvec3, kvec123) * (W3/b1_fid) ) * Z1_1 * Z1_2 * Z1_3;

	return Z3 - (1.0/3.0) * R_A + (1.0/6.0) * R_B;
}

double SN_bispectrum(double * kvec, double * los, double b1, double f, double Pshot, double Bshot) {

	double mu = MU(kvec, los);

	return Bshot * b1 + 2 * f * mu * mu * (1 + Pshot);

}

//
//double V1S1(double * kvec1, double * kvec2, double * los, double f, double b1, double b1_fid, double R) {
//	double k1 = NORM(kvec1);
//	double k2 = NORM(kvec2);
//    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
//	double W2 = exp( - pow(k2 * R, 2) / 2.0);
//	double kvec12[3] = PLUS(kvec1, kvec2);
//	double kn  = DOT(kvec12, los); 
//	double result = (-1.0/b1_fid) * LV1(kvec1, kvec12) * W1 * Z1_Bias(kvec1, los, f, b1) * (f) * (kn) * LV1(kvec1,los) / 2.0
//                      + (-1.0/b1_fid) * LV1(kvec2, kvec12) * W2 * Z1_Bias(kvec2, los, f, b1) * (f) * (kn) * LV1(kvec2,los) / 2.0;
//	return result;
//
//}
//
//double D1S1(double * kvec1, double * kvec2, double * los, double f, double b1, double b1_fid, double R) {
//	double k1 = NORM(kvec1);
//	double k2 = NORM(kvec2);
//    	double W1 = exp( - pow(k1 * R, 2) / 2.0);
//	double W2 = exp( - pow(k2 * R, 2) / 2.0);
//	double kvec12[3] = PLUS(kvec1, kvec2);
//	double result = (-1.0/b1_fid) * LV1(kvec1, kvec12) * W1 * Z1_Bias(kvec1, los, f, b1) / 2.0
//                      + (-1.0/b1_fid) * LV1(kvec2, kvec12) * W2 * Z1_Bias(kvec2, los, f, b1) / 2.0;
//	return result;
//}
//
//double Z2_Bias_Reconstructed(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2, double b1_fid, double R) {
//	return Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) 
//	     + V1S1(kvec1, kvec2, los, f, b1, b1_fid, R)
//	     + b1 * D1S1(kvec1, kvec2, los, f, b1, b1_fid, R);
//} 

double ExpDamping(double * kvec, double * los, double sigma2_perp, double sigma2_para ) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);
	double mu2 = mu * mu;
	double lnD = - k * k * ( (1.0 - mu2) * sigma2_perp + mu2 * sigma2_para ) / 2.0;
	double D = exp(lnD);
	return D;

}

double IntegrandDampingBAO(double k, double rbao) {
	
	double result = f_pk_no_wiggle(k) * (1 - 3 * sin(k * rbao) / (k * rbao) + 6  * sin(k * rbao) / pow((k * rbao), 3) - 6 * cos(k * rbao) / pow((k * rbao), 2));

	return result;

}

double IntegrandDampingBAO_add(double k, double rbao) {

	double result = f_pk_no_wiggle(k) * (- sin(k * rbao) / (k * rbao) + 3  * sin(k * rbao) / pow((k * rbao), 3) - 3 * cos(k * rbao) / pow((k * rbao), 2));

        return result;

}

double Sig2(double rbao, double ks) {

        double cst = 1 / (6 * M_PI * M_PI);
        double Nintegration = 1000;
        double step = ks / Nintegration;
        double integral = 0.;
        for (int index=0; index<Nintegration; index++){
		double eps = .00001;
		double ki = step * index + eps;
		double kf = ki + step;
		double ri = IntegrandDampingBAO(ki, rbao);
		double rf = IntegrandDampingBAO(kf, rbao);
		integral += step * (ri + rf) / 2;
	}
        double result = cst * integral;
	return result;
}

double dSig2(double rbao, double ks) {

        double cst = 1 / (2 * M_PI * M_PI);
        double Nintegration = 1000;
        double step = ks / Nintegration;
        double integral = 0.;
        for (int index=0; index<Nintegration; index++){
        	double eps = .00001;
                double ki = step * index + eps;
                double kf = ki + step;
                double ri = IntegrandDampingBAO_add(ki, rbao);
                double rf = IntegrandDampingBAO_add(kf, rbao);
                integral += step * (ri + rf) / 2;
	}
        double result = cst * integral;
        return result;

}

double ExpDamping_Ivanov(double * kvec, double * los, double f, double Sigma2, double dSigma2) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);
	double mu2 = mu * mu;
	double lnD = - Sigma2 * k * k * (1 + f * mu2 * (2 + f)) - dSigma2 * k * k * f * f * mu2 * (mu2 -1);
	double D = exp(lnD);
	return D;

}

double ExpDamping_RealSpace(double * kvec, double sigma2_perp) {

	double k = NORM(kvec);
	double lnD = - k * k * ( sigma2_perp ) / 2.0;
	double D = exp(lnD);
	return D;

}

double Powerspectrum_Tree(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double s8, double f, double b1) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	
	double P = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * f_pk(k);
	return P / alpha3;
}

double Powerspectrum_Tree(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	
	double P = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * f_pk(k);
	return P / alpha3;
}

double Powerspectrum_Tree_Damping(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1, 
									double Sigma2, double dSigma2) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);

    
	double BAO = f_pk(k) - f_pk_no_wiggle(k);
	double D = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);
	
	double P = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * (f_pk_no_wiggle(k) + BAO * D);
	return P / alpha3;
}

double Powerspectrum_Tree_Damping_for_1loop(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1, 
									double Sigma2, double dSigma2) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);

	double BAO = f_pk(k) - f_pk_no_wiggle(k);
	double D = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);
	
	double P = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * (f_pk_no_wiggle(k) + BAO * D * ( 1 - log(D) ) );
	return P / alpha3;
}

double Powerspectrum_Counterterm(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1, 
											double c0, double c1, double c2, double ch, double knl) {
	
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double mu = MU(kvec, los);
	double k = NORM(kvec);

	double ct = ( c0 + c1 * f * mu * mu + c2 * f * f * pow(mu,4) ) * k * k;
	
	double P = - 2 * Z1_Bias(kvec, los, f, b1) * f_pk(k) * ct;
    
    double Ph = - ch * pow(f * mu * k , 4) * Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * f_pk(k);
    
    double Ptot = P + Ph;
    
	return Ptot / alpha3;
}

double Powerspectrum_Counterterm_Damping(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1, 
											double c0, double c1, double c2, double ch, double knl, double Sigma2, double dSigma2) {
	
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double mu = MU(kvec, los);
	double k = NORM(kvec);

	double BAO = f_pk(k) - f_pk_no_wiggle(k);
	double D = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);

	double ct = ( c0 + c1 * f * mu * mu + c2 * f * f * pow(mu,4) ) * k * k;
	
	double P = - 2 * Z1_Bias(kvec, los, f, b1) * (f_pk_no_wiggle(k) + BAO * D) * ct;
    
    double Ph = - ch * pow(f * mu * k , 4) * Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * (f_pk_no_wiggle(k) + BAO * D);
    
    double Ptot = P + Ph;
    
    
	return Ptot / alpha3;
}

//double Powerspectrum_1loop_Damping( double * kvec_in, double * pvec,
//		               double * los, double alpha_perp, double alpha_parallel,
//			       double f, double b1, double b2, double bG2, double b3,
//			       double bG3, double bG2d, double bGamma3, double Sigma2, double dSigma2) {
//
//    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
//	double kvec[3] = { 0.0, 0.0, 0.0 };
//	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
//	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
//	double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
//	double k_M_p = NORM(kvec_M_pvec);
//	double p = NORM(pvec);
//	double k = NORM(kvec);
//
//	double BAOk = f_pk(k) - f_pk_no_wiggle(k);
//	double BAOp = f_pk(p) - f_pk_no_wiggle(p);
//	double BAOkMp = f_pk(k_M_p) - f_pk_no_wiggle(k_M_p);
//
//	double Dk = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);
//	double Dp = ExpDamping_Ivanov(pvec, los, f, Sigma2, dSigma2);
//	double DkMp = ExpDamping_Ivanov(kvec_M_pvec, los, f, Sigma2, dSigma2);
//
//	double P22 = 2.0 * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2)
//				* (f_pk_no_wiggle(k_M_p) + BAOkMp * DkMp) * (f_pk_no_wiggle(p) + BAOp * Dp)
//				- 2.0 / 4.0 * b2 * b2 * (f_pk_no_wiggle(p) + BAOp * Dp) * (f_pk_no_wiggle(p) + BAOp * Dp);
//
//	double P13 = 6.0 * Z1_Bias(kvec, los, f, b1)
//	           	* Z3_Bias_G(pvec, M_pvec, kvec, los, f, b1, b2, bG2, b3, bG3, bG2d, bGamma3)
//		   		* (f_pk_no_wiggle(p) + BAOp * Dp) * (f_pk_no_wiggle(k) + BAOk * Dk);
//
//	double P1loop = P22 + P13;
//
//	return P1loop / alpha3;
//}

double Powerspectrum_1loop( double * kvec_in, double * pvec,
		               double * los, double alpha_perp, double alpha_parallel,
			       double f, double b1, double b2, double bG2, double bGamma3) {

   double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
	double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
	double k_M_p = NORM(kvec_M_pvec);
	double kvec_P_pvec[3] = { kvec[0] + pvec[0], kvec[1] + pvec[1], kvec[2] + pvec[2] };
	double k_P_p = NORM(kvec_P_pvec);
	double p = NORM(pvec);
	double k = NORM(kvec);

	double P22 = 0.0;

	if (k_M_p > p) {
	    P22 += 2.0 * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2)
				* f_pk(k_M_p) * f_pk(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk(p) * f_pk(p);
	}

	if (k_P_p > p) {
	    P22 += 2.0 * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2)
				* f_pk(k_P_p) * f_pk(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk(p) * f_pk(p);
	}

	double P13 = 6.0 * Z1_Bias(kvec, los, f, b1)
	           	* Z3_Bias_G(pvec, M_pvec, kvec, los, f, b1, b2, bG2, bGamma3)
		   		* f_pk(p) * f_pk(k);

	double P1loop = P22 + P13;

	return P1loop / alpha3;
}

//double Powerspectrum_1loop_Damping( double * kvec_in, double * pvec,
//		               double * los, double alpha_perp, double alpha_parallel,
//			       double f, double b1, double b2, double bG2, double bGamma3, double Sigma2, double dSigma2) {

//   double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
//	double kvec[3] = { 0.0, 0.0, 0.0 };
//	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
//	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
//	double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
//	double k_M_p = NORM(kvec_M_pvec);
//	double kvec_P_pvec[3] = { kvec[0] + pvec[0], kvec[1] + pvec[1], kvec[2] + pvec[2] };
//	double k_P_p = NORM(kvec_P_pvec);
//	double p = NORM(pvec);
//	double k = NORM(kvec);

//	double BAOk = f_pk(k) - f_pk_no_wiggle(k);
//	double BAOp = f_pk(p) - f_pk_no_wiggle(p);
//	double BAOkMp = f_pk(k_M_p) - f_pk_no_wiggle(k_M_p);
//	double BAOkPp = f_pk(k_P_p) - f_pk_no_wiggle(k_P_p);

//	double Dk = ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2);
//	double Dp = ExpDamping_Ivanov(pvec, los, f, Sigma2, dSigma2);
//	double DkMp = ExpDamping_Ivanov(kvec_M_pvec, los, f, Sigma2, dSigma2);
//	double DkPp = ExpDamping_Ivanov(kvec_P_pvec, los, f, Sigma2, dSigma2);

//	double P22 = 0.0;

//	if (k_M_p > p) {
//	    P22 += 2.0 * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2)
//				* (f_pk_no_wiggle(k_M_p) + BAOkMp * DkMp) * (f_pk_no_wiggle(p) + BAOp * Dp)
//				- 2.0 / 4.0 * b2 * b2 * (f_pk_no_wiggle(p) + BAOp * Dp) * (f_pk_no_wiggle(p) + BAOp * Dp);
//	}

//	if (k_P_p > p) {
//	    P22 += 2.0 * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2)
//				* (f_pk_no_wiggle(k_P_p) + BAOkPp * DkPp) * (f_pk_no_wiggle(p) + BAOp * Dp)
//				- 2.0 / 4.0 * b2 * b2 * (f_pk_no_wiggle(p) + BAOp * Dp) * (f_pk_no_wiggle(p) + BAOp * Dp);
//	}

//	double P13 = 6.0 * Z1_Bias(kvec, los, f, b1)
//	           	* Z3_Bias_G(pvec, M_pvec, kvec, los, f, b1, b2, bG2, bGamma3)
//		   		* (f_pk_no_wiggle(p) + BAOp * Dp) * (f_pk_no_wiggle(k) + BAOk * Dk);

//	double P1loop = P22 + P13;

//	return P1loop / alpha3;
//}

double Powerspectrum_1loop_Damping( double * kvec_in, double * pvec,
		               double * los, double alpha_perp, double alpha_parallel,
			       double f, double b1, double b2, double bG2, double bGamma3, double Sigma2, double dSigma2) {

   double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
	double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
	double k_M_p = NORM(kvec_M_pvec);
	double kvec_P_pvec[3] = { kvec[0] + pvec[0], kvec[1] + pvec[1], kvec[2] + pvec[2] };
	double k_P_p = NORM(kvec_P_pvec);
	double p = NORM(pvec);
	double k = NORM(kvec);

	double P22_nw = 0.0;
	double P22_tot = 0.0;

	if (k_M_p > p) {
	    P22_nw += 2.0 * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2)
				* f_pk_no_wiggle(k_M_p) * f_pk_no_wiggle(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk_no_wiggle(p) * f_pk_no_wiggle(p);
	    P22_tot += 2.0 * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2)
				* f_pk(k_M_p) * f_pk(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk(p) * f_pk(p);
	}

	if (k_P_p > p) {
	    P22_nw += 2.0 * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2)
				* f_pk_no_wiggle(k_P_p) * f_pk_no_wiggle(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk_no_wiggle(p) * f_pk_no_wiggle(p);
	    P22_tot += 2.0 * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2)
				* f_pk(k_P_p) * f_pk(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk(p) * f_pk(p);
	}

	double P13_nw = 6.0 * Z1_Bias(kvec, los, f, b1)
	           	* Z3_Bias_G(pvec, M_pvec, kvec, los, f, b1, b2, bG2, bGamma3)
		   		* f_pk_no_wiggle(p) * f_pk_no_wiggle(k);
    double P13_tot = 6.0 * Z1_Bias(kvec, los, f, b1)
	           	* Z3_Bias_G(pvec, M_pvec, kvec, los, f, b1, b2, bG2, bGamma3)
		   		* f_pk(p) * f_pk(k);

	double P1loop_nw = P22_nw + P13_nw;
	double P1loop_tot = P22_tot + P13_tot;
    double P1loop = P1loop_nw + ExpDamping_Ivanov(kvec, los, f, Sigma2, dSigma2) * (P1loop_tot - P1loop_nw);

	return P1loop / alpha3;
}

double Powerspectrum_1loop_Damping_NW( double * kvec_in, double * pvec,
		               double * los, double alpha_perp, double alpha_parallel,
			       double f, double b1, double b2, double bG2,double bGamma3) {

   double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
	double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
	double k_M_p = NORM(kvec_M_pvec);
	double kvec_P_pvec[3] = { kvec[0] + pvec[0], kvec[1] + pvec[1], kvec[2] + pvec[2] };
	double k_P_p = NORM(kvec_P_pvec);
	double p = NORM(pvec);
	double k = NORM(kvec);

	double P22 = 0.0;

	if (k_M_p > p) {
	    P22 += 2.0 * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2)
				* f_pk_no_wiggle(k_M_p) * f_pk_no_wiggle(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk_no_wiggle(p) * f_pk_no_wiggle(p);
	}

	if (k_P_p > p) {
	    P22 += 2.0 * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2)
				* f_pk_no_wiggle(k_P_p) * f_pk_no_wiggle(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk_no_wiggle(p) * f_pk_no_wiggle(p);
	}

	double P13 = 6.0 * Z1_Bias(kvec, los, f, b1)
	           	* Z3_Bias_G(pvec, M_pvec, kvec, los, f, b1, b2, bG2, bGamma3)
		   		* f_pk_no_wiggle(p) * f_pk_no_wiggle(k);

	double P1loop = P22 + P13;

	return P1loop / alpha3;
}

double Powerspectrum_1loop_Damping_tot

( double * kvec_in, double * pvec,
		               double * los, double alpha_perp, double alpha_parallel,
			       double f, double b1, double b2, double bG2, double bGamma3) {

   double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
	double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
	double k_M_p = NORM(kvec_M_pvec);
	double kvec_P_pvec[3] = { kvec[0] + pvec[0], kvec[1] + pvec[1], kvec[2] + pvec[2] };
	double k_P_p = NORM(kvec_P_pvec);
	double p = NORM(pvec);
	double k = NORM(kvec);

	double P22 = 0.0;

	if (k_M_p > p) {
	    P22 += 2.0 * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(pvec, kvec_M_pvec, los, f, b1, b2, bG2)
				* f_pk(k_M_p) * f_pk(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk(p) * f_pk(p);
	}

	if (k_P_p > p) {
	    P22 += 2.0 * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2) * Z2_Bias_G(M_pvec, kvec_P_pvec, los, f, b1, b2, bG2)
				* f_pk(k_P_p) * f_pk(p)
				- 2.0 / 4.0 * b2 * b2 * f_pk(p) * f_pk(p);
	}

	double P13 = 6.0 * Z1_Bias(kvec, los, f, b1)
	           	* Z3_Bias_G(pvec, M_pvec, kvec, los, f, b1, b2, bG2, bGamma3)
		   		* f_pk(p) * f_pk(k);

	double P1loop = P22 + P13;

	return P1loop / alpha3;
}

double Powerspectrum_NonLinearFitting(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double sigma2_perp, double sigma2_para) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
    double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
    double BAO = f_pk(k) - f_pk_no_wiggle(k);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
    double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);
	double P_lin = (G+MC);

	double k_in = NORM(kvec_in);
	if ( (k_in > kmin_nl) && (k_in < kmax_nl) ) {
	    double mu = MU(kvec, los);
	    double P_nl = f_pk_nl_0(k) + f_pk_nl_2(k) * (3.0 * mu * mu - 1.0) / 2.0;
	    return P_nl / alpha3;
	} else {
	    return P_lin / alpha3;
	}
}


double Powerspectrum_NonLinearFitting_Window(double * kvec_in, double * pvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double sigma2_perp, double sigma2_para, double volume) {

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double pvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(pvec_in, los, alpha_perp, alpha_parallel, pvec);
	double kvec[3] = {kvec_in[0], kvec_in[1], kvec_in[2]};
	double kvec_M_pvec[3] = {kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2]};
	double kp = NORM(kvec_M_pvec);

	double Pnl = Powerspectrum_NonLinearFitting(pvec, los, alpha_perp, alpha_parallel,
	             sigma8, f, b1, sigma2_perp, sigma2_para);
	
	return volume * W2(kp, R) * Pnl;
}



std::complex<double> Powerspectrum_Tree_Aniso(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double g20, double g21_real, double g21_imag, double g22_real, double g22_imag) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);

	std::complex<double> g2M[5];
	g2M[0].real(g22_real);
	g2M[0].imag(-1.0 * g22_imag);

	g2M[1].real(-1.0 * g21_real);
	g2M[1].imag(g21_imag);

	g2M[2].real(g20);
	g2M[2].imag(0.0);

	g2M[3].real(g21_real);
	g2M[3].imag(g21_imag);

	g2M[4].real(g22_real);
	g2M[4].imag(g22_imag);

	std::complex<double> LL = 0.0;
	int iM = 0;
	for(int M = -2; M <=2; M++) {
	    std::complex<double> Ylm = calcSphericalHarmonics(2, M, kvec);
	    LL += (g2M[iM]) * (Ylm);
	    iM++;
	}
	std::complex<double> P = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2) * f_pk(k) * (1.0 + LL);
	return P / alpha3;
}

double Powerspectrum_Tree_BAO(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
	                      double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
     	double BAO = f_pk(k) - f_pk_no_wiggle(k);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
     	double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);
	return (G+MC) / alpha3;
}

double Powerspectrum_Tree_BAO_Fitting(
        double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
	    double sigma2_perp, double sigma2_para, 
        double A20, double A11, double A02, 
        double A30, double A21, double A12, double A03,
        double A40, double A31, double A22, double A13, double A04) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	
    double k = NORM(kvec);
    double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
    double BAO = f_pk(k) - f_pk_no_wiggle(k);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
    double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);

	double mu = MU(kvec, los);
    double mu2 = mu * mu;

    double Fit2 =  pow(k,2) * ( A20 * pow(1.0-mu2, 2) + A11 * pow(1.0-mu2, 1) * pow(mu2, 1) + A02 * pow(mu2, 2) );
    double Fit3 =  pow(k,3) * ( A30 * pow(1.0-mu2, 3) + A21 * pow(1.0-mu2, 2) * pow(mu2, 1) + A12 * pow(1.0-mu2, 1) * pow(mu2, 2) + A03 * pow(mu2, 3) );
    double Fit4 =  pow(k,4) * ( A40 * pow(1.0-mu2, 4) + A31 * pow(1.0-mu2, 3) * pow(mu2, 1) + A22 * pow(1.0-mu2, 2) * pow(mu2, 2) + A13 * pow(1.0-mu2, 1) * pow(mu2, 3) + A04 * pow(mu2, 4) );
    
    double MC_Fit = D * D * MC * (Fit2 + Fit3 + Fit4);

	return ( G + MC + MC_Fit ) / alpha3;

}


double Powerspectrum_SPT1loop( double * kvec_in, double * pvec,
		               double * los, double alpha_perp, double alpha_parallel,
			       double sigma8, double f, double b1, double b2, double b3, double bK2,
			       double bK3, double bDK, double bO, double sigma2_perp, double sigma2_para) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double M_pvec[3] = { - pvec[0], - pvec[1], - pvec[2] };
	double kvec_M_pvec[3] = { kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2] };
	double k_M_p = NORM(kvec_M_pvec);
	double p = NORM(pvec);
	double k = NORM(kvec);

	double P22 = 2.0 * Z2_Bias(kvec_M_pvec, pvec, los, f, b1, b2, bK2)
		   * Z2_Bias(kvec_M_pvec, pvec, los, f, b1, b2, bK2) * pow(sigma8, 4)
		   * f_pk_no_wiggle(p) * f_pk_no_wiggle(k_M_p);

//	double mu = MU(kvec, los);
//	double mu2 = mu * mu;
//	double lnD2 = k * k * ( (1.0 - mu2) * sigma2_perp + mu2 * sigma2_para );
//	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
//     	double P22 =  Kaiser * lnD2 * f_pk_no_wiggle(k);

	double P13 = 6.0 * Z1_Bias(kvec, los, f, b1)
	           * Z3_Bias(kvec, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) * pow(sigma8, 4)
		   * f_pk_no_wiggle(p) * f_pk_no_wiggle(k);

     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);

	double P1loop = D * D * (P22 + P13);

//	double mu = MU(kvec, los);
//	double mu2 = mu * mu;
//	double lnD2 = k * k * ( (1.0 - mu2) * sigma2_perp + mu2 * sigma2_para );
//     	double BAO = f_pk(k) - f_pk_no_wiggle(k);
//	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
//     	double PP = D * D * Kaiser * BAO * lnD2;

//	P1loop = P1loop + PP;
//	P1loop = PP;


	return P1loop / alpha3;
}


double Powerspectrum_SPT2loop( double * kvec_in, double * pvec1, double * pvec2,
		               double * los, double alpha_perp, double alpha_parallel,
			       double sigma8, double f, double b1, double b2, double b3, double bK2,
			       double bK3, double bDK, double bO, double sigma2_perp, double sigma2_para) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double M_pvec1[3] = { - pvec1[0], - pvec1[1], - pvec1[2] };
	double M_pvec2[3] = { - pvec2[0], - pvec2[1], - pvec2[2] };

	double kvec_M_pvec1[3] = { kvec[0] - pvec1[0], kvec[1] - pvec1[1], kvec[2] - pvec1[2] };
	double k_M_p1 = NORM(kvec_M_pvec1);

	double kvec_M_pvec12[3] = { kvec[0] - pvec1[0] - pvec2[0], kvec[1] - pvec1[1] - pvec2[1], kvec[2] - pvec1[2] - pvec2[2] };
	double k_M_p12 = NORM(kvec_M_pvec12);


	double p1 = NORM(pvec1);
	double p2 = NORM(pvec2);

	double k = NORM(kvec);


//	double P22 = 2.0 * Z2_Bias(kvec_M_pvec1, pvec1, los, f, b1, b2, bK2)
//		   * Z2_Bias(kvec_M_pvec1, pvec1, los, f, b1, b2, bK2) * pow(sigma8, 4)
//		   * f_pk_no_wiggle(p1) * f_pk_no_wiggle(k_M_p1);
//	
//	double P13 = 6.0 * Z1_Bias(kvec, los, f, b1)
//	           * Z3_Bias(kvec, pvec1, M_pvec1, los, f, b1, b2, b3, bK2, bK3, bDK, bO) * pow(sigma8, 4)
//		   * f_pk_no_wiggle(p1) * f_pk_no_wiggle(k);

	double P15 = 30.0 * Z1_Bias(kvec, los, f, b1)
	           * Z5_Bias(kvec, pvec1, M_pvec1, pvec2, M_pvec2, los, f, b1, b2, b3, bK2, bK3, bDK, bO) * pow(sigma8, 6)
		   * f_pk_no_wiggle(p1)* f_pk_no_wiggle(p2) * f_pk_no_wiggle(k);

	double P24 = 24.0 * Z2_Bias(kvec_M_pvec1, pvec1, los, f, b1, b2, bK2)
	           * Z4_Bias(kvec_M_pvec1, pvec1, pvec2, M_pvec2, los, f, b1, b2, b3, bK2, bK3, bDK, bO) * pow(sigma8, 6)
		   * f_pk_no_wiggle(k_M_p1)* f_pk_no_wiggle(p2) * f_pk_no_wiggle(p1);

	double P33_b = 6.0 * Z3_Bias(kvec_M_pvec12, pvec1, pvec2, los, f, b1, b2, b3, bK2, bK3, bDK, bO)
	             * Z3_Bias(kvec_M_pvec12, pvec1, pvec2, los, f, b1, b2, b3, bK2, bK3, bDK, bO) * pow(sigma8, 6)
		     * f_pk_no_wiggle(p1)* f_pk_no_wiggle(p2) * f_pk_no_wiggle(k_M_p12);

	double P33_a = 9.0 
	           * Z3_Bias(kvec, pvec1, M_pvec1, los, f, b1, b2, b3, bK2, bK3, bDK, bO)
	           * Z3_Bias(kvec, pvec2, M_pvec2, los, f, b1, b2, b3, bK2, bK3, bDK, bO) * pow(sigma8, 6)
		   * f_pk_no_wiggle(p1) *f_pk_no_wiggle(p2) *  f_pk_no_wiggle(k);

     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);

	double P2loop = D * D * (P15+P24+P33_b+P33_a);

	return P2loop / alpha3;
}

double Powerspectrum_Tree_A_perp(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
	                         double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double mu = MU(kvec, los);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
     	double A_perp = D * D * Kaiser * f_pk_no_wiggle(k) * (k*k) * (1.0 - mu * mu) * sigma2_perp / 2.0;
	return A_perp / alpha3;
}

double Powerspectrum_Tree_A_para(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
	                         double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double mu = MU(kvec, los);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
     	double A_para = D * D * Kaiser * f_pk_no_wiggle(k) * (k*k) * (mu * mu) * sigma2_para / 2.0;
	return A_para / alpha3;
}

double Powerspectrum_Tree_B_perp(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
	                         double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double mu = MU(kvec, los);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
     	double A_perp = D * D * Kaiser * f_pk_no_wiggle(k)
	              * pow((k*k) * (1.0 - mu * mu) * sigma2_perp / 2.0,2) /  2.0;
	return A_perp / alpha3;
}

double Powerspectrum_Tree_B_para(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
	                         double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double mu = MU(kvec, los);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
     	double A_para = D * D * Kaiser * f_pk_no_wiggle(k) * pow((k*k) * (mu * mu) * sigma2_para / 2.0, 2) / 2.0;
	return A_para / alpha3;
}

double Powerspectrum_Tree_B_perp_para(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
	                         double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double mu = MU(kvec, los);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
     	double A_para = D * D * Kaiser * f_pk_no_wiggle(k) * (k*k) * (1.0 - mu * mu) * (k*k) * (mu * mu) * sigma2_para * sigma2_perp / 4.0 / 2.0;
	return A_para / alpha3;
}

double Powerspectrum_Tree_NoWiggle(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double Kaiser = Z1_Bias(kvec, los, f, b1) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2);
	double P = Kaiser * f_pk_no_wiggle(k);
	return P / alpha3;
}

double Powerspectrum_Tree_KSZ(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double aH_tau_T0_over_c) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);

	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double P = 2.0 * aH_tau_T0_over_c * f * mu * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2) * f_pk(k) / k;
	return P / alpha3;
}

double Powerspectrum_Tree_KSZ_VV(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double aH_tau_T0_over_c1, double aH_tau_T0_over_c2) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);

	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double f2 = f * f;
	double mu2 = mu * mu;
	double k2 = k * k;

	b1 = 0.0;
	double P = aH_tau_T0_over_c1 * aH_tau_T0_over_c2 * f2 * mu2 * pow(sigma8, 2) * f_pk(k) / k2;
	return P / alpha3;
}

double Powerspectrum_Tree_Window(double * kvec_in, double * pvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double volume) {
	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double pvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(pvec_in, los, alpha_perp, alpha_parallel, pvec);
	double kvec[3] = {kvec_in[0], kvec_in[1], kvec_in[2]};
	double kvec_M_pvec[3] = {kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2]};
	double p = NORM(pvec);
	double kp = NORM(kvec_M_pvec);
	double P = Z1_Bias(pvec, los, f, b1) * Z1_Bias(pvec, los, f, b1) * pow(sigma8, 2) * f_pk(p) / alpha3;
	return volume * W2(kp, R) * P;
}

double Powerspectrum_Tree_Window_IC(double * kvec_in, double * pvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double volume) {
	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double pvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(pvec_in, los, alpha_perp, alpha_parallel, pvec);
	double kvec[3] = {kvec_in[0], kvec_in[1], kvec_in[2]};
	double kvec_M_pvec[3] = {kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2]};
	double k = NORM(kvec);
	double p = NORM(pvec);
	double kp = NORM(kvec_M_pvec);
	double P = Z1_Bias(pvec, los, f, b1) * Z1_Bias(pvec, los, f, b1) * pow(sigma8, 2) * f_pk(p) / alpha3;
//	return volume * ( W2(kp, R) - 2.0 * W1(k,R) * W1(kp, R) * W1(p, R) + W2(k,R) * W2(p,R) ) * P;
	return volume * ( - 2.0 * W1(k,R) * W1(kp, R) * W1(p, R) + W2(k,R) * W2(p,R) ) * P;
}

double Powerspectrum_Tree_Window_IC_Approx(double * kvec_in, double * pvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double volume) {
	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double pvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(pvec_in, los, alpha_perp, alpha_parallel, pvec);
	double kvec[3] = {kvec_in[0], kvec_in[1], kvec_in[2]};
//	double kvec_M_pvec[3] = {kvec[0] - pvec[0], kvec[1] - pvec[1], kvec[2] - pvec[2]};
	double k = NORM(kvec);
	double p = NORM(pvec);
//	double kp = NORM(kvec_M_pvec);
	double P = Z1_Bias(pvec, los, f, b1) * Z1_Bias(pvec, los, f, b1) * pow(sigma8, 2) * f_pk(p) / alpha3;
//	return volume * ( W2(kp, R) - W2(k,R) * W2(p,R) ) * P;
	return volume * ( - W2(k,R) * W2(p,R) ) * P;
}


double Bispectrum_Tree(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double bK2) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_FoG(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double bK2,
			double c1, double c2, double knl) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
        double alpha6 = alpha3 * alpha3;

        double kvec1[3] = { 0.0, 0.0, 0.0 };
        double kvec2[3] = { 0.0, 0.0, 0.0 };
        double kvec3[3] = { 0.0, 0.0, 0.0 };

        calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
        calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
        calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

        double k1 = NORM(kvec1);
        double k2 = NORM(kvec2);
        double k3 = NORM(kvec3);

        double K12 = 2.0 * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) * Z1_Bias_FoG(kvec1, los, f, b1, c1, c2, knl) * Z1_Bias_FoG(kvec2, los, f, b1, c1, c2, knl) * pow(sigma8, 4);
        double K13 = 2.0 * Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2) * Z1_Bias_FoG(kvec1, los, f, b1, c1, c2, knl) * Z1_Bias_FoG(kvec3, los, f, b1, c1, c2, knl) * pow(sigma8, 4);
        double K23 = 2.0 * Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2) * Z1_Bias_FoG(kvec2, los, f, b1, c1, c2, knl) * Z1_Bias_FoG(kvec3, los, f, b1, c1, c2, knl) * pow(sigma8, 4);

        double B = K12 * f_pk(k1) * f_pk(k2)
                 + K13 * f_pk(k1) * f_pk(k3)
                 + K23 * f_pk(k2) * f_pk(k3);

        return B / alpha6;

}

double Bispectrum_Tree_Ivanov_Damping(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double bK2, double Sigma2, double dSigma2) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
        double alpha6 = alpha3 * alpha3;

        double kvec1[3] = { 0.0, 0.0, 0.0 };
        double kvec2[3] = { 0.0, 0.0, 0.0 };
        double kvec3[3] = { 0.0, 0.0, 0.0 };

        calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
        calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
        calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

        double k1 = NORM(kvec1);
        double k2 = NORM(kvec2);
        double k3 = NORM(kvec3);

	double D1 = ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2);
	double D2 = ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2);
	double D3 = ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2);

        double K12 = 2.0 * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
        double K13 = 2.0 * Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
        double K23 = 2.0 * Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
    	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
    	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);

	double GG = D1 * D2 * K12 * BAO1 * BAO2
              + D1 * D3 * K13 * BAO1 * BAO3
              + D2 * D3 * K23 * BAO2 * BAO3;

	double GM_MG = D1 * K12 * BAO1 * f_pk_no_wiggle(k2)
                 + D2 * K12 * f_pk_no_wiggle(k1) * BAO2

                 + D1 * K13 * BAO1 * f_pk_no_wiggle(k3)
                 + D3 * K13 * f_pk_no_wiggle(k1) * BAO3

                 + D2 * K23 * BAO2 * f_pk_no_wiggle(k3)
                 + D3 * K23 * f_pk_no_wiggle(k2) * BAO3;

	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
              + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
              + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

        double B = (GG + GM_MG + NW);

        return B / alpha6;

}

double Bispectrum_Tree_FoG_Damping(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1, double b2, double bG2,
			double c1, double c2, double knl, double Sigma2, double dSigma2) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
        double alpha6 = alpha3 * alpha3;

        double kvec1[3] = { 0.0, 0.0, 0.0 };
        double kvec2[3] = { 0.0, 0.0, 0.0 };
        double kvec3[3] = { 0.0, 0.0, 0.0 };

        calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
        calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
        calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

        double k1 = NORM(kvec1);
        double k2 = NORM(kvec2);
        double k3 = NORM(kvec3);
    
        double kvec12[3] = PLUS(kvec1, kvec2);
        double kvec13[3] = PLUS(kvec1, kvec3);
        double kvec23[3] = PLUS(kvec2, kvec3);

		double D1 = ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2);
		double D2 = ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2);
		double D3 = ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2);
    
        double sqD1 = sqrt(D1);
        double sqD2 = sqrt(D2);
        double sqD3 = sqrt(D3);
    
        double sqD12 = sqrt(ExpDamping_Ivanov(kvec12, los, f, Sigma2, dSigma2));
		double sqD13 = sqrt(ExpDamping_Ivanov(kvec13, los, f, Sigma2, dSigma2));
		double sqD23 = sqrt(ExpDamping_Ivanov(kvec23, los, f, Sigma2, dSigma2));
    
    
    
        double K12 = 2.0 * Z2_Bias_G(kvec1, kvec2, los, f, b1, b2, bG2) * Z1_Bias_FoG(kvec1, los, f, b1, c1, c2, knl) * Z1_Bias_FoG(kvec2, los, f, b1, c1, c2, knl);
        double K13 = 2.0 * Z2_Bias_G(kvec1, kvec3, los, f, b1, b2, bG2) * Z1_Bias_FoG(kvec1, los, f, b1, c1, c2, knl) * Z1_Bias_FoG(kvec3, los, f, b1, c1, c2, knl);
        double K23 = 2.0 * Z2_Bias_G(kvec2, kvec3, los, f, b1, b2, bG2) * Z1_Bias_FoG(kvec2, los, f, b1, c1, c2, knl) * Z1_Bias_FoG(kvec3, los, f, b1, c1, c2, knl);

        double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
    	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
    	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);

		double GG = D1 * D2 * K12 * BAO1 * BAO2 //sqD1 * sqD2 * sqD12 * K12 * BAO1 * BAO2
              + D1 * D3 * K13 * BAO1 * BAO3 //sqD1 * sqD3 * sqD13 * K13 * BAO1 * BAO3
              + D2 * D3 * K23 * BAO2 * BAO3; //sqD2 * sqD3 * sqD23 * K23 * BAO2 * BAO3;

		double GM_MG = D1 * K12 * BAO1 * f_pk_no_wiggle(k2)
                 + D2 * K12 * f_pk_no_wiggle(k1) * BAO2

                 + D1 * K13 * BAO1 * f_pk_no_wiggle(k3)
                 + D3 * K13 * f_pk_no_wiggle(k1) * BAO3

                 + D2 * K23 * BAO2 * f_pk_no_wiggle(k3)
                 + D3 * K23 * f_pk_no_wiggle(k2) * BAO3;

		double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
              + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
              + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

        double B = (GG + GM_MG + NW);

        return B / alpha6;

}

double Bispectrum_SN_FoG_Damping(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1,
			double c1, double c2, double knl,  double Pshot, double Bshot, double Sigma2, double dSigma2) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
        double alpha6 = alpha3 * alpha3;

        double kvec1[3] = { 0.0, 0.0, 0.0 };
        double kvec2[3] = { 0.0, 0.0, 0.0 };
        double kvec3[3] = { 0.0, 0.0, 0.0 };

        calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
        calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
        calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

        double k1 = NORM(kvec1);
        double k2 = NORM(kvec2);
        double k3 = NORM(kvec3);

		double D1 = ExpDamping_Ivanov(kvec1, los, f, Sigma2, dSigma2);
		double D2 = ExpDamping_Ivanov(kvec2, los, f, Sigma2, dSigma2);
		double D3 = ExpDamping_Ivanov(kvec3, los, f, Sigma2, dSigma2);

        double K1 = SN_bispectrum(kvec1, los, b1, f, Pshot, Bshot) * Z1_Bias_FoG(kvec1, los, f, b1, c1, c2, knl);
        double K2 = SN_bispectrum(kvec2, los, b1, f, Pshot, Bshot) * Z1_Bias_FoG(kvec2, los, f, b1, c1, c2, knl);
        double K3 = SN_bispectrum(kvec3, los, b1, f, Pshot, Bshot) * Z1_Bias_FoG(kvec3, los, f, b1, c1, c2, knl);

        double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
    	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
    	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);

		double B = K1 * (f_pk_no_wiggle(k1) + D1 * BAO1) 
				+ K2 * (f_pk_no_wiggle(k2) + D2 * BAO2)
				+ K3 * (f_pk_no_wiggle(k3) + D3 * BAO3);

        return B / alpha6;

}

double Primordial_Matter_Bispectrum_Local_IRresum(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double f, double Sigma2, double dSigma2) {

        double k1 = NORM(kvec1_in);
        double k2 = NORM(kvec2_in);
        double k3 = NORM(kvec3_in);

        double M1 = 0.0;
        double M2 = 0.0;
        double M3 = 0.0;

        double pk1 = f_pk_no_wiggle(k1) + ( f_pk(k1)-f_pk_no_wiggle(k1) ) * ExpDamping_Ivanov(kvec1_in, los, f, Sigma2, dSigma2);
        double pk2 = f_pk_no_wiggle(k2) + ( f_pk(k2)-f_pk_no_wiggle(k2) ) * ExpDamping_Ivanov(kvec2_in, los, f, Sigma2, dSigma2);
        double pk3 = f_pk_no_wiggle(k3) + ( f_pk(k3)-f_pk_no_wiggle(k3) ) * ExpDamping_Ivanov(kvec3_in, los, f, Sigma2, dSigma2);

        if (f_pk_pri(k1) < 1.0e-20) {
            M1 = 0.0;
        } else {
            M1 = sqrt( pk1 / f_pk_pri(k1) );
        }
        if (f_pk_pri(k2) < 1.0e-20) {
            M2 = 0.0;
        } else {
            M2 = sqrt( pk2 / f_pk_pri(k2) );
        }
        if (f_pk_pri(k3) < 1.0e-20) {
            M3 = 0.0;
        } else {
            M3 = sqrt( pk3 / f_pk_pri(k3) );
        }

        double B = 2.0 * ( f_pk_pri(k1) * f_pk_pri(k2) + f_pk_pri(k1) * f_pk_pri(k3) + f_pk_pri(k2) * f_pk_pri(k3) );

        return M1 * M2 * M3 * B;
}

double Primordial_Matter_Bispectrum_Equilateral_IRresum(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double f, double Sigma2, double dSigma2) {

        double k1 = NORM(kvec1_in);
        double k2 = NORM(kvec2_in);
        double k3 = NORM(kvec3_in);

        double pk1 = f_pk_no_wiggle(k1) + ( f_pk(k1)-f_pk_no_wiggle(k1) ) * ExpDamping_Ivanov(kvec1_in, los, f, Sigma2, dSigma2);
        double pk2 = f_pk_no_wiggle(k2) + ( f_pk(k2)-f_pk_no_wiggle(k2) ) * ExpDamping_Ivanov(kvec2_in, los, f, Sigma2, dSigma2);
        double pk3 = f_pk_no_wiggle(k3) + ( f_pk(k3)-f_pk_no_wiggle(k3) ) * ExpDamping_Ivanov(kvec3_in, los, f, Sigma2, dSigma2);

        double M1 = 0.0;
        double M2 = 0.0;
        double M3 = 0.0;

        if (f_pk_pri(k1) < 1.0e-20) {
            M1 = 0.0;
        } else {
            M1 = sqrt( pk1 / f_pk_pri(k1) );
        }
        if (f_pk_pri(k2) < 1.0e-20) {
            M2 = 0.0;
        } else {
            M2 = sqrt( pk2 / f_pk_pri(k2) );
        }
        if (f_pk_pri(k3) < 1.0e-20) {
            M3 = 0.0;
        } else {
            M3 = sqrt( pk3 / f_pk_pri(k3) );
        }

        double B1 = - ( f_pk_pri(k1) * f_pk_pri(k2) + f_pk_pri(k1) * f_pk_pri(k3) + f_pk_pri(k2) * f_pk_pri(k3) );
        double B2 = - 2.0 * pow( f_pk_pri(k1) * f_pk_pri(k2) * f_pk_pri(k3), 2.0/3.0 );
        double B3 = pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k3)
                  + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k3)
                  + pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k2)
                  + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k2)
                  + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k1)
                  + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k1);
        double B = 6.0 * (B1 + B2 + B3);

        return M1 * M2 * M3 * B;
}

double Primordial_Matter_Bispectrum_Orthogonal_IRresum(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double f, double Sigma2, double dSigma2) {

        double k1 = NORM(kvec1_in);
        double k2 = NORM(kvec2_in);
        double k3 = NORM(kvec3_in);

        double pk1 = f_pk_no_wiggle(k1) + ( f_pk(k1)-f_pk_no_wiggle(k1) ) * ExpDamping_Ivanov(kvec1_in, los, f, Sigma2, dSigma2);
        double pk2 = f_pk_no_wiggle(k2) + ( f_pk(k2)-f_pk_no_wiggle(k2) ) * ExpDamping_Ivanov(kvec2_in, los, f, Sigma2, dSigma2);
        double pk3 = f_pk_no_wiggle(k3) + ( f_pk(k3)-f_pk_no_wiggle(k3) ) * ExpDamping_Ivanov(kvec3_in, los, f, Sigma2, dSigma2);

        double M1 = 0.0;
        double M2 = 0.0;
        double M3 = 0.0;

        if (f_pk_pri(k1) < 1.0e-20) {
            M1 = 0.0;
        } else {
            M1 = sqrt( pk1 / f_pk_pri(k1) );
        }
        if (f_pk_pri(k2) < 1.0e-20) {
            M2 = 0.0;
        } else {
            M2 = sqrt( pk2 / f_pk_pri(k2) );
        }
        if (f_pk_pri(k3) < 1.0e-20) {
            M3 = 0.0;
        } else {
            M3 = sqrt( pk3 / f_pk_pri(k3) );
        }

        double B1 = - 3.0 * ( f_pk_pri(k1) * f_pk_pri(k2) + f_pk_pri(k1) * f_pk_pri(k3) + f_pk_pri(k2) * f_pk_pri(k3) );
        double B2 = - 8.0 * pow( f_pk_pri(k1) * f_pk_pri(k2) * f_pk_pri(k3), 2.0/3.0 );
        double B3 = 3.0 * ( pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k3)
                              + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k3)
                              + pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k2)
                              + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k2)
                              + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k1)
                              + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k1) );
        double B = 6.0 * (B1 + B2 + B3);

        return M1 * M2 * M3 * B;
}

double Primordial_Matter_Bispectrum_Orthogonal_LSS_IRresum(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double f, double Sigma2, double dSigma2) {

        double k1 = NORM(kvec1_in);
        double k2 = NORM(kvec2_in);
        double k3 = NORM(kvec3_in);

        double pk1 = f_pk_no_wiggle(k1) + ( f_pk(k1)-f_pk_no_wiggle(k1) ) * ExpDamping_Ivanov(kvec1_in, los, f, Sigma2, dSigma2);
        double pk2 = f_pk_no_wiggle(k2) + ( f_pk(k2)-f_pk_no_wiggle(k2) ) * ExpDamping_Ivanov(kvec2_in, los, f, Sigma2, dSigma2);
        double pk3 = f_pk_no_wiggle(k3) + ( f_pk(k3)-f_pk_no_wiggle(k3) ) * ExpDamping_Ivanov(kvec3_in, los, f, Sigma2, dSigma2);

        double M1 = 0.0;
        double M2 = 0.0;
        double M3 = 0.0;

        if (f_pk_pri(k1) < 1.0e-20) {
            M1 = 0.0;
        } else {
            M1 = sqrt( pk1 / f_pk_pri(k1) );
        }
        if (f_pk_pri(k2) < 1.0e-20) {
            M2 = 0.0;
        } else {
            M2 = sqrt( pk2 / f_pk_pri(k2) );
        }
        if (f_pk_pri(k3) < 1.0e-20) {
            M3 = 0.0;
        } else {
            M3 = sqrt( pk3 / f_pk_pri(k3) );
        }

        double k1sk2 = k1 / k2;
        double k1sk3 = k1 / k3;
        double k2sk3 = k2 / k3;

        double k1_2sk2k3 = k1sk2 * k1sk3;
        double k2_2sk1k3 = k2sk3 / k1sk2;
        double k3_2sk1k2 = 1 / k1sk3 / k2sk3;

        double p = 1 / (-21 + 743 / (7 * (20 * M_PI * M_PI - 193)));

        double cst = - (2 + 60 * p);

        double K2perm = k1_2sk2k3 * (- 1 - 9 * p + k1_2sk2k3 * p) - 1 / k1_2sk2k3 * 20 * p
        				+ k2_2sk1k3 * (- 1 - 9 * p + k2_2sk1k3 * p) - 1 / k2_2sk1k3 * 20 * p
        				+ k3_2sk1k2 * (- 1 - 9 * p + k3_2sk1k2 * p) - 1 / k3_2sk1k2 * 20 * p;

        double K5perm = k1sk2 * ( 1 + 15 * p - 6 * p * k1sk3 * k1sk3 + 15 * p * k1sk2)
        				+ k1sk3 * ( 1 + 15 * p - 6 * p * k1sk2 * k1sk2 + 15 * p * k1sk3)
        				+ k2sk3 * ( 1 + 15 * p - 6 * p / k1sk2 / k1sk2 + 15 * p * k2sk3)
        				+ 1 / k1sk2 * ( 1 + 15 * p - 6 * p * k2sk3 * k2sk3 + 15 * p / k1sk2)
        				+ 1 / k1sk3 * ( 1 + 15 * p - 6 * p / k2sk3 / k2sk3 + 15 * p / k1sk3)
        				+ 1 / k2sk3 * ( 1 + 15 * p - 6 * p / k1sk3 / k1sk3 + 15 * p / k2sk3);

        double B = 6 * pow( f_pk_pri(k1) * f_pk_pri(k2) * f_pk_pri(k3), 2.0/3.0 ) * (cst + K2perm + K5perm);

        return M1 * M2 * M3 * B;
}

double Bispectrum_PNG_Damping(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double f, double b1,
			double fnlloc, double fnlequi, double fnlortho, double Sigma2, double dSigma2) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
        double alpha6 = alpha3 * alpha3;

        double kvec1[3] = { 0.0, 0.0, 0.0 };
        double kvec2[3] = { 0.0, 0.0, 0.0 };
        double kvec3[3] = { 0.0, 0.0, 0.0 };

        calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
        calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
        calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

        double Bloc = Primordial_Matter_Bispectrum_Local_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
		double Bequi = Primordial_Matter_Bispectrum_Equilateral_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
		double Bortho = Primordial_Matter_Bispectrum_Orthogonal_LSS_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

		double B = Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * (fnlloc * Bloc + fnlequi * Bequi + fnlortho * Bortho);

        return B / alpha6;

}

double Bispectrum_Kernel_b1_b1_b1(double * kvec1, double * kvec2) {
	
	double K = 2.0 * F2(kvec1, kvec2) * D1() * D1();

    return K;

}

double Bispectrum_Kernel_b1_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * F2(kvec1, kvec2) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1()) + 2.0 * V2(kvec1, kvec2, los) * D1() * D1();

	K *= f;

    return K;

}

double Bispectrum_Kernel_b1_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * F2(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los) + 2.0 * V2(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );

	K *= f * f;
    
    return K;

}

/*********/
double Bispectrum_Kernel_b1_b1_b2(double * kvec1, double * kvec2) {

	double K = 2.0 * (D1D1(kvec1, kvec2)/2.0) * D1() * D1();

    return K;

}

double Bispectrum_Kernel_b1_b2_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (D1D1(kvec1, kvec2)/2.0) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1());

	K *= f;

    return K;

}

double Bispectrum_Kernel_b2_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (D1D1(kvec1, kvec2)/2.0) * V1(kvec1, los) * V1(kvec2, los);

	K *= f * f;

    return K;

}

/********/

double Bispectrum_Kernel_b1_b1_bG2(double * kvec1, double * kvec2) {

	double K = 2.0 * G1G1(kvec1, kvec2) * D1() * D1();

    return K;

}

double Bispectrum_Kernel_b1_bG2_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * G1G1(kvec1, kvec2) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );

	K *= f;

    return K;

}

double Bispectrum_Kernel_bG2_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * G1G1(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los);

	K *= f * f;

    return K;

}

/************/

double Bispectrum_Kernel_b1_b1_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * D1V1(kvec1, kvec2, los) * D1() * D1();

	K *= f;

    return K;

}

double Bispectrum_Kernel_b1_b1_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * D1V1(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );

	K += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * D1() * D1();

	K *= f * f;

    return K;

}

double Bispectrum_Kernel_b1_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * D1V1(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);

	K += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );

	K *= f * f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_f_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * V1(kvec1, los) * V1(kvec2, los);

	K *= f * f * f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * V2(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);

	K *= f * f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_b1_b1(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * (KFOG(kvec1, los) + KFOG(kvec2, los)) * F2(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c1_b1_b2(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * (KFOG(kvec1, los) + KFOG(kvec2, los)) * D1D1(kvec1, kvec2)/2.;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_b1_bG2(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * (KFOG(kvec1, los) + KFOG(kvec2, los)) * G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c1_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (KFOG(kvec1, los) * V1(kvec2, los) + KFOG(kvec2, los) * V1(kvec1, los)) * F2(kvec1,kvec2);

	K += 2.0 * (KFOG(kvec1, los) + KFOG(kvec2, los)) * V2(kvec1, kvec2, los);

	K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_b1_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (KFOG(kvec1, los) + KFOG(kvec2, los)) * D1V1(kvec1, kvec2, los);

	K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_b1_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (KFOG(kvec1, los) + KFOG(kvec2, los)) * (V1V1(kvec1, kvec2, los)/2.0);

	K += 2.0 * (KFOG(kvec1, los) * V1(kvec2, los) + KFOG(kvec2, los) * V1(kvec1, los)) * D1V1(kvec1,kvec2,los);

	K *= f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_b2_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (KFOG(kvec1, los) * V1(kvec2, los) + KFOG(kvec2, los) * V1(kvec1, los)) / 2.0;

	K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_bG2_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (KFOG(kvec1, los) * V1(kvec2, los) + KFOG(kvec2, los) * V1(kvec1, los)) * G1G1(kvec1, kvec2);

	K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (KFOG(kvec1, los) * V1(kvec2, los) + KFOG(kvec2, los) * V1(kvec1, los)) * V2(kvec1, kvec2, los);

	K *= f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * (KFOG(kvec1, los) * V1(kvec2, los) + KFOG(kvec2, los) * V1(kvec1, los)) * V1V1(kvec1, kvec2, los)/2;

	K *= f * f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c1_b1(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * KFOG(kvec1, los) * KFOG(kvec2, los) * F2(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c1_b2(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * KFOG(kvec1, los) * KFOG(kvec2, los) / 2.0;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c1_bG2(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * KFOG(kvec1, los) * KFOG(kvec2, los) * G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c1_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * KFOG(kvec1, los) * KFOG(kvec2, los) * V2(kvec1, kvec2, los);

	K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c1_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * KFOG(kvec1, los) * KFOG(kvec2, los) * D1V1(kvec1, kvec2, los);

	K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c1_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * KFOG(kvec1, los) * KFOG(kvec2, los) * V1V1(kvec1, kvec2, los)/2;

	K *= f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_c2_b1_b1(double * kvec1, double * kvec2, double * los) {

        double K = 2.0 * (KFOG2(kvec1, los) + KFOG2(kvec2, los)) * F2(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c2_b1_b2(double * kvec1, double * kvec2, double * los) {

        double K = 2.0 * (KFOG2(kvec1, los) + KFOG2(kvec2, los)) * D1D1(kvec1, kvec2)/2.;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_b1_bG2(double * kvec1, double * kvec2, double * los) {

        double K = 2.0 * (KFOG2(kvec1, los) + KFOG2(kvec2, los)) * G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c2_b1_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG2(kvec1, los) * V1(kvec2, los) + KFOG2(kvec2, los) * V1(kvec1, los)) * F2(kvec1,kvec2);

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_b1_b1_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG2(kvec1, los) + KFOG2(kvec2, los)) * D1V1(kvec1, kvec2, los);

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_b1_f_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG2(kvec1, los) + KFOG2(kvec2, los)) * (V1V1(kvec1, kvec2, los)/2.0);

        K += 2.0 * (KFOG2(kvec1, los) * V1(kvec2, los) + KFOG2(kvec2, los) * V1(kvec1, los)) * D1V1(kvec1,kvec2,los);

		K *= f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_b2_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG2(kvec1, los) * V1(kvec2, los) + KFOG2(kvec2, los) * V1(kvec1, los)) / 2.0;

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_bG2_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG2(kvec1, los) * V1(kvec2, los) + KFOG2(kvec2, los) * V1(kvec1, los)) * G1G1(kvec1, kvec2);

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_f_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG2(kvec1, los) * V1(kvec2, los) + KFOG2(kvec2, los) * V1(kvec1, los)) * V2(kvec1, kvec2, los);

		K *= f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG2(kvec1, los) * V1(kvec2, los) + KFOG2(kvec2, los) * V1(kvec1, los)) * V1V1(kvec1, kvec2, los)/2;

		K *= f * f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c2_b1(double * kvec1, double * kvec2, double * los) {

        double K = 2.0 * (KFOG(kvec1, los) * KFOG2(kvec2, los) + KFOG2(kvec1, los) * KFOG(kvec2, los)) * F2(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c2_b2(double * kvec1, double * kvec2, double * los) {

        double K = 2.0 * (KFOG(kvec1, los) * KFOG2(kvec2, los) + KFOG2(kvec1, los) * KFOG(kvec2, los)) / 2.0;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c2_bG2(double * kvec1, double * kvec2, double * los) {

        double K = 2.0 * (KFOG(kvec1, los) * KFOG2(kvec2, los) + KFOG2(kvec1, los) * KFOG(kvec2, los)) * G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c2_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG(kvec1, los) * KFOG2(kvec2, los) + KFOG2(kvec1, los) * KFOG(kvec2, los)) * V2(kvec1, kvec2, los);

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c2_b1_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG(kvec1, los) * KFOG2(kvec2, los) + KFOG2(kvec1, los) * KFOG(kvec2, los)) * D1V1(kvec1, kvec2, los);

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c1_c2_f_f(double * kvec1, double * kvec2, double * los, double f) {

        double K = 2.0 * (KFOG(kvec1, los) * KFOG2(kvec2, los) + KFOG2(kvec1, los) * KFOG(kvec2, los)) * V1V1(kvec1, kvec2, los)/2;

		K *= f * f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_c2_b1(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * KFOG2(kvec1, los) * KFOG2(kvec2, los) * F2(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c2_c2_b2(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * KFOG2(kvec1, los) * KFOG2(kvec2, los) / 2.0;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_c2_bG2(double * kvec1, double * kvec2, double * los) {

	double K = 2.0 * KFOG2(kvec1, los) * KFOG2(kvec2, los) * G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Bispectrum_Kernel_c2_c2_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * KFOG2(kvec1, los) * KFOG2(kvec2, los) * V2(kvec1, kvec2, los);

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_c2_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * KFOG2(kvec1, los) * KFOG2(kvec2, los) * D1V1(kvec1, kvec2, los);

		K *= f;

    return K;

}

/*********/

double Bispectrum_Kernel_c2_c2_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double K = 2.0 * KFOG2(kvec1, los) * KFOG2(kvec2, los) * V1V1(kvec1, kvec2, los)/2;

		K *= f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_Bshot_b1_b1() {

	double K = 1.0;
	
    return K;

}

/*********/

double Bispectrum_Kernel_Bshot_b1_f(double * kvec1, double * los, double f) {

	double K = V1(kvec1, los);

	K *= f;
	
    return K;

}
/*********/

double Bispectrum_Kernel_Bshot_b1_c1(double * kvec1, double * los) {

	double K = KFOG(kvec1, los);
	
    return K;

}
/*********/

double Bispectrum_Kernel_Bshot_b1_c2(double * kvec1, double * los) {

	double K = KFOG2(kvec1, los);
	
    return K;

}
/*********/

double Bispectrum_Kernel_Pshot_f_b1(double * kvec1, double * los, double f) {

	double K = 2.0 * V1(kvec1, los);

	K *= f;
	
    return K;

}
/*********/

double Bispectrum_Kernel_Pshot_f_f(double * kvec1, double * los, double f) {

	double K = 2.0 * V1(kvec1, los) * V1(kvec1, los);

	K *= f * f;
	
    return K;

}
/*********/

double Bispectrum_Kernel_Pshot_f_c1(double * kvec1, double * los, double f) {

	double K = 2.0 * V1(kvec1, los) * KFOG(kvec1, los);

	K *= f;
	
    return K;

}

double Bispectrum_Kernel_Pshot_f_c2(double * kvec1, double * los, double f) {

	double K = 2.0 * V1(kvec1, los) * KFOG2(kvec1, los);

	K *= f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlloc_b1_b1_b1(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = Primordial_Matter_Bispectrum_Local_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlloc_b1_b1_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = (V1(kvec1, los) + V1(kvec2, los) + V1(kvec3, los)) * Primordial_Matter_Bispectrum_Local_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlloc_b1_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = (V1(kvec1, los) * V1(kvec2, los) + V1(kvec1, los) * V1(kvec3, los) + V1(kvec2, los) * V1(kvec3, los)) * Primordial_Matter_Bispectrum_Local_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlloc_f_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = V1(kvec1, los) * V1(kvec2, los) * V1(kvec3, los) * Primordial_Matter_Bispectrum_Local_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlequi_b1_b1_b1(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = - Primordial_Matter_Bispectrum_Equilateral_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlequi_b1_b1_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = - (V1(kvec1, los) + V1(kvec2, los) + V1(kvec3, los)) * Primordial_Matter_Bispectrum_Equilateral_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlequi_b1_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = - (V1(kvec1, los) * V1(kvec2, los) + V1(kvec1, los) * V1(kvec3, los) + V1(kvec2, los) * V1(kvec3, los)) * Primordial_Matter_Bispectrum_Equilateral_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlequi_f_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = - V1(kvec1, los) * V1(kvec2, los) * V1(kvec3, los) * Primordial_Matter_Bispectrum_Equilateral_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_b1_b1_b1(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = Primordial_Matter_Bispectrum_Orthogonal_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_b1_b1_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = (V1(kvec1, los) + V1(kvec2, los) + V1(kvec3, los)) * Primordial_Matter_Bispectrum_Orthogonal_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_b1_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = (V1(kvec1, los) * V1(kvec2, los) + V1(kvec1, los) * V1(kvec3, los) + V1(kvec2, los) * V1(kvec3, los)) * Primordial_Matter_Bispectrum_Orthogonal_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_f_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = V1(kvec1, los) * V1(kvec2, los) * V1(kvec3, los) * Primordial_Matter_Bispectrum_Orthogonal_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_LSS_b1_b1_b1(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = Primordial_Matter_Bispectrum_Orthogonal_LSS_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_LSS_b1_b1_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = (V1(kvec1, los) + V1(kvec2, los) + V1(kvec3, los)) * Primordial_Matter_Bispectrum_Orthogonal_LSS_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_LSS_b1_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = (V1(kvec1, los) * V1(kvec2, los) + V1(kvec1, los) * V1(kvec3, los) + V1(kvec2, los) * V1(kvec3, los)) * Primordial_Matter_Bispectrum_Orthogonal_LSS_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f;
	
    return K;

}

/*********/

double Bispectrum_Kernel_PNG_fnlortho_LSS_f_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f, double Sigma2, double dSigma2) {

	double K = V1(kvec1, los) * V1(kvec2, los) * V1(kvec3, los) * Primordial_Matter_Bispectrum_Orthogonal_LSS_IRresum(kvec1, kvec2, kvec3, los, f, Sigma2, dSigma2);

	K *= f * f * f;
	
    return K;

}

/*********/
/*********/

double Power_Spectrum_Kernel_Tree_b1_b1() {

	double K = 1.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_Tree_b1_f(double * kvec1, double * los, double f) {

	double K = 2.0 * V1(kvec1, los);

	K *= f;
	
    return K;

}

/*********/

double Power_Spectrum_Kernel_Tree_f_f(double * kvec1, double * los, double f) {

	double K = V1(kvec1, los) * V1(kvec1, los);

	K *= f * f;
	
    return K;

}

/*********/
/*********/

double Power_Spectrum_Kernel_Counterterm_c0_b1(double * kvec) {

	double k = NORM(kvec);

	double K = k * k / 0.3 / 0.3;

    return K;

}

/*********/

double Power_Spectrum_Kernel_Counterterm_c0_f(double * kvec, double * los, double f) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double K = mu * mu * k * k / 0.3 / 0.3;

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_Counterterm_c1_b1(double * kvec, double * los) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double K = mu * mu * k * k / 0.3 / 0.3;

    return K;

}

/*********/

double Power_Spectrum_Kernel_Counterterm_c1_f(double * kvec, double * los, double f) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double K = pow(mu, 4) * k * k / 0.3 / 0.3;

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_Counterterm_c2_b1(double * kvec, double * los) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double K = pow(mu, 4) * k * k / 0.3 / 0.3;

    return K;

}

/*********/

double Power_Spectrum_Kernel_Counterterm_c2_f(double * kvec, double * los, double f) {

	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double K = pow(mu, 6) * k * k / 0.3 / 0.3;

	K *= f;

    return K;

}

/*********/
/*********/

double Power_Spectrum_Kernel_1loop_22_b2_b2() {

	double K = 1.0 / 4.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b2_bG2(double * kvec1, double * kvec2) {

	double K = G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_b2(double * kvec1, double * kvec2) {

	double K = F2(kvec1, kvec2);

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b2_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);

	double mu = MU(kvec, los);

	double K = mu * mu * G2(kvec1, kvec2);

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_b2_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = 1.0/ 2.0 * mu * k * (mu1 / k1 + mu2 / k2);

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b2_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = 1.0/ 2.0 * mu * k * (mu1 / k1 * mu2 * mu2 + mu2 / k2 * mu1 * mu1);

	K *= f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_bG2_bG2(double * kvec1, double * kvec2) {

	double K = G1G1(kvec1, kvec2) * G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_bG2(double * kvec1, double * kvec2) {

	double K = 2.0 * F2(kvec1, kvec2) * G1G1(kvec1, kvec2);

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_bG2_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);

	double mu = MU(kvec, los);

	double K = 2.0 * mu * mu * G2(kvec1, kvec2) * G1G1(kvec1, kvec2);

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_bG2_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = mu * k * (mu1 / k1 + mu2 / k2) * G1G1(kvec1, kvec2);

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_bG2_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = mu * k * (mu1 / k1 * mu2 * mu2 + mu2 / k2 * mu1 * mu1) * G1G1(kvec1, kvec2);

	K *= f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_b1(double * kvec1, double * kvec2) {

	double K = F2(kvec1, kvec2) * F2(kvec1, kvec2);

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);

	double mu = MU(kvec, los);

	double K = 2.0 * F2(kvec1, kvec2) * mu * mu * G2(kvec1, kvec2);

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_b1_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = mu * k * (mu1 / k1 + mu2 / k2) * F2(kvec1, kvec2);

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = mu * k * (mu1 / k1 * mu2 * mu2 + mu2 / k2 * mu1 * mu1) * F2(kvec1, kvec2);

	K += mu * k * (mu1 / k1+ mu2 / k2) * G2(kvec1, kvec2);

	K *= f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);

	double mu = MU(kvec, los);

	double K = pow(mu, 4) * G2(kvec1, kvec2) * G2(kvec1, kvec2);

	K *= f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = k * pow(mu, 3) * (mu1 / k1 * mu2 * mu2 + mu2 / k2 * mu1 * mu1) * G2(kvec1, kvec2);

	K *= pow(f, 3);

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_b1_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = mu * mu * k * k * (mu1 / k1 + mu2 / k2) * (mu1 / k1 + mu2 / k2);

	K *= f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_b1_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = mu * mu * k * k * (mu1 / k1 + mu2 / k2) * (mu1 / k1 * k2 * k2 + mu2 / k2 * k1 * k1);

	K *= f * f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_22_f_f_f_f(double * kvec1, double * kvec2, double * los, double f) {

	double kvec[3] = PLUS(kvec1,kvec2);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);

	double K = mu * mu * k * k * (mu1 / k1 * k2 * k2 + mu2 / k2 * k1 * k1) * (mu1 / k1 * k2 * k2 + mu2 / k2 * k1 * k1);

	K *= f * f * f * f;

    return K;

}

/*********/
/*********/

double Power_Spectrum_Kernel_1loop_13_b1_b3() {

	double K = 1.0 / 6.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_bG3(double * kvec1, double * kvec2, double * kvec3) {

	double K = G1G1G1(kvec1, kvec2, kvec3) / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_bG2d(double * kvec1, double * kvec2, double * kvec3) {

	double K = G1G1(kvec1, kvec2) + G1G1(kvec1, kvec3) + G1G1(kvec2, kvec3);

    return K / 3.0;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_bGamma3(double * kvec1, double * kvec2, double * kvec3) {

	double K = 2.0 * (G_F(kvec1,kvec2,kvec3) - G_G(kvec1,kvec2,kvec3)) / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_b1(double * kvec1, double * kvec2, double * kvec3) {

	double K = F3(kvec1,kvec2,kvec3);

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double mu = MU(kvec, los);
	double mu3 = MU(kvec3, los);

	double K = mu * mu * G3(kvec1, kvec2, kvec3);

	K += mu3 * mu3 * F3(kvec1,kvec2,kvec3);

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_b1_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double K = mu * mu * k * k / 2.0 * (DVV(kvec2,kvec3,los) + DVV(kvec1,kvec3,los) + DVV(kvec1,kvec2,los)) / 3.0;

	K *= f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_f_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double mu = MU(kvec, los);
	double mu3 = MU(kvec3, los);

	double K = mu * mu * k * k / 2.0 * (VVV(kvec1,kvec2,kvec3,los) + VVV(kvec2,kvec1,kvec3,los) + VVV(kvec3,kvec1,kvec2,los));

	K += mu3 * mu3 * mu * mu * k * k / 2.0 * (DVV(kvec2,kvec3,los) + DVV(kvec1,kvec3,los) + DVV(kvec1,kvec2,los));

	K *= f * f * f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_b1_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double mu = MU(kvec, los);

	double K = mu * k * DV2(kvec1,kvec2,kvec3,los);

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double mu = MU(kvec, los);
	double mu3 = MU(kvec3, los);

	double K = mu * k * VV2(kvec1,kvec2,kvec3,los);

	K += mu3 * mu3 * mu * k * DV2(kvec1,kvec2,kvec3,los);

	K *= f * f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_b2(double * kvec1, double * kvec2, double * kvec3) {

	double K = F2(kvec1,kvec2) + F2(kvec1,kvec3) + F2(kvec2,kvec3);

    return K / 3.0;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_bG2(double * kvec1, double * kvec2, double * kvec3) {

	double K = 2.0 * G_F(kvec1,kvec2,kvec3) / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_b2_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);
	double mu3 = MU(kvec3, los);

	double K = mu * k / 2.0 * (mu1 / k1 + mu2 / k2 + mu3 / k3);

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b1_bG2_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);
	double mu3 = MU(kvec3, los);

	double K = mu * k * (mu1 / k1 * G1G1(kvec1, kvec2) + mu2 / k2 * G1G1(kvec1, kvec3) + mu3 / k3 * G1G1(kvec2, kvec3));

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b3_f(double * kvec, double * los, double f) {

	double mu = MU(kvec, los);

	double K = mu * mu / 6.0;

	K *= f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_bG3_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 * G1G1G1(kvec1, kvec2, kvec3);

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_bG2d_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 * (G1G1(kvec1, kvec2) + G1G1(kvec1, kvec3) + G1G1(kvec2, kvec3));

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_bGamma3_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 * 2.0 * (G_F(kvec1,kvec2,kvec3) - G_G(kvec1,kvec2,kvec3));

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double mu3 = MU(kvec3, los);
	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double mu = MU(kvec, los);

	double K = mu3 * mu3 * mu * mu * G3(kvec1,kvec2,kvec3);

	K *= f * f;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_f_f_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double mu = MU(kvec, los);
	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 * mu * mu * k * k / 2.0 * (VVV(kvec1,kvec2,kvec3,los) + VVV(kvec2,kvec1,kvec3,los) + VVV(kvec3,kvec1,kvec2,los));

	K *= pow(f,4) / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_f_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double mu = MU(kvec, los);
	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 * mu * k * VV2(kvec1,kvec2,kvec3,los);

	K *= pow(f,3) / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b2_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 * (F2(kvec1,kvec2) + F2(kvec1,kvec3) + F2(kvec2,kvec3));

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_bG2_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double mu3 = MU(kvec3, los);

	double K = 2.0 * mu3 * mu3 * G_F(kvec1,kvec2,kvec3);

	K *= f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_b2_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);
	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 * mu * k / 2.0 * (mu1 / k1 + mu2 / k2 + mu3 / k3);

	K *= f * f / 3.0;

    return K;

}

/*********/

double Power_Spectrum_Kernel_1loop_13_bG2_f_f(double * kvec1, double * kvec2, double * kvec3, double * los, double f) {

	double kvec[3] = PLUS3(kvec1,kvec2,kvec3);
	double k = NORM(kvec);
	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double mu = MU(kvec, los);
	double mu1 = MU(kvec1, los);
	double mu2 = MU(kvec2, los);
	double mu3 = MU(kvec3, los);

	double K = mu3 * mu3 *mu * k * (mu1 / k1 * G1G1(kvec1, kvec2) + mu2 / k2 * G1G1(kvec1, kvec3) + mu3 / k3 * G1G1(kvec2, kvec3));

	K *= f * f / 3.0;

    return K;

}

/*********/
/*********/

double Powerspectrum_LocalMean(double * kvec_in, double * evec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f,  
                               double b1, double b2, double bK2, 
                               double sigma2_perp, double sigma2_para, double nmean, double volume) {

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double kvec[3] = { kvec_in[0], kvec_in[1],  kvec_in[2] };
	double evec[3] = { evec_in[0], evec_in[1],  evec_in[2] };
	double M_kvec_M_evec[3] = { - kvec_in[0] - evec_in[0], - kvec_in[1] - evec_in[1], - kvec_in[2] - evec_in[2] };

    double B = Bispectrum_Tree(kvec, M_kvec_M_evec, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);
    double BN = B + (1.0/nmean) * Powerspectrum_NonLinearFitting(M_kvec_M_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
                  + (1.0/nmean) * Powerspectrum_NonLinearFitting(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para);

//                  + (1.0/nmean) * Powerspectrum_Tree(M_kvec_M_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
//                  + (1.0/nmean) * Powerspectrum_Tree(evec, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double emag = NORM(evec);
	return - 2.0 * W2(emag, R) * BN;

}



double Bispectrum_Tree_b2(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * ( D1D1(kvec1, kvec2) / 2.0 ) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * ( D1D1(kvec1, kvec3) / 2.0 ) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * ( D1D1(kvec2, kvec3) / 2.0 ) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_bK2(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * K1K1(kvec1, kvec3) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * K1K1(kvec2, kvec3) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_RealSpace_DarkMatter_Growth(double * kvec1, double * kvec2, double * kvec3, double sigma8) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2_Growth(kvec1, kvec2) * pow(sigma8, 4);
	double K13 = 2.0 * F2_Growth(kvec1, kvec3) * pow(sigma8, 4);
	double K23 = 2.0 * F2_Growth(kvec2, kvec3) * pow(sigma8, 4);

    double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
             + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
             + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
    return B;

}


double Bispectrum_Tree_NoWiggle_RealSpace_DarkMatter_Shift(double * kvec1, double * kvec2, double * kvec3, double sigma8) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2_Shift(kvec1, kvec2) * pow(sigma8, 4);
	double K13 = 2.0 * F2_Shift(kvec1, kvec3) * pow(sigma8, 4);
	double K23 = 2.0 * F2_Shift(kvec2, kvec3) * pow(sigma8, 4);

    double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
             + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
             + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
    return B;

}

double Bispectrum_Tree_NoWiggle_RealSpace_DarkMatter_Tidal(double * kvec1, double * kvec2, double * kvec3, double sigma8) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2_Tidal(kvec1, kvec2) * pow(sigma8, 4);
	double K13 = 2.0 * F2_Tidal(kvec1, kvec3) * pow(sigma8, 4);
	double K23 = 2.0 * F2_Tidal(kvec2, kvec3) * pow(sigma8, 4);

    double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
             + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
             + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
    return B;

}

double Bispectrum_Tree_BAO_RealSpace_DarkMatter_Growth(double * kvec1, double * kvec2, double * kvec3, double sigma8, double sigma2_perp) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

    double D1 = ExpDamping_RealSpace(kvec1, sigma2_perp);
    double D2 = ExpDamping_RealSpace(kvec2, sigma2_perp);
    double D3 = ExpDamping_RealSpace(kvec3, sigma2_perp);

	double K12 = 2.0 * F2_Growth(kvec1, kvec2) * pow(sigma8, 4);
	double K13 = 2.0 * F2_Growth(kvec1, kvec3) * pow(sigma8, 4);
	double K23 = 2.0 * F2_Growth(kvec2, kvec3) * pow(sigma8, 4);

    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
    double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
    double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
    
    double GG = D1 * D2 * D3 * K12 * BAO1 * BAO2
    	      + D1 * D2 * D3 * K13 * BAO1 * BAO3
    	      + D1 * D2 * D3 * K23 * BAO2 * BAO3;
    
    double GM_MG = D1 * D1 * K12 * BAO1 * f_pk_no_wiggle(k2)
    	         + D2 * D2 * K12 * f_pk_no_wiggle(k1) * BAO2
    
    	         + D1 * D1 * K13 * BAO1 * f_pk_no_wiggle(k3)
    	         + D3 * D3 * K13 * f_pk_no_wiggle(k1) * BAO3
    
    	         + D2 * D2 * K23 * BAO2 * f_pk_no_wiggle(k3)
    	         + D3 * D3 * K23 * f_pk_no_wiggle(k2) * BAO3;
    
    double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
              + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
              + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
    
    return (GG + GM_MG + NW);

}

double Bispectrum_Tree_BAO_RealSpace_DarkMatter_Shift(double * kvec1, double * kvec2, double * kvec3, double sigma8, double sigma2_perp) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

    double D1 = ExpDamping_RealSpace(kvec1, sigma2_perp);
    double D2 = ExpDamping_RealSpace(kvec2, sigma2_perp);
    double D3 = ExpDamping_RealSpace(kvec3, sigma2_perp);

	double K12 = 2.0 * F2_Shift(kvec1, kvec2) * pow(sigma8, 4);
	double K13 = 2.0 * F2_Shift(kvec1, kvec3) * pow(sigma8, 4);
	double K23 = 2.0 * F2_Shift(kvec2, kvec3) * pow(sigma8, 4);

    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
    double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
    double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
    
    double GG = D1 * D2 * D3 * K12 * BAO1 * BAO2
    	      + D1 * D2 * D3 * K13 * BAO1 * BAO3
    	      + D1 * D2 * D3 * K23 * BAO2 * BAO3;
    
    double GM_MG = D1 * D1 * K12 * BAO1 * f_pk_no_wiggle(k2)
    	         + D2 * D2 * K12 * f_pk_no_wiggle(k1) * BAO2
    
    	         + D1 * D1 * K13 * BAO1 * f_pk_no_wiggle(k3)
    	         + D3 * D3 * K13 * f_pk_no_wiggle(k1) * BAO3
    
    	         + D2 * D2 * K23 * BAO2 * f_pk_no_wiggle(k3)
    	         + D3 * D3 * K23 * f_pk_no_wiggle(k2) * BAO3;
    
    double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
              + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
              + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
    
    return (GG + GM_MG + NW);

}

double Bispectrum_Tree_BAO_RealSpace_DarkMatter_Tidal(double * kvec1, double * kvec2, double * kvec3, double sigma8, double sigma2_perp) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

    double D1 = ExpDamping_RealSpace(kvec1, sigma2_perp);
    double D2 = ExpDamping_RealSpace(kvec2, sigma2_perp);
    double D3 = ExpDamping_RealSpace(kvec3, sigma2_perp);

	double K12 = 2.0 * F2_Tidal(kvec1, kvec2) * pow(sigma8, 4);
	double K13 = 2.0 * F2_Tidal(kvec1, kvec3) * pow(sigma8, 4);
	double K23 = 2.0 * F2_Tidal(kvec2, kvec3) * pow(sigma8, 4);

    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
    double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
    double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
    
    double GG = D1 * D2 * D3 * K12 * BAO1 * BAO2
    	      + D1 * D2 * D3 * K13 * BAO1 * BAO3
    	      + D1 * D2 * D3 * K23 * BAO2 * BAO3;
    
    double GM_MG = D1 * D1 * K12 * BAO1 * f_pk_no_wiggle(k2)
    	         + D2 * D2 * K12 * f_pk_no_wiggle(k1) * BAO2
    
    	         + D1 * D1 * K13 * BAO1 * f_pk_no_wiggle(k3)
    	         + D3 * D3 * K13 * f_pk_no_wiggle(k1) * BAO3
    
    	         + D2 * D2 * K23 * BAO2 * f_pk_no_wiggle(k3)
    	         + D3 * D3 * K23 * f_pk_no_wiggle(k2) * BAO3;
    
    double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
              + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
              + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
    
    return (GG + GM_MG + NW);

}

double Bispectrum_Tree_NoWiggle(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double bK2) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

     	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return B / alpha6;

}


double Bispectrum_Tree_BAO(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double bK2,
	                   double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1 = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2 = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3 = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
     	double GG = D1 * D2 * D3 * K12 * BAO1 * BAO2
     		  + D1 * D2 * D3 * K13 * BAO1 * BAO3
     		  + D1 * D2 * D3 * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1 * D1 * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2 * D2 * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1 * D1 * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3 * D3 * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2 * D2 * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3 * D3 * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_Reconstructed(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, 
	                             double f, double b1, double b2, double bK2, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z2_Bias_Reconstructed(kvec1, kvec2, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * Z2_Bias_Reconstructed(kvec1, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * Z2_Bias_Reconstructed(kvec2, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_Reconstructed(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, 
	                                      double f, double b1, double b2, double bK2, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z2_Bias_Reconstructed(kvec1, kvec2, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * Z2_Bias_Reconstructed(kvec1, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * Z2_Bias_Reconstructed(kvec2, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

    	return B / alpha6;

}

double Primordial_Matter_Bispectrum_Local(double * kvec1_in, double * kvec2_in, double * kvec3_in) {

	double k1 = NORM(kvec1_in);
	double k2 = NORM(kvec2_in);
	double k3 = NORM(kvec3_in);

	double M1 = 0.0;
	double M2 = 0.0;
	double M3 = 0.0;

	if (f_pk_pri(k1) < 1.0e-20) {
	    M1 = 0.0;
	} else {
	    M1 = sqrt( f_pk(k1) / f_pk_pri(k1) );
	}
	if (f_pk_pri(k2) < 1.0e-20) {
	    M2 = 0.0;
	} else {
	    M2 = sqrt( f_pk(k2) / f_pk_pri(k2) );
	}
	if (f_pk_pri(k3) < 1.0e-20) {
	    M3 = 0.0;
	} else {
	    M3 = sqrt( f_pk(k3) / f_pk_pri(k3) );
	}
	
	double B = 2.0 * ( f_pk_pri(k1) * f_pk_pri(k2) + f_pk_pri(k1) * f_pk_pri(k3) + f_pk_pri(k2) * f_pk_pri(k3) );
	
	return M1 * M2 * M3 * B;
}

double Primordial_Matter_Bispectrum_Local_R(double * kvec1_in, double * kvec2_in, double * kvec3_in, double R) {

	double k1 = NORM(kvec1_in);
	double k2 = NORM(kvec2_in);
	double k3 = NORM(kvec3_in);

	double M1 = 0.0;
	double M2 = 0.0;
	double M3 = 0.0;

	if (f_pk_pri(k1) < 1.0e-20) {
	    M1 = 0.0;
	} else {
	    M1 = sqrt( f_pk(k1) / f_pk_pri(k1) );
	}
	if (f_pk_pri(k2) < 1.0e-20) {
	    M2 = 0.0;
	} else {
	    M2 = sqrt( f_pk(k2) / f_pk_pri(k2) );
	}
	if (f_pk_pri(k3) < 1.0e-20) {
	    M3 = 0.0;
	} else {
	    M3 = sqrt( f_pk(k3) / f_pk_pri(k3) );
	}
	
	double B = 2.0 * ( f_pk_pri(k1) * f_pk_pri(k2) + f_pk_pri(k1) * f_pk_pri(k3) + f_pk_pri(k2) * f_pk_pri(k3) );

	return M1 * M2 * M3 * B * W1(k1, R) * W1(k2, R) * W1(k3, R);
}


double Primordial_Matter_Bispectrum_Equilateral(double * kvec1_in, double * kvec2_in, double * kvec3_in) {

	double k1 = NORM(kvec1_in);
	double k2 = NORM(kvec2_in);
	double k3 = NORM(kvec3_in);

	double M1 = 0.0;
	double M2 = 0.0;
	double M3 = 0.0;

	if (f_pk_pri(k1) < 1.0e-20) {
	    M1 = 0.0;
	} else {
	    M1 = sqrt( f_pk(k1) / f_pk_pri(k1) );
	}
	if (f_pk_pri(k2) < 1.0e-20) {
	    M2 = 0.0;
	} else {
	    M2 = sqrt( f_pk(k2) / f_pk_pri(k2) );
	}
	if (f_pk_pri(k3) < 1.0e-20) {
	    M3 = 0.0;
	} else {
	    M3 = sqrt( f_pk(k3) / f_pk_pri(k3) );
	}

	double B1 = - ( f_pk_pri(k1) * f_pk_pri(k2) + f_pk_pri(k1) * f_pk_pri(k3) + f_pk_pri(k2) * f_pk_pri(k3) );
	double B2 = - 2.0 * pow( f_pk_pri(k1) * f_pk_pri(k2) * f_pk_pri(k3), 2.0/3.0 );
	double B3 = pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k3)
		  + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k3)
		  + pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k2)
		  + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k2)
		  + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k1)
		  + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k1);
	double B = 6.0 * (B1 + B2 + B3);

	return M1 * M2 * M3 * B;
}

double Primordial_Matter_Bispectrum_Orthogonal(double * kvec1_in, double * kvec2_in, double * kvec3_in) {

	double k1 = NORM(kvec1_in);
	double k2 = NORM(kvec2_in);
	double k3 = NORM(kvec3_in);

	double M1 = 0.0;
	double M2 = 0.0;
	double M3 = 0.0;

	if (f_pk_pri(k1) < 1.0e-20) {
	    M1 = 0.0;
	} else {
	    M1 = sqrt( f_pk(k1) / f_pk_pri(k1) );
	}
	if (f_pk_pri(k2) < 1.0e-20) {
	    M2 = 0.0;
	} else {
	    M2 = sqrt( f_pk(k2) / f_pk_pri(k2) );
	}
	if (f_pk_pri(k3) < 1.0e-20) {
	    M3 = 0.0;
	} else {
	    M3 = sqrt( f_pk(k3) / f_pk_pri(k3) );
	}

	double B1 = - 3.0 * ( f_pk_pri(k1) * f_pk_pri(k2) + f_pk_pri(k1) * f_pk_pri(k3) + f_pk_pri(k2) * f_pk_pri(k3) );
	double B2 = - 8.0 * pow( f_pk_pri(k1) * f_pk_pri(k2) * f_pk_pri(k3), 2.0/3.0 );
	double B3 = 3.0 * ( pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k3)
		              + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k3)
		              + pow(f_pk_pri(k1), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k2)
		              + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k1), 2.0/3.0) * f_pk_pri(k2)
		              + pow(f_pk_pri(k2), 1.0/3.0) * pow(f_pk_pri(k3), 2.0/3.0) * f_pk_pri(k1)
		              + pow(f_pk_pri(k3), 1.0/3.0) * pow(f_pk_pri(k2), 2.0/3.0) * f_pk_pri(k1) );
	double B = 6.0 * (B1 + B2 + B3);

	return M1 * M2 * M3 * B;
}

double Bispectrum_NonGaussian_Local(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, 
	                                double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1);
	double B = Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3);
	B = pow(sigma8, 3) * B; 
	return factor * B / alpha6;

}

double Bispectrum_NonGaussian_Equilateral(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los,
                                   	      double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1);
	double B = Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3);
	B = pow(sigma8, 3) * B; 
	return factor * B / alpha6;

}

double Bispectrum_NonGaussian_Orthogonal(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los,
                                   	     double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1);
	double B = Primordial_Matter_Bispectrum_Orthogonal(kvec1, kvec2, kvec3);
	B = pow(sigma8, 3) * B; 
	return factor * B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                    double alpha_perp, double alpha_parallel, double sigma8, double f, 
    					    double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO
					    ) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2) * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2) * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * Z2_Bias(kvec3, kvec2, los, f, b1, b2, bK2) * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * Z2_Bias(kvec2, kvec1, los, f, b1, b2, bK2) * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * Z2_Bias(kvec3, kvec1, los, f, b1, b2, bK2) * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B_113_A = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias(kvec2_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias(kvec1_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec3_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec1_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k3)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec3_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec2_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k3);

	double B_113_B = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias(kvec3, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec2, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec1, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3) * f_pk(p);

	double B = B_122_A + B_122_B + B_113_A + B_113_B;
	B = pow(sigma8, 5) * B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Equilateral(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                    double alpha_perp, double alpha_parallel, double sigma8, double f, 
    					    double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO
					    ) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z2_Bias(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * Z2_Bias(kvec2, kvec3, los, f, b1, b2, bK2) * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * Z2_Bias(kvec1, kvec3, los, f, b1, b2, bK2) * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * Z2_Bias(kvec1, kvec2, los, f, b1, b2, bK2) * Z2_Bias(kvec1_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * Z2_Bias(kvec3, kvec2, los, f, b1, b2, bK2) * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * Z2_Bias(kvec2, kvec1, los, f, b1, b2, bK2) * Z2_Bias(kvec2_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * Z2_Bias(kvec3, kvec1, los, f, b1, b2, bK2) * Z2_Bias(kvec3_M_pvec, pvec, los, f, b1, b2, bK2) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B_113_A = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias(kvec2_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias(kvec1_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec3_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec1_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k3)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec3_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec2_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k3);

	double B_113_B = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias(kvec3, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec2, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias(kvec1, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3) * f_pk(p);

	double B = B_122_A + B_122_B + B_113_A + B_113_B;
	B = pow(sigma8, 5) * B;
	return B / alpha6;

}



double Bispectrum_NonGaussian_From_PB_Local_Reconstructed(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                    double alpha_perp, double alpha_parallel, double sigma8, double f, 
    					    double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
	               * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * Z2_Bias_Reconstructed(kvec2, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * Z2_Bias_Reconstructed(kvec1, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * Z2_Bias_Reconstructed(kvec1, kvec2, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * Z2_Bias_Reconstructed(kvec3, kvec2, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * Z2_Bias_Reconstructed(kvec2, kvec1, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * Z2_Bias_Reconstructed(kvec3, kvec1, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B_113_A = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias_Reconstructed(kvec2_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias_Reconstructed(kvec1_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec3_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec1_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k3)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec3_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec2_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k3);

	double B_113_B = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias_Reconstructed(kvec3, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec2, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec1, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3) * f_pk(p);

	double B = B_122_A + B_122_B + B_113_A + B_113_B;
	B = pow(sigma8, 5) * B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Equilateral_Reconstructed(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                    double alpha_perp, double alpha_parallel, double sigma8, double f, 
    					    double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
	               * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec3, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R) 
		       * Z1_Bias(kvec2, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R)
		       * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_P_pvec, M_pvec, los, f, b1, b2, bK2, b1_fid, R)
		       * Z1_Bias(kvec1, los, f, b1) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * Z2_Bias_Reconstructed(kvec2, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * Z2_Bias_Reconstructed(kvec1, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec3, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * Z2_Bias_Reconstructed(kvec1, kvec2, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec1_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * Z2_Bias_Reconstructed(kvec3, kvec2, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec2, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * Z2_Bias_Reconstructed(kvec2, kvec1, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec2_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * Z2_Bias_Reconstructed(kvec3, kvec1, los, f, b1, b2, bK2, b1_fid, R) * Z2_Bias_Reconstructed(kvec3_M_pvec, pvec, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec1, los, f, b1)
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B_113_A = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias_Reconstructed(kvec2_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias_Reconstructed(kvec1_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec3_P_pvec, kvec1, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k1)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec1_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1_P_pvec, M_kvec1, M_pvec) * f_pk(k3)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec3_P_pvec, kvec2, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec3_P_pvec, M_kvec3, M_pvec) * f_pk(k2)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec2_P_pvec, kvec3, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec2_P_pvec, M_kvec2, M_pvec) * f_pk(k3);

	double B_113_B = 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z3_Bias_Reconstructed(kvec3, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec2, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3) * f_pk(p)

		       + 3.0 * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z3_Bias_Reconstructed(kvec1, pvec, M_pvec, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
	               * Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3) * f_pk(p);

	double B = B_122_A + B_122_B + B_113_A + B_113_B;
	B = pow(sigma8, 5) * B;
	return B / alpha6;

}

std::complex<double> Bispectrum_NonGaussian_Aniso(double * kvec1_in, double * kvec2_in, double * kvec3_in,
                              	                  double * los, 
	                            double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
				    double f_NL_local, double fNL_lambda20, 
				    double fNL_lambda21_real, double fNL_lambda21_imag,
				    double fNL_lambda22_real, double fNL_lambda22_imag
				    ) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double M1 = 0.0;
	double M2 = 0.0;
	double M3 = 0.0;

	if (f_pk_pri(k1) < 1.0e-20) {
	    M1 = 0.0;
	} else {
	    M1 = sigma8 * sqrt( f_pk(k1) / f_pk_pri(k1) );
	}
	if (f_pk_pri(k2) < 1.0e-20) {
	    M2 = 0.0;
	} else {
	    M2 = sigma8 * sqrt( f_pk(k2) / f_pk_pri(k2) );
	}
	if (f_pk_pri(k3) < 1.0e-20) {
	    M3 = 0.0;
	} else {
	    M3 = sigma8 * sqrt( f_pk(k3) / f_pk_pri(k3) );
	}

	double factor = Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
	              * M1 * M2 * M3;
	
	std::complex<double> fNL_lambda2M[5];

	fNL_lambda2M[0].real(fNL_lambda22_real);
	fNL_lambda2M[0].imag(-1.0 * fNL_lambda22_imag);

	fNL_lambda2M[1].real(-1.0 * fNL_lambda21_real);
	fNL_lambda2M[1].imag(fNL_lambda21_imag);

	fNL_lambda2M[2].real(fNL_lambda20);
	fNL_lambda2M[2].imag(0.0);

	fNL_lambda2M[3].real(fNL_lambda21_real);
	fNL_lambda2M[3].imag(fNL_lambda21_imag);

	fNL_lambda2M[4].real(fNL_lambda22_real);
	fNL_lambda2M[4].imag(fNL_lambda22_imag);

	std::complex<double> L12 = 0.0;
	std::complex<double> L13 = 0.0;
	std::complex<double> L23 = 0.0;
	std::complex<double> _I_(0.0,1.0);
	int iM = 0;
	for(int M = -2; M <=2; M++) {
	    std::complex<double> Ylm1 = calcSphericalHarmonics(2, M, kvec1);
	    std::complex<double> Ylm2 = calcSphericalHarmonics(2, M, kvec2);
	    std::complex<double> Ylm3 = calcSphericalHarmonics(2, M, kvec3);

	    L12 += (fNL_lambda2M[iM]) * (Ylm1 + Ylm2);
	    L13 += (fNL_lambda2M[iM]) * (Ylm1 + Ylm3);
	    L23 += (fNL_lambda2M[iM]) * (Ylm2 + Ylm3);

	    iM++;
	}

	std::complex<double> B = 2.0 * (f_NL_local + L12 ) * f_pk_pri(k1) * f_pk_pri(k2)
	                       + 2.0 * (f_NL_local + L13 ) * f_pk_pri(k1) * f_pk_pri(k3)
	                       + 2.0 * (f_NL_local + L23 ) * f_pk_pri(k2) * f_pk_pri(k3);
	
	return factor * B / alpha6;

}

double T_PartOfT2211(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double sigma8, double f, double b1, double b2, double bK2) {

	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);

	double M_kvec1[3] = {- kvec1[0], - kvec1[1], - kvec1[2]};
	double M_kvec2[3] = {- kvec2[0], - kvec2[1], - kvec2[2]};

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k24 = NORM(kvec24);
	double k23 = NORM(kvec23);
	double k14 = NORM(kvec14);
	double k13 = NORM(kvec13);

	double T2211 = Z2_Bias(M_kvec1, kvec14, los, f, b1, b2, bK2) * Z2_Bias(M_kvec2, kvec23, los, f, b1, b2, bK2)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0
		     + Z2_Bias(M_kvec1, kvec13, los, f, b1, b2, bK2) * Z2_Bias(M_kvec2, kvec24, los, f, b1, b2, bK2)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0;

	return T2211;
}

double T_PartOfT2211_b2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double sigma8, double f, double b1) {

	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);

	double M_kvec1[3] = {- kvec1[0], - kvec1[1], - kvec1[2]};
	double M_kvec2[3] = {- kvec2[0], - kvec2[1], - kvec2[2]};

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k24 = NORM(kvec24);
	double k23 = NORM(kvec23);
	double k14 = NORM(kvec14);
	double k13 = NORM(kvec13);

	double T2211 = Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0) * ( D1D1(M_kvec2, kvec23) / 2.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + ( D1D1(M_kvec1, kvec14) / 2.0 ) * Z2_Bias(M_kvec2, kvec23, los, f, b1, 0.0, 0.0)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + Z2_Bias(M_kvec1, kvec13, los, f, b1, 0.0, 0.0) * ( D1D1(M_kvec2, kvec24) / 2.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0

		     + ( D1D1(M_kvec1, kvec13) / 2.0 ) * Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0;

	return T2211;
}

double T_PartOfT2211_bK2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double sigma8, double f, double b1) {

	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);

	double M_kvec1[3] = {- kvec1[0], - kvec1[1], - kvec1[2]};
	double M_kvec2[3] = {- kvec2[0], - kvec2[1], - kvec2[2]};

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k24 = NORM(kvec24);
	double k23 = NORM(kvec23);
	double k14 = NORM(kvec14);
	double k13 = NORM(kvec13);

	double T2211 = Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0) * K1K1(M_kvec2, kvec23)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + K1K1(M_kvec1, kvec14) * Z2_Bias(M_kvec2, kvec23, los, f, b1, 0.0, 0.0)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + Z2_Bias(M_kvec1, kvec13, los, f, b1, 0.0, 0.0) * K1K1(M_kvec2, kvec24)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0

		     + K1K1(M_kvec1, kvec13) * Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0;

	return T2211;
}

double T_PartOfT2211_b2_b2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double sigma8, double f, double b1) {

	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);

	double M_kvec1[3] = {- kvec1[0], - kvec1[1], - kvec1[2]};
	double M_kvec2[3] = {- kvec2[0], - kvec2[1], - kvec2[2]};

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k24 = NORM(kvec24);
	double k23 = NORM(kvec23);
	double k14 = NORM(kvec14);
	double k13 = NORM(kvec13);

	double T2211 = ( D1D1(M_kvec1, kvec14) / 2.0 ) * ( D1D1(M_kvec2, kvec23) / 2.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) 
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + ( D1D1(M_kvec1, kvec13) / 2.0 ) * ( D1D1(M_kvec2, kvec24) / 2.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) 
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0;

	return T2211;
}

double T_PartOfT2211_bK2_bK2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double sigma8, double f, double b1) {

	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);

	double M_kvec1[3] = {- kvec1[0], - kvec1[1], - kvec1[2]};
	double M_kvec2[3] = {- kvec2[0], - kvec2[1], - kvec2[2]};

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k24 = NORM(kvec24);
	double k23 = NORM(kvec23);
	double k14 = NORM(kvec14);
	double k13 = NORM(kvec13);

	double T2211 = K1K1(M_kvec1, kvec14) * K1K1(M_kvec2, kvec23)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) 
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + K1K1(M_kvec1, kvec13) * K1K1(M_kvec2, kvec24)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) 
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0;

	return T2211;
}

double T_PartOfT2211_b2_bK2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double sigma8, double f, double b1) {

	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);

	double M_kvec1[3] = {- kvec1[0], - kvec1[1], - kvec1[2]};
	double M_kvec2[3] = {- kvec2[0], - kvec2[1], - kvec2[2]};

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k24 = NORM(kvec24);
	double k23 = NORM(kvec23);
	double k14 = NORM(kvec14);
	double k13 = NORM(kvec13);

	double T2211 = K1K1(M_kvec1, kvec14) * ( D1D1(M_kvec2, kvec23) / 2.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1)
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + ( D1D1(M_kvec1, kvec14) / 2.0 ) * K1K1(M_kvec2, kvec23)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) 
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0

		     + K1K1(M_kvec1, kvec13) * ( D1D1(M_kvec2, kvec24) / 2.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) 
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0

		     + ( D1D1(M_kvec1, kvec13) / 2.0 ) * K1K1(M_kvec2, kvec24)
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1)
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0;

	return T2211;
}


double Trispectrum_Tree(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		        double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T2211 = T_PartOfT2211(kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1, b2, bK2)
		     + T_PartOfT2211(kvec1, kvec3, kvec2, kvec4, los, sigma8, f, b1, b2, bK2)
		     + T_PartOfT2211(kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1, b2, bK2)
		     + T_PartOfT2211(kvec2, kvec3, kvec1, kvec4, los, sigma8, f, b1, b2, bK2)
		     + T_PartOfT2211(kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1, b2, bK2)
		     + T_PartOfT2211(kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1, b2, bK2);

	double T3111 = Z3_Bias(kvec1, kvec2, kvec3, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + Z3_Bias(kvec1, kvec2, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + Z3_Bias(kvec1, kvec3, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + Z3_Bias(kvec2, kvec3, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO) 
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 4.0 * T2211 + 6.0 * T3111;  
	return T / alpha9;
}


double Trispectrum_Tree_b2(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T2211 = T_PartOfT2211_b2(kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2(kvec1, kvec3, kvec2, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2(kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_b2(kvec2, kvec3, kvec1, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2(kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_b2(kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1);

	double T3111 = Z3_Bias_b2(kvec1, kvec2, kvec3, los, f) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + Z3_Bias_b2(kvec1, kvec2, kvec4, los, f) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + Z3_Bias_b2(kvec1, kvec3, kvec4, los, f) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + Z3_Bias_b2(kvec2, kvec3, kvec4, los, f) 
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 4.0 * T2211 + 6.0 * T3111;  
	return T / alpha9;
}

double Trispectrum_Tree_bK2(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T2211 = T_PartOfT2211_bK2(kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2(kvec1, kvec3, kvec2, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2(kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2(kvec2, kvec3, kvec1, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2(kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2(kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1);

	double T3111 = Z3_Bias_bK2(kvec1, kvec2, kvec3, los, f) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + Z3_Bias_bK2(kvec1, kvec2, kvec4, los, f) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + Z3_Bias_bK2(kvec1, kvec3, kvec4, los, f) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + Z3_Bias_bK2(kvec2, kvec3, kvec4, los, f) 
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 4.0 * T2211 + 6.0 * T3111;  
	return T / alpha9;
}

double Trispectrum_Tree_b2_b2(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double T2211 = T_PartOfT2211_b2_b2(kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_b2(kvec1, kvec3, kvec2, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_b2(kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_b2(kvec2, kvec3, kvec1, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_b2(kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_b2(kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1);

	double T = 4.0 * T2211;

	return T / alpha9;
}

double Trispectrum_Tree_b2_bK2(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double T2211 = T_PartOfT2211_b2_bK2(kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_bK2(kvec1, kvec3, kvec2, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_bK2(kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_bK2(kvec2, kvec3, kvec1, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_bK2(kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_b2_bK2(kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1);

	double T = 4.0 * T2211;

	return T / alpha9;
}

double Trispectrum_Tree_bK2_bK2(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double T2211 = T_PartOfT2211_bK2_bK2(kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2_bK2(kvec1, kvec3, kvec2, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2_bK2(kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2_bK2(kvec2, kvec3, kvec1, kvec4, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2_bK2(kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1)
		     + T_PartOfT2211_bK2_bK2(kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1);

	double T = 4.0 * T2211;

	return T / alpha9;
}


double Trispectrum_Tree_b3(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T3111 = ( D1D1D1(kvec1, kvec2, kvec3) / 6.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + ( D1D1D1(kvec1, kvec2, kvec4) / 6.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + ( D1D1D1(kvec1, kvec3, kvec4) / 6.0 )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + ( D1D1D1(kvec2, kvec3, kvec4) / 6.0 )
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 6.0 * T3111;  
	return T / alpha9;
}

double Trispectrum_Tree_bK3(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T3111 = ( K1K1K1(kvec1, kvec2, kvec3) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + ( K1K1K1(kvec1, kvec2, kvec4) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + ( K1K1K1(kvec1, kvec3, kvec4) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + ( K1K1K1(kvec2, kvec3, kvec4) )
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 6.0 * T3111;  
	return T / alpha9;
}

double Trispectrum_Tree_bDK(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T3111 = ( D1K1K1(kvec1, kvec2, kvec3) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + ( D1K1K1(kvec1, kvec2, kvec4) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + ( D1K1K1(kvec1, kvec3, kvec4) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + ( D1K1K1(kvec2, kvec3, kvec4) )
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 6.0 * T3111;  
	return T / alpha9;
}


double Trispectrum_Tree_bO(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		           double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T3111 = ( O3(kvec1, kvec2, kvec3) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + ( O3(kvec1, kvec2, kvec4) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + ( O3(kvec1, kvec3, kvec4) )
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + ( O3(kvec2, kvec3, kvec4) )
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 6.0 * T3111;  
	return T / alpha9;
}



double T_PartOfT2211_Reconstructed(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * los, double sigma8, double f, double b1, double b2, double bK2, double b1_fid, double R) {

	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);

	double M_kvec1[3] = {- kvec1[0], - kvec1[1], - kvec1[2]};
	double M_kvec2[3] = {- kvec2[0], - kvec2[1], - kvec2[2]};

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k24 = NORM(kvec24);
	double k23 = NORM(kvec23);
	double k14 = NORM(kvec14);
	double k13 = NORM(kvec13);

	double T2211 = Z2_Bias_Reconstructed(M_kvec1, kvec14, los, f, b1, b2, bK2, b1_fid, R) 
	             * Z2_Bias_Reconstructed(M_kvec2, kvec23, los, f, b1, b2, bK2, b1_fid, R)
		     * Z1_Bias(kvec1, los, f, b1)
		     * Z1_Bias(kvec2, los, f, b1)
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k14) + f_pk(k23) ) / 2.0
		     + Z2_Bias_Reconstructed(M_kvec1, kvec13, los, f, b1, b2, bK2, b1_fid, R) 
		     * Z2_Bias_Reconstructed(M_kvec2, kvec24, los, f, b1, b2, bK2, b1_fid, R)
		     * Z1_Bias(kvec1, los, f, b1)
		     * Z1_Bias(kvec2, los, f, b1)
		     * pow(sigma8,6) * f_pk(k1) * f_pk(k2) * ( f_pk(k13) + f_pk(k24) ) / 2.0;

	return T2211;
}


double Trispectrum_Tree_Reconstructed(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in,
		                      double * los, double alpha_perp, double alpha_parallel, double sigma8, double f,
		                      double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha9 = alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double T2211 = T_PartOfT2211_Reconstructed(kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1, b2, bK2, b1_fid, R)
		     + T_PartOfT2211_Reconstructed(kvec1, kvec3, kvec2, kvec4, los, sigma8, f, b1, b2, bK2, b1_fid, R)
		     + T_PartOfT2211_Reconstructed(kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1, b2, bK2, b1_fid, R)
		     + T_PartOfT2211_Reconstructed(kvec2, kvec3, kvec1, kvec4, los, sigma8, f, b1, b2, bK2, b1_fid, R)
		     + T_PartOfT2211_Reconstructed(kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1, b2, bK2, b1_fid, R)
		     + T_PartOfT2211_Reconstructed(kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1, b2, bK2, b1_fid, R);

	double T3111 = Z3_Bias_Reconstructed(kvec1, kvec2, kvec3, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		     + Z3_Bias_Reconstructed(kvec1, kvec2, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k2) * f_pk(k4)
		     + Z3_Bias_Reconstructed(kvec1, kvec3, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
		     * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k1) * f_pk(k3) * f_pk(k4)
		     + Z3_Bias_Reconstructed(kvec2, kvec3, kvec4, los, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R) 
		     * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * pow(sigma8, 6) * f_pk(k2) * f_pk(k3) * f_pk(k4);

	double T = 4.0 * T2211 + 6.0 * T3111;  
	return T / alpha9;
}

double Bispectrum_Tree_ReconForCovariance(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, 
	                                  double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double bK2,
					  double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z2_Bias_Reconstructed(kvec1, kvec2, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias_D(kvec1, los, f, b1, b1_fid, R) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 4);
	double K13 = 2.0 * Z2_Bias_Reconstructed(kvec1, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias_D(kvec1, los, f, b1, b1_fid, R) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);
	double K23 = 2.0 * Z2_Bias_D(kvec2, kvec3, los, f, b1, b2, bK2, b1_fid, R) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 4);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Powerspectrum_Tree_ReconForCovariance(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b1_fid, double R) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);

	double k = NORM(kvec);
	double P = Z1_Bias_D(kvec, los, f, b1, b1_fid, R) * Z1_Bias_D(kvec, los, f, b1, b1_fid, R) * pow(sigma8, 2) * f_pk(k);
	return P / alpha3;
}

double Powerspectrum_Tree_ReconForPBCrossCovariance(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b1_fid, double R) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);

	double k = NORM(kvec);
	double P = Z1_Bias_D(kvec, los, f, b1, b1_fid, R) * Z1_Bias(kvec, los, f, b1) * pow(sigma8, 2) * f_pk(k);
	return P / alpha3;
}


double Q_PartOfT22211(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double sigma8, double f, double b1) {

	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);
	double kvec15[3] = PLUS(kvec1, kvec5);

	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec25[3] = PLUS(kvec2, kvec5);

	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double k13 = NORM(kvec13);
	double k14 = NORM(kvec14);
	double k15 = NORM(kvec15);

	double k23 = NORM(kvec23);
	double k24 = NORM(kvec24);
	double k25 = NORM(kvec25);

	double T22211 = Z2_Bias(M_kvec1, kvec13, los, f, b1, 0.0, 0.0) * Z2_Bias(kvec13,  kvec25, los, f, b1, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec25, los, f, b1, 0.0, 0.0)
	              * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 8) * f_pk(k13) * f_pk(k25) * f_pk(k1) * f_pk(k2)
		      + Z2_Bias(M_kvec1, kvec13, los, f, b1, 0.0, 0.0) * Z2_Bias(kvec13,  kvec24, los, f, b1, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0)
	              * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 8) * f_pk(k13) * f_pk(k24) * f_pk(k1) * f_pk(k2)
		      + Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0) * Z2_Bias(kvec14,  kvec25, los, f, b1, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec25, los, f, b1, 0.0, 0.0)
	              * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 8) * f_pk(k14) * f_pk(k25) * f_pk(k1) * f_pk(k2)
		      + Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0) * Z2_Bias(kvec14,  kvec23, los, f, b1, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec23, los, f, b1, 0.0, 0.0)
	              * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 8) * f_pk(k14) * f_pk(k23) * f_pk(k1) * f_pk(k2)
		      + Z2_Bias(M_kvec1, kvec15, los, f, b1, 0.0, 0.0) * Z2_Bias(kvec15,  kvec23, los, f, b1, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec23, los, f, b1, 0.0, 0.0)
	              * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 8) * f_pk(k15) * f_pk(k23) * f_pk(k1) * f_pk(k2)
		      + Z2_Bias(M_kvec1, kvec15, los, f, b1, 0.0, 0.0) * Z2_Bias(kvec15,  kvec24, los, f, b1, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0)
	              * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * pow(sigma8, 8) * f_pk(k15) * f_pk(k24) * f_pk(k1) * f_pk(k2);
	return T22211;
}

double Q_PartOfT32111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double sigma8, double f, double b1) {

	double kvec14[3] = PLUS(kvec1, kvec4);

	double kvec24[3] = PLUS(kvec2, kvec4);

	double kvec34[3] = PLUS(kvec3, kvec4);

	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double k14 = NORM(kvec14);

	double k24 = NORM(kvec24);

	double k34 = NORM(kvec34);

	double T32111 = Z2_Bias(M_kvec3, kvec34, los, f, b1, 0.0, 0.0) * Z3_Bias(kvec1, kvec2, kvec34, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		      * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 8) * f_pk(k34) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		      + Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0) * Z3_Bias(kvec1, kvec3, kvec24, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		      * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 8) * f_pk(k24) * f_pk(k1) * f_pk(k2) * f_pk(k3)
		      + Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0) * Z3_Bias(kvec2, kvec3, kvec14, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		      * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * pow(sigma8, 8) * f_pk(k14) * f_pk(k1) * f_pk(k2) * f_pk(k3);
	
	return T32111;
}

double Q_PartOfT41111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * los, double sigma8, double f, double b1) {


	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);
	double k5 = NORM(kvec5);

	double T41111 = Z4_Bias(kvec1, kvec2, kvec3, kvec4, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		      * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) 
		      * pow(sigma8, 8) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4)
		      + Z4_Bias(kvec1, kvec2, kvec3, kvec5, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		      * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec5, los, f, b1) 
		      * pow(sigma8, 8) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k5)
		      + Z4_Bias(kvec1, kvec2, kvec4, kvec5, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		      * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec5, los, f, b1) 
		      * pow(sigma8, 8) * f_pk(k1) * f_pk(k2) * f_pk(k4) * f_pk(k5)
		      + Z4_Bias(kvec1, kvec3, kvec4, kvec5, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		      * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec5, los, f, b1) 
		      * pow(sigma8, 8) * f_pk(k1) * f_pk(k3) * f_pk(k4) * f_pk(k5)
		      + Z4_Bias(kvec2, kvec3, kvec4, kvec5, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		      * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec5, los, f, b1) 
		      * pow(sigma8, 8) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k5);

	return T41111;

}

double P5spectrum_Tree(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in, double * kvec5_in,
		       double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha12 = alpha3 * alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };
	double kvec5[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);
	calcTrueWavevector(kvec5_in, los, alpha_perp, alpha_parallel, kvec5);

	double T22211 = Q_PartOfT22211(kvec4, kvec5, kvec1, kvec2, kvec3, los, sigma8, f, b1) 
		     
		      + Q_PartOfT22211(kvec4, kvec1, kvec5, kvec2, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT22211(kvec4, kvec2, kvec5, kvec1, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT22211(kvec4, kvec3, kvec5, kvec1, kvec2, los, sigma8, f, b1) 
		      + Q_PartOfT22211(kvec5, kvec1, kvec4, kvec2, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT22211(kvec5, kvec2, kvec4, kvec1, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT22211(kvec5, kvec3, kvec4, kvec1, kvec2, los, sigma8, f, b1) 
	
		      + Q_PartOfT22211(kvec1, kvec2, kvec4, kvec5, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT22211(kvec1, kvec3, kvec4, kvec5, kvec2, los, sigma8, f, b1) 
		      + Q_PartOfT22211(kvec2, kvec3, kvec4, kvec5, kvec1, los, sigma8, f, b1);

	/************/
	double T32111 = Q_PartOfT32111(kvec1, kvec2, kvec3, kvec4, kvec5, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec2, kvec3, kvec5, kvec4, los, sigma8, f, b1) 

		      + Q_PartOfT32111(kvec1, kvec2, kvec4, kvec3, kvec5, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec2, kvec4, kvec5, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec3, kvec4, kvec2, kvec5, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec3, kvec4, kvec5, kvec2, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec2, kvec3, kvec4, kvec1, kvec5, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec2, kvec3, kvec4, kvec5, kvec1, los, sigma8, f, b1) 

		      + Q_PartOfT32111(kvec1, kvec2, kvec5, kvec3, kvec4, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec2, kvec5, kvec4, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec3, kvec5, kvec2, kvec4, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec3, kvec5, kvec4, kvec2, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec2, kvec3, kvec5, kvec1, kvec4, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec2, kvec3, kvec5, kvec4, kvec1, los, sigma8, f, b1) 

		      + Q_PartOfT32111(kvec1, kvec5, kvec4, kvec2, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec1, kvec5, kvec4, kvec3, kvec2, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec2, kvec5, kvec4, kvec1, kvec3, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec2, kvec5, kvec4, kvec3, kvec1, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec3, kvec5, kvec4, kvec1, kvec2, los, sigma8, f, b1) 
		      + Q_PartOfT32111(kvec3, kvec5, kvec4, kvec2, kvec1, los, sigma8, f, b1);
		/*******/
	
	double T41111 = Q_PartOfT41111(kvec1, kvec5, kvec4, kvec2, kvec3, los, sigma8, f, b1);

	double P5 = 8.0 * T22211 + 12.0 * T32111 + 24.0 * T41111;
	return P5 / alpha12;
}


double T6_PartOfT222211_A(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);

	double M_kvec1[3] = MINUS_XVEC(kvec1);
	double M_kvec2[3] = MINUS_XVEC(kvec2);

	double kvec13[3] = PLUS(kvec1, kvec3);
	double kvec14[3] = PLUS(kvec1, kvec4);
	double kvec23[3] = PLUS(kvec2, kvec3);
	double kvec24[3] = PLUS(kvec2, kvec4);

	double k13 = NORM(kvec13);
	double k14 = NORM(kvec14);
	double k23 = NORM(kvec23);
	double k24 = NORM(kvec24);

	double kvec135[3] = PLUS3(kvec1, kvec3, kvec5);
	double kvec246[3] = PLUS3(kvec2, kvec4, kvec6);
	double kvec136[3] = PLUS3(kvec1, kvec3, kvec6);
	double kvec245[3] = PLUS3(kvec2, kvec4, kvec5);
	double kvec146[3] = PLUS3(kvec1, kvec4, kvec6);
	double kvec235[3] = PLUS3(kvec2, kvec3, kvec5);
	double kvec145[3] = PLUS3(kvec1, kvec4, kvec5);
	double kvec236[3] = PLUS3(kvec2, kvec3, kvec6);

	double k135 = NORM(kvec135);
	double k246 = NORM(kvec246);
	double k136 = NORM(kvec136);
	double k245 = NORM(kvec245);
	double k146 = NORM(kvec146);
	double k235 = NORM(kvec235);
	double k145 = NORM(kvec145);
	double k236 = NORM(kvec236);

	double T222211 = Z2_Bias( kvec13, kvec246, los, f, b1, 0.0, 0.0)
		       * Z2_Bias( kvec24, kvec135, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec1,  kvec13, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec2,  kvec24, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) 
		       * Z1_Bias(kvec2, los, f, b1) 
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k13) * f_pk(k24) * ( f_pk(k246) + f_pk(k135) ) / 2.0

		       + Z2_Bias( kvec13, kvec245, los, f, b1, 0.0, 0.0)
		       * Z2_Bias( kvec24, kvec136, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec1,  kvec13, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec2,  kvec24, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) 
		       * Z1_Bias(kvec2, los, f, b1) 
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k13) * f_pk(k24) * ( f_pk(k245) + f_pk(k136) ) / 2.0

		       + Z2_Bias( kvec23, kvec146, los, f, b1, 0.0, 0.0)
		       * Z2_Bias( kvec14, kvec235, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec1,  kvec14, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec2,  kvec23, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) 
		       * Z1_Bias(kvec2, los, f, b1) 
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k14) * f_pk(k23) * ( f_pk(k146) + f_pk(k235) ) / 2.0

		       + Z2_Bias( kvec23, kvec145, los, f, b1, 0.0, 0.0)
		       * Z2_Bias( kvec14, kvec236, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec1,  kvec14, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec2,  kvec23, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) 
		       * Z1_Bias(kvec2, los, f, b1) 
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k14) * f_pk(k23) * ( f_pk(k145) + f_pk(k236) ) / 2.0;

	return T222211;

}


double T6_PartOfT222211_B(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {
	
	double T222211 = T6_PartOfT222211_A(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_A(kvec1, kvec2, kvec3, kvec5, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_A(kvec1, kvec2, kvec3, kvec6, kvec4, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT222211_A(kvec1, kvec2, kvec4, kvec5, kvec3, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_A(kvec1, kvec2, kvec4, kvec6, kvec3, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT222211_A(kvec1, kvec2, kvec5, kvec6, kvec3, kvec4, los, sigma8, f, b1);

	return T222211;
}

double T6_222211(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double T222211 = T6_PartOfT222211_B(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec1, kvec3, kvec2, kvec4, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec1, kvec4, kvec2, kvec3, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec1, kvec5, kvec2, kvec3, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec1, kvec6, kvec2, kvec3, kvec4, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec2, kvec3, kvec1, kvec4, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec2, kvec4, kvec1, kvec3, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec2, kvec5, kvec1, kvec3, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec2, kvec6, kvec1, kvec3, kvec4, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec3, kvec4, kvec1, kvec2, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec3, kvec5, kvec1, kvec2, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec3, kvec6, kvec1, kvec2, kvec4, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec4, kvec5, kvec1, kvec2, kvec3, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec4, kvec6, kvec1, kvec2, kvec3, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT222211_B(kvec5, kvec6, kvec1, kvec2, kvec3, kvec4, los, sigma8, f, b1);

	return T222211;
}

double T6_PartOfT322111_A_1(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double M_kvec1[3] = MINUS_XVEC(kvec1);
	double M_kvec2[3] = MINUS_XVEC(kvec2);
	double M_kvec3[3] = MINUS_XVEC(kvec3);

	double kvec34[3] = PLUS(kvec3, kvec4);
	double kvec35[3] = PLUS(kvec3, kvec5);
	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec25[3] = PLUS(kvec2, kvec5);
	double kvec14[3] = PLUS(kvec1, kvec4);
	double kvec15[3] = PLUS(kvec1, kvec5);

	double k34 = NORM(kvec34);
	double k35 = NORM(kvec35);
	double k24 = NORM(kvec24);
	double k25 = NORM(kvec25);
	double k14 = NORM(kvec14);
	double k15 = NORM(kvec15);

	double kvec126[3] = PLUS3(kvec1, kvec2, kvec6);
	double kvec136[3] = PLUS3(kvec1, kvec3, kvec6);
	double kvec236[3] = PLUS3(kvec2, kvec3, kvec6);

	double k126 = NORM(kvec126);
	double k136 = NORM(kvec136);
	double k236 = NORM(kvec236);

	double T322111 = Z3_Bias(M_kvec1, M_kvec2, kvec126, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec3, kvec34, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(kvec34, kvec126, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k34) * f_pk(k126)
		
		       + Z3_Bias(M_kvec1, M_kvec2, kvec126, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec3, kvec35, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(kvec35, kvec126, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k35) * f_pk(k126)

		       + Z3_Bias(M_kvec1, M_kvec3, kvec136, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(kvec24, kvec136, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k24) * f_pk(k136)

		       + Z3_Bias(M_kvec1, M_kvec3, kvec136, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec2, kvec25, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(kvec25, kvec136, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k25) * f_pk(k136)

		       + Z3_Bias(M_kvec2, M_kvec3, kvec236, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(kvec14, kvec236, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k14) * f_pk(k236)

		       + Z3_Bias(M_kvec2, M_kvec3, kvec236, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec1, kvec15, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(kvec15, kvec236, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k15) * f_pk(k236);

	return T322111;

}

double T6_PartOfT322111_A_2(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double M_kvec1[3] = MINUS_XVEC(kvec1);
	double M_kvec2[3] = MINUS_XVEC(kvec2);
	double M_kvec3[3] = MINUS_XVEC(kvec3);

	double kvec14[3] = PLUS(kvec1, kvec4);
	double kvec15[3] = PLUS(kvec1, kvec5);
	double kvec24[3] = PLUS(kvec2, kvec4);
	double kvec25[3] = PLUS(kvec2, kvec5);
	double kvec34[3] = PLUS(kvec3, kvec4);
	double kvec35[3] = PLUS(kvec3, kvec5);

	double k14 = NORM(kvec14);
	double k15 = NORM(kvec15);
	double k24 = NORM(kvec24);
	double k25 = NORM(kvec25);
	double k34 = NORM(kvec34);
	double k35 = NORM(kvec35);

	double T322111 = Z3_Bias(kvec1, kvec24, kvec35, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec3, kvec35, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k24) * f_pk(k35)		

		       + Z3_Bias(kvec1, kvec25, kvec34, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec2, kvec25, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec3, kvec34, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3)* f_pk(k25) * f_pk(k34) 

		       + Z3_Bias(kvec2, kvec14, kvec35, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec3, kvec35, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3)* f_pk(k14) * f_pk(k35) 

		       + Z3_Bias(kvec2, kvec15, kvec34, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec1, kvec15, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec3, kvec34, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3)* f_pk(k15) * f_pk(k34) 

		       + Z3_Bias(kvec3, kvec14, kvec25, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec1, kvec14, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec2, kvec25, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3)* f_pk(k14) * f_pk(k25) 

		       + Z3_Bias(kvec3, kvec15, kvec24, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z2_Bias(M_kvec1, kvec15, los, f, b1, 0.0, 0.0)
		       * Z2_Bias(M_kvec2, kvec24, los, f, b1, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3)* f_pk(k15) * f_pk(k24);

	return T322111;

}


double T6_PartOfT322111_A(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	return T6_PartOfT322111_A_1(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1) 
	     + T6_PartOfT322111_A_2(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

}

double T6_PartOfT322111_B(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double T322111 = T6_PartOfT322111_A(kvec3, kvec4, kvec5, kvec1, kvec2, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec2, kvec4, kvec5, kvec1, kvec3, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec2, kvec3, kvec5, kvec1, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec2, kvec3, kvec4, kvec1, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec1, kvec4, kvec5, kvec2, kvec3, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec1, kvec3, kvec5, kvec2, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec1, kvec3, kvec4, kvec2, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec1, kvec2, kvec5, kvec3, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec1, kvec2, kvec4, kvec3, kvec5, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT322111_A(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	return T322111;

}


double T6_322111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double T322111 = T6_PartOfT322111_B(kvec2, kvec3, kvec4, kvec5, kvec6, kvec1, los, sigma8, f, b1)
		       + T6_PartOfT322111_B(kvec1, kvec3, kvec4, kvec5, kvec6, kvec2, los, sigma8, f, b1)
		       + T6_PartOfT322111_B(kvec1, kvec2, kvec4, kvec5, kvec6, kvec3, los, sigma8, f, b1)
		       + T6_PartOfT322111_B(kvec1, kvec2, kvec3, kvec5, kvec6, kvec4, los, sigma8, f, b1)
		       + T6_PartOfT322111_B(kvec1, kvec2, kvec3, kvec4, kvec6, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT322111_B(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	return T322111;

}

double T6_PartOfT331111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double M_kvec1[3] = MINUS_XVEC(kvec1);
	double M_kvec2[3] = MINUS_XVEC(kvec2);
	double M_kvec3[3] = MINUS_XVEC(kvec3);
	double M_kvec4[3] = MINUS_XVEC(kvec4);

	double kvec125[3] = PLUS3(kvec1, kvec2, kvec5);
	double kvec346[3] = PLUS3(kvec3, kvec4, kvec6);
	double kvec135[3] = PLUS3(kvec1, kvec3, kvec5);
	double kvec246[3] = PLUS3(kvec2, kvec4, kvec6);
	double kvec145[3] = PLUS3(kvec1, kvec4, kvec5);
	double kvec236[3] = PLUS3(kvec2, kvec3, kvec6);
	double kvec235[3] = PLUS3(kvec2, kvec3, kvec5);
	double kvec146[3] = PLUS3(kvec1, kvec4, kvec6);
	double kvec245[3] = PLUS3(kvec2, kvec4, kvec5);
	double kvec136[3] = PLUS3(kvec1, kvec3, kvec6);
	double kvec345[3] = PLUS3(kvec3, kvec4, kvec5);
	double kvec126[3] = PLUS3(kvec1, kvec2, kvec6);

	double k125 = NORM(kvec125);
	double k346 = NORM(kvec346);
	double k135 = NORM(kvec135);
	double k246 = NORM(kvec246);
	double k145 = NORM(kvec145);
	double k236 = NORM(kvec236);
	double k235 = NORM(kvec235);
	double k146 = NORM(kvec146);
	double k245 = NORM(kvec245);
	double k136 = NORM(kvec136);
	double k345 = NORM(kvec345);
	double k126 = NORM(kvec126);

	double T331111 = Z3_Bias(M_kvec1, M_kvec2, kvec125, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z3_Bias(M_kvec3, M_kvec4, kvec346, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * ( f_pk(k125) + f_pk(k346) ) / 2.0
		       
		       + Z3_Bias(M_kvec1, M_kvec3, kvec135, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z3_Bias(M_kvec2, M_kvec4, kvec246, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * ( f_pk(k135) + f_pk(k246) ) / 2.0

		       + Z3_Bias(M_kvec1, M_kvec4, kvec145, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z3_Bias(M_kvec2, M_kvec3, kvec236, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * ( f_pk(k145) + f_pk(k236) ) / 2.0

		       + Z3_Bias(M_kvec2, M_kvec3, kvec235, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z3_Bias(M_kvec1, M_kvec4, kvec146, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * ( f_pk(k235) + f_pk(k146) ) / 2.0

		       + Z3_Bias(M_kvec2, M_kvec4, kvec245, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z3_Bias(M_kvec1, M_kvec3, kvec136, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * ( f_pk(k245) + f_pk(k136) ) / 2.0

		       + Z3_Bias(M_kvec3, M_kvec4, kvec345, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z3_Bias(M_kvec1, M_kvec2, kvec126, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * ( f_pk(k345) + f_pk(k126) ) / 2.0;

	return T331111;

}

double T6_331111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double T331111 = T6_PartOfT331111(kvec3, kvec4, kvec5, kvec6, kvec1, kvec2, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec2, kvec4, kvec5, kvec6, kvec1, kvec3, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec2, kvec3, kvec5, kvec6, kvec1, kvec4, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec2, kvec3, kvec4, kvec6, kvec1, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec2, kvec3, kvec4, kvec5, kvec1, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec4, kvec5, kvec6, kvec2, kvec3, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec3, kvec5, kvec6, kvec2, kvec4, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec3, kvec4, kvec6, kvec2, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec3, kvec4, kvec5, kvec2, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec2, kvec5, kvec6, kvec3, kvec4, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec2, kvec4, kvec6, kvec3, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec2, kvec4, kvec5, kvec3, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec2, kvec3, kvec6, kvec4, kvec5, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec2, kvec3, kvec5, kvec4, kvec6, los, sigma8, f, b1)
		       + T6_PartOfT331111(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	return T331111;

}

double T6_PartOfT421111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);

	double M_kvec1[3] = MINUS_XVEC(kvec1);
	double M_kvec2[3] = MINUS_XVEC(kvec2);
	double M_kvec3[3] = MINUS_XVEC(kvec3);
	double M_kvec4[3] = MINUS_XVEC(kvec4);

	double kvec15[3] = PLUS(kvec1, kvec5);
	double kvec25[3] = PLUS(kvec2, kvec5);
	double kvec35[3] = PLUS(kvec3, kvec5);
	double kvec45[3] = PLUS(kvec4, kvec5);

	double kvec16[3] = PLUS(kvec1, kvec6);
	double kvec26[3] = PLUS(kvec2, kvec6);
	double kvec36[3] = PLUS(kvec3, kvec6);
	double kvec46[3] = PLUS(kvec4, kvec6);

	double k15 = NORM(kvec15);
	double k25 = NORM(kvec25);
	double k35 = NORM(kvec35);
	double k45 = NORM(kvec45);

	double k16 = NORM(kvec16);
	double k26 = NORM(kvec26);
	double k36 = NORM(kvec36);
	double k46 = NORM(kvec46);

	double T421111_A = Z4_Bias(kvec1, kvec2, kvec3, kvec45, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec4, kvec45, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k45)
		         + Z4_Bias(kvec1, kvec2, kvec4, kvec35, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec3, kvec35, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k35)
		         + Z4_Bias(kvec1, kvec3, kvec4, kvec25, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec25, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k25)
		         + Z4_Bias(kvec2, kvec3, kvec4, kvec15, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec1, kvec15, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k15); 


	double T421111_B = Z4_Bias(kvec1, kvec2, kvec3, kvec46, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec4, kvec46, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k46)
		         + Z4_Bias(kvec1, kvec2, kvec4, kvec36, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec3, kvec36, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k36)
		         + Z4_Bias(kvec1, kvec3, kvec4, kvec26, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec2, kvec26, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k26)
		         + Z4_Bias(kvec2, kvec3, kvec4, kvec16, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) * Z2_Bias(M_kvec1, kvec16, los, f, b1, 0.0, 0.0)
		         * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1)
		         * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k16); 

	return T421111_A + T421111_B;

}


double T6_421111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double T421111 = T6_PartOfT421111(kvec3, kvec4, kvec5, kvec6, kvec1, kvec2, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec2, kvec4, kvec5, kvec6, kvec1, kvec3, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec2, kvec3, kvec5, kvec6, kvec1, kvec4, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec2, kvec3, kvec4, kvec6, kvec1, kvec5, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec2, kvec3, kvec4, kvec5, kvec1, kvec6, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec4, kvec5, kvec6, kvec2, kvec3, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec3, kvec5, kvec6, kvec2, kvec4, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec3, kvec4, kvec6, kvec2, kvec5, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec3, kvec4, kvec5, kvec2, kvec6, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec2, kvec5, kvec6, kvec3, kvec4, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec2, kvec4, kvec6, kvec3, kvec5, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec2, kvec4, kvec5, kvec3, kvec6, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec2, kvec3, kvec6, kvec4, kvec5, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec2, kvec3, kvec5, kvec4, kvec6, los, sigma8, f, b1)
	               + T6_PartOfT421111(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	return T421111;

}

double T6_511111(double * kvec1, double * kvec2, double * kvec3, double * kvec4, double * kvec5, double * kvec6, double * los, double sigma8, double f, double b1) {

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double k4 = NORM(kvec4);
	double k5 = NORM(kvec5);
	double k6 = NORM(kvec6);

	double T511111 = Z5_Bias(kvec1, kvec2, kvec3, kvec4, kvec5, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1)
		       * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec5, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k5)
		       + Z5_Bias(kvec1, kvec2, kvec3, kvec4, kvec6, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1)
		       * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec6, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k6)
		       + Z5_Bias(kvec1, kvec2, kvec3, kvec5, kvec6, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1)
		       * Z1_Bias(kvec3, los, f, b1) * Z1_Bias(kvec5, los, f, b1) * Z1_Bias(kvec6, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k3) * f_pk(k5) * f_pk(k6)
		       + Z5_Bias(kvec1, kvec2, kvec4, kvec5, kvec6, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec2, los, f, b1)
		       * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec5, los, f, b1) * Z1_Bias(kvec6, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k2) * f_pk(k4) * f_pk(k5) * f_pk(k6)
		       + Z5_Bias(kvec1, kvec3, kvec4, kvec5, kvec6, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		       * Z1_Bias(kvec1, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec5, los, f, b1) * Z1_Bias(kvec6, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k1) * f_pk(k3) * f_pk(k4) * f_pk(k5) * f_pk(k6)
		       + Z5_Bias(kvec2, kvec3, kvec4, kvec5, kvec6, los, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		       * Z1_Bias(kvec2, los, f, b1) * Z1_Bias(kvec3, los, f, b1)
		       * Z1_Bias(kvec4, los, f, b1) * Z1_Bias(kvec5, los, f, b1) * Z1_Bias(kvec6, los, f, b1)
		       * pow(sigma8, 10) * f_pk(k2) * f_pk(k3) * f_pk(k4) * f_pk(k5) * f_pk(k6);


	return T511111;

}

double P6spectrum_Tree(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * kvec4_in, double * kvec5_in, double * kvec6_in,
		       double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha15 = alpha3 * alpha3 * alpha3 * alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };
	double kvec4[3] = { 0.0, 0.0, 0.0 };
	double kvec5[3] = { 0.0, 0.0, 0.0 };
	double kvec6[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
	calcTrueWavevector(kvec4_in, los, alpha_perp, alpha_parallel, kvec4);
	calcTrueWavevector(kvec5_in, los, alpha_perp, alpha_parallel, kvec5);
	calcTrueWavevector(kvec6_in, los, alpha_perp, alpha_parallel, kvec6);

	double T222211 = T6_222211(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);
	
	double T322111 = T6_322111(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	double T331111 = T6_331111(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	double T421111 = T6_421111(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	double T511111 = T6_511111(kvec1, kvec2, kvec3, kvec4, kvec5, kvec6, los, sigma8, f, b1);

	double P6 = 120.0 * T511111 + 48.0 * T421111 + 36.0 * T331111 + 24.0 * T322111 + 16.0 * T222211;
	return P6 / alpha15;

}


double Cov_PP_G(double * kvec, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double P = Powerspectrum_Tree(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (1.0/nmean);
	return ( P * P  ) / volume;
}

double Cov_PP_G_NL(double * kvec, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double sigma2_perp, double sigma2_para, double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double P = Powerspectrum_NonLinearFitting(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) + (1.0/nmean);
	return ( P * P  ) / volume;

}


double Cov_PP_KSZ_G(double * kvec, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double nmean, double volume,
              	    double aH_tau_T0_over_c1, double aH_tau_T0_over_c2, double sigma2_vv, double R2_N) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Pdv1 = Powerspectrum_Tree_KSZ(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1, aH_tau_T0_over_c1);
	double Pdv2 = Powerspectrum_Tree_KSZ(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1, aH_tau_T0_over_c2);
	double Pvv  = Powerspectrum_Tree_KSZ_VV(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1, aH_tau_T0_over_c1, aH_tau_T0_over_c2)
	            + aH_tau_T0_over_c1 * aH_tau_T0_over_c2 * (1.0 + R2_N) * (sigma2_vv/nmean);
	double P = Powerspectrum_Tree(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (1.0/nmean);

	return 2.0 * ( Pdv1 * Pdv2 + Pvv * P ) / volume;
}

double Cov_PP_GALAXY_KSZ_G(double * kvec, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double nmean, double volume, double aH_tau_T0_over_c) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Pdv = Powerspectrum_Tree_KSZ(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1, aH_tau_T0_over_c);
	double P = Powerspectrum_Tree(kvec, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (1.0/nmean);

	return 2.0 * ( Pdv * P ) / volume;
}


double Cov_PP_NG(double * kvec1, double * kvec2,
	         double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, 
		 double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};
	double kvec12[3] = PLUS(kvec1, kvec2);
	double Minus_kvec12[3] = { - kvec12[0], - kvec12[1], - kvec12[2]};

	double k1_2[3] = { kvec1[0] - kvec2[0], kvec1[1] - kvec2[1], kvec1[2] - kvec2[2] };
	double k2_1[3] = { kvec2[0] - kvec1[0], kvec2[1] - kvec1[1], kvec2[2] - kvec1[2] };

	double cov =  Trispectrum_Tree(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, b3, bK2, bK3, bDK, bO)
	           + (1.0/nmean) * ( Bispectrum_Tree(kvec1, kvec2, Minus_kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2)
	             	           + Bispectrum_Tree(Minus_kvec1, Minus_kvec2, kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
	             	           + Bispectrum_Tree(kvec1, Minus_kvec2, k2_1, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
	             	           + Bispectrum_Tree(Minus_kvec1, kvec2, k1_2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) )
	           + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(k1_2, los, alpha_perp, alpha_parallel, sigma8, f, b1) );

	return alpha3 * cov / volume;
}

double Cov_PP_NG_b2(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};
	double kvec12[3] = PLUS(kvec1, kvec2);
	double Minus_kvec12[3] = { - kvec12[0], - kvec12[1], - kvec12[2]};

	double k1_2[3] = { kvec1[0] - kvec2[0], kvec1[1] - kvec2[1], kvec1[2] - kvec2[2] };
	double k2_1[3] = { kvec2[0] - kvec1[0], kvec2[1] - kvec1[1], kvec2[2] - kvec1[2] };

	double cov =  Trispectrum_Tree_b2(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1)
	           + (1.0/nmean) * ( Bispectrum_Tree_b2(kvec1, kvec2, Minus_kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1)
	             	           + Bispectrum_Tree_b2(Minus_kvec1, Minus_kvec2, kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
	             	           + Bispectrum_Tree_b2(kvec1, Minus_kvec2, k2_1, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
	             	           + Bispectrum_Tree_b2(Minus_kvec1, kvec2, k1_2, los, alpha_perp, alpha_parallel, sigma8, f, b1) );

	return alpha3 * cov / volume;
}

double Cov_PP_NG_bK2(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};
	double kvec12[3] = PLUS(kvec1, kvec2);
	double Minus_kvec12[3] = { - kvec12[0], - kvec12[1], - kvec12[2]};

	double k1_2[3] = { kvec1[0] - kvec2[0], kvec1[1] - kvec2[1], kvec1[2] - kvec2[2] };
	double k2_1[3] = { kvec2[0] - kvec1[0], kvec2[1] - kvec1[1], kvec2[2] - kvec1[2] };

	double cov =  Trispectrum_Tree_bK2(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1)
	           + (1.0/nmean) * ( Bispectrum_Tree_bK2(kvec1, kvec2, Minus_kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1)
	             	           + Bispectrum_Tree_bK2(Minus_kvec1, Minus_kvec2, kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
	             	           + Bispectrum_Tree_bK2(kvec1, Minus_kvec2, k2_1, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
	             	           + Bispectrum_Tree_bK2(Minus_kvec1, kvec2, k1_2, los, alpha_perp, alpha_parallel, sigma8, f, b1) );

	return alpha3 * cov / volume;
}

double Cov_PP_NG_b2_b2(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double cov =  Trispectrum_Tree_b2_b2(kvec1, Minus_kvec1, kvec2, Minus_kvec2,
	             	                     los, alpha_perp, alpha_parallel, sigma8, f, b1);

	return alpha3 * cov / volume;
}

double Cov_PP_NG_b2_bK2(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double cov =  Trispectrum_Tree_b2_bK2(kvec1, Minus_kvec1, kvec2, Minus_kvec2,
	             	                     los, alpha_perp, alpha_parallel, sigma8, f, b1);

	return alpha3 * cov / volume;
}

double Cov_PP_NG_bK2_bK2(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double cov =  Trispectrum_Tree_bK2_bK2(kvec1, Minus_kvec1, kvec2, Minus_kvec2,
	             	                     los, alpha_perp, alpha_parallel, sigma8, f, b1);

	return alpha3 * cov / volume;
}


double Cov_PP_NG_b3(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double cov =  Trispectrum_Tree_b3(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	return alpha3 * cov / volume;
}

double Cov_PP_NG_bK3(
        double * kvec1, double * kvec2, double * los, 
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		double nmean, double volume
        ) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double cov =  Trispectrum_Tree_bK3(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	return alpha3 * cov / volume;
}


double Cov_PP_NG_bDK(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double cov =  Trispectrum_Tree_bDK(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	return alpha3 * cov / volume;
}

double Cov_PP_NG_bO(double * kvec1, double * kvec2,
	            double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double cov =  Trispectrum_Tree_bO(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	return alpha3 * cov / volume;
}

double Cov_PP_NG_Reconstructed(double * kvec1, double * kvec2,
	                       double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, 
		               double nmean, double volume, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};
	double kvec12[3] = PLUS(kvec1, kvec2);
	double Minus_kvec12[3] = { - kvec12[0], - kvec12[1], - kvec12[2]};

	double k1_2[3] = { kvec1[0] - kvec2[0], kvec1[1] - kvec2[1], kvec1[2] - kvec2[2] };
	double k2_1[3] = { kvec2[0] - kvec1[0], kvec2[1] - kvec1[1], kvec2[2] - kvec1[2] };

	double cov =  Trispectrum_Tree_Reconstructed(kvec1, Minus_kvec1, kvec2, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, b3, bK2, bK3, bDK, bO, b1_fid, R)

	           + (1.0/nmean) * ( Bispectrum_Tree_ReconForCovariance(Minus_kvec12,kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2, b1_fid, R)
	             	           + Bispectrum_Tree_ReconForCovariance(kvec12, Minus_kvec1, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2, b1_fid, R) 
	             	           + Bispectrum_Tree_ReconForCovariance(k2_1, kvec1, Minus_kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2, b1_fid, R) 
	             	           + Bispectrum_Tree_ReconForCovariance(k1_2, Minus_kvec1, kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2, b1_fid, R) )

	           + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree_ReconForCovariance(kvec12, los, alpha_perp, alpha_parallel, sigma8, f, b1, b1_fid, R) 
			                  + Powerspectrum_Tree_ReconForCovariance(k1_2, los, alpha_perp, alpha_parallel, sigma8, f, b1, b1_fid, R) );

	return alpha3 * cov / volume;
}

double Cov_PP_NG_BeatCoupling(
        double * kvec1, double * kvec3, double * evec, double * los, 
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double nmean, double volume) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double M_k13[3] = { - kvec1[0] - kvec3[0], - kvec1[1] - kvec3[1], - kvec1[2] - kvec3[2] };
	double M_k14[3] = { - kvec1[0] - kvec4[0], - kvec1[1] - kvec4[1], - kvec1[2] - kvec4[2] };
	double M_k23[3] = { - kvec2[0] - kvec3[0], - kvec2[1] - kvec3[1], - kvec2[2] - kvec3[2] };
	double M_k24[3] = { - kvec2[0] - kvec4[0], - kvec2[1] - kvec4[1], - kvec2[2] - kvec4[2] };

	double k13[3] = { kvec1[0] + kvec3[0], kvec1[1] + kvec3[1], kvec1[2] + kvec3[2] };
	double k14[3] = { kvec1[0] + kvec4[0], kvec1[1] + kvec4[1], kvec1[2] + kvec4[2] };

	double T =  Trispectrum_Tree(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, b3, bK2, bK3, bDK, bO)
	         + (1.0/nmean) * ( Bispectrum_Tree(kvec1, kvec3, M_k13, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
                             + Bispectrum_Tree(kvec1, kvec4, M_k14, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
                             + Bispectrum_Tree(kvec2, kvec3, M_k23, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
                             + Bispectrum_Tree(kvec2, kvec4, M_k24, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
				             )
	         + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(k13, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(k14, los, alpha_perp, alpha_parallel, sigma8, f, b1) );

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}

double Cov_PP_NG_BeatCoupling_b2(
        double * kvec1, double * kvec3, double * evec, double * los, 
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
        double nmean, double volume) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    nmean *= alpha3;
    
    double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
    double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};
    
    double M_k13[3] = { - kvec1[0] - kvec3[0], - kvec1[1] - kvec3[1], - kvec1[2] - kvec3[2] };
    double M_k14[3] = { - kvec1[0] - kvec4[0], - kvec1[1] - kvec4[1], - kvec1[2] - kvec4[2] };
    double M_k23[3] = { - kvec2[0] - kvec3[0], - kvec2[1] - kvec3[1], - kvec2[2] - kvec3[2] };
    double M_k24[3] = { - kvec2[0] - kvec4[0], - kvec2[1] - kvec4[1], - kvec2[2] - kvec4[2] };
    
    double T =  Trispectrum_Tree_b2(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1)
             + (1.0/nmean) * ( Bispectrum_Tree_b2(kvec1, kvec3, M_k13, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
                             + Bispectrum_Tree_b2(kvec1, kvec4, M_k14, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
                             + Bispectrum_Tree_b2(kvec2, kvec3, M_k23, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
                             + Bispectrum_Tree_b2(kvec2, kvec4, M_k24, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
    			             );
    
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double emag = NORM(evec);
    return alpha3 * W2(emag, R) * T;
}

double Cov_PP_NG_BeatCoupling_bK2(
        double * kvec1, double * kvec3, double * evec, double * los, 
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
        double nmean, double volume) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    nmean *= alpha3;
    
    double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
    double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};
    
    double M_k13[3] = { - kvec1[0] - kvec3[0], - kvec1[1] - kvec3[1], - kvec1[2] - kvec3[2] };
    double M_k14[3] = { - kvec1[0] - kvec4[0], - kvec1[1] - kvec4[1], - kvec1[2] - kvec4[2] };
    double M_k23[3] = { - kvec2[0] - kvec3[0], - kvec2[1] - kvec3[1], - kvec2[2] - kvec3[2] };
    double M_k24[3] = { - kvec2[0] - kvec4[0], - kvec2[1] - kvec4[1], - kvec2[2] - kvec4[2] };
    
    double T =  Trispectrum_Tree_bK2(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1)
             + (1.0/nmean) * ( Bispectrum_Tree_bK2(kvec1, kvec3, M_k13, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
                             + Bispectrum_Tree_bK2(kvec1, kvec4, M_k14, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
                             + Bispectrum_Tree_bK2(kvec2, kvec3, M_k23, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
                             + Bispectrum_Tree_bK2(kvec2, kvec4, M_k24, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
    			             );
    
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double emag = NORM(evec);
    return alpha3 * W2(emag, R) * T;
}


double Cov_PP_NG_BeatCoupling_b2_b2(
        double * kvec1, double * kvec3, double * evec, double * los,
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
        double nmean, double volume) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    nmean *= alpha3;
    
    double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
    double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};
    
    double T = Trispectrum_Tree_b2_b2(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1);
    
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double emag = NORM(evec);
    return alpha3 * W2(emag, R) * T;
}



double Cov_PP_NG_BeatCoupling_b2_bK2(
        double * kvec1, double * kvec3, double * evec, double * los,
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
        double nmean, double volume) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    nmean *= alpha3;
    
    double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
    double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};
    
    double T =  Trispectrum_Tree_b2_bK2(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1);
    
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double emag = NORM(evec);
    return alpha3 * W2(emag, R) * T;

}



double Cov_PP_NG_BeatCoupling_bK2_bK2(double * kvec1, double * kvec3, double * evec, 
		              double * los, double alpha_perp, double alpha_parallel, 
			      double sigma8, double f, double b1, double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double T =  Trispectrum_Tree_bK2_bK2(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}



double Cov_PP_NG_BeatCoupling_b3(double * kvec1, double * kvec3, double * evec, 
		              double * los, double alpha_perp, double alpha_parallel, 
			      double sigma8, double f, double b1, double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double T =  Trispectrum_Tree_b3(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}

double Cov_PP_NG_BeatCoupling_bK3(double * kvec1, double * kvec3, double * evec, 
		              double * los, double alpha_perp, double alpha_parallel, 
			      double sigma8, double f, double b1, double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double T =  Trispectrum_Tree_bK3(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}

double Cov_PP_NG_BeatCoupling_bDK(double * kvec1, double * kvec3, double * evec, 
		              double * los, double alpha_perp, double alpha_parallel, 
			      double sigma8, double f, double b1, double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double T =  Trispectrum_Tree_bDK(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}

double Cov_PP_NG_BeatCoupling_bO(double * kvec1, double * kvec3, double * evec, 
		              double * los, double alpha_perp, double alpha_parallel, 
			      double sigma8, double f, double b1, double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double T =  Trispectrum_Tree_bO(kvec1, kvec2, kvec3, kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}

double Cov_PP_NG_LocalMean(
            double * kvec1, double * kvec2, double * evec, double * los, 
            double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, 
            double b2, double b3, double bK2, double bK3, double bDK, double bO,
            double nmean, double volume) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    nmean *= alpha3;
    
    double kvec1_Minus_evec[3] = {kvec1[0] - evec[0], kvec1[1] - evec[1], kvec1[2] - evec[2]};
    double kvec2_Minus_evec[3] = {kvec2[0] - evec[0], kvec2[1] - evec[1], kvec2[2] - evec[2]};
    double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
    double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};
    
    double B1 = 
            Bispectrum_Tree(kvec1_Minus_evec, Minus_kvec1, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
            + (1.0/pow(nmean,1)) * ( Powerspectrum_Tree(kvec1_Minus_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1) );
    double B2 =
            Bispectrum_Tree(kvec2_Minus_evec, Minus_kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
            + (1.0/pow(nmean,1)) * ( Powerspectrum_Tree(kvec2_Minus_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1) );
    
    double T = 
            (-2.0) * B1 * Powerspectrum_Tree(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1)
            + (-2.0) * B2 * Powerspectrum_Tree(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1)
            + (+4.0) * ( Powerspectrum_Tree(evec, los, alpha_perp, alpha_parallel, sigma8, f, b1) + 1.0/nmean ) 
            * Powerspectrum_Tree(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
            * Powerspectrum_Tree(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1);
    
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double emag = NORM(evec);
    return alpha3 * W2(emag, R) * T;
    
}

double Cov_PP_NG_LocalMean_NL(
        double * kvec1, double * kvec2, double * evec, double * los,
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double sigma2_perp, double sigma2_para,
        double nmean, double volume) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    nmean *= alpha3;
    
    double kvec1_Minus_evec[3] = {kvec1[0] - evec[0], kvec1[1] - evec[1], kvec1[2] - evec[2]};
    double kvec2_Minus_evec[3] = {kvec2[0] - evec[0], kvec2[1] - evec[1], kvec2[2] - evec[2]};
    double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
    double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};
    
    double B1 =
            Bispectrum_Tree(kvec1_Minus_evec, Minus_kvec1, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
            + (1.0/pow(nmean,1)) * ( Powerspectrum_NonLinearFitting(kvec1_Minus_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) 
            + Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) );

    double B2
            = Bispectrum_Tree(kvec2_Minus_evec, Minus_kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
            + (1.0/pow(nmean,1)) * ( Powerspectrum_NonLinearFitting(kvec2_Minus_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
            + Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) );
    
    double T
            = (-2.0) * B1 * Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
            + (-2.0) * B2 * Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
            + (+4.0) * ( Powerspectrum_NonLinearFitting(evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) + 1.0/nmean ) 
            * Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) 
            * Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para);
    
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double emag = NORM(evec);
    return alpha3 * W2(emag, R) * T;
    
}

double Cov_PP_NG_LocalMean_NL_Sigma2B(
        double * kvec1, double * kvec2, double * evec, double * los,
        double alpha_perp, double alpha_parallel, double sigma8, double f, double b1,
        double b2, double b3, double bK2, double bK3, double bDK, double bO, 
        double sigma2_perp, double sigma2_para,
        double nmean, double volume, double sigma2_b) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    nmean *= alpha3;
    
    double kvec1_Minus_evec[3] = {kvec1[0] - evec[0], kvec1[1] - evec[1], kvec1[2] - evec[2]};
    double kvec2_Minus_evec[3] = {kvec2[0] - evec[0], kvec2[1] - evec[1], kvec2[2] - evec[2]};
    double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
    double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};
    
    double B1 =
            Bispectrum_Tree(kvec1_Minus_evec, Minus_kvec1, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
            + (1.0/pow(nmean,1)) * ( Powerspectrum_NonLinearFitting(kvec1_Minus_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) 
            + Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) );

    double B2
            = Bispectrum_Tree(kvec2_Minus_evec, Minus_kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
            + (1.0/pow(nmean,1)) * ( Powerspectrum_NonLinearFitting(kvec2_Minus_evec, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
            + Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) );
    
    double T
            = (-2.0) * B1 * Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
            + (-2.0) * B2 * Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
            + (+4.0) * ( sigma2_b * volume ) 
            * Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para) 
            * Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para);
    
    double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
    double emag = NORM(evec);
    return alpha3 * W2(emag, R) * T;
    
}


double Cov_PP_NG_LocalMean_NL_b2(double * kvec1, double * kvec2, double * evec, 
		                 double * los, double alpha_perp, double alpha_parallel,
				 double sigma8, double f, double b1, double sigma2_perp, double sigma2_para,
				 double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec1_Minus_evec[3] = {kvec1[0] - evec[0], kvec1[1] - evec[1], kvec1[2] - evec[2]};
	double kvec2_Minus_evec[3] = {kvec2[0] - evec[0], kvec2[1] - evec[1], kvec2[2] - evec[2]};
	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double B1 = Bispectrum_Tree_b2(kvec1_Minus_evec, Minus_kvec1, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double B2 = Bispectrum_Tree_b2(kvec2_Minus_evec, Minus_kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double T = (-2.0) * B1 * Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
	         + (-2.0) * B2 * Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para);

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;

}

double Cov_PP_NG_LocalMean_NL_bK2(double * kvec1, double * kvec2, double * evec, 
		                 double * los, double alpha_perp, double alpha_parallel,
				 double sigma8, double f, double b1, double sigma2_perp, double sigma2_para,
				 double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec1_Minus_evec[3] = {kvec1[0] - evec[0], kvec1[1] - evec[1], kvec1[2] - evec[2]};
	double kvec2_Minus_evec[3] = {kvec2[0] - evec[0], kvec2[1] - evec[1], kvec2[2] - evec[2]};
	double Minus_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2]};
	double Minus_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2]};

	double B1 = Bispectrum_Tree_bK2(kvec1_Minus_evec, Minus_kvec1, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double B2 = Bispectrum_Tree_bK2(kvec2_Minus_evec, Minus_kvec2, evec, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double T = (-2.0) * B1 * Powerspectrum_NonLinearFitting(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para)
	         + (-2.0) * B2 * Powerspectrum_NonLinearFitting(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1, sigma2_perp, sigma2_para);

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;

}

double Cov_PP_SN(double * kvec1, double * kvec3, double * evec, 
		 double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, 
		 double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double M_k12[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	double M_k34[3] = { - kvec3[0] - kvec4[0], - kvec3[1] - kvec4[1], - kvec3[2] - kvec4[2] };

	double k12[3] = { kvec1[0] + kvec2[0], kvec1[1] + kvec2[1], kvec1[2] + kvec2[2] };

	double T =  (1.0/nmean) * ( Bispectrum_Tree(kvec1, kvec2, M_k12, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2)
	           	          + Bispectrum_Tree(kvec3, kvec4, M_k34, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
			 	  )
	         + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1) )
	         + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(kvec3, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1) )
	         + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(k12, los, alpha_perp, alpha_parallel, sigma8, f, b1) )
		 + (1.0/pow(nmean,3));

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}


double Cov_PP_SN_LocalMean(double * kvec1, double * kvec3, double * evec, 
		    double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO, 
		    double nmean, double volume) {

        double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec2[3] = { - kvec1[0] + evec[0], - kvec1[1] + evec[1], - kvec1[2] + evec[2]};
	double kvec4[3] = { - kvec3[0] - evec[0], - kvec3[1] - evec[1], - kvec3[2] - evec[2]};

	double M_k12[3] = { - kvec1[0] - kvec2[0], - kvec1[1] - kvec2[1], - kvec1[2] - kvec2[2] };
	double M_k34[3] = { - kvec3[0] - kvec4[0], - kvec3[1] - kvec4[1], - kvec3[2] - kvec4[2] };

	double k12[3] = { kvec1[0] + kvec2[0], kvec1[1] + kvec2[1], kvec1[2] + kvec2[2] };

	double T_GM =  (1.0/nmean) * ( Bispectrum_Tree(kvec1, kvec2, M_k12, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2)
	           	          + Bispectrum_Tree(kvec3, kvec4, M_k34, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2) 
			 	  )
	            + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1) )
	            + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(kvec3, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(kvec4, los, alpha_perp, alpha_parallel, sigma8, f, b1) )
	            + (1.0/pow(nmean,2)) * ( Powerspectrum_Tree(k12, los, alpha_perp, alpha_parallel, sigma8, f, b1) )
		    + (1.0/pow(nmean,3));

	double T_LM = (2.0 / nmean) * ( Powerspectrum_Tree(evec, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (1.0/nmean) ) 
	            * ( Powerspectrum_Tree(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1) + Powerspectrum_Tree(kvec3, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (2.0/nmean) );

	double T = - T_GM + T_LM;

	double R = pow( 3.0 * volume / (4.0 * M_PI), 1.0/3.0);
	double emag = NORM(evec);
	return alpha3 * W2(emag, R) * T;
}

double Cov_PB_NG_PB_Reconstructed(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f,
              	                  double b1, double nmean, double volume, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double P1 = Powerspectrum_Tree(kvec1_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double P2 = Powerspectrum_Tree_ReconForPBCrossCovariance(kvec1_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b1_fid, R);
	double P3 = Powerspectrum_Tree_ReconForPBCrossCovariance(kvec1_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b1_fid, R);

	double P_N = P1 + (1.0/nmean);
	double B_N = (1.0/nmean) * (P2 + P3);

	return ( 2.0 * P_N * B_N  ) / volume;
}


double Cov_PB_NG_PB(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f,
              	    double b1, double b2, double bK2,  double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double P1 = Powerspectrum_Tree(kvec1_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double P2 = Powerspectrum_Tree(kvec2_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double P3 = Powerspectrum_Tree(kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double B = Bispectrum_Tree(kvec1_in, kvec2_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);

	double P_N = P1 + (1.0/nmean);
	double B_N = B + (1.0/nmean) * (P2 + P3);

	return ( 2.0 * P_N * B_N  ) / volume;
}

double Cov_PB_NG_P5(double * kvec_in, double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f,
              	    double b1, double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kk1_in[3] = PLUS(kvec_in, kvec1_in);
	double kk2_in[3] = PLUS(kvec_in, kvec2_in);
	double kk3_in[3] = PLUS(kvec_in, kvec3_in);

	double Mkvec_in[3] = { - kvec_in[0], - kvec_in[1], - kvec_in[2] };

	double Mkk1_in[3] = PLUS(Mkvec_in, kvec1_in);
	double Mkk2_in[3] = PLUS(Mkvec_in, kvec2_in);
	double Mkk3_in[3] = PLUS(Mkvec_in, kvec3_in);

	double P5 = P5spectrum_Tree(kvec_in, Mkvec_in, kvec1_in, kvec2_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double T4 = Trispectrum_Tree(kk1_in, Mkvec_in, kvec2_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		  + Trispectrum_Tree(kk2_in, Mkvec_in, kvec1_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		  + Trispectrum_Tree(kk3_in, Mkvec_in, kvec1_in, kvec2_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		  + Trispectrum_Tree(Mkk1_in, kvec_in, kvec2_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		  + Trispectrum_Tree(Mkk2_in, kvec_in, kvec1_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		  + Trispectrum_Tree(Mkk3_in, kvec_in, kvec1_in, kvec2_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

	double B = Bispectrum_Tree(kk1_in, Mkk2_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0) 
		 + Bispectrum_Tree(kk1_in, Mkk3_in, kvec2_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0) 
		 + Bispectrum_Tree(kk2_in, Mkk1_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0)
		 + Bispectrum_Tree(kk2_in, Mkk3_in, kvec1_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0)
		 + Bispectrum_Tree(kk3_in, Mkk1_in, kvec2_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0)
		 + Bispectrum_Tree(kk3_in, Mkk2_in, kvec1_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0);

	double cov = P5 + T4 / nmean + B / pow(nmean, 2); 
	
	return alpha3 * cov / volume;

}

double Cov_BB_G_PPP(double * kvec1, double * kvec2, double * kvec3, double * los, double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double P1 = Powerspectrum_Tree(kvec1, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (1.0/nmean);
	double P2 = Powerspectrum_Tree(kvec2, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (1.0/nmean);
	double P3 = Powerspectrum_Tree(kvec3, los, alpha_perp, alpha_parallel, sigma8, f, b1) + (1.0/nmean);

	return ( P1 * P2 * P3 ) / volume / alpha3;
}

double Cov_BB_NG_BB(double * kvec1_in, double * kvec2_in, double * kvec3_in, 
	            double * kvec1_dash_in, double * kvec2_dash_in, double * kvec3_dash_in, double * los, 
                    double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double b2, double bK2,
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double B = Bispectrum_Tree(kvec1_in, kvec2_in, kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);
	double P2 = Powerspectrum_Tree(kvec2_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double P3 = Powerspectrum_Tree(kvec3_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double B_N = B + (1.0/nmean) * (P2 + P3);

	double B_dash = Bispectrum_Tree(kvec1_dash_in, kvec2_dash_in, kvec3_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);
	double P2_dash = Powerspectrum_Tree(kvec2_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double P3_dash = Powerspectrum_Tree(kvec3_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double B_N_dash = B_dash + (1.0/nmean) * (P2_dash + P3_dash);

	return ( B_N * B_N_dash ) / volume;

}

double Cov_BB_NG_PT(double * kvec1_in, double * kvec2_in, double * kvec3_in, 
	            double * kvec2_dash_in, double * kvec3_dash_in, double * los, 
                    double alpha_perp, double alpha_parallel, double sigma8, double f, 
		    double b1, double b2, double b3, double bK2, double bK3, double bDK, double bO,
		    double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double P = Powerspectrum_Tree(kvec1_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double P_N = P + (1.0/nmean);

	double T = Trispectrum_Tree(kvec2_in, kvec3_in, kvec2_dash_in, kvec3_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, b3, bK2, bK3, bDK, bO);
	
	double kvec22_dash_in[3] = { kvec2_in[0] + kvec2_dash_in[0], kvec2_in[1] + kvec2_dash_in[1], kvec2_in[2] + kvec2_dash_in[2] }; 
	double kvec23_dash_in[3] = { kvec2_in[0] + kvec3_dash_in[0], kvec2_in[1] + kvec3_dash_in[1], kvec2_in[2] + kvec3_dash_in[2] }; 
	double kvec32_dash_in[3] = { kvec3_in[0] + kvec2_dash_in[0], kvec3_in[1] + kvec2_dash_in[1], kvec3_in[2] + kvec2_dash_in[2] }; 
	double kvec33_dash_in[3] = { kvec3_in[0] + kvec3_dash_in[0], kvec3_in[1] + kvec3_dash_in[1], kvec3_in[2] + kvec3_dash_in[2] }; 

	double M_kvec2_in[3] = MINUS_XVEC(kvec2_in);
	double M_kvec3_in[3] = MINUS_XVEC(kvec3_in);
	double M_kvec2_dash_in[3] = MINUS_XVEC(kvec2_dash_in);
	double M_kvec3_dash_in[3] = MINUS_XVEC(kvec3_dash_in);

	double Ba = Bispectrum_Tree(kvec22_dash_in, M_kvec2_in, M_kvec2_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);
	double Bb = Bispectrum_Tree(kvec23_dash_in, M_kvec2_in, M_kvec3_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);
	double Bc = Bispectrum_Tree(kvec32_dash_in, M_kvec3_in, M_kvec2_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);
	double Bd = Bispectrum_Tree(kvec33_dash_in, M_kvec3_in, M_kvec3_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, b2, bK2);

	double Pa = Powerspectrum_Tree(kvec22_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);
	double Pb = Powerspectrum_Tree(kvec23_dash_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double T_N = T + (1.0/nmean) * (Ba + Bb + Bc +Bd) + (1.0/pow(nmean,2)) * (Pa + Pb);

	return ( P_N * T_N ) / volume;

}

double Cov_BB_NG_P6(double * kvec1_in, double * kvec2_in, double * kvec3_in, 
	            double * kvec4_in, double * kvec5_in, double * kvec6_in, double * los, 
                    double alpha_perp, double alpha_parallel, double sigma8, double f, double b1, double nmean, double volume) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	nmean *= alpha3;

	double kvec14_in[3] = PLUS(kvec1_in, kvec4_in);
	double kvec15_in[3] = PLUS(kvec1_in, kvec5_in);
	double kvec16_in[3] = PLUS(kvec1_in, kvec6_in);

	double kvec24_in[3] = PLUS(kvec2_in, kvec4_in);
	double kvec25_in[3] = PLUS(kvec2_in, kvec5_in);
	double kvec26_in[3] = PLUS(kvec2_in, kvec6_in);

	double kvec34_in[3] = PLUS(kvec3_in, kvec4_in);
	double kvec35_in[3] = PLUS(kvec3_in, kvec5_in);
	double kvec36_in[3] = PLUS(kvec3_in, kvec6_in);

	double P6 = P6spectrum_Tree(kvec1_in, kvec2_in, kvec3_in, kvec4_in, kvec5_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double P5 = P5spectrum_Tree(kvec14_in, kvec2_in, kvec3_in, kvec5_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1)
	          + P5spectrum_Tree(kvec15_in, kvec2_in, kvec3_in, kvec4_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
		  + P5spectrum_Tree(kvec16_in, kvec2_in, kvec3_in, kvec4_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
		  + P5spectrum_Tree(kvec24_in, kvec1_in, kvec3_in, kvec5_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
		  + P5spectrum_Tree(kvec25_in, kvec1_in, kvec3_in, kvec4_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
		  + P5spectrum_Tree(kvec26_in, kvec1_in, kvec3_in, kvec4_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1)
		  + P5spectrum_Tree(kvec34_in, kvec1_in, kvec2_in, kvec5_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1) 
		  + P5spectrum_Tree(kvec35_in, kvec1_in, kvec2_in, kvec4_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1)
		  + P5spectrum_Tree(kvec36_in, kvec1_in, kvec2_in, kvec4_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1);

	double T = Trispectrum_Tree(kvec14_in, kvec25_in, kvec3_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec14_in, kvec36_in, kvec2_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec14_in, kvec26_in, kvec3_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec14_in, kvec35_in, kvec2_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec15_in, kvec24_in, kvec3_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		 + Trispectrum_Tree(kvec15_in, kvec36_in, kvec2_in, kvec4_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec15_in, kvec26_in, kvec3_in, kvec4_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		 + Trispectrum_Tree(kvec15_in, kvec34_in, kvec2_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec16_in, kvec25_in, kvec3_in, kvec4_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		 + Trispectrum_Tree(kvec16_in, kvec34_in, kvec2_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec16_in, kvec24_in, kvec3_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		 + Trispectrum_Tree(kvec16_in, kvec35_in, kvec2_in, kvec4_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec25_in, kvec36_in, kvec1_in, kvec4_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		 + Trispectrum_Tree(kvec26_in, kvec35_in, kvec1_in, kvec4_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec24_in, kvec36_in, kvec1_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		 + Trispectrum_Tree(kvec26_in, kvec34_in, kvec1_in, kvec5_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
		 + Trispectrum_Tree(kvec24_in, kvec35_in, kvec1_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) 
		 + Trispectrum_Tree(kvec25_in, kvec34_in, kvec1_in, kvec6_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

	double B = Bispectrum_Tree(kvec14_in, kvec25_in, kvec36_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0)
		 + Bispectrum_Tree(kvec14_in, kvec26_in, kvec35_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0)
		 + Bispectrum_Tree(kvec15_in, kvec24_in, kvec36_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0) 
		 + Bispectrum_Tree(kvec15_in, kvec26_in, kvec34_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0)
		 + Bispectrum_Tree(kvec16_in, kvec24_in, kvec35_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0) 
		 + Bispectrum_Tree(kvec16_in, kvec25_in, kvec34_in, los, alpha_perp, alpha_parallel, sigma8, f, b1, 0.0, 0.0);

	double cov = P6 + P5 / nmean + T / pow(nmean,2) + B / pow(nmean,3);

	return alpha3 * cov / volume;

}

/*****************************/
/* decomposed power spectrum */
/*****************************/

//double Z1_Bias(double * kvec1, double * los, double f, double b1) {
//	double mu = MU(kvec1, los);
//	return (b1) * D1()
//	     + (f)  * V1(kvec1, los);
//}
//
//
//double Z2_Bias(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2) {
//	return  (f) * V2(kvec1, kvec2, los)
//	      + (f * f) * ( V1V1(kvec1, kvec2, los) / 2.0)
//	      + (b1 * f) * D1V1(kvec1, kvec2, los) 
//	      
//	      +  (b1)  * F2(kvec1, kvec2) 
//	      +  (b2)  * (D1D1(kvec1, kvec2) / 2.0)
//	      +  (bK2) * K1K1(kvec1, kvec2);
//}

//double Z2_Bias_Reconstructed(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2, double b1_fid, double R) {
//	return  (f) * V2(kvec1, kvec2, los)
//	      + (f * f) * ( V1V1(kvec1, kvec2, los) / 2.0)
//	      + (b1 * f) * D1V1(kvec1, kvec2, los) 
//	      
//	      +  (b1)  * F2(kvec1, kvec2) 
//	      +  (b2)  * (D1D1(kvec1, kvec2) / 2.0)
//	      +  (bK2) * K1K1(kvec1, kvec2)
//
//	      + (b1 * b1) * Z1S1_b1_b1(kvec1, kvec2, los, b1_fid, R) 
//	      + (b1 * f ) * Z1S1_b1_f( kvec1, kvec2, los, b1_fid, R) 
//	      + (f  * f ) * Z1S1_f_f(  kvec1, kvec2, los, b1_fid, R) 
//}


/********/
/* Tree */
/********/

double Powerspectrum_Tree_b1_b1(double * kvec_in, double * los, double alpha_perp, double alpha_parallel) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double P = D1() * D1()  * f_pk(k);
	return P / alpha3;
}

double Powerspectrum_Tree_b1_f(double * kvec_in, double * los, double alpha_perp, double alpha_parallel) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double P = 2.0 * D1() * V1(kvec, los)  * f_pk(k);
	return P / alpha3;
}

double Powerspectrum_Tree_f_f(double * kvec_in, double * los, double alpha_perp, double alpha_parallel) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double P = V1(kvec, los) * V1(kvec, los) * f_pk(k);
	return P / alpha3;
}


/*******/
/* BAO */
/*******/

double Powerspectrum_Tree_BAO_b1_b1(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
    double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
     	
	double BAO = f_pk(k) - f_pk_no_wiggle(k);

	double Kaiser = D1() * D1();

    double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);
	return (G+MC) / alpha3;

}

double Powerspectrum_Tree_BAO_b1_f(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
     	
	double BAO = f_pk(k) - f_pk_no_wiggle(k);

	double Kaiser = 2.0 * D1() * V1(kvec, los);

    double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);
	return (G+MC) / alpha3;
}

double Powerspectrum_Tree_BAO_f_f(double * kvec_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
     	double D = ExpDamping(kvec, los, sigma2_perp, sigma2_para);
     	
	double BAO = f_pk(k) - f_pk_no_wiggle(k);

	double Kaiser = V1(kvec, los) * V1(kvec, los);

     	double G = D * D * Kaiser * BAO;
	double MC = Kaiser * f_pk_no_wiggle(k);
	return (G+MC) / alpha3;
}

/*****************/
/* Tree_NoWiggle */
/*****************/

double Powerspectrum_Tree_NoWiggle_b1_b1(double * kvec_in, double * los, double alpha_perp, double alpha_parallel) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double P = D1() * D1()  * f_pk_no_wiggle(k);
	return P / alpha3;
}

double Powerspectrum_Tree_NoWiggle_b1_f(double * kvec_in, double * los, double alpha_perp, double alpha_parallel) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double P = 2.0 * D1() * V1(kvec, los)  * f_pk_no_wiggle(k);
	return P / alpha3;
}

double Powerspectrum_Tree_NoWiggle_f_f(double * kvec_in, double * los, double alpha_perp, double alpha_parallel) {
	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double kvec[3] = { 0.0, 0.0, 0.0 };
	calcTrueWavevector(kvec_in, los, alpha_perp, alpha_parallel, kvec);
	double k = NORM(kvec);
	double P = V1(kvec, los) * V1(kvec, los) * f_pk_no_wiggle(k);
	return P / alpha3;
}

/*************************/
/* decomposed bispectrum */
/*************************/

/********/
/* Tree */
/********/

/*****************************************/
/* 6 * 2 * 2 = 24 terms - 10 = 14 terms. */
/*****************************************/

double Bispectrum_Tree_b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2(kvec1, kvec2) * D1() * D1();
	double K13 = 2.0 * F2(kvec1, kvec3) * D1() * D1();
	double K23 = 2.0 * F2(kvec2, kvec3) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2(kvec1, kvec2) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1()) + 2.0 * V2(kvec1, kvec2, los) * D1() * D1();
	double K13 = 2.0 * F2(kvec1, kvec3) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1()) + 2.0 * V2(kvec1, kvec3, los) * D1() * D1();
	double K23 = 2.0 * F2(kvec2, kvec3) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1()) + 2.0 * V2(kvec2, kvec3, los) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los) + 2.0 * V2(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * F2(kvec1, kvec3) * V1(kvec1, los) * V1(kvec3, los) + 2.0 * V2(kvec1, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * F2(kvec2, kvec3) * V1(kvec2, los) * V1(kvec3, los) + 2.0 * V2(kvec2, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

/*********/
double Bispectrum_Tree_b2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * D1() * D1();
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * D1() * D1();
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_b2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1());
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1());
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1());

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_b2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

/********/

double Bispectrum_Tree_bK2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * D1() * D1();
	double K13 = 2.0 * K1K1(kvec1, kvec3) * D1() * D1();
	double K23 = 2.0 * K1K1(kvec2, kvec3) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_bK2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * K1K1(kvec1, kvec3) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * K1K1(kvec2, kvec3) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_bK2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * K1K1(kvec1, kvec3) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * K1K1(kvec2, kvec3) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

/************/

double Bispectrum_Tree_b1f_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * D1() * D1();
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * D1() * D1();
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_b1f_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	K12 += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * D1() * D1();
	K13 += 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * D1() * D1();
	K23 += 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_b1f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * V1(kvec2, los) * V1(kvec3, los);

	K12 += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	K13 += 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	K23 += 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

/*********/

double Bispectrum_Tree_ff_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

/*********/

double Bispectrum_Tree_f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * V2(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * V2(kvec1, kvec3, los) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * V2(kvec2, kvec3, los) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_Reconstructed_b1b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_b1(kvec1, kvec2, b1_fid, R) * D1() * D1();
	double K13 = 2.0 * Z1S1_b1_b1(kvec1, kvec3, b1_fid, R) * D1() * D1();
	double K23 = 2.0 * Z1S1_b1_b1(kvec2, kvec3, b1_fid, R) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_Reconstructed_b1b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_b1(kvec1, kvec2, b1_fid, R) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1()) 
	           + 2.0 * Z1S1_b1_f(kvec1, kvec2, los, b1_fid, R) * D1() * D1() ;
	double K13 = 2.0 * Z1S1_b1_b1(kvec1, kvec3, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1())
	           + 2.0 * Z1S1_b1_f(kvec1, kvec3, los, b1_fid, R) * D1() * D1() ;
	double K23 = 2.0 * Z1S1_b1_b1(kvec2, kvec3, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1()) 
	           + 2.0 * Z1S1_b1_f(kvec2, kvec3, los, b1_fid, R) * D1() * D1() ;

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_Reconstructed_b1b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_b1(kvec1, kvec2, b1_fid, R) * V1(kvec1, los) * V1(kvec2, los)
	           + 2.0 * Z1S1_b1_f(kvec1, kvec2, los, b1_fid, R) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1())
	           + 2.0 * Z1S1_f_f(kvec1, kvec2, los, b1_fid, R) * D1() * D1();
	double K13 = 2.0 * Z1S1_b1_b1(kvec1, kvec3, b1_fid, R) * V1(kvec1, los) * V1(kvec3, los) 
	           + 2.0 * Z1S1_b1_f(kvec1, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1())
	           + 2.0 * Z1S1_f_f(kvec1, kvec3, los, b1_fid, R) * D1() * D1();
	double K23 = 2.0 * Z1S1_b1_b1(kvec2, kvec3, b1_fid, R) * V1(kvec2, los) * V1(kvec3, los) 
	           + 2.0 * Z1S1_b1_f(kvec2, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1())
		   + 2.0 * Z1S1_f_f(kvec2, kvec3, los, b1_fid, R) * D1() * D1();

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_Reconstructed_b1f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_f(kvec1, kvec2, los, b1_fid, R) * V1(kvec1, los) * V1(kvec2, los)
	           + 2.0 * Z1S1_f_f(kvec1, kvec2, los, b1_fid, R) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1());
	double K13 = 2.0 * Z1S1_b1_f(kvec1, kvec3, los, b1_fid, R) * V1(kvec1, los) * V1(kvec3, los)
	           + 2.0 * Z1S1_f_f(kvec1, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1());
	double K23 = 2.0 * Z1S1_b1_f(kvec2, kvec3, los, b1_fid, R) * V1(kvec2, los) * V1(kvec3, los)
		   + 2.0 * Z1S1_f_f(kvec2, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1());

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}

double Bispectrum_Tree_Reconstructed_ff_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_f_f(kvec1, kvec2, los, b1_fid, R) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * Z1S1_f_f(kvec1, kvec3, los, b1_fid, R) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * Z1S1_f_f(kvec2, kvec3, los, b1_fid, R) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk(k1) * f_pk(k2)
	         + K13 * f_pk(k1) * f_pk(k3)
	         + K23 * f_pk(k2) * f_pk(k3);

	return B / alpha6;

}


/****************/
/*** NoWiggle ***/
/****************/
double Bispectrum_Tree_NoWiggle_b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2(kvec1, kvec2) * D1() * D1();
	double K13 = 2.0 * F2(kvec1, kvec3) * D1() * D1();
	double K23 = 2.0 * F2(kvec2, kvec3) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2(kvec1, kvec2) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1()) + 2.0 * V2(kvec1, kvec2, los) * D1() * D1();
	double K13 = 2.0 * F2(kvec1, kvec3) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1()) + 2.0 * V2(kvec1, kvec3, los) * D1() * D1();
	double K23 = 2.0 * F2(kvec2, kvec3) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1()) + 2.0 * V2(kvec2, kvec3, los) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * F2(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los) + 2.0 * V2(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * F2(kvec1, kvec3) * V1(kvec1, los) * V1(kvec3, los) + 2.0 * V2(kvec1, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * F2(kvec2, kvec3) * V1(kvec2, los) * V1(kvec3, los) + 2.0 * V2(kvec2, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

/*********/
double Bispectrum_Tree_NoWiggle_b2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * D1() * D1();
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * D1() * D1();
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_b2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1());
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1());
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1());

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_b2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

/********/

double Bispectrum_Tree_NoWiggle_bK2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * D1() * D1();
	double K13 = 2.0 * K1K1(kvec1, kvec3) * D1() * D1();
	double K23 = 2.0 * K1K1(kvec2, kvec3) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_bK2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * K1K1(kvec1, kvec3) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * K1K1(kvec2, kvec3) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_bK2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * K1K1(kvec1, kvec3) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * K1K1(kvec2, kvec3) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

/************/

double Bispectrum_Tree_NoWiggle_b1f_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * D1() * D1();
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * D1() * D1();
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_b1f_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	K12 += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * D1() * D1();
	K13 += 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * D1() * D1();
	K23 += 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_b1f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * V1(kvec2, los) * V1(kvec3, los);

	K12 += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	K13 += 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	K23 += 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

/*********/

double Bispectrum_Tree_NoWiggle_ff_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

/*********/

double Bispectrum_Tree_NoWiggle_f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * V2(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * V2(kvec1, kvec3, los) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * V2(kvec2, kvec3, los) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_Reconstructed_b1b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_b1(kvec1, kvec2, b1_fid, R) * D1() * D1();
	double K13 = 2.0 * Z1S1_b1_b1(kvec1, kvec3, b1_fid, R) * D1() * D1();
	double K23 = 2.0 * Z1S1_b1_b1(kvec2, kvec3, b1_fid, R) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_Reconstructed_b1b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_b1(kvec1, kvec2, b1_fid, R) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1()) 
	           + 2.0 * Z1S1_b1_f(kvec1, kvec2, los, b1_fid, R) * D1() * D1() ;
	double K13 = 2.0 * Z1S1_b1_b1(kvec1, kvec3, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1())
	           + 2.0 * Z1S1_b1_f(kvec1, kvec3, los, b1_fid, R) * D1() * D1() ;
	double K23 = 2.0 * Z1S1_b1_b1(kvec2, kvec3, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1()) 
	           + 2.0 * Z1S1_b1_f(kvec2, kvec3, los, b1_fid, R) * D1() * D1() ;

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_Reconstructed_b1b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_b1(kvec1, kvec2, b1_fid, R) * V1(kvec1, los) * V1(kvec2, los)
	           + 2.0 * Z1S1_b1_f(kvec1, kvec2, los, b1_fid, R) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1())
	           + 2.0 * Z1S1_f_f(kvec1, kvec2, los, b1_fid, R) * D1() * D1();
	double K13 = 2.0 * Z1S1_b1_b1(kvec1, kvec3, b1_fid, R) * V1(kvec1, los) * V1(kvec3, los) 
	           + 2.0 * Z1S1_b1_f(kvec1, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1())
	           + 2.0 * Z1S1_f_f(kvec1, kvec3, los, b1_fid, R) * D1() * D1();
	double K23 = 2.0 * Z1S1_b1_b1(kvec2, kvec3, b1_fid, R) * V1(kvec2, los) * V1(kvec3, los) 
	           + 2.0 * Z1S1_b1_f(kvec2, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1())
		   + 2.0 * Z1S1_f_f(kvec2, kvec3, los, b1_fid, R) * D1() * D1();

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_Reconstructed_b1f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_b1_f(kvec1, kvec2, los, b1_fid, R) * V1(kvec1, los) * V1(kvec2, los)
	           + 2.0 * Z1S1_f_f(kvec1, kvec2, los, b1_fid, R) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1());
	double K13 = 2.0 * Z1S1_b1_f(kvec1, kvec3, los, b1_fid, R) * V1(kvec1, los) * V1(kvec3, los)
	           + 2.0 * Z1S1_f_f(kvec1, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1());
	double K23 = 2.0 * Z1S1_b1_f(kvec2, kvec3, los, b1_fid, R) * V1(kvec2, los) * V1(kvec3, los)
		   + 2.0 * Z1S1_f_f(kvec2, kvec3, los, b1_fid, R) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1());

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}

double Bispectrum_Tree_NoWiggle_Reconstructed_ff_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double b1_fid, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);

	double K12 = 2.0 * Z1S1_f_f(kvec1, kvec2, los, b1_fid, R) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * Z1S1_f_f(kvec1, kvec3, los, b1_fid, R) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * Z1S1_f_f(kvec2, kvec3, los, b1_fid, R) * V1(kvec2, los) * V1(kvec3, los);

	double B = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
	         + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
	         + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);

	return B / alpha6;

}


/***************/
/***** BAO *****/
/***************/

double Bispectrum_Tree_BAO_b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

    double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
    double alpha6 = alpha3 * alpha3;
    
    double kvec1[3] = { 0.0, 0.0, 0.0 };
    double kvec2[3] = { 0.0, 0.0, 0.0 };
    double kvec3[3] = { 0.0, 0.0, 0.0 };
    
    calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
    calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
    calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);
    
    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
    double k3 = NORM(kvec3);
    double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
    double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
    double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);
    
    double K12 = 2.0 * F2(kvec1, kvec2) * D1() * D1();
    double K13 = 2.0 * F2(kvec1, kvec3) * D1() * D1();
    double K23 = 2.0 * F2(kvec2, kvec3) * D1() * D1();
    
    double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
    double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
    double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
    double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
              + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
              + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
    
    double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
                 + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
    
                 + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
                 + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
                 
                 + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
                 + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
    
    double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
              + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
              + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
    
    return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * F2(kvec1, kvec2) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1()) + 2.0 * V2(kvec1, kvec2, los) * D1() * D1();
	double K13 = 2.0 * F2(kvec1, kvec3) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1()) + 2.0 * V2(kvec1, kvec3, los) * D1() * D1();
	double K23 = 2.0 * F2(kvec2, kvec3) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1()) + 2.0 * V2(kvec2, kvec3, los) * D1() * D1();

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
        double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
 
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * F2(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los) + 2.0 * V2(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * F2(kvec1, kvec3) * V1(kvec1, los) * V1(kvec3, los) + 2.0 * V2(kvec1, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * F2(kvec2, kvec3) * V1(kvec2, los) * V1(kvec3, los) + 2.0 * V2(kvec2, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}


double Bispectrum_Tree_BAO_b2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * D1() * D1();
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * D1() * D1();
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * D1() * D1();

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     

     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_b2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * (D1() * V1(kvec2, los) + V1(kvec1, los) * D1());
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * (D1() * V1(kvec3, los) + V1(kvec1, los) * D1());
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * (D1() * V1(kvec3, los) + V1(kvec2, los) * D1());

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
     	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     

     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_b2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * (D1D1(kvec1, kvec2)/2.0) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * (D1D1(kvec1, kvec3)/2.0) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * (D1D1(kvec2, kvec3)/2.0) * V1(kvec2, los) * V1(kvec3, los);

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
      	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}


double Bispectrum_Tree_BAO_bK2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * D1() * D1();
	double K13 = 2.0 * K1K1(kvec1, kvec3) * D1() * D1();
	double K23 = 2.0 * K1K1(kvec2, kvec3) * D1() * D1();

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
       	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_bK2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * K1K1(kvec1, kvec3) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * K1K1(kvec2, kvec3) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
       	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_bK2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * K1K1(kvec1, kvec2) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * K1K1(kvec1, kvec3) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * K1K1(kvec2, kvec3) * V1(kvec2, los) * V1(kvec3, los);

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
       	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}



double Bispectrum_Tree_BAO_b1f_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * D1() * D1();
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * D1() * D1();
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * D1() * D1();

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
       	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_b1f_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

	K12 += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * D1() * D1();
	K13 += 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * D1() * D1();
	K23 += 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * D1() * D1();

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
       	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_b1f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * D1V1(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * D1V1(kvec1, kvec3, los) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * D1V1(kvec2, kvec3, los) * V1(kvec2, los) * V1(kvec3, los);

	K12 += 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * ( D1() * V1(kvec2, los) + V1(kvec1, los) * D1() );
	K13 += 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * ( D1() * V1(kvec3, los) + V1(kvec1, los) * D1() );
	K23 += 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * ( D1() * V1(kvec3, los) + V1(kvec2, los) * D1() );

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
       	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}


double Bispectrum_Tree_BAO_ff_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * V1(kvec2, los) * V1(kvec3, los);

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
       	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

double Bispectrum_Tree_BAO_f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
     	double D1_damp = ExpDamping(kvec1, los, sigma2_perp, sigma2_para);
     	double D2_damp = ExpDamping(kvec2, los, sigma2_perp, sigma2_para);
     	double D3_damp = ExpDamping(kvec3, los, sigma2_perp, sigma2_para);

	double K12 = 2.0 * V2(kvec1, kvec2, los) * V1(kvec1, los) * V1(kvec2, los);
	double K13 = 2.0 * V2(kvec1, kvec3, los) * V1(kvec1, los) * V1(kvec3, los);
	double K23 = 2.0 * V2(kvec2, kvec3, los) * V1(kvec2, los) * V1(kvec3, los);

     	double BAO1 = f_pk(k1) - f_pk_no_wiggle(k1);
     	double BAO2 = f_pk(k2) - f_pk_no_wiggle(k2);
     	double BAO3 = f_pk(k3) - f_pk_no_wiggle(k3);
     
	double GG = D1_damp * D2_damp * D3_damp * K12 * BAO1 * BAO2
     		  + D1_damp * D2_damp * D3_damp * K13 * BAO1 * BAO3
     		  + D1_damp * D2_damp * D3_damp * K23 * BAO2 * BAO3;
     
     	double GM_MG = D1_damp * D1_damp * K12 * BAO1 * f_pk_no_wiggle(k2)
     		     + D2_damp * D2_damp * K12 * f_pk_no_wiggle(k1) * BAO2
     
     		     + D1_damp * D1_damp * K13 * BAO1 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K13 * f_pk_no_wiggle(k1) * BAO3
     
     		     + D2_damp * D2_damp * K23 * BAO2 * f_pk_no_wiggle(k3)
     		     + D3_damp * D3_damp * K23 * f_pk_no_wiggle(k2) * BAO3;
     
     	double NW = K12 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k2)
     	          + K13 * f_pk_no_wiggle(k1) * f_pk_no_wiggle(k3)
     	          + K23 * f_pk_no_wiggle(k2) * f_pk_no_wiggle(k3);
     
     	return (GG + GM_MG + NW) / alpha6;

}

/***************/
/* nonGaussian */
/***************/
// Z1 = (b1) * D1() + (f) * V1(kvec1, los);

double Bispectrum_NonGaussian_Local_b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = D1() * D1() * D1();
	double B = Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}

double Bispectrum_NonGaussian_Local_b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = D1() * D1() * ( V1(kvec1, los) + V1(kvec2, los) + V1(kvec3, los) );
	double B = Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}


double Bispectrum_NonGaussian_Local_b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = D1() * ( V1(kvec1, los) * V1(kvec2, los) + V1(kvec1, los) * V1(kvec3, los) + V1(kvec2, los) * V1(kvec3, los) );
	double B = Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}

double Bispectrum_NonGaussian_Local_f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = ( V1(kvec1, los) * V1(kvec2, los) * V1(kvec3, los) );
	double B = Primordial_Matter_Bispectrum_Local(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}


double Bispectrum_NonGaussian_Equilateral_b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = D1() * D1() * D1();
	double B = Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}

double Bispectrum_NonGaussian_Equilateral_b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = D1() * D1() * ( V1(kvec1, los) + V1(kvec2, los) + V1(kvec3, los) );
	double B = Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}


double Bispectrum_NonGaussian_Equilateral_b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = D1() * ( V1(kvec1, los) * V1(kvec2, los) + V1(kvec1, los) * V1(kvec3, los) + V1(kvec2, los) * V1(kvec3, los) );
	double B = Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}

double Bispectrum_NonGaussian_Equilateral_f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double factor = ( V1(kvec1, los) * V1(kvec2, los) * V1(kvec3, los) );
	double B = Primordial_Matter_Bispectrum_Equilateral(kvec1, kvec2, kvec3);
	return factor * B / alpha6;
}


/* NonGaussian from PB */


double Bispectrum_NonGaussian_From_PB_Local_b1_b1_b1_R(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                               double alpha_perp, double alpha_parallel, double R) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p  = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec1_M_pvec, kvec2_P_pvec, kvec3, R) * f_pk(p) * W1(p, R) * W1(p, R)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec2_M_pvec, kvec1_P_pvec, kvec3, R) * f_pk(p) * W1(p, R) * W1(p, R)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec1_M_pvec, kvec3_P_pvec, kvec2, R) * f_pk(p) * W1(p, R) * W1(p, R)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec3_M_pvec, kvec1_P_pvec, kvec2, R) * f_pk(p) * W1(p, R) * W1(p, R)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec2_M_pvec, kvec3_P_pvec, kvec1, R) * f_pk(p) * W1(p, R) * W1(p, R)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec3_M_pvec, kvec2_P_pvec, kvec1, R) * f_pk(p) * W1(p, R) * W1(p, R);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec2_M_pvec, M_kvec2, pvec, R) * f_pk(k3) * W1(k3, R) * W1(k3, R)
		       
		       + 2.0 * F2(kvec1, kvec3) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec1_M_pvec, M_kvec1, pvec, R) * f_pk(k3) * W1(k3, R) * W1(k3, R)

		       + 2.0 * F2(kvec1, kvec2) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec1_M_pvec, M_kvec1, pvec, R) * f_pk(k2) * W1(k2, R) * W1(k2, R)
		       
		       + 2.0 * F2(kvec3, kvec2) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec3_M_pvec, M_kvec3, pvec, R) * f_pk(k2) * W1(k2, R) * W1(k2, R)

		       + 2.0 * F2(kvec2, kvec1) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec2_M_pvec, M_kvec2, pvec, R) * f_pk(k1) * W1(k1, R) * W1(k1, R)
		       
		       + 2.0 * F2(kvec3, kvec1) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local_R(kvec3_M_pvec, M_kvec3, pvec, R) * f_pk(k1) * W1(k1, R) * W1(k1, R);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b1_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p  = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * F2(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * F2(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * F2(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * F2(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * F2(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * F2(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_b2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_b2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_bK2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_bK2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * K1K1(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * K1K1(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * K1K1(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * K1K1(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * K1K1(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * K1K1(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b1_b1f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_b1f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}




double Bispectrum_NonGaussian_From_PB_Local_b1_ff_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_ff_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * F2(kvec1_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * F2(kvec1_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec3_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * F2(kvec2_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * F2(kvec3_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * F2(kvec2, kvec3) * V2(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * F2(kvec1, kvec3) * V2(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * F2(kvec1, kvec2) * V2(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * F2(kvec3, kvec2) * V2(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * F2(kvec2, kvec1) * V2(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * F2(kvec3, kvec1) * V2(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * F2(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * F2(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * F2(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * F2(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * F2(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * F2(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * F2(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * F2(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * F2(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * F2(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * F2(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * F2(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b2_b2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b2_b2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b2_bK2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b2_bK2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * K1K1(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * K1K1(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * K1K1(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * K1K1(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * K1K1(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * K1K1(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b2_b1f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b2_b1f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b2_ff_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b2_ff_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b2_f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (D1D1(kvec1_M_pvec, pvec)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (D1D1(kvec2_M_pvec, pvec)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (D1D1(kvec3_M_pvec, pvec)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (D1D1(kvec2, kvec3)/2.0) * V2(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (D1D1(kvec1, kvec3)/2.0) * V2(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (D1D1(kvec1, kvec2)/2.0) * V2(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (D1D1(kvec3, kvec2)/2.0) * V2(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (D1D1(kvec2, kvec1)/2.0) * V2(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (D1D1(kvec3, kvec1)/2.0) * V2(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p  = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * F2(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * F2(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * F2(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * F2(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * F2(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * F2(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * F2(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * F2(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * F2(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_b2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_b2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_bK2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_bK2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * K1K1(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * K1K1(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * K1K1(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * K1K1(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * K1K1(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * K1K1(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_bK2_b1f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                               double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_b1f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}




double Bispectrum_NonGaussian_From_PB_Local_bK2_ff_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_ff_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_bK2_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * K1K1(kvec1_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * K1K1(kvec1_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec3_M_pvec, pvec) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * K1K1(kvec2_M_pvec, pvec) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * K1K1(kvec3_M_pvec, pvec) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * K1K1(kvec2, kvec3) * V2(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * K1K1(kvec1, kvec3) * V2(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * K1K1(kvec1, kvec2) * V2(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * K1K1(kvec3, kvec2) * V2(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * K1K1(kvec2, kvec1) * V2(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * K1K1(kvec3, kvec1) * V2(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p  = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * F2(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * F2(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * F2(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * F2(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * F2(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * F2(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_b2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_b2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_bK2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_bK2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * K1K1(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * K1K1(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * K1K1(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * K1K1(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * K1K1(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * K1K1(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_b1f_b1f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_b1f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}




double Bispectrum_NonGaussian_From_PB_Local_b1f_ff_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_ff_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_b1f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * D1V1(kvec1_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * D1V1(kvec1_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * D1V1(kvec2_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * D1V1(kvec3_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * D1V1(kvec2, kvec3, los) * V2(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * D1V1(kvec1, kvec3, los) * V2(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * D1V1(kvec1, kvec2, los) * V2(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * D1V1(kvec3, kvec2, los) * V2(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * D1V1(kvec2, kvec1, los) * V2(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * D1V1(kvec3, kvec1, los) * V2(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p  = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * F2(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * F2(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * F2(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * F2(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * F2(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * F2(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * F2(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * F2(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * F2(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * F2(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * F2(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * F2(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_b2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_b2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_bK2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                              double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_bK2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * K1K1(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * K1K1(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * K1K1(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * K1K1(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * K1K1(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * K1K1(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_ff_b1f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                              double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_b1f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_ff_ff_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_ff_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_ff_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * (V1V1(kvec2, kvec3, los)/2.0) * V2(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * (V1V1(kvec1, kvec3, los)/2.0) * V2(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * (V1V1(kvec1, kvec2, los)/2.0) * V2(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * (V1V1(kvec3, kvec2, los)/2.0) * V2(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * (V1V1(kvec2, kvec1, los)/2.0) * V2(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * (V1V1(kvec3, kvec1, los)/2.0) * V2(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

/*NAONAO*/
//double Z2_Bias(double * kvec1, double * kvec2, double * los, double f, double b1, double b2, double bK2) {
//	return   (f) *  V2(kvec1, kvec2, los)
//	      +  (f * f) * (V1V1(kvec1, kvec2, los) / 2.0)
//	      +  b1 * f * D1V1(kvec1, kvec2, los) 
//	      +  b1 * F2(kvec1, kvec2) 
//	      +  b2 * (D1D1(kvec1, kvec2) /2.0)
//	      +  bK2 * K1K1(kvec1, kvec2);
//}

/* HITOMI */

double Bispectrum_NonGaussian_From_PB_Local_f_b1_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p  = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * F2(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * F2(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * F2(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_b1_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * F2(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * F2(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * F2(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * F2(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * F2(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * F2(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * F2(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * F2(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * F2(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_b2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_b2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (D1D1(kvec1_P_pvec, M_pvec)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (D1D1(kvec3_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (D1D1(kvec2_P_pvec, M_pvec)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * (D1D1(kvec1_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * (D1D1(kvec2_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * (D1D1(kvec3_M_pvec, pvec)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_bK2_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * K1K1(kvec1_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * K1K1(kvec2_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * K1K1(kvec3_M_pvec, pvec) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_bK2_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * K1K1(kvec1_P_pvec, M_pvec) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * K1K1(kvec3_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * K1K1(kvec2_P_pvec, M_pvec) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * K1K1(kvec2_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * K1K1(kvec1_M_pvec, pvec) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * K1K1(kvec1_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * K1K1(kvec3_M_pvec, pvec) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * K1K1(kvec2_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * K1K1(kvec3_M_pvec, pvec) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}


double Bispectrum_NonGaussian_From_PB_Local_f_b1f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * D1V1(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * D1V1(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * D1V1(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_b1f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * D1V1(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * D1V1(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * D1V1(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * D1V1(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * D1V1(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * D1V1(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}




double Bispectrum_NonGaussian_From_PB_Local_f_ff_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_ff_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (V1V1(kvec1_P_pvec, M_pvec, los)/2.0) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * (V1V1(kvec3_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * (V1V1(kvec2_P_pvec, M_pvec, los)/2.0) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * (V1V1(kvec1_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * (V1V1(kvec2_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * (V1V1(kvec3_M_pvec, pvec, los)/2.0) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_f_b1(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                             double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * V2(kvec1_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * V2(kvec2_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * V2(kvec3_M_pvec, pvec, los) * D1()
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	return B / alpha6;

}

double Bispectrum_NonGaussian_From_PB_Local_f_f_f(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * pvec, double * los, 
	                                            double alpha_perp, double alpha_parallel) {

	double alpha3 = alpha_perp * alpha_perp * alpha_parallel;
	double alpha6 = alpha3 * alpha3;

	double kvec1[3] = { 0.0, 0.0, 0.0 };
	double kvec2[3] = { 0.0, 0.0, 0.0 };
	double kvec3[3] = { 0.0, 0.0, 0.0 };

	calcTrueWavevector(kvec1_in, los, alpha_perp, alpha_parallel, kvec1);
	calcTrueWavevector(kvec2_in, los, alpha_perp, alpha_parallel, kvec2);
	calcTrueWavevector(kvec3_in, los, alpha_perp, alpha_parallel, kvec3);

	double kvec1_M_pvec[3] = { kvec1[0] - pvec[0], kvec1[1] - pvec[1], kvec1[2] - pvec[2] };
	double kvec2_M_pvec[3] = { kvec2[0] - pvec[0], kvec2[1] - pvec[1], kvec2[2] - pvec[2] };
	double kvec3_M_pvec[3] = { kvec3[0] - pvec[0], kvec3[1] - pvec[1], kvec3[2] - pvec[2] };

	double kvec1_P_pvec[3] = { kvec1[0] + pvec[0], kvec1[1] + pvec[1], kvec1[2] + pvec[2] };
	double kvec2_P_pvec[3] = { kvec2[0] + pvec[0], kvec2[1] + pvec[1], kvec2[2] + pvec[2] };
	double kvec3_P_pvec[3] = { kvec3[0] + pvec[0], kvec3[1] + pvec[1], kvec3[2] + pvec[2] };
	
	double M_kvec1[3] = { - kvec1[0], - kvec1[1], - kvec1[2] };
	double M_kvec2[3] = { - kvec2[0], - kvec2[1], - kvec2[2] };
	double M_kvec3[3] = { - kvec3[0], - kvec3[1], - kvec3[2] };
	double M_pvec[3]  = { - pvec[0], - pvec[1], - pvec[2] };

	double k1 = NORM(kvec1);
	double k2 = NORM(kvec2);
	double k3 = NORM(kvec3);
	double p = NORM(pvec);
	
	double B_122_A = 2.0 * V2(kvec1_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec2_P_pvec, kvec3) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec3, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec1_P_pvec, kvec3) * f_pk(p)
	
		       + 2.0 * V2(kvec1_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, kvec3_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec3_M_pvec, pvec, los) * V2(kvec1_P_pvec, M_pvec, los) * V1(kvec2, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec1_P_pvec, kvec2) * f_pk(p)

		       + 2.0 * V2(kvec2_M_pvec, pvec, los) * V2(kvec3_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, kvec3_P_pvec, kvec1) * f_pk(p)

	    	       + 2.0 * V2(kvec3_M_pvec, pvec, los) * V2(kvec2_P_pvec, M_pvec, los) * V1(kvec1, los) 
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, kvec2_P_pvec, kvec1) * f_pk(p);

	double B_122_B = 2.0 * V2(kvec2, kvec3, los) * V2(kvec2_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k3)
		       
		       + 2.0 * V2(kvec1, kvec3, los) * V2(kvec1_M_pvec, pvec, los) * V1(kvec3, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k3)

		       + 2.0 * V2(kvec1, kvec2, los) * V2(kvec1_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec1_M_pvec, M_kvec1, pvec) * f_pk(k2)
		       
		       + 2.0 * V2(kvec3, kvec2, los) * V2(kvec3_M_pvec, pvec, los) * V1(kvec2, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k2)

		       + 2.0 * V2(kvec2, kvec1, los) * V2(kvec2_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec2_M_pvec, M_kvec2, pvec) * f_pk(k1)
		       
		       + 2.0 * V2(kvec3, kvec1, los) * V2(kvec3_M_pvec, pvec, los) * V1(kvec1, los)
	               * Primordial_Matter_Bispectrum_Local(kvec3_M_pvec, M_kvec3, pvec) * f_pk(k1);

	double B = B_122_A + B_122_B;
	B = B;
	return B / alpha6;

}



#endif

