#ifndef __kernel2__
#define __kernel2__

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

/**********************************************************************************/
/* second order */

double Growth(double * kvec1, double * kvec2) {
	return 1.0;
}

double Shift(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
    double result = 0.0;
	if( (k1 > pk_kmin) && (k2 > pk_kmin)) {
	    result = 0.5 * mu * (k1/k2 + k2/k1);
    }
    return result;
}

double Tidal(double * kvec1, double * kvec2) {
	double mu = MU(kvec1, kvec2);
	return (mu * mu - 1.0/3.0);
}

double DV(double * kvec1, double * kvec2, double * los) {
    double kvec12[3] = PLUS(kvec1, kvec2);
    double kn = DOT(kvec12, los);
    double k1n = DOT(kvec1, los);
    double k2n = DOT(kvec2, los);

    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
    double result = 0.0;
	if( (k1 > pk_kmin) && (k2 > pk_kmin)) {
        result = 0.5 * kn * ( k1n / pow(k1,2) + k2n / pow(k2,2) );
    }   
    return result;
}

double V11(double * kvec1, double * kvec2, double * los) {

    double kvec12[3] = PLUS(kvec1, kvec2);
    double kn = DOT(kvec12, los);
    double k1n = DOT(kvec1, los);
    double k2n = DOT(kvec2, los);

    double k1 = NORM(kvec1);
    double k2 = NORM(kvec2);
    double result = 0.0;
	if( (k1 > pk_kmin) && (k2 > pk_kmin)) {
        result = 0.5 * kn * kn * ( k1n / pow(k1,2) ) * ( k2n / pow(k2,2) );
    }   
    return result;
}


/**********************************************************************************/

double LOS_MU12(double * kvec1, double * kvec2, double * los) {
    double kvec12[3] = PLUS(kvec1, kvec2);
    double mu = MU(kvec12, los);
    return (mu * mu);
}

double LOS_MU1_PLUS_MU2(double * kvec1, double * kvec2, double * los) {
    double mu1 = MU(kvec1, los);
    double mu2 = MU(kvec2, los);
    return (mu1 * mu1) + (mu2 * mu2);
}

double LOS_MU1MU2(double * kvec1, double * kvec2, double * los) {
    double mu1 = MU(kvec1, los);
    double mu2 = MU(kvec2, los);
    return (mu1 * mu1) * (mu2 * mu2);
}

/**********************************************************************************/

/**********************************************************************************/

/********/
double FG_b3_f0(double * kvec1, double * kvec2, double * los) {
    return Growth(kvec1, kvec2);
}
double FS_b3_f0(double * kvec1, double * kvec2, double * los) {
    return Shift(kvec1, kvec2);
}
double FT_b3_f0(double * kvec1, double * kvec2, double * los) {
    return Tidal(kvec1, kvec2);
}
/********/
double FG_b2_f1(double * kvec1, double * kvec2, double * los) {
    return Growth(kvec1, kvec2) * LOS_MU1_PLUS_MU2(kvec1, kvec2, los);
}
double FS_b2_f1(double * kvec1, double * kvec2, double * los) {
    return Shift(kvec1, kvec2)  * LOS_MU1_PLUS_MU2(kvec1, kvec2, los);
}
double FT_b2_f1(double * kvec1, double * kvec2, double * los) {
    return Tidal(kvec1, kvec2)  * LOS_MU1_PLUS_MU2(kvec1, kvec2, los);
}
/********/
double FG_b1_f2(double * kvec1, double * kvec2, double * los) {
    return Growth(kvec1, kvec2) * LOS_MU1MU2(kvec1, kvec2, los);
}
double FS_b1_f2(double * kvec1, double * kvec2, double * los) {
    return Shift(kvec1, kvec2)  * LOS_MU1MU2(kvec1, kvec2, los);
}
double FT_b1_f2(double * kvec1, double * kvec2, double * los) {
    return Tidal(kvec1, kvec2)  * LOS_MU1MU2(kvec1, kvec2, los);
}
/********/
/********/
double GG_b2_f1(double * kvec1, double * kvec2, double * los) {
    return Growth(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los);
}
double GS_b2_f1(double * kvec1, double * kvec2, double * los) {
    return Shift(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los);
}
double GT_b2_f1(double * kvec1, double * kvec2, double * los) {
    return Tidal(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los);
}
/********/
double GG_b1_f2(double * kvec1, double * kvec2, double * los) {
    return Growth(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los) * LOS_MU1_PLUS_MU2(kvec1, kvec2, los);
}
double GS_b1_f2(double * kvec1, double * kvec2, double * los) {
    return Shift(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los)  * LOS_MU1_PLUS_MU2(kvec1, kvec2, los);
}
double GT_b1_f2(double * kvec1, double * kvec2, double * los) {
    return Tidal(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los)  * LOS_MU1_PLUS_MU2(kvec1, kvec2, los);
}
/********/
double GG_b0_f3(double * kvec1, double * kvec2, double * los) {
    return Growth(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los) * LOS_MU1MU2(kvec1, kvec2, los);
}
double GS_b0_f3(double * kvec1, double * kvec2, double * los) {
    return Shift(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los)  * LOS_MU1MU2(kvec1, kvec2, los);
}
double GT_b0_f3(double * kvec1, double * kvec2, double * los) {
    return Tidal(kvec1, kvec2) * LOS_MU12(kvec1, kvec2, los)  * LOS_MU1MU2(kvec1, kvec2, los);
}
/********/
double b3_f1(double * kvec1, double * kvec2, double * los) {
    return DV(kvec1, kvec2, los);
}

double b2_f2(double * kvec1, double * kvec2, double * los) {
    return DV(kvec1, kvec2, los) * LOS_MU1_PLUS_MU2(kvec1, kvec2, los) + V11(kvec1, kvec2, los);
}

double b1_f3(double * kvec1, double * kvec2, double * los) {
    return DV(kvec1, kvec2, los) * LOS_MU1MU2(kvec1, kvec2, los) + V11(kvec1, kvec2, los) * LOS_MU1_PLUS_MU2(kvec1, kvec2, los);
}

double b0_f4(double * kvec1, double * kvec2, double * los) {
    return V11(kvec1, kvec2, los) * LOS_MU1MU2(kvec1, kvec2, los);
}
/********/
double KernelFunctionTemplate(double * kvec1, double * kvec2, double * los, char * parameters_in) {
    
    std::string parameters = parameters_in;

    if(0) {

    } else if (parameters == "FG_b3_f0") {
        return FG_b3_f0(kvec1, kvec2, los);
    } else if (parameters == "FS_b3_f0") {
        return FS_b3_f0(kvec1, kvec2, los);
    } else if (parameters == "FT_b3_f0") {
        return FT_b3_f0(kvec1, kvec2, los);
 
    } else if (parameters == "FG_b2_f1") {
        return FG_b2_f1(kvec1, kvec2, los);
    } else if (parameters == "FS_b2_f1") {
        return FS_b2_f1(kvec1, kvec2, los);
    } else if (parameters == "FT_b2_f1") {
        return FT_b2_f1(kvec1, kvec2, los);
  
    } else if (parameters == "FG_b1_f2") {
        return FG_b1_f2(kvec1, kvec2, los);
    } else if (parameters == "FS_b1_f2") {
        return FS_b1_f2(kvec1, kvec2, los);
    } else if (parameters == "FT_b1_f2") {
        return FT_b1_f2(kvec1, kvec2, los);
   
    } else if (parameters == "GG_b2_f1") {
        return GG_b2_f1(kvec1, kvec2, los);
    } else if (parameters == "GS_b2_f1") {
        return GS_b2_f1(kvec1, kvec2, los);
    } else if (parameters == "GT_b2_f1") {
        return GT_b2_f1(kvec1, kvec2, los);
    
    } else if (parameters == "GG_b1_f2") {
        return GG_b1_f2(kvec1, kvec2, los);
    } else if (parameters == "GS_b1_f2") {
        return GS_b1_f2(kvec1, kvec2, los);
    } else if (parameters == "GT_b1_f2") {
        return GT_b1_f2(kvec1, kvec2, los);

    } else if (parameters == "GG_b0_f3") {
        return GG_b0_f3(kvec1, kvec2, los);
    } else if (parameters == "GS_b0_f3") {
        return GS_b0_f3(kvec1, kvec2, los);
    } else if (parameters == "GT_b0_f3") {
        return GT_b0_f3(kvec1, kvec2, los);

    } else if (parameters == "b3_f1") {
        return b3_f1(kvec1, kvec2, los);
    } else if (parameters == "b2_f2") {
        return b2_f2(kvec1, kvec2, los);
    } else if (parameters == "b1_f3") {
        return b1_f3(kvec1, kvec2, los);
    } else if (parameters == "b0_f4") {
        return b0_f4(kvec1, kvec2, los);
    } else {
        return 0.0;
    }
}

/********/

double Bispectrum_Tree_BAO_Template(
        double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, 
        double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para,
        char * parameters) {

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
    
    double K12 = 2.0 * KernelFunctionTemplate(kvec1, kvec2, los, parameters);
    double K13 = 2.0 * KernelFunctionTemplate(kvec1, kvec3, los, parameters);
    double K23 = 2.0 * KernelFunctionTemplate(kvec2, kvec3, los, parameters);
    
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

double Bispectrum_Tree_BAO_FS_b3_f0(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

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
    
    double K12 = 2.0 * FS_b3_f0(kvec1, kvec2, los);
    double K13 = 2.0 * FS_b3_f0(kvec1, kvec3, los);
    double K23 = 2.0 * FS_b3_f0(kvec2, kvec3, los);
    
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

double Bispectrum_Tree_BAO_FT_b3_f0(double * kvec1_in, double * kvec2_in, double * kvec3_in, double * los, double alpha_perp, double alpha_parallel, double sigma2_perp, double sigma2_para) {

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
    
    double K12 = 2.0 * FT_b3_f0(kvec1, kvec2, los);
    double K13 = 2.0 * FT_b3_f0(kvec1, kvec3, los);
    double K23 = 2.0 * FT_b3_f0(kvec2, kvec3, los);
    
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







#endif

