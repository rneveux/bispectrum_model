import os
import argparse
from classy import Class
import numpy as np
from powerbispectrum import ComputePowerBiSpectrum

def create_parser():
    """
    Creates a command line parser
    """
    parser = argparse.ArgumentParser(description='Computation of power and bispectrum kernels')
    parser.add_argument('-redshift', type=float, help='effective redhshift', required=True)
    parser.add_argument('-estimator', type=str, help='pk or bk', required=True)
    parser.add_argument('-spectrum_part', type=str, default='tree', help='tree, SN, PNG for bk or tree, 1loop, counterterm for pk', required=False)
    parser.add_argument('-kernel', type=int, help='number of kernel to emulate in list', required=True)
    parser.add_argument('-ells', type=str, help='multipoles of power or bispectrum - (ell1 ell2) ELL', required=True)
    parser.add_argument('-diag', type=str, default='Yes', help='file with full kernel or diag only for bk', required=False)
    parser.add_argument('-directory', type=str, default='/home/rneveux/kernels/bk/', help='name of directory', required=False)
    parser.add_argument('-AP', type=str, help='include AP effect', default='no', required=False)
    parser.add_argument('-f_fit', type=str, help='f in kernels', default='no', required=False)
    return parser

def main(args):
    """
    Main function that runs the computation
    """

    k = np.array([0.01159058, 0.02081665, 0.03048062, 0.04036437, 0.05036595,
       0.06034213, 0.0702755 , 0.08021871, 0.09019841, 0.10017723,
       0.11014627, 0.12013268, 0.13012506, 0.14012919, 0.15012158,
       0.16013557, 0.1701332 , 0.18010434, 0.19008708, 0.20007673])
    z = args.redshift

    # create a dictionary mapping for all_kernels
    all_kernels_mapping = {
        'bk': {
            'tree': ['b1_b1_b1', 'b1_b1_b2', 'b1_b1_bG2', 'b1_b1_f', 'b1_b1_b1_f', 'b1_b1_f_f', 
                     'b1_b2_f', 'b1_bG2_f', 'b1_f_f', 
                     'b1_f_f_f', 'b2_f_f', 'bG2_f_f', 'f_f_f', 'f_f_f_f', 
                     'c1_b1_b1', 'c1_b1_b2', 'c1_b1_bG2', 'c1_b1_f', 'c1_b1_b1_f', 'c1_b1_f_f', 'c1_b2_f', 
                     'c1_bG2_f', 'c1_f_f', 'c1_f_f_f', 'c1_c1_b1', 'c1_c1_b2', 'c1_c1_bG2', 'c1_c1_f', 
                     'c1_c1_b1_f', 'c1_c1_f_f', 'c2_b1_b1', 'c2_b1_b2', 'c2_b1_bG2', 'c2_b1_f', 'c2_b1_b1_f', 
                     'c2_b1_f_f', 'c2_b2_f', 'c2_bG2_f', 'c2_f_f', 'c2_f_f_f', 'c2_c1_b1', 'c2_c1_b2', 
                     'c2_c1_bG2', 'c2_c1_f', 'c2_c1_b1_f', 'c2_c1_f_f', 'c2_c2_b1', 'c2_c2_b2', 'c2_c2_bG2', 
                     'c2_c2_f', 'c2_c2_b1_f', 'c2_c2_f_f'],
            'SN': ['Bshot_b1_b1', 'Bshot_b1_f', 'Bshot_b1_c1', 'Bshot_b1_c2', 
                   'Pshot_f_b1', 'Pshot_f_f', 'Pshot_f_c1', 'Pshot_f_c2'],
            'PNG': ['fnlloc_b1_b1_b1', 'fnlloc_b1_b1_f', 'fnlloc_b1_f_f', 'fnlloc_f_f_f', 
                    'fnlequi_b1_b1_b1', 'fnlequi_b1_b1_f', 'fnlequi_b1_f_f', 'fnlequi_f_f_f', 
                    'fnlortho_b1_b1_b1', 'fnlortho_b1_b1_f', 'fnlortho_b1_f_f', 'fnlortho_f_f_f', 
                    'fnlortho_LSS_b1_b1_b1', 'fnlortho_LSS_b1_b1_f', 'fnlortho_LSS_b1_f_f', 'fnlortho_LSS_f_f_f']
        },
        'pk': {
            'tree': ["b1_b1", "b1_f", "f_f"],
            '1loop': ["b2_b2", "b2_bG2", "b1_b2", "b2_f", "b1_b2_f", "b2_f_f", "bG2_bG2", "b1_bG2", "bG2_f", 
                      "b1_bG2_f", "bG2_f_f", "b1_b1", "b1_f", "b1_b1_f", "b1_f_f", "f_f", "f_f_f", "b1_b1_f_f", "b1_f_f_f", "f_f_f_f", 
                      "b1_b3", "b1_bG3", "b1_bG2d", "b1_bGamma3", "b1_b1", "b1_f", "b1_b1_f_f", "b1_f_f_f", "b1_b1_f", 
                      "b1_f_f", "b1_b2", "b1_bG2", "b1_b2_f", "b1_bG2_f", "b3_f", "bG3_f", "bG2d_f", "bGamma3_f", "f_f", 
                      "f_f_f_f", "f_f_f", "b2_f", "bG2_f", "b2_f_f", "bG2_f_f"],
            'counterterm': ["c0_b1", "c0_f", "c1_b1", "c1_f", "c2_b1", "c2_f"]
        }
    }

    try:
        all_kernels = all_kernels_mapping[args.estimator][args.spectrum_part]
    except KeyError:
        print(f"Invalid estimator: {args.estimator} or spectrum_part: {args.spectrum_part}")
        return

    h = .6777
    omega_b = 0.048*h**2
    omega_cdm = 0.307115*h**2-omega_b
    lnAs = 3.0364
    ns = .96

    '''params_cosmo_s8 = {
        'output': 'mPk',
        'z_max_pk': 3.,
        'P_k_max_h/Mpc': 50.,
        'omega_cdm': omega_cdm,
        'omega_b': omega_b,
        'h': h,
        'sigma8': sigma8_input,
        'n_s': ns
    }

    c_s8 = Class()
    c_s8.set(params_cosmo_s8)
    c_s8.set({'output':'mPk',
           })

    c_s8.compute()

    As_rescaled = c_s8.A_s()/c_s8.sigma8()**2*sigma8_input**2'''

    params_cosmo = {
            'z_pk':z,
            'output': 'mPk',
            'z_max_pk': 3.,
            'P_k_max_h/Mpc': 50.,
            'omega_cdm': omega_cdm,
            'omega_b': omega_b,
            'h': h,
            'ln10^{10}A_s': lnAs,
            'n_s': ns,
            'N_ur':2.0328, 
            'N_ncdm':1, 
            'omega_ncdm':0.0006442
        }

    cosmo = Class()
    cosmo.set(params_cosmo)
    cosmo.compute()
    
    if args.AP!='no':

        #Abacus cosmology
        h_fid = .6736
        omega_cdm_fid = .12
        omega_b_fid = .02237
        As_fid = 2.0830e-9
        ns_fid = .9649

        D = cosmo.angular_distance(z)
        H = cosmo.Hubble(z)

        params_cosmo_fid = {
        'omega_cdm': omega_cdm_fid,
        'omega_b': omega_b_fid,
        'h': h_fid,
        'A_s': As_fid,
        'n_s': ns_fid,
        'N_ur':2.0328, 
        'N_ncdm':1, 
        'omega_ncdm':0.0006442
        }

        cosmo_fid = Class()
        cosmo_fid.set(params_cosmo_fid)
        cosmo_fid.compute()

        Dfid = cosmo_fid.angular_distance(z)
        Hfid = cosmo_fid.Hubble(z)

        aperp = D/Dfid
        apar = Hfid/H

    else:

        aperp = 0
        apar = 0
        
    f_fit = False
    if args.f_fit!='no':
        f_fit = True

    c = ComputePowerBiSpectrum(params_cosmo, z)
    c.initial_power_spectrum(cosmo, aperp=aperp, apar=apar)

    name_file = os.path.join(args.directory, args.ells, args.spectrum_part, all_kernels[args.kernel])

    if args.estimator == 'bk':
        ell1, ell2, ELL = map(int, args.ells)
        c.kernel_computation_Bk(all_kernels[args.kernel], k * h, ell1, ell2, ELL, integrand=args.spectrum_part, to_save=name_file, f_fit=f_fit)
    elif args.estimator == 'pk':
        ELL = int(args.ells)
        c.kernel_computation_Pk(all_kernels[args.kernel], k * h, ELL, integrand=args.spectrum_part, to_save=name_file)

if __name__ == "__main__":
    parser = create_parser()
    cmdline = parser.parse_args()
    main(cmdline)
