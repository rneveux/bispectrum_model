import os
import argparse
import numpy as np
from classy import Class
from powerbispectrum import ComputePowerBiSpectrum

def create_parser():
    """
    Creates a command line parser
    """
    parser = argparse.ArgumentParser(description='Computation of bispectrum kernels')
    parser.add_argument('-cosmo', type=int, help='number of test cosmology', required=True)
    parser.add_argument('-redshift', type=float, help='effective redhshift', required=True)
    parser.add_argument('-ells', type=str, help='multipoles of bispectrum - ell1 ell2 ELL', required=True)
    parser.add_argument('-diag', type=str, default='Yes', help='file with full kernel or diag only', required=False)
    parser.add_argument('-directory', type=str, default='/home/rneveux/test/', help='name of directory', required=False)
    parser.add_argument('-AP', type=str, help='include AP effect', default='no', required=False)
    parser.add_argument('-k_list', type=str, help='file name of list of wavenumber to compute bk', required=True)
    return parser

def main(args):
    """
    Main function that runs the computation
    """

    k = np.loadtxt(args.k_list)
    z = args.redshift
    ells = args.ells

    # create a dictionary mapping for all_kernels
    part_bispectrum = ['tree','SN','PNG']
    all_kernels_mapping = {
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
    }
    
    os.makedirs(os.path.join(args.directory, ells), exist_ok=True)
    name_file = os.path.join(args.directory, ells, f'bk_kernels_{args.cosmo}.txt')

    if not os.path.exists(name_file):

        cosmo = np.load('/home/rneveux/bispectrum/theory/cosmologies/lnAs/eft_cosmologies_noDQ1.npy')[args.cosmo]

        h = cosmo[2]
        omega_cdm = cosmo[0]
        omega_b = cosmo[1]
        lnAs = cosmo[3]
        ns = cosmo[4]
        
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

        c = ComputePowerBiSpectrum(params_cosmo, z)
        c.initial_power_spectrum(cosmo, aperp=aperp, apar=apar)
    
        ell1, ell2, ELL = map(int, ells)
        bk_kernels = []
        
        for part in part_bispectrum:
            for ker in all_kernels_mapping[part]:
                c.kernel_computation_Bk(ker, k, ell1, ell2, ELL, integrand=part)
                bk_kernels.append(c.BK['K'])
        bk_kernels = np.array(bk_kernels)
        
        np.savetxt(name_file, bk_kernels)

if __name__ == "__main__":
    parser = create_parser()
    cmdline = parser.parse_args()
    main(cmdline)
