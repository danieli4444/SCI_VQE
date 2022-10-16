import psi4
import numpy as np
import os 
import json
import re
from helper_CI import Determinant, HamiltonianGenerator
from householder import get_tridiagonal_mat

PSI4_LOG_FILE = 'psi4_output_log.dat' # the general output log of all psi4 computations
# CI_CONFIGS_FILE = "CI_configs.txt" # holds the fci/cisd results
CI_RESULT_FILE = "CI_result.json"
CI_DIR_NAME = "CI_matrices"   


class CIResult:
    def __init__(self,molecule,basis,nuclear_repulsion_energy,ci_type,
        ci_ground_state_energy,ci_total_num_configs) -> None:

        self.molecule = molecule
        self.basis_set = basis
        self.nuclear_repulsion_energy = nuclear_repulsion_energy
        self.ci_type = ci_type
        self.ci_ground_state_energy = ci_ground_state_energy
        self.ci_total_num_configs = ci_total_num_configs
    
    def toJson(self,filename):
        with open(filename,'w+') as f:
             json.dump(self.__dict__,f,indent=4)
            


class SimpleDet:
    def __init__(self,rank,coefficent,ci_sym_repr,occu_repr) -> None:
        self.rank = rank
        self.coeffcient = coefficent
        self.c1_sym_repr = ci_sym_repr
        self.occu_repr = occu_repr
 
class CIMatrixGenerator:

    def __init__(self, molecule,mol_name,basis_set):
        """_summary_

        Args:
            molecule (psi4 molecule definition string ): contains also the basis set.
            num_congifs (_type_): number of desired configurations to include in the ci matrix.
        """
        self.molecule = molecule
        self.basis = basis_set
        self.mol_name = mol_name
        self.mol_dir_path ,self.mol_dir_name,self.mol_ci_dir = self.create_molecule_dir()
        

        self.num_ci_configurations = 0
        self.ci_matrix = 0 
        self.ci_ground_energy = 0
        self.ci_ground_state = 0
        self.ci_type = '' # whether to use fci, cisd or else in the psi4 CI energy calculation. 
        # will serve later for extracting the most significant configs.

        self.ci_most_important_configs_list = 0
        self.ci_total_num_configs = 0
        self.ci_filepath = ''

    def create_molecule_dir(self):
        mol_dir_name = self.mol_name
        path = os.getcwd() + '/' + mol_dir_name
        if not os.path.exists(path):
            os.mkdir(path) 

        ci_path = path + '/' + CI_DIR_NAME
        if not os.path.exists(ci_path):
            os.mkdir(ci_path) 

        return path,mol_dir_name,ci_path

    @staticmethod
    def parse_psi4_outputfile(filename):
        """parses psi4 ouputfile and returns:
            * a sorted list of most important determinants in the orbitals occupation representation
            * total number of determinants in CI calculation

        Args:
            filename (_type_, optional): _description_. Defaults to PSI4_LOG_FILE.
        """

        print("parsing psi4 outputfile...")

        num_total_dets = 0 
        num_important_dets = 0
        most_important_det_list = []
        with open(file=filename) as raw_output_file:

            output_string = raw_output_file.read()
            lines = output_string.splitlines()

            for id,line in enumerate(lines):
                # get total number of determinants
                if "requires" in line:
                    x = line.split(' ')
                    i = x.index("requires")
                    num_total_dets = int(x[i+1])
                
                # get all determinants
                if "important determinants" in line:
                    det_idx = id
                    x = line.split(' ')
                    i = x.index('most')
                    num_important_dets = int(x[i-1])
                    
            if num_total_dets == 0:
                raise("File Parsing Error")

            # parse all important dets and save to a list:
            relevant_lines = lines[ det_idx + 2 : det_idx + 2 + num_important_dets]
            for line in relevant_lines:
                words = line.split(' ')
                for w in words:
                    if w.isdigit():
                        rank = int(w)
                        break
                nums = [float(s) for s in words if re.match(r'^-?\d+(?:\.\d+)$', s) is not None]
                coeff = nums[0]

                z = line.index(")")
                config_part = line[z+3:]
                orbitals = config_part.split(' ')
                alphas = []
                betas = []
                for orbital in orbitals[:-1]:
                    orbit_num_str = re.findall(r'\d+', orbital)[0]
                    orbit_num = int(orbit_num_str) - 1
                    spin_num = orbital.replace(orbit_num_str,'')[1]
                    if spin_num == 'A':
                        alphas.append(orbit_num)
                    if spin_num == 'B':
                        betas.append(orbit_num)
                    if spin_num == 'X':
                        alphas.append(orbit_num)
                        betas.append(orbit_num)
                occupation_representation = (alphas,betas)
                det = SimpleDet(rank,coeff,ci_sym_repr=orbitals,occu_repr=occupation_representation)
                most_important_det_list.append(det)

        return most_important_det_list,num_total_dets


    def run_psi4_ci(self, configs):
        # run CI calculation 
        
        # psi4.set_memory('500 MB')

        # h2o = psi4.geometry("""
        # O
        # H 1 0.955
        # H 1 0.955 2 105
        # symmetry c1

        # """)
        ci_type = self.ci_type
        mol = psi4.geometry(self.molecule)
        self.shift = mol.nuclear_repulsion_energy()
        output_file_name = self.mol_dir_name + '/' + PSI4_LOG_FILE
        psi4.core.set_output_file(output_file_name, False)
        psi4.core.prepare_options_for_module('DETCI')

        psi4.set_options({'basis': self.basis,
                        'scf_type': 'pk',
                        'R_CONVERGENCE': 1e-4,
                        'E_CONVERGENCE': 1e-4,
                        'NUM_DETS_PRINT': configs})

        print("starting fci computation!!!!")
        (ci_ground_energy,scf_wfn) = psi4.energy(ci_type,return_wfn=True)
        psi4.core.clean()

        best_configs_list,ci_total_dets = CIMatrixGenerator.parse_psi4_outputfile(output_file_name)
        
        return (best_configs_list, ci_ground_energy, ci_total_dets)


    def save_ci_results(self):
        filename = self.ci_filepath
        output = CIResult(self.molecule,self.basis,self.shift,
            self.ci_type,self.ci_ground_energy,self.ci_total_num_configs)
        output.toJson(filename)


    def generate_ci_results(self,ci_type, max_configs=128):
        """ runs psi4 ci for the molecule and saves all relevant info

        Args:
            max_configs (_int_): the number of configs needs to be set less than total number possible configurations.
        """
        self.ci_type = ci_type
        # get a list of most imporant configurations ,ci energy and total number of configs in CI
        (best_configs, ci_energy, ci_total_dets) = self.run_psi4_ci(max_configs)


        self.ci_ground_energy = ci_energy
        self.ci_most_important_configs_list = best_configs
        self.ci_total_num_configs = ci_total_dets

        self.ci_filepath = self.mol_dir_name + '/' + CI_RESULT_FILE 
        self.save_ci_results()
        print("Computed and save CI congifs in {0}".format(self.ci_filepath))


    def get_reduced_cimatrix(self,num_desired_configs):
        """ returns a ci matrix

        Args:
            num_desired_configs (_type_): top congiration from ci computation from which to construct a ci matrix

        """
        
        mol = psi4.geometry(self.molecule)

        psi4.set_options({'basis': self.basis,
                        'scf_type': 'pk'})
        # First compute SCF energy using Psi4
        scf_e, wfn = psi4.energy('SCF', return_wfn=True)
        psi4.core.clean()

        mints = psi4.core.MintsHelper(wfn.basisset())

        # Grab data from wavfunction class
        C = wfn.Ca()
        MO = np.asarray(mints.mo_spin_eri(C, C))

        H = np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())

        # Update H, transform to MO basis and tile for alpha/beta spin
        H = np.einsum('uj,vi,uv', C, C, H)
        H = np.repeat(H, 2, axis=0)
        H = np.repeat(H, 2, axis=1)

        # Make H block diagonal
        spin_ind = np.arange(H.shape[0], dtype=int) % 2
        H *= (spin_ind.reshape(-1, 1) == spin_ind)

        Hamiltonian_generator = HamiltonianGenerator(H, MO)

        new_det_list = []
        for i in range(num_desired_configs):
            det = self.ci_most_important_configs_list[i]
            (a,b) = det.occu_repr
            new_det_list.append(Determinant(alphaObtList=a, betaObtList=b))

        my_Hamiltonian_matrix = Hamiltonian_generator.generateMatrix(new_det_list)
        # tridiagonal_H = get_tridiagonal_mat(my_Hamiltonian_matrix)

        hamiltonian_matrix_file_name = self.mol_ci_dir + '/' + self.mol_name + '_cimat_' + '_' + str(num_desired_configs) + '.out'
        # save ci matrix
        np.savetxt(hamiltonian_matrix_file_name,my_Hamiltonian_matrix)
        print("computed and saved CI matrix to file {0}".format(hamiltonian_matrix_file_name))
                

def generate_h20_ci_matrices():
    print("generating ci matrices!")
    h20_mol = """
        O
        H 1 0.955
        H 1 0.955 2 105
        symmetry c1
        """
    basis = 'sto-3g'
    # number of wanted configurations (currently tested only powers of 2 ):
    num_configs = 130
    mol_name = 'H2O_tridiagonal'
    ci1 = CIMatrixGenerator(h20_mol,mol_name,basis)
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type) # use full ci for finding the most significant configs
    
    for i in range(1,7):
        num_configs = 2**i
        ci1.get_reduced_cimatrix(num_configs)    


def generate_h2_ci_matrix(dist):
    print("generating ci matrices!")
    h2_mol = """
            H
            H 1 {0}
            symmetry c1
            """.format(dist)

    basis = 'sto-3g'
    # number of wanted configurations (currently tested only powers of 2 ):
    num_configs = 2
    mol_name='H2_' + '{0}'.format(round(dist,1))
    ci1 = CIMatrixGenerator(h2_mol,mol_name,basis) 
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type,max_configs=2) # use full ci for finding the most significant configs
    
    ci1.get_reduced_cimatrix(num_configs)    


def generate_LiH_ci_matrix(dist):
    print("generating ci matrices!")
    h2_mol = """
            Li
            H 1 {0}
            symmetry c1
            """.format(dist)

    basis = 'sto-3g'
    # number of wanted configurations (currently tested only powers of 2 ):
    num_configs = 8
    mol_name='LiH_' + '{0}'.format(round(dist,1))
    ci1 = CIMatrixGenerator(h2_mol,mol_name,basis) 
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type,max_configs=32) # use full ci for finding the most significant configs

    ci1.get_reduced_cimatrix(num_configs) 
    

def generate_H2_vs_distance():
    distance = np.arange(0.3,2.1,0.1)
    for d in distance:
        generate_h2_ci_matrix(d)


def generate_LiH_vs_distance():
    distance = np.arange(0.5,4.2,0.2)
    distance = [round(d,2) for d in distance]
    for d in distance:
        generate_LiH_ci_matrix(d)



def generate_ammonia_ci_matrix():
    print("generating ci matrices!")
    NH3_mol =  """

        N 0.0 0.0 0.2851841
        H -0.4697490 0.8136292 -0.0950614
        H -0.4697490 -0.8136292 -0.0950614
        H 0.9394980 0.0 -0.0950614
        symmetry c1
        """
    basis = 'sto-3g'
    # number of wanted configurations (currently tested only powers of 2 ):
    num_configs = 2
    mol_name = 'NH3'
    ci1 = CIMatrixGenerator(NH3_mol,mol_name,basis)
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type) # use full ci for finding the most significant configs

def generate_NH3_eq_ci_matrix():
    # equilbrium geometry taken from HF sto-3g geometry- https://cccbdb.nist.gov/energy3x.asp?method=1&basis=20&charge=0
    print("generating ci matrices!")
    NH3_mol =  """

        N 0.0    0.0     0.128
        H 0.0    0.941   -0.298
        H 0.815 -0.470  -0.298
        H -0.815 -0.470  -0.298
        symmetry c1
        """
    basis = 'sto-3g'
    # number of wanted configurations (currently tested only powers of 2 ):
    mol_name = 'NH3_NIST'
    ci1 = CIMatrixGenerator(NH3_mol,mol_name,basis)
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type,max_configs=128) # use full ci for finding the most significant configs
    
    for i in range(1,8):
        num_configs = 2**i
        ci1.get_reduced_cimatrix(num_configs)    


def generate_benzene_ci_matrix():
    benzene = """
    C          0.710500000000    0.000000000000   -1.230622098778
    C          1.421000000000    0.000000000000    0.000000000000
    C          0.710500000000    0.000000000000    1.230622098778
    C         -0.710500000000    0.000000000000    1.230622098778
    C         -0.710500000000    0.000000000000   -1.230622098778
    C         -1.421000000000    0.000000000000    0.000000000000
    H          1.254500000000    0.000000000000   -2.172857738095
    H         -1.254500000000    0.000000000000    2.172857738095
    H          2.509000000000    0.000000000000    0.000000000000
    H          1.254500000000    0.000000000000    2.172857738095
    H         -1.254500000000    0.000000000000   -2.172857738095
    H         -2.509000000000    0.000000000000    0.000000000000
    symmetry c1

    """
    """
        C        0.0000      1.3862 	0.0000 
        C 	     1.2005 	 0.6931 	0.0000 
        C 	     1.2005 	-0.6931 	0.0000
        C 	     0.0000 	-1.3862 	0.0000
        C 	    -1.2005 	-0.6931 	0.0000
        C 	    -1.2005 	 0.6931 	0.0000
        H 	     0.0000 	 2.4617 	0.0000
        H 	     2.1319 	 1.2309 	0.0000
        H 	     2.1319 	-1.2309 	0.0000
        H 	     0.0000 	-2.4617 	0.0000
        H 	    -2.1319 	-1.2309 	0.0000
        H       -2.1319 	 1.2309 	0.0000
        symmetry c1
    """

    # benzene sto3g optimized from NIST by MP2/STO-3G :
    benzene_nist = """
        C    0.000 	1.416 	0.000
        C 	1.226 	0.708 	0.000
        C 	1.226 	-0.708 	0.000
        C 	0.000 	-1.416 	0.000
        C 	-1.226 	-0.708 	0.000
        C 	-1.226 	0.708 	0.000
        H 	0.000 	2.517 	0.000
        H 	2.180 	1.259 	0.000
        H 	2.180 	-1.259 	0.000
        H 	0.000 	-2.517 	0.000
        H 	-2.180 	-1.259 	0.000
        H 	-2.180 	1.259 	0.000
        symmetry c1
    """
    basis = 'sto-3g'
    mol_name = 'benzene_nist_MP2_32768'
    ci1 = CIMatrixGenerator(benzene_nist,mol_name,basis)
    ci_type = 'cisdt'
    ci1.generate_ci_results(ci_type,max_configs=32768) # use full ci for finding the most significant configs
    
    for i in range(14,16):
        num_configs = 2**i
        ci1.get_reduced_cimatrix(num_configs)    



def generate_C2H4_ci_matrix():
    #equilbrium geometry taken from HF sto-3g geometry- https://cccbdb.nist.gov/energy3x.asp?method=1&basis=20&charge=0
    
    C2H4 = """
    C          0.0   0.0      0.653
    C          0.0   0.0     -0.653
    H          0.0   0.916    1.229
    H          0.0  -0.916    1.229
    H          0.0  -0.916   -1.229
    H          0.0   0.916   -1.229
    
    symmetry c1

    """

    basis = 'sto-3g'
    mol_name = 'C2H4'
    ci1 = CIMatrixGenerator(C2H4,mol_name,basis)
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type,max_configs=16384) # use full ci for finding the most significant configs
    
    for i in range(1,14):
        num_configs = 2**i
        ci1.get_reduced_cimatrix(num_configs)    


def generate_BeH2_matrices():
    print("generating ci matrices!")
    BeH2_mol =  """

        0 1
        Be 0 0 0 
        H 0 0 1.3
        H 0 0 -1.3
        symmetry c1
        """
    basis = 'sto-3g'
    # number of wanted configurations (currently tested only powers of 2 ):
    num_configs = 2
    mol_name = 'BeH2_1.3'
    ci1 = CIMatrixGenerator(BeH2_mol,mol_name,basis)
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type,max_configs=128) # use full ci for finding the most significant configs
    
     
    for i in range(1,8):
        num_configs = 2**i
        ci1.get_reduced_cimatrix(num_configs)    

def generate_BeH2_ci_matrix(dist):
    print("generating ci matrices!")
    BeH2_mol =  """

        0 1
        Be 0 0 0 
        H 0 0 {0}
        H 0 0 -{0}
        symmetry c1
        """.format(dist)
    basis = 'sto-3g'
    # number of wanted configurations (currently tested only powers of 2 ):
    num_configs = 16
    mol_name='BeH2_' + '{0}'.format(round(dist,1))
    ci1 = CIMatrixGenerator(BeH2_mol,mol_name,basis) 
    ci_type = 'fci'
    ci1.generate_ci_results(ci_type,max_configs=16) # use full ci for finding the most significant configs

    ci1.get_reduced_cimatrix(num_configs) 


def generate_BeH2_vs_distance():
    distance = np.arange(0.5,4.2,0.2)
    distance = [round(d,2) for d in distance]
    for d in distance:
        generate_BeH2_ci_matrix(d)

    
if __name__ == "__main__":
    # generate_BeH2_matrices()
    # generate_BeH2_vs_distance()
    # generate_h20_ci_matrices()
    # generate_h2_ci_matrix()
    # generate_ammonia_ci_matrix()
    # generate_benzene_ci_matrix()
    # generate_H2_vs_distance()
    # generate_LiH_vs_distance()
    # generate_C2H4_ci_matrix()
    generate_NH3_eq_ci_matrix()
