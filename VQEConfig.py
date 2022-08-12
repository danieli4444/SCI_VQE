import numpy as np
import json
import os
from time import time,ctime
from CIMatrixGenerator import CI_RESULT_FILE


VQECONFIG_FILE = 'vqe_config.json'

class VQEConfig:
    """ 
    the aim of this class is to create a wll defined configuration for the VQE runner
    so the VQERunner results can be perfectly recreated.
    for this we store the VQE relevant parameters in the vqe_config.json file in the VQEresult dir under molecule dir.
    it is important to notice that a molecule is a specific molecule geomtry.
    """
    def __init__(self,mol_dir,ci_matrix_path,
            output_name = '',
            backend_type='simulator',
            backend_name='None',
            provider = 'ibm-q/open/main',
            encoding_type='efficient',
            simulator='state_vector',
            ansatz_layers=1,
            entanglement_type='linear',
            max_iterations=100,
            shots=1000,
            use_noise=False,
            use_error_mit=False,
            cals_matrix_refresh_period=30) -> None:
        
        self.output_name = output_name
        self.mol_dir = mol_dir # path to the directory where the ci files are located. the VQE results will be saved here.
        self.ci_matrix_filepath = ci_matrix_path
        self.backend_type = backend_type # can be simulator or real hardware
        self.backend_name = backend_name # if its a simulator with noise should be enterted or if backend_type is real hardware.else leave empty.
        self.provider = provider
        self.encoding_type = encoding_type
        self.ci_params_dict = self.get_ci_params() # contains the nuclear repulsion,basis set, fci energy and other params
        self.shift = self.ci_params_dict['nuclear_repulsion_energy']
        self.simulator = simulator
        
        self.ansatz_layers = ansatz_layers
        self.entanglement_type = entanglement_type
        self.max_iterations = max_iterations
        self.shots = shots
        self.use_noise = use_noise
        self.use_error_mit = use_error_mit
        self.cals_matrix_refresh_period = cals_matrix_refresh_period

        self.output_dir = self.create_output_dir() 

    def get_ci_params(self):
        """
        the Nuclear - Nuclear potential energy for the molecule which is computed in by the CIMatrixGenerator.
        """
        filename = self.mol_dir + '/' + CI_RESULT_FILE
        with open(filename,'r') as f:
            j = json.load(f)
        return j

    def toJson(self):
        filename = self.output_dir + '/' + VQECONFIG_FILE
        with open(filename,'w') as f:
             json.dump(self.__dict__,f,indent=4)

    def create_output_dir(self):
        dirname = self.mol_dir + '/' + self.output_name + '_VQE_result_' + ctime()
        path = os.getcwd() + '/' + dirname
        print("\nVQERunner output dir path:")
        print(path)
        print("\n")
        if not os.path.exists(path):
            os.mkdir(path)            
        return dirname
