"""

Result:
    data for running VQE for H2 with chosen distances and different VQE configurations.

"""

import json
from VQEConfig import *
from VQERunner import *
import os
import shutil
from matplotlib import pyplot as plt
from time import ctime

running_options = ['statevector','noisy_simulator','IBM_real']



def clear_vqe_results(molecules_dir,distances,dir_to_delete):
    for d in distances:
        mol_name= molecules_dir+ '/LiH_' + '{0}'.format(round(d,1))
        for x in os.walk(mol_name):
            subdir_name = x[0]
            if dir_to_delete in subdir_name:
                # os.system("rm -f {0}".format(subdir_name))
                shutil.rmtree(subdir_name)
                print("deleted {0}".format(subdir_name))


def generate_vqe_results(distances,molecules_dir,running_options):
    # running_options = ['statevector','noisy_simulator','IBM_real']

    for run_opt in running_options:
        if run_opt == 'statevector':
            output_name = 'statevector'        
            simulator = 'statevector_simulator'      
            for d in distances:
                mol_name = '/LiH_' + '{0}'.format(d)
                mol_dir= molecules_dir+ mol_name
                ci_matrix_file = '{0}/CI_matrices/{1}_cimat__8.out'.format(mol_dir,mol_name)

                vqeconfig = VQEConfig(
                    mol_dir=mol_dir,
                    ci_matrix_path=ci_matrix_file,
                    output_name = output_name,
                    backend_type='simulator',
                    backend_name='None',
                    encoding_type='efficient',
                    simulator= simulator,
                    ansatz_layers=2,
                    entanglement_type='circular',
                    max_iterations=100,
                    shots=10000,
                    use_noise=False,
                    use_error_mit=False,
                    cals_matrix_refresh_period=30
                )
                new_runner = VQERunner(vqeconfig)
                new_runner.run_VQE_qasm()
            
        if run_opt == 'noisy_simulator':
            output_name = 'noisy_simulator'        
            simulator = 'qasm_simulator'      
            for d in distances:
                mol_name = '/LiH_' + '{0}'.format(d)
                mol_dir= molecules_dir+ mol_name
                ci_matrix_file = '{0}/CI_matrices/{1}_cimat__8.out'.format(mol_dir,mol_name)

                vqeconfig = VQEConfig(
                    mol_dir=mol_dir,
                    ci_matrix_path=ci_matrix_file,
                    output_name = output_name,
                    backend_type='simulator',
                    backend_name='ibmq_santiago',
                    encoding_type='efficient',
                    simulator='qasm_simulator',
                    ansatz_layers=2,
                    entanglement_type='circular',
                    max_iterations=100,
                    shots=20000,
                    use_noise=True,
                    use_error_mit=True,
                    cals_matrix_refresh_period=30
                )
                new_runner = VQERunner(vqeconfig)
                new_runner.run_VQE_qasm()

        if run_opt == 'IBM_real':
            output_name = 'IBM_real'        
            IBM_distances = [2.1,4.1]
            # distances = []
            for d in IBM_distances:
                mol_name = '/LiH_' + '{0}'.format(d)
                mol_dir= molecules_dir+ mol_name
                ci_matrix_file = '{0}/CI_matrices/{1}_cimat__8.out'.format(mol_dir,mol_name)

                vqeconfig = VQEConfig(
                    mol_dir=mol_dir,
                    ci_matrix_path=ci_matrix_file,
                    output_name = output_name,
                    backend_type='IBMQ',
                    backend_name='ibmq_quito',
                    encoding_type='efficient',
                    simulator='',
                    ansatz_layers=2,
                    entanglement_type='circular',
                    max_iterations=100,
                    shots=20000,
                    use_noise=True,
                    use_error_mit=True,
                    cals_matrix_refresh_period=30
                )
                new_runner = VQERunner(vqeconfig)
                new_runner.run_VQE_qasm()

    
def extract_all_results(output_file,distances,molecules_dir,running_options):
    
    total_results = {}
    output_filename = molecules_dir + '/' + output_file +'_' + ctime() + '.json'
    for run_opt in running_options:
        exact_energy_list = []
        vqe_optimal_energy_list = []
        vqe_std_list = []
        vqe_minimal_energy_list = []
        vqe_averaged_energy_list = []
        vqe_averaged_energy_std_list = []
        for d in distances:
            mol_name = '/LiH_' + '{0}'.format(d)
            mol_dir= molecules_dir+ mol_name
            results_dir_list = [x[0] for x in os.walk(mol_dir)]
            for result_dir in results_dir_list:
                if run_opt in result_dir:
                    vqe_result_file = result_dir + '/VQE_output.json'
                    vqe_config_file = result_dir + '/vqe_config.json'
                    with open(vqe_result_file) as json_file:
                        vqe_result = json.load(json_file)
                        vqe_optimal_energy_list.append(float(vqe_result['optimal_energy']))
                        exact_energy_list.append(float(vqe_result['exact_ground_energy']))
                        vqe_minimal_energy_list.append(float(vqe_result['minimal_energy']))
                        vqe_std_list.append(float(vqe_result['optimal_energy_std']))
                        vqe_averaged_energy_list.append(float(vqe_result['averaged_energy']))
                        vqe_averaged_energy_std_list.append(float(vqe_result['averaged_energy_std']))
                        


        current_opt_results = {'distance': distances,
            'vqe_optimal_energy_list': vqe_optimal_energy_list,
            'vqe_minimal_energy_list': vqe_minimal_energy_list,
            'exact_energy_list': exact_energy_list,
            'vqe_std_list': vqe_std_list,
            'vqe_averaged_energy_list': vqe_averaged_energy_list,
            'vqe_averaged_energy_std_list': vqe_averaged_energy_std_list 
            }

        total_results['{0}'.format(run_opt)] = current_opt_results
                    
    with open(output_filename,'w') as f:
        json.dump(total_results,f,indent=4)
        
def extract_all_results_2(output_file,distances,molecules_dir,running_options):
    total_results = {}
    output_filename = molecules_dir + '/' + output_file +'_' + ctime() + '.json'
    for run_opt in running_options:
        exact_energy_list = []
        vqe_optimal_energy_list = []
        vqe_std_list = []
        vqe_minimal_energy_list = []
        vqe_averaged_energy_list = []
        vqe_averaged_energy_std_list = []
        for d in distances:
            mol_name = '/LiH_' + '{0}'.format(d)
            mol_dir= molecules_dir+ mol_name
            results_dir_list = [x[0] for x in os.walk(mol_dir)]
            for result_dir in results_dir_list:
                if run_opt in result_dir:

                    vqe_result_file = result_dir + '/VQE_output.json'
                    vqe_config_file = result_dir + '/vqe_config.json'
                    with open(vqe_result_file) as json_file:
                        vqe_result = json.load(json_file)
                        vqe_optimal_energy_list.append(float(vqe_result['optimal_energy']))
                        exact_energy_list.append(float(vqe_result['exact_ground_energy']))
                        vqe_minimal_energy_list.append(float(vqe_result['minimal_energy']))
                        vqe_std_list.append(float(vqe_result['optimal_energy_std']))
                        vqe_averaged_energy_list.append(float(vqe_result['averaged_energy']))
                       
                        with open(vqe_config_file) as config_file:
                            vqe_config = json.load(config_file)
                            shift = float(vqe_config['shift'])

                        vqe_values = vqe_result['vqe_result_params']['values'][-11:-1]
                        vqe_values = [val+shift for val in vqe_values]
                        average_energy = np.average(vqe_values)
                        average_energy_std = np.std(vqe_values)
                        vqe_averaged_energy_std_list.append(average_energy_std)
                        


        current_opt_results = {'distance': distances,
            'vqe_optimal_energy_list': vqe_optimal_energy_list,
            'vqe_minimal_energy_list': vqe_minimal_energy_list,
            'exact_energy_list': exact_energy_list,
            'vqe_std_list': vqe_std_list,
            'vqe_averaged_energy_list': vqe_averaged_energy_list,
            'vqe_averaged_energy_std_list': vqe_averaged_energy_std_list 
            }

        total_results['{0}'.format(run_opt)] = current_opt_results
                    
    with open(output_filename,'w') as f:
        json.dump(total_results,f,indent=4)


if __name__ == "__main__":
    molecules_dir = 'LiH_dist_8configs'
    # distances = np.arange(0.5,4.2,0.2)
    # distances = [round(d,2) for d in distances]
    distances = []
    IBM_distances = [2.1,4.1]
    # clear_vqe_results(molecules_dir,distances,'statevector')
    # clear_vqe_results(molecules_dir,distances,'noisy_simulator')
    # clear_vqe_results(molecules_dir,IBM_distances,'IBM_real')

    generate_vqe_results(distances,molecules_dir ,running_options)

    output_file = 'LiH_dist_8configs'
    extract_all_results(output_file,distances, molecules_dir,running_options)
    # extract_all_results_2(output_file,distances, molecules_dir,running_options)
