"""
the purpose is to generate the data to output a graph similar to the H2O ion trap paper -
https://arxiv.org/abs/1902.10171s
"""

"""

Result:
    data for running VQE for H2 with chosen distances and different VQE configurations.

"""

import json
from VQEConfig import *
from VQERunner import *
import os
import shutil
from time import ctime

running_options = ['exact','statevector','noisy_simulator','IBM_real']

final_energy_list = []
final_averaged_energy = []
final_averaged_std = []

def clear_vqe_results(molecules_dir,configs,dir_to_delete):
    for x in os.walk(molecules_dir):
        subdir_name = x[0]
        if dir_to_delete in subdir_name:
            # os.system("rm -f {0}".format(subdir_name))
            shutil.rmtree(subdir_name)
            print("deleted {0}".format(subdir_name))


def generate_vqe_results(configs,noisy_configs,IBM_configs,molecules_dir,running_options):

    for run_opt in running_options:
        if run_opt == 'statevector':
            output_name = 'statevector'        
            simulator = 'statevector_simulator'  
            # configs =   [2, 4, 8, 16, 32, 64, 128]    
            layers =    [0, 1, 2, 3, 7, 11, 18]
            for idx,c in enumerate(configs):
                ci_matrix_file = '{0}/CI_matrices/{0}_cimat__{1}.out'.format(molecules_dir,c)
                output_name2 = output_name + '_{0}'.format(c)
                if c > 16:
                    iterations = 5000
                else:
                    iterations = 500
                vqeconfig = VQEConfig(
                    mol_dir=molecules_dir,
                    ci_matrix_path=ci_matrix_file,
                    output_name = output_name2,
                    backend_type='simulator',
                    backend_name='None',
                    encoding_type='efficient',
                    simulator= simulator,
                    ansatz_layers=layers[idx],
                    entanglement_type='circular',
                    max_iterations=iterations,
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
            layers =    [0, 1, 2, 3, 7]
            iterations = [300,300,300,300,300]

            for idx,c in enumerate(noisy_configs):
                final_averaged_energy_list = []
                ci_matrix_file = '{0}/CI_matrices/{0}_cimat__{1}.out'.format(molecules_dir,c)
                output_name2 = output_name + '_{0}'.format(c)
                vqeconfig = VQEConfig(
                    mol_dir=molecules_dir,
                    ci_matrix_path=ci_matrix_file,
                    output_name = output_name2,
                    backend_type='simulator',
                    backend_name='ibmq_santiago',
                    encoding_type='efficient',
                    simulator='qasm_simulator',
                    ansatz_layers=layers[idx],
                    entanglement_type='circular',
                    max_iterations=iterations[idx],
                    shots=100000,
                    use_noise=True,
                    use_error_mit=True,
                    cals_matrix_refresh_period=30
                )
                new_runner = VQERunner(vqeconfig)
                for i in range(1,2):
                    print("i = ",i)
                    e = new_runner.run_VQE_qasm()
                    final_averaged_energy_list.append(e)

                final_energy = np.average(final_averaged_energy_list)
                final_std = np.std(final_averaged_energy_list)
                global final_averaged_energy,final_averaged_std
                final_averaged_energy.append(final_energy)
                final_averaged_std.append(final_std)
                global final_energy_list
                final_energy_list.append(final_averaged_energy_list.copy())
                print(final_averaged_energy_list)
                

        if run_opt == 'IBM_real':
            output_name = 'IBM_real'     
            # layers =    [0, 1, 2, 3, 4, 5 ]
            layers = [0,1,2]
            # iterations = [30,80,100,80,100,100]
            iterations = [100,100,100,100]
            for idx,c in enumerate(IBM_configs):
                ci_matrix_file = '{0}/CI_matrices/{0}_cimat__{1}.out'.format(molecules_dir,c)
                output_name2 = output_name + '_{0}'.format(c)
                vqeconfig = VQEConfig(
                    mol_dir=molecules_dir,
                    ci_matrix_path=ci_matrix_file,
                    output_name = output_name2,
                    backend_type='IBMQ',
                    backend_name='ibm_nairobi',
                    encoding_type='efficient',
                    simulator='',
                    ansatz_layers=layers[idx],
                    entanglement_type='circular',
                    max_iterations=iterations[idx],
                    shots=20000,
                    use_noise=True,
                    use_error_mit=True,
                    cals_matrix_refresh_period=30
                )
                new_runner = VQERunner(vqeconfig)
                new_runner.run_VQE_qasm()


    
def extract_all_results(output_file,configs,noisy_configs,IBM_configs,molecules_dir,running_options):
    
    total_results = {}
    output_filename = molecules_dir + '/' + output_file +'_' + ctime() + '.json'
    for run_opt in running_options:
        exact_energy_list = []
        vqe_optimal_energy_list = []
        vqe_std_list = []
        vqe_minimal_energy_list = []
        vqe_averaged_energy_list = []
        vqe_averaged_energy_std_list = []
        vqe_layers_used = []
        # running_options = ['statevector','noisy_simulator','IBM_real']
        if run_opt == 'noisy_simulator':
            configs = noisy_configs
        elif run_opt == 'IBM_real':
            configs = IBM_configs
             
        if run_opt == 'exact':
            exact_ci_energy_list = []
            exact_configs_list = []
            num_qubits_list = []
            ci_mat_dir = molecules_dir + '/CI_matrices'
            filename = molecules_dir+ '/' + 'CI_result.json'
            with open(filename,'r') as f:
                j = json.load(f)
                nuclear_repulsion = j['nuclear_repulsion_energy']

            for ci_mat_file in os.listdir(ci_mat_dir):
                print(ci_mat_file)
                configs_num = ci_mat_file.replace('.out','')
                configs_num = [int(s) for s in configs_num.split('_') if s.isdigit()][-1]

                ci_mat_path = ci_mat_dir + '/' + ci_mat_file
                ci_mat = np.loadtxt(ci_mat_path)
                eigenvals, eigenstates = np.linalg.eigh(ci_mat)
                energy = np.real(eigenvals[0]) + nuclear_repulsion
                exact_configs_list.append(configs_num)
                num_qubits_list.append(np.log2(configs_num))
                exact_ci_energy_list.append(energy)

            exact_configs_list, exact_ci_energy_list = zip(*sorted(zip(exact_configs_list, exact_ci_energy_list)))

            print(exact_configs_list,num_qubits_list,exact_ci_energy_list)
            
            current_opt_results = {'configs': exact_configs_list,
            'exact_energy_list': exact_ci_energy_list,
            }
        
        else:
            for c in configs:
                results_dir_list = [x[0] for x in os.walk(molecules_dir)]
                for result_dir in results_dir_list:
                    current_dir = run_opt + '_{0}'.format(c) 
                    if current_dir in result_dir:
                        vqe_result_file = result_dir + '/VQE_output.json'
                        with open(vqe_result_file) as json_file:
                            global final_averaged_std, final_energy_list,final_averaged_energy
                            vqe_result = json.load(json_file)
                            vqe_optimal_energy_list.append(float(vqe_result['optimal_energy']))
                            exact_energy_list.append(float(vqe_result['exact_ground_energy']))
                            vqe_minimal_energy_list.append(float(vqe_result['minimal_energy']))
                            vqe_std_list.append(float(vqe_result['optimal_energy_std']))
                            vqe_averaged_energy_list.append(float(vqe_result['averaged_energy']))
                            vqe_averaged_energy_std_list.append(float(vqe_result['averaged_energy_std']))
                            vqe_layers_used.append(float(vqe_result['ansatz_layers']))
       


            global final_averaged_std, final_energy_list,final_averaged_energy
            current_opt_results = {'configs': configs,
                'vqe_optimal_energy_list': vqe_optimal_energy_list,
                'vqe_minimal_energy_list': vqe_minimal_energy_list,
                'exact_energy_list': exact_energy_list,
                'vqe_std_list': vqe_std_list,
                'vqe_averaged_energy_list': vqe_averaged_energy_list,
                'vqe_averaged_energy_std_list': vqe_averaged_energy_std_list,
                'vqe_layers_used': vqe_layers_used,
                'final_averaged_energy': final_averaged_energy,
                'final_averaged_std': final_averaged_std,
                'final_energy_list':final_energy_list}


        total_results['{0}'.format(run_opt)] = current_opt_results
                    
    with open(output_filename,'w') as f:
        json.dump(total_results,f,indent=4)
        


if __name__ == "__main__":
    molecules_dir = 'C2H4'
    # configs = [2, 4, 8, 16, 32, 64]
    noiseless_configs = [2,4,8,16,32,64,128]
    noisy_configs = []
    # noisy_configs = []
    # IBM_configs = [2,4,8]
    IBM_configs = []
    clear_vqe_results(molecules_dir,noiseless_configs,'statevector')
    # clear_vqe_results(molecules_dir,noiseless_configs,'noisy_simulator')
    # clear_vqe_results(molecules_dir,configs,'IBM_real')

    # exact_daig_configs = [2,4,8,16,32,128,256,512,1024,2048,4096]
    generate_vqe_results(noiseless_configs,noisy_configs,IBM_configs,molecules_dir ,running_options)

    output_file = 'C2H4'
    extract_all_results(output_file,noiseless_configs,noisy_configs, IBM_configs, molecules_dir,running_options)
