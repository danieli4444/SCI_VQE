import json
from turtle import distance
import numpy as np
from matplotlib import pyplot as plt

running_options = ['statevector','noisy_simulator','IBM_real']

def plot_results(results_file_path,output_file,error_output_file):
    with open(results_file_path,'r') as f:
        total_results = json.load(f)
    
    # state_vector = total_results['statevector']
    # noisy_simulator = total_results['noisy_simulator']
    # IBM_real = total_results['IBM_real']
    chemical_accuracy = 0.0016
    plt.figure()
    plt.xlabel("number of configurations ")
    plt.ylabel("Energy [Ha]")
    text = "H2O - Energy vs number of configurations"
    plt.title(text)

    configs = total_results['statevector']['configs']
    exact_energies = np.array(total_results['statevector']['exact_energy_list'])
    l = len(exact_energies)
    chemical_accuracy_line = np.array(l*[chemical_accuracy])

    # exact energies (diagonalization)
    plt.plot(configs,exact_energies,color='green',label='exact',linestyle='dashed')
    plt.fill_between(configs,exact_energies - chemical_accuracy_line, exact_energies + chemical_accuracy_line,color='grey')
    
    statevector = total_results['statevector']
    plt.errorbar(configs,statevector['vqe_averaged_energy_list'],label='Noisless',color='blue',fmt='o')
    plt.xticks(configs)


    
    noisy_simulator = total_results['noisy_simulator']
    noisy_configs = noisy_simulator['configs'][0:-1]
    plt.errorbar(noisy_configs,noisy_simulator['vqe_averaged_energy_list'][0:-1],yerr=noisy_simulator['vqe_averaged_energy_std_list'][0:-1],
        ecolor='purple' ,capsize=2,label='Noisy',color='purple',fmt='o')
    
    # plt.errorbar(noisy_configs,noisy_simulator['final_averaged_energy'],yerr=noisy_simulator['final_averaged_std'],
    #     ecolor='green' ,capsize=2,label='Noisy_final_95%',color='green',fmt='o')
    


    # std2 = noisy_simulator['vqe_averaged_energy_std_list']
    # std2 = [2* z for z in std2]
    # plt.errorbar(noisy_configs,noisy_simulator['vqe_averaged_energy_list'],yerr=std2,
    #     ecolor='pink' ,capsize=2,label='Noisy_95%',color='green',fmt='o')

    IBM_real = total_results['IBM_real']
    IBM_configs= IBM_real['configs'] 
    IBM_configs = [2,4]
    IBM_VQE = IBM_real['vqe_averaged_energy_list'][0:2]
    plt.errorbar(IBM_configs,IBM_VQE,yerr=IBM_real['vqe_averaged_energy_std_list'][0:2],
        ecolor = 'red',capsize=2,label='IBMQ Lima VQE',color='red',fmt='s')
    plt.legend()
    plt.savefig(output_file)


    FCI_energy = -75.01156040647396,

    plt.show()


   
if __name__ == "__main__":
    results_file_name = "H2O_dev_Wed Jun 15 11:35:49 2022"
    molecule_dir = "H2O_dev"
    results_file_path = molecule_dir+ '/' + results_file_name + '.json'
    output_file = molecule_dir+ '/' + results_file_name +'.png'
    error_output_file = molecule_dir + '/error_' + results_file_name + '.png'
    plot_results(results_file_path,output_file,error_output_file)