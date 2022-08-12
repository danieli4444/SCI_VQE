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
    plt.xlabel("# Determinants ",fontsize=15)
    plt.ylabel("Energy [Ha]",fontsize=15)
    text = "H2O - Energy vs number of configurations"

    configs = total_results['statevector']['configs']
    exact_energies = np.array(total_results['statevector']['exact_energy_list'])
    l = len(exact_energies)
    chemical_accuracy_line = np.array(l*[chemical_accuracy])

    # exact energies (diagonalization)
    exact_plot, = plt.plot(configs,exact_energies,'-gd',color='black',label='Exact diagonalization',linestyle='dotted')
    plt.fill_between(configs,exact_energies - chemical_accuracy_line, exact_energies + chemical_accuracy_line,color='lightgrey')
    
    statevector = total_results['statevector']
    noiseless_plot = plt.errorbar(configs,statevector['vqe_averaged_energy_list'],label='Noisless statevector simulation',color='blue',fmt='o')
    plt.xticks(configs)
    y_axis = np.arange(-228.255, max(exact_energies)+0.1, 0.01)
    plt.yticks(y_axis)

    
    noisy_simulator = total_results['noisy_simulator']
    noisy_configs = noisy_simulator['configs']
    vqe10iterationsavg = noisy_simulator['final_averaged_energy']
    vqe10iterationsstd = noisy_simulator['final_averaged_std']
    noisy_plot = plt.errorbar(noisy_configs,vqe10iterationsavg,yerr=vqe10iterationsstd,
        ecolor='purple' ,capsize=2,label='Noisy simulation (IBMQ Santiago)',color='purple',fmt='^')
    

    IBM_real = total_results['IBM_real']
    IBM_configs= IBM_real['configs']
    IBM_VQE = IBM_real['vqe_averaged_energy_list']
    IBM_plot = plt.errorbar(IBM_configs,IBM_VQE,yerr=IBM_real['vqe_averaged_energy_std_list'],
        ecolor = 'red',capsize=2,label='Real Hardware (IBMQ Lima)',color='red',fmt='s')


    # FCI_energy = -228.2551930614793
    # fci = 5*[FCI_energy]
    # x = [2,4,8,16,32]
    # fci_plot = plt.plot(x,fci,color='green',label='FCI (sto-3g)',linestyle='dashed')
    # chemical_accuracy_line = np.array(5*[chemical_accuracy])
    # plt.fill_between(x,fci - chemical_accuracy_line, fci + chemical_accuracy_line,color='palegreen')
    
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_file)
    plt.show()


   
if __name__ == "__main__":
    results_file_name = "benzene_lima_Sat Jul 16 10:55:45 2022"
    molecule_dir = "benzene_nist_MP2_32768"
    results_file_path = molecule_dir+ '/' + results_file_name + '.json'
    output_file = molecule_dir+ '/' + results_file_name +'.png'
    error_output_file = molecule_dir + '/error_' + 'new' + results_file_name + '.png'
    plot_results(results_file_path,output_file,error_output_file)