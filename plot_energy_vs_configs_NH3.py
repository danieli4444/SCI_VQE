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
    fig = plt.figure(figsize=(5.7,5.2),dpi=100)
    plt.xlabel("#Determinants [#Qubits] ",fontsize=15)
    plt.ylabel("Energy [Ha]",fontsize=15)

    configs = total_results['statevector']['configs']
    exact_energies = np.array(total_results['statevector']['exact_energy_list'])
    l = len(exact_energies)
    chemical_accuracy_line = np.array(l*[chemical_accuracy])

    FCI_energy = -55.524589065582546
    fci = 6*[FCI_energy]
    x = [2,4,8,16,32,64]
    fci_plot = plt.plot(x,fci,color='green',label='FCI (sto-3g)',linestyle='dashed')
    chemical_accuracy_line = np.array(6*[chemical_accuracy])
    plt.fill_between(x,fci - chemical_accuracy_line, fci + chemical_accuracy_line,color='palegreen',zorder=2)

    # exact energies (diagonalization)
    exact_plot, = plt.plot(configs,exact_energies,'-gd',color='grey',label='SCI exact diagonalization',linestyle='dotted')
    plt.fill_between(configs,exact_energies - chemical_accuracy_line, exact_energies + chemical_accuracy_line,color='lightgrey')
    
    statevector = total_results['statevector']
    noiseless_plot = plt.errorbar(configs,statevector['vqe_averaged_energy_list'],label='Noiseless statevector simulation',color='blue',fmt='o')
    xticks_q = ['2[1]','4[2]','8[3]','16[4]','32[5]','64[6]']
    plt.xticks([2, 4, 8, 16, 32, 64], ["2\n[1]", "4\n[2]", "8\n[3]", "16\n[4]", "32\n[5]","64\n[6]"])
    # plt.xticks(configs,labels=xticks_q)
    y_axis = np.arange(-55.52, max(exact_energies)+0.1, 0.005)
    plt.yticks(y_axis)

    
    noisy_simulator = total_results['noisy_simulator']
    noisy_configs = noisy_simulator['configs']
    vqe10iterationsavg = noisy_simulator['vqe_averaged_energy_list']
    vqe10iterationsstd = noisy_simulator['vqe_averaged_energy_std_list']
    noisy_plot = plt.errorbar(noisy_configs,vqe10iterationsavg,yerr=vqe10iterationsstd,
        ecolor='purple' ,capsize=2,label='Noisy simulation (IBMQ Santiago)',color='purple',fmt='^')
    

    IBM_real = total_results['IBM_real']
    IBM_configs= IBM_real['configs']
    IBM_VQE = IBM_real['vqe_averaged_energy_list']
    IBM_plot = plt.errorbar(IBM_configs,IBM_VQE,yerr=IBM_real['vqe_averaged_energy_std_list'],
        ecolor = 'red',capsize=2,label='Real Hardware (IBMQ Nairobi)',color='red',fmt='s')


   
    
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_file)
    plt.show()


   
if __name__ == "__main__":
    results_file_name = "NH3_NIST_Wed Oct 12 10:18:05 2022"
    molecule_dir = "NH3_NIST"
    results_file_path = molecule_dir+ '/' + results_file_name + '.json'
    output_file = molecule_dir+ '/' + results_file_name +'.png'
    error_output_file = molecule_dir + '/error_' + results_file_name + '.png'
    plot_results(results_file_path,output_file,error_output_file)