import json
from turtle import distance
import numpy as np
from matplotlib import pyplot as plt

running_options = ['statevector','noisy_simulator','IBM_real']

def plot_results(results_file_path,output_file,error_output_file):

    FCI_energy = -77.23211357330942

    with open(results_file_path,'r') as f:
        total_results = json.load(f)
    
    # state_vector = total_results['statevector']
    # noisy_simulator = total_results['noisy_simulator']
    # IBM_real = total_results['IBM_real']
    chemical_accuracy = 0.0016

    fig = plt.figure(figsize=(5.7,5.2),dpi=100)
    ax = fig.add_subplot(111)


    configs = total_results['statevector']['configs']

    exact_configs = total_results['exact']['configs']
    exact_energies = np.array(total_results['exact']['exact_energy_list'])


    l = len(exact_energies)
    chemical_accuracy_line = np.array(l*[chemical_accuracy])

    fci = len(exact_configs)*[FCI_energy]
    fci_plot = plt.plot(exact_configs,fci,color='green',label='FCI (sto-3g)',linestyle='dashed')
    chemical_accuracy_line = np.array(len(exact_configs)*[chemical_accuracy])
    plt.fill_between(exact_configs,fci - chemical_accuracy_line, fci + chemical_accuracy_line,color='palegreen',zorder=2)
    

    # exact energies (diagonalization)
    plt.plot(exact_configs,exact_energies,'-gd',color='grey',label='SCI exact diagonalization',linestyle='dotted')
    plt.fill_between(exact_configs,exact_energies - chemical_accuracy_line, exact_energies + chemical_accuracy_line,color='lightgrey')
    
    statevector = total_results['statevector']
    noiseless_plot = plt.errorbar(configs,statevector['vqe_averaged_energy_list'],yerr=statevector['vqe_averaged_energy_std_list'],
        label='Noisless statevector simulation',color='blue',fmt='o')
    # plt.xticks(exact_configs)
    ax.set_xscale('log', basex=2)

    plt.xticks([2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 ], 
        ["2\n[1]", "4\n[2]", "8\n[3]", "16\n[4]", "32\n[5]", "64\n[6]", "128\n[7]",
        "256\n[8]", "","1024\n[10]", "", "4096\n[12]" ])
    y_axis = np.arange(-77.232, max(exact_energies)+0.1, 0.010)
    plt.yticks(y_axis)

    
    noisy_simulator = total_results['noisy_simulator']
    noisy_configs = noisy_simulator['configs'] 
    noisy_configs = noisy_configs[:-1]
    vqe10iterationsavg = noisy_simulator['vqe_averaged_energy_list']
    vqe10iterationsstd = noisy_simulator['vqe_averaged_energy_std_list']
    vqe10iterationsavg = vqe10iterationsavg[:-1]
    vqe10iterationsstd = vqe10iterationsstd[:-1]
    noisy_plot = plt.errorbar(noisy_configs,vqe10iterationsavg,yerr=vqe10iterationsstd,
        ecolor='purple' ,capsize=2,label='Noisy simulation (IBMQ Santiago)',color='purple',fmt='^')
    

    IBM_real = total_results['IBM_real']
    IBM_configs= IBM_real['configs']
    IBM_VQE = IBM_real['vqe_averaged_energy_list']
    IBM_plot = plt.errorbar(IBM_configs,IBM_VQE,yerr=IBM_real['vqe_averaged_energy_std_list'],
        ecolor = 'red',capsize=2,label='Real Hardware (IBMQ Nairobi)',color='red',fmt='s')


    
    
    
    
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.xlabel("#Determinants [#Qubits]",fontsize=15)
    plt.ylabel("Energy [Ha]",fontsize=15)

    plt.legend(loc='upper right',fontsize='small')
    plt.tight_layout()

    plt.savefig(output_file)

    plt.show()


   
if __name__ == "__main__":
    results_file_name = "C2H4_Wed Sep 21 12:10:23 2022"
    molecule_dir = "C2H4"
    results_file_path = molecule_dir+ '/' + results_file_name + '.json'
    output_file = molecule_dir+ '/' + results_file_name +'.png'
    error_output_file = molecule_dir + '/error_' + results_file_name + '.png'
    plot_results(results_file_path,output_file,error_output_file)