import json
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
    # fig = plt.figure(figsize=(5.7,5.2))
    # ax = fig.add_subplot(111)
    fig = plt.figure(figsize=(5.7,5.2),dpi=100)
   
   
    distances = total_results['statevector']['distance']
  
    exact_energies = np.array(total_results['statevector']['exact_energy_list'])
    l = len(exact_energies)
    chemical_accuracy_line = np.array(l*[chemical_accuracy])

    FCI_energies = total_results['FCI']['FCI_energy']
    plt.plot(distances,FCI_energies,color='green',label='FCI (sto-3g)',linestyle="dashed")
    plt.fill_between(distances,FCI_energies - chemical_accuracy_line, FCI_energies + chemical_accuracy_line,color='palegreen',zorder=2)
    

    # exact energies (diagonalization)
    plt.plot(distances,exact_energies,'-gd',color='grey',label='SCI exact diagonalization',linestyle="dotted")
    
    
    statevector = total_results['statevector']
    plt.errorbar(distances,statevector['vqe_averaged_energy_list'],label='Noisless statevector simulation',color='blue',fmt='o')

    noisy_simulator = total_results['noisy_simulator']
    plt.errorbar(distances,noisy_simulator['vqe_averaged_energy_list'],label='Noisy simulation (IBMQ Santiago)',color='purple',fmt='^',
        yerr = noisy_simulator['vqe_averaged_energy_std_list'],ecolor = 'purple',capsize=2)
    
    
    IBM_real = total_results['IBM_real']
    # IBM_distances = IBM_real['distance']
    IBM_distances = [0.5, 0.9, 1.3, 2.1, 3.1]
    
    plt.errorbar(IBM_distances,IBM_real['vqe_averaged_energy_list'],yerr=IBM_real['vqe_averaged_energy_std_list'],
        ecolor = 'red',capsize=2,label='Real hardware (IBMQ Oslo)',color='red',fmt='s')
        # yerr = error_bar,linestyle="dotted",ecolor = 'red',label='std',capsize=3)
    single_qubit_dist = [1.3]
    single_qubit_energy = -15.561835591016527
    single_qubit_std = 0.006538740464803411
    plt.errorbar(single_qubit_dist,single_qubit_energy,yerr=single_qubit_std,
        ecolor = 'orange',capsize=2,label='Real hardware-single qubit (IBMQ Oslo)',color='orange',fmt='*',markersize=8)
    
       # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    plt.xlabel("Interatomic distance [Angstrom]",fontsize=15)
    plt.ylabel("Energy [Ha]",fontsize=15)
    plt.xticks(distances[::2])

    plt.legend(fontsize = 10)
   
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()
    
    
    fig = plt.figure(figsize=(5.7,5.2),dpi=100)
    ax = fig.add_subplot(111)
    # fig = plt.figure(figsize=(5.7,5.2),dpi=100)
    
    z = np.array(len(distances)*[0])
    plt.plot(distances,len(distances)*[0],color='green',label='FCI (sto-3g)',linestyle='dashed')
    plt.fill_between(distances,z - chemical_accuracy_line, z + chemical_accuracy_line,color='palegreen',zorder=2)


    diff_diag = []
    for idx,d in enumerate(distances):
        p = exact_energies[idx] - FCI_energies[idx]
        diff_diag.append(p)
    
    plt.plot(distances,diff_diag,'-gd',label='SCI exact diagonalization',
        color='grey',linestyle='dotted')
    plt.fill_between(distances,diff_diag - chemical_accuracy_line, diff_diag + chemical_accuracy_line,color='lightgrey')


    diff_statevector = []
    for idx,d in enumerate(distances):
        p = statevector['vqe_averaged_energy_list'][idx] - FCI_energies[idx]
        diff_statevector.append(p)
    
    plt.errorbar(distances,diff_statevector,yerr=statevector['vqe_averaged_energy_std_list'],
        ecolor = 'blue',capsize=2,label='Noiseless statevector simulation',color='blue',fmt='o')


    
    # diff_real = []
    # for idx,d in enumerate(IBM_distances):
    #     z = IBM_real['vqe_averaged_energy_list'][idx] - IBM_real['exact_energy_list'][idx]
    #     diff_real.append(z)
    
    # plt.errorbar(IBM_distances,diff_real,yerr=IBM_real['vqe_averaged_energy_std_list'],
    #     ecolor = 'red',capsize=2,label='Real hardware (IBMQ Oslo)',color='red',fmt='s'

  

    diff_sim = []
    for idx,d in enumerate(distances):
        p = noisy_simulator['vqe_averaged_energy_list'][idx] - FCI_energies[idx]
        diff_sim.append(p)
    
    plt.errorbar(distances,diff_sim,yerr=noisy_simulator['vqe_averaged_energy_std_list'],
        ecolor = 'purple',capsize=2,label='Noisy simulation (IBMQ Santiago)',color='purple',fmt='^')

    single_qubit_dist = [1.3]
    single_qubit_energy_diff = single_qubit_energy + 15.595047080804816
    single_qubit_std = 0.006538740464803411
    plt.errorbar(single_qubit_dist,single_qubit_energy_diff,yerr=single_qubit_std,
        ecolor = 'orange',capsize=2,label='Real hardware-single qubit (IBMQ Oslo)',color='orange',fmt='*',markersize=8)


    # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.xlabel("Interatomic distance [Angstrom]",fontsize=15)
    plt.ylabel("Energy error [Ha]",fontsize=15)
    plt.xticks(distances[::2])
    plt.legend(loc='upper left',fontsize = 10)

    # plt.tight_layout()
    plt.savefig(error_output_file)
    plt.show()

if __name__ == "__main__":
    # generate_vqe_results()
    # results_file_path = 'H2_dist/H2_E_vs_dist.json'
    molecules_dir = 'BeH2_dist'
    results_file_name = "BeH2_E_vs_dist_Sun Oct  9 01:01:03 2022"
    results_file_path = molecules_dir+ '/' + results_file_name + '.json'
    output_file = molecules_dir+ '/' + results_file_name +'.png'
    error_output_file = molecules_dir + '/error_' + results_file_name + '.png'
    plot_results(results_file_path,output_file,error_output_file)