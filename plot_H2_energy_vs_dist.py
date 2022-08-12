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
    
    fig = plt.figure(figsize=(5.5,5))
    ax = fig.add_subplot(111)
    plt.xlabel("Interatomic distance [Angstrom]",fontsize=15)
    plt.ylabel("Energy [Ha]",fontsize=15)
    
    distances = total_results['statevector']['distance']
    plt.xticks(distances[::2])

    exact_energies = np.array(total_results['statevector']['exact_energy_list'])
    l = len(exact_energies)
    chemical_accuracy_line = np.array(l*[chemical_accuracy])

    # exact energies (diagonalization)
    plt.plot(distances,exact_energies,'-gd',color='black',label='Exact diagonalization',linestyle="dotted")
    plt.fill_between(distances,exact_energies - chemical_accuracy_line, exact_energies + chemical_accuracy_line,color='grey')
    
    statevector = total_results['statevector']
    plt.errorbar(distances,statevector['vqe_averaged_energy_list'],label='Noisless statevector simulation',color='blue',fmt='o')

    noisy_simulator = total_results['noisy_simulator']
    plt.errorbar(distances,noisy_simulator['vqe_averaged_energy_list'],label='Noisy simulation (IBMQ Lima)',color='purple',fmt='^',
        yerr = noisy_simulator['vqe_averaged_energy_std_list'],ecolor = 'purple',capsize=2)
    
    
    IBM_real = total_results['IBM_real']
    IBM_distances = [0.3, 0.5, 0.7, 0.8,1.2, 2.0]
    # IBM_distances = []
    plt.errorbar(IBM_distances,IBM_real['vqe_averaged_energy_list'],yerr=IBM_real['vqe_averaged_energy_std_list'],
        ecolor = 'red',capsize=2,label='Real hardware (IBMQ Lima)',color='red',fmt='s')
        # yerr = error_bar,linestyle="dotted",ecolor = 'red',label='std',capsize=3)
    
    
    # ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.legend()
    plt.tight_layout()

    plt.savefig(output_file)
    plt.show()


    fig = plt.figure(figsize=(5.5,4.5))
    ax = fig.add_subplot(111)

  

    diff_sim = []
    for idx,d in enumerate(distances):
        p = noisy_simulator['vqe_averaged_energy_list'][idx] - exact_energies[idx]
        diff_sim.append(p)
    
    plt.errorbar(distances,diff_sim,yerr=noisy_simulator['vqe_averaged_energy_std_list'],
        ecolor = 'purple',capsize=2,label='Noisy simulation (IBMQ Lima)',color='purple',fmt='^')

    diff_real = []
    for idx,d in enumerate(IBM_distances):
        i = distances.index(d)
        z = IBM_real['vqe_averaged_energy_list'][idx] - exact_energies[i]
        diff_real.append(z)
    
    plt.errorbar(IBM_distances,diff_real,yerr=IBM_real['vqe_averaged_energy_std_list'],
        ecolor = 'red',capsize=2,label='Real hardware (IBMQ Lima)',color='red',fmt='s')
  
    z = np.array(len(distances)*[0])
    plt.errorbar(distances,len(distances)*[0],color='green',linestyle='dashed')
    plt.fill_between(distances,z - chemical_accuracy_line, z + chemical_accuracy_line,color='palegreen')

    text = "Hydrogen Molecule - Energy error (VQE - Exact)"
    # y_axis = np.arange(-0.006, 0.007, 0.001)
    # plt.yticks(y_axis)
    # ratio = 0.007011596814686329
    # asp = 1.0/ax.get_data_ratio()
    # ax.set_aspect('auto')
    plt.xlabel("Interatomic distance [Angstrom]",fontsize=15)
    plt.ylabel("Energy error [Ha]",fontsize=15)
    plt.xticks(distances[::2])
    plt.legend(loc='upper left')

    plt.tight_layout()
    plt.savefig(error_output_file)
    plt.show()

if __name__ == "__main__":
    # generate_vqe_results()
    # results_file_path = 'H2_dist/H2_E_vs_dist.json'
    molecules_dir = 'H2_dist'
    results_file_name = "H2_E_vs_dist_Tue Jun 14 20:25:13 2022"
    results_file_path = molecules_dir+ '/' + results_file_name + '.json'
    output_file = molecules_dir+ '/' + results_file_name +'.png'
    error_output_file = molecules_dir + '/error_' + results_file_name + '.png'
    plot_results(results_file_path,output_file,error_output_file)