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
    fig = plt.figure(figsize=(5.7,5.2))
    ax = fig.add_subplot(111)
   
    distances = total_results['statevector']['distance']
  

    exact_energies = np.array(total_results['statevector']['exact_energy_list'])
    l = len(exact_energies)
    chemical_accuracy_line = np.array(l*[chemical_accuracy])

    # exact energies (diagonalization)
    plt.plot(distances,exact_energies,'-gd',color='black',label='Exact diagonalization',linestyle="dotted")
    plt.fill_between(distances,exact_energies - chemical_accuracy_line, exact_energies + chemical_accuracy_line,color='grey')
    
    statevector = total_results['statevector']
    plt.errorbar(distances,statevector['vqe_averaged_energy_list'],label='Noisless statevector simulation',color='blue',fmt='o')

    noisy_simulator = total_results['noisy_simulator']
    plt.errorbar(distances,noisy_simulator['vqe_averaged_energy_list'],label='Noisy simulation (IBMQ Santiago)',color='purple',fmt='^',
        yerr = noisy_simulator['vqe_averaged_energy_std_list'],ecolor = 'purple',capsize=2)
    
    
    IBM_real = total_results['IBM_real']
    IBM_distances = IBM_real['distance']

    plt.errorbar(IBM_distances,IBM_real['vqe_averaged_energy_list'],yerr=IBM_real['vqe_averaged_energy_std_list'],
        ecolor = 'red',capsize=2,label='Real hardware (IBMQ Quito)',color='red',fmt='s')
        # yerr = error_bar,linestyle="dotted",ecolor = 'red',label='std',capsize=3)
   
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    plt.xlabel("Interatomic distance [Angstrom]",fontsize=15)
    plt.ylabel("Energy [Ha]",fontsize=15)
    plt.xticks(distances[::2])

    plt.legend()
   
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()    
    
    fig = plt.figure(figsize=(5.7,5.2))
    ax = fig.add_subplot(111)
    
    
    diff_real = []
    for idx,d in enumerate(IBM_distances):
        z = IBM_real['vqe_averaged_energy_list'][idx] - IBM_real['exact_energy_list'][idx]
        diff_real.append(z)
    
    # plt.errorbar(IBM_distances,diff_real,yerr=IBM_real['vqe_averaged_energy_std_list'],
    #     ecolor = 'red',capsize=2,label='Real hardware (IBMQ Lima)',color='red',fmt='s')

    diff_sim = []
    for idx,d in enumerate(distances):
        p = noisy_simulator['vqe_averaged_energy_list'][idx] - exact_energies[idx]
        diff_sim.append(p)
    
    plt.errorbar(distances,diff_sim,yerr=noisy_simulator['vqe_averaged_energy_std_list'],
        ecolor = 'purple',capsize=2,label='Noisy simulation (IBMQ Santiago)',color='purple',fmt='^')

    z = np.array(len(distances)*[0])
    plt.errorbar(distances,len(distances)*[0],color='green',linestyle='dashed')
    plt.fill_between(distances,z - chemical_accuracy_line, z + chemical_accuracy_line,color='palegreen')

    text = "Hydrogen Molecule - Energy error (VQE - Exact)"
    print(ax.get_data_ratio())
    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    plt.xlabel("Interatomic distance [Angstrom]",fontsize=15)
    plt.ylabel("Energy error [Ha]",fontsize=15)
    plt.xticks(distances[::2])
    plt.legend(loc='upper left')
    figure_size = plt.gcf().get_size_inches()
    print(figure_size)
    # plt.tight_layout()
    plt.savefig(error_output_file)
    plt.show()

if __name__ == "__main__":
    # generate_vqe_results()
    # results_file_path = 'H2_dist/H2_E_vs_dist.json'
    molecules_dir = 'LiH_dist_8configs'
    results_file_name = "LiH_dist_8configs_Tue Jul 19 11:12:40 2022"
    results_file_path = molecules_dir+ '/' + results_file_name + '.json'
    output_file = molecules_dir+ '/' + results_file_name +'.png'
    error_output_file = molecules_dir + '/error_' + results_file_name + '.png'
    plot_results(results_file_path,output_file,error_output_file)