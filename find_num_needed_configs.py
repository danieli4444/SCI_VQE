"""
for a specific molecule, take the available CI matrices and plot a chart that compares the CI diagonlization (with numpy)
result with the FCI result/ CIDSD (or else).
this helps to show what is the minimal number of configurations needed (or qubits) to reach chemical accuracy (from the FCI solution)
"""

import os
import numpy as np
from matplotlib import pyplot as plt

CHEMICAL_ACCURACY = 0.0016 # hartree

def create_chart(ci_mat_dir,compared_energy,chart_name):
    ci_energy_list = []
    num_configs_list = []
    num_qubits_list = []
    x = os.listdir(ci_mat_dir)
    for ci_mat_file in os.listdir(ci_mat_dir):
        configs_num = ci_mat_file.replace('.out','')
        x = configs_num.split('_')
        configs_num = [int(s) for s in configs_num.split('_') if s.isdigit()][0]

        ci_mat_path = ci_mat_dir + '/' + ci_mat_file
        ci_mat = np.loadtxt(ci_mat_path)
        eigenvals, eigenstates = np.linalg.eigh(ci_mat)
        energy = np.real(eigenvals[0])
        num_configs_list.append(configs_num)
        num_qubits_list.append(np.log2(configs_num))
        ci_energy_list.append(energy)
    print(num_configs_list,num_qubits_list,ci_energy_list)
    
    

    fig = plt.figure(figsize=(8,5.2))
    # plt.errorbar(counts,exact,yerr=self.chemical_accuracy,label='exact',linestyle="dashed",color='green',ecolor = 'black',capsize=3)
    # plt.plot(num_configs,exact,label='exact',linestyle="",marker="o")
    # plt.plot(num_configs,fci,label='fci',linestyle="dotted")
    # plt.xticks(counts)
    plt.xlabel("Number of configurations")
    plt.ylabel("Energy [Ha]")
    
    plt.xticks(num_configs_list,num_configs_list)
    plt.plot(num_configs_list,ci_energy_list,'o',label='Selected CI')
    plt.errorbar(num_configs_list,len(num_configs_list)*[compared_energy], yerr = CHEMICAL_ACCURACY,ecolor = 'green',label='Full CI',capsize=3)
    plt.legend()
    # text = "exact ci E=" + str(round(ci_ground_energy,8)) + ' vqe best=' + str(round(self.optimal_energy ,8))
    
    filename = ci_mat_dir + '/' + chart_name + '_energy_vs_configs.png'
    plt.savefig(filename)
    plt.show()


if __name__ == "__main__":
    # ci_mat_dir = 'benzene_nist_opt_1024/CI_matrices'
    # cisdt_energy = -228.25005814809177 - 204.8096752115309
    # chart_name = 'benzene'
    # create_chart(ci_mat_dir,cisdt_energy,chart_name)

    # ci_mat_dir = 'H2O_dev/CI_matrices'
    # fci_energy = -75.01156040647366 - 9.215017810268293
    # chart_name = 'H2O'
    # create_chart(ci_mat_dir,fci_energy,chart_name)


    ci_mat_dir = 'BeH2_1.3/CI_matrices'
    fci_energy = -15.595047080804692 - 3.4600048389961535
    chart_name = 'BeH2_1.3'
    create_chart(ci_mat_dir,fci_energy,chart_name)
