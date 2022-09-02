"""
this module receives given qubitOp and runs a VQE hardware efficient algorithm to find the lowest 
EigenValue and EigenState

"""

from fileinput import filename
from qiskit.aqua.algorithms import VQE
from qiskit import Aer
from qiskit.circuit.library import RealAmplitudes,TwoLocal
from qiskit.aqua import QuantumInstance
from qiskit.aqua.components.optimizers import COBYLA,SPSA,SLSQP,L_BFGS_B
from qiskit.providers.aer.noise import NoiseModel
from qiskit.aqua.operators.legacy.op_converter import to_tpb_grouped_weighted_pauli_operator,WeightedPauliOperator

from qiskit import IBMQ
from qiskit.providers.ibmq import least_busy

from CIMatrixGenerator import CI_RESULT_FILE
from MatrixEncoder import MatrixEncoder
from ExactSolver import diagonalize_qubitOp_Numpy
from VQEConfig import *

import numpy as np
import json
import os
import sys
from time import time,ctime
from matplotlib import pyplot as plt

class VQERunner:
    def __init__(self,vqeconfig) -> VQEConfig:
        
        # CI matrix:
        # self.ci_matrix_filepath = ci_matrix_filepath
        # self.encoding_type = encoding_type
        # repulsion energy for final energy output
        # self.shift = shift

        self._vqe_config = vqeconfig
        self._vqe_config.toJson()
        self.ci_matrix = np.loadtxt(self._vqe_config.ci_matrix_filepath,dtype=complex)
        # Matrix Encoder:
        matrix_encoder = MatrixEncoder(self.ci_matrix)
        
        self.qubitOp = matrix_encoder.encoded_matrix_operator
        self.raw_matrix_sparsity = matrix_encoder.raw_matrix_sparsity
        # compute the TPB grouped form for analyzing the minimum number of measurements possible:
        from qiskit.aqua.operators.legacy import TPBGroupedWeightedPauliOperator
        TP = TPBGroupedWeightedPauliOperator(self.qubitOp.paulis,self.qubitOp.basis)
        x = TP.sorted_grouping(self.qubitOp)
        self.TPB_grouped_form = x.__str__()

        if self._vqe_config.encoding_type == 'efficient':
            self.paulis_list_len = len(self.qubitOp.paulis)
            self.paulis_list = self.qubitOp.paulis
        elif self._vqe_config.encoding_type == 'TPB':
            raise("TPB grouping Not supported in VQE...")

            
            # H_Operator = to_tpb_grouped_weighted_pauli_operator(H_MatrixOperator, 
                                                    # TPBGroupedWeightedPauliOperator.sorted_grouping)
            # print(H_Operator.print_details())
            # self.paulis_list = self.qubitOp.paulis

        self.paulis_list_len = len(self.qubitOp.paulis)
        self.num_qubits = matrix_encoder.num_qubits

        # exact diagonalization results:
        self.exact_ground_energy, self.exact_ground_state = self.exact_diag() 
        self.exact_ground_energy += self._vqe_config.shift
        # qubitOp exact diagonalization results:
        self.exact_qubitOp_energy, self.exact_qubitOp_g_state = self.exact_qubitOp_diag()

       

        # VQE params:
        # self.ansatz_layers = 0
        # self.entanglement_type = ''
        self.ansatz_ops = ''
        self.ansatz_layers = self._vqe_config.ansatz_layers
        self.ansatz_depth = 0
        self.ansatz_initial_state = 0
        self.optimizer_name = ''
        # self.max_iterations = 0
        # self.shots = 0
        self.backend_properties = ''
        self.backend_config = ''
        self.noise_model = ''
        # self.use_noise = True
        self.result = 0
        self.optimal_energy = 0
        self.optimal_energy_std = 0
        self.minimal_energy = 0
        self.averaged_energy = 0 # the average energy of last 10 vqe iterations 
        self.averaged_energy_std = 0
        self.optimal_state = 0
        self.vqe_result_params = []
        self.actual_iterations = 0
        self.number_params = 0

       
        self.start_time = 0
        self.end_time = 0
        self.total_runtime = 0
        self.error = 0 # in terms of chemical accuracy
        self.chemical_accuracy = 0.0016
        
    @property
    def vqe_config(self):
        return self._vqe_config

    def toJson(self):

        filename = self._vqe_config.output_dir + '/' + 'VQE_output.json'
        new_dict = {}
        for key in self.__dict__.keys():
            if key== "vqe_result_params":
                new_dict[key] = self.__dict__[key]
            else:
                new_dict[key] = str(self.__dict__[key])
        with open(filename,'w') as f:
             json.dump(new_dict,f,indent=4)


    def generate_plots(self):
        plots_data_output_file = self._vqe_config.output_dir + '/energy_plot_data.txt'
        counts = self.vqe_result_params['counts']
        energies = self.vqe_result_params['values']    
        energies = list(map(lambda x:x+self._vqe_config.shift, energies))
        values_std = self.vqe_result_params['values_std'] 
        ci_ground_energy = self.exact_ground_energy 

        with open(plots_data_output_file,'w') as f:
            f.write("counts:")
            f.write(counts.__str__())
            f.write("\n energies:")
            f.write(energies.__str__())
            f.write("\n std_values:")
            f.write(values_std.__str__())
            f.write("\n ci_ground energy:")
            f.write(ci_ground_energy.__str__())
            
        # ignore_the_first_vals
        ig = 30
        if len(counts) < ig:
            ig = len(counts)

        error_bar = values_std
        exact = len(counts) * [0]

        plt.figure()
        plt.errorbar(counts,exact,yerr=self.chemical_accuracy,label='exact',linestyle="dashed",color='green',ecolor = 'black',capsize=3)
        # plt.plot(num_configs,exact,label='exact',linestyle="",marker="o")
        # plt.plot(num_configs,fci,label='fci',linestyle="dotted")
        # plt.xticks(counts)
        plt.xlabel("vqe iteration")

        plt.ylabel("Energy [Ha]")
        plt.errorbar(counts,energies-ci_ground_energy, yerr = error_bar,linestyle="dotted",ecolor = 'red',label='std',capsize=3)
        
        text = "exact ci E=" + str(round(ci_ground_energy,8)) + ' vqe best=' + str(round(self.optimal_energy ,8))
        plt.title('vqe energy error (e - exact) |'+ text)
        plt.legend()

        filename = self._vqe_config.output_dir + '/energy_vs_iter.png'
        plt.savefig(filename)
        # plt.show()

   


    def exact_diag(self):
        # eigenvals, eigenstates = np.linalg.eigh(self.ci_matrix)
        import scipy
        eigenvals, eigenstates = scipy.linalg.eig(self.ci_matrix)
        idx = eigenvals.argsort()
        eigenvals = eigenvals[idx]
        eigenstates = eigenstates[:,idx]
        print(eigenvals)
        return np.real(eigenvals[0]),eigenstates[0]

    def exact_qubitOp_diag(self):
        (ground_energy,ground_state) = diagonalize_qubitOp_Numpy(self.qubitOp)
        return ground_energy,ground_state


    def run_VQE_qasm(self):
        """ Run VQE 
        """


        ##### hyper - params initialization:

        qubitOp = self.qubitOp
        num_qubits = self.num_qubits
        layers = self._vqe_config.ansatz_layers 
        simulator = self._vqe_config.simulator
        
        
        print("\nstarting VQE!\n")
        print("num_qubits:",num_qubits)
        print("num_layers: ",layers)
        print("max iterations:",self._vqe_config.max_iterations )

        # set backend:
        def get_job_status(job_id, job_status, queue_position, job):
            # print("job_id")
            # print(job_id)
            # print("job_status")
            # print(job_status)
            # print("queue_position")
            # print(queue_position)
            # print(job)
            pass


        from qiskit.ignis.mitigation import CompleteMeasFitter

        if self._vqe_config.use_noise:
            from qiskit.test.mock import FakeLima,FakeSantiago
            
            if self._vqe_config.backend_type =='IBMQ': # in case we want to run on real hardware:
                
                IBMQ.save_account("81ab8c4babe3018d036e12c049e5a00461be47586925fd6d322e7348e338a5d1896beec9282e033c7bfcc2899da1832ac8a75ce080b558e1d55f2c39cfcddb1f",
                    hub='ibm-q-research-2',group='bar-ilan-uni-1', project='main')
                provider = IBMQ.load_account()
                print(provider.backends())
                # z = IBMQ.providers() 
                # provider2 = IBMQ.get_provider(hub='ibm-q-research-2',group='bar-ilan-uni-1', project='main')
                # y = provider.backends(simulator=False, operational=True)

                backend_name = self._vqe_config.backend_name
                # backend = provider2.get_backend(backend_name)
                backend = provider.get_backend(backend_name)
                self.backend_properties = backend.properties()
                self.backend_config = backend.configuration().__dict__
                self.noise_model = simulator + backend_name
                noise_model = None


            elif self._vqe_config.backend_type =='simulator':

                backend_name = self._vqe_config.backend_name

                if backend_name == 'ibmq_santiago':
                    fakebackend = FakeSantiago()                
                # if backend_name == 'ibmq_lima':
                #     fakebackend = FakeLima()
                else:
                    provider = IBMQ.load_account()
                    print(provider.backends())
                    fakebackend = provider.get_backend(backend_name)
                # backend = provider.get_backend(backend_name)
                noise_model = NoiseModel.from_backend(fakebackend)

                backend = Aer.get_backend(simulator)
                self.backend_properties = backend.properties()
                self.backend_config = backend.configuration().__dict__
                self.noise_model = simulator + backend_name

            if self._vqe_config.use_error_mit:
                qi = QuantumInstance(backend=backend,shots=self._vqe_config.shots,job_callback=get_job_status,
                noise_model=noise_model, measurement_error_mitigation_cls=CompleteMeasFitter,
                cals_matrix_refresh_period=self._vqe_config.cals_matrix_refresh_period)
            else:
                qi = QuantumInstance(backend=backend,shots=self._vqe_config.shots,job_callback=get_job_status,
                    noise_model=noise_model)
                    

        else:
            backend = Aer.get_backend(simulator)
            self.backend_properties = backend.properties()
            self.backend_config = backend.configuration().__dict__
            self.noise_model = ''
            qi = QuantumInstance(backend=backend,shots=self._vqe_config.shots,job_callback=get_job_status)
 
        # set per-vqe-iteration info for further analysis:
        counts = []
        values = []
        values_std = []
        parameters_list = []
        def store_intermediate_result(eval_count, parameters, mean, std):
            print("\n")
            print("Finished VQE iteration {0}".format(eval_count))
            print("energy = ",self._vqe_config.shift + mean)
            print("std=",std)
            print("\n")

            counts.append(eval_count)
            values.append(mean) 
            values_std.append(std)
            parameters_list.append(parameters.tolist())

        # set a hardware efficient Ansatz
        ansatz = RealAmplitudes(qubitOp.num_qubits, entanglement=self._vqe_config.entanglement_type, reps=layers)
        # ansatz =TwoLocal(qubitOp.num_qubits,['ry','rz'], 'cx',entanglement=self._vqe_config.entanglement_type,reps=layers)
        self.ansatz_initial_state = ansatz_initial_point = [0.]*ansatz.num_parameters
        # self.ansatz_draw_file = self._vqe_config.output_dir + "/ansatz_circuit.txt"
        # 
        x = ansatz.decompose()
        self.ansatz_ops = x.count_ops()
        self.ansatz_depth = x.depth()
        x.draw(output='mpl', filename=self._vqe_config.output_dir + '/ansatz_circuit.png')

        # set an optimizer       
        optimizer = COBYLA(maxiter=self._vqe_config.max_iterations)
        # optimizer = SPSA(maxiter=self._vqe_config.max_iterations)
        # optimizer = SLSQP(maxiter=self._vqe_config.max_iterations)
        # optimizer = L_BFGS_B(maxiter=self._vqe_config.max_iterations)

        # self.optimizer_name = cobyla.print_options()
        print(self.optimizer_name)

        # create a VQE algorithm obj and run
        vqe = VQE(var_form=ansatz, optimizer=optimizer, callback=store_intermediate_result, 
                quantum_instance=qi, initial_point=ansatz_initial_point)
        
        t1 = time()
        self.start_time = ctime()
        result = vqe.compute_minimum_eigenvalue(qubitOp)
        self.end_time = ctime()
        self.total_runtime = round(time() - t1,3)

        # save results:
        self.result = result.__dict__
        self.optimal_state = result.optimal_point
        # self.optimal_energy = np.real(result['eigenvalue'])+ shift
        self.optimal_energy = vqe.get_optimal_cost() + self._vqe_config.shift
        self.optimal_energy_std = values_std[-1]
        self.minimal_energy = min(values) + self._vqe_config.shift
        # self.optimal_state = vqe.get_optimal_vector()
        # self.optimal_circuit = vqe.get_optimal_circuit()
        f = self._vqe_config.output_dir + '/' + 'optimal_circuit.png'
        optimal_circuit = vqe.get_optimal_circuit()
        optimal_circuit = optimal_circuit.decompose()
        optimal_circuit.draw("mpl",filename=f)
        print("VQE Result:\n")
        print("Final Energy = ",self.optimal_energy)
        print("minimal Energy = ",self.minimal_energy)
        print("the exact ci energy = ",self.exact_ground_energy)
        print("Optimal state measurement:")
        print(self.optimal_state)
        self.number_params = len(parameters_list[0])
        print("Number of optimization params:",len(parameters_list[0]))
        print("total runtime = ",self.total_runtime)
        self.vqe_result_params = {'counts':counts,'values':values,
            'values_std':values_std,'last_param_list':parameters_list[-1]}

        self.actual_iterations = counts[-1]
        if self.actual_iterations > 14:
            x = values[-11:-1]
            x = [(y + self._vqe_config.shift) for y in x]
            self.averaged_energy = np.average(x)
            self.averaged_energy_std = np.std(x)
            print("averaged energy = ",self.averaged_energy)
            
        Error = np.abs(self.exact_ground_energy - self.optimal_energy)
        accu_ratio = round(Error/self.chemical_accuracy,3)
        
        self.error = accu_ratio
        print("the Error = {0} Hartree or {1} X chemical accuracy \n\n\n".format(round(Error,8),accu_ratio))

        self.toJson()
        self.generate_plots()
        print("Finished!")
        return self.averaged_energy


    #     IBMQ.save_account("81ab8c4babe3018d036e12c049e5a00461be47586925fd6d322e7348e338a5d1896beec9282e033c7bfcc2899da1832ac8a75ce080b558e1d55f2c39cfcddb1f",
    #         hub='ibm-q-research-2',group='bar-ilan-uni-1', project='main')

    #     IBMQ.load_account()



def run_vqe_qasm():

    vqeconfig = VQEConfig(
            mol_dir='NH3/',
            ci_matrix_path='NH3/CI_matrices/NH3_cimat__4.out',
            output_name='config64',
            backend_type='simulator',
            backend_name='ibmq_santiago',
            encoding_type='efficient',
            simulator='statevector_simulator',
            ansatz_layers=1,
            entanglement_type='circular',
            max_iterations=100,
            shots=100000,
            use_noise=False,
            use_error_mit=False,
            cals_matrix_refresh_period=30
            
    )

    # real hardware config:
    # vqeconfig = VQEConfig(
    #         mol_dir='H2',
    #         ci_matrix_path='H2/CI_matrices/H2_cimat__2.out',
    #         backend_type='IBMQ',
    #         backend_name='ibmq_lima',
    #         encoding_type='efficient',
    #         simulator='',
    #         ansatz_layers=0,
    #         entanglement_type='circular',
    #         max_iterations=30,
    #         shots=20000,
    #         use_noise=True,
    #         use_error_mit=True,
    #         cals_matrix_refresh_period=30
            
    # )

    new_runner = VQERunner(vqeconfig)
    e = new_runner.run_VQE_qasm()
    return e


if __name__ == "__main__":
    vqe_energy_list = []
    for i in range(1,2):
        energy = run_vqe_qasm()
        vqe_energy_list.append(energy)
    print(vqe_energy_list)
    print(np.average(vqe_energy_list))
    print(np.std(vqe_energy_list))
    
    # run_vqe_ibmq()



    