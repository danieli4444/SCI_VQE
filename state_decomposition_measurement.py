""" 
implementation of decomposition the CI matrix which is Real and symmetrical.
We assume the quantumcircuit doesnt rotate the qubits to complex plain and therefore 
decompose all of the needed state parameters for the energy computation.

implemented on 1 qubits system only

"""

import numpy as np
from qiskit import QuantumCircuit, assemble, Aer
from qiskit.visualization import plot_histogram, plot_bloch_vector
from math import sqrt, pi
from matplotlib import pyplot as plt


# 1 qubit dimension


def get_diagonal(qc,shots):
    sim = Aer.get_backend('qasm_simulator')  # Tell Qiskit how to simulate our circuit

    qc.measure_all() 
    qobj = assemble(qc)
    result = sim.run(qobj,shots=shots).result()
    counts = result.get_counts()
    # qc.draw(output='mpl')
    # plot_histogram(counts_clean)
    # plt.show()
    diag = []
    for key in counts:
        diag.append(counts[key]/shots)
    return diag


def get_pauli_expectation(op,qc,shots):

    from qiskit.aqua import QuantumInstance
    from qiskit.aqua.operators import PauliExpectation, CircuitSampler, StateFn,CircuitStateFn

    psi = CircuitStateFn(qc)
    # define your backend or quantum instance (where you can add settings)
    backend = Aer.get_backend('qasm_simulator') 
    q_instance = QuantumInstance(backend, shots=shots)

    # define the state to sample
    measurable_expression = StateFn(op, is_measurement=True).compose(psi) 

    # convert to expectation value
    expectation = PauliExpectation().convert(measurable_expression)  

    # get state sampler (you can also pass the backend directly)
    sampler = CircuitSampler(q_instance).convert(expectation) 

    expectation_val = sampler.eval().real
    # evaluate
    # print('Sampled:', expectation_val)  

    # qc.draw(output='mpl')
    # plt.show()
    return expectation_val

# create Real Amplitude state
alpha = 0.5
beta = np.sqrt(3)/2
initial_state = [alpha, beta]
print("created state")
print(initial_state)
shots = 1024

qc = QuantumCircuit(1) # Must redefine qc
qc.initialize(initial_state, 0) # Initialize the 0th qubit in the state `initial_state`

print("measuring diagonal")
diag = get_diagonal(qc,shots)
a_2 = diag[0]
b_2 = diag[1]

print("measuring the X part")
from qiskit.aqua.operators import X, Y, Z, I
qc_x = QuantumCircuit(1) # Must redefine qc
qc_x.initialize(initial_state, 0) # Initialize the 0th qubit in the state `initial_state`

print("\n Successfully computed needed state parameters using only 2 measurements instead of 4!")
x_val = get_pauli_expectation(X,qc_x,shots)
print("computed all needed expecation values alpha^2 beta^2 2*alpha*beta:")
print(a_2)
print(b_2)
print(x_val)

print("\n original alpha^2 beta^2 2*alpha*beta are:")
print(alpha**2)
print(beta**2)
print(2*alpha*beta)


print("\n\n")
print("calculating M energy via decomposition:")
print("M = ")
a = -8.417745856086088452e+01
b =  1.522088975946591649e-01
c =  1.522088975946591649e-01
d = -8.233178498775147602e+01


M = np.array([[ a,b],[ c,d]])
print(M)
print("exact energy:")
state = np.array(initial_state)
exact_energy = np.dot(state,np.dot(M,state)) # calc <psi|M|psi>
print(exact_energy)

state_energy = a* a_2 + d*b_2 + b*x_val
print("state energy:")
print(state_energy)
# plt.show()
