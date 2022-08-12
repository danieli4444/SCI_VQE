from qiskit.aqua.algorithms import NumPyEigensolver
import numpy as np

def diagonalize_qubitOp_Numpy(qubitOp):
    """return lowest eigenvalue and eigenstate

    Args:
        qubitOp (type - qiskit.aqua.operators.WeightedPauliOperator): a weighted pauli operator
    """
    exact_result = NumPyEigensolver(qubitOp).run()
    exact_energy = np.real(exact_result.eigenvalues)[0]
    return (exact_energy,exact_result.eigenstates[0])

