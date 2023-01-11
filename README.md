the project goal is to find the ground state energy of a molecule using log2(#qubits) in comparison to standard VQE implementations.
This is achieved using:
* Building CI matrix with pre-optimized set of configurations from fci computation (performed via psi4 package). - CIMatrixGenerator.py.
* Enconding the matrix to a qubit efficient qubit operator - MatrixEncoder.py.
* finding the ground state energy using hardware-efficient VQE - VQERunner.py

important notes:
- the code assumes low qubit number and therefore relatively small ci matrices - O(64X64). 
  for larger size the code needs to be optimized from memory and computation perspectives.



prerequisites:
* run the code with psi4 conda python intepreter (needed to be set manualy in the VScode editor) 
    psi4 conds installation - https://psicode.org/psi4manual/master/conda.html (used Psi4 1.5 Conda version)
* install qiskit old version (0.18.3 for example) containing qiskit aqua and qiskit chemistry 
