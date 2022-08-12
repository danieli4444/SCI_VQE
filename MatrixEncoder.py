import numpy as np
from qiskit.aqua.operators import WeightedPauliOperator
from qiskit.quantum_info import Pauli
from itertools import product
from functools import reduce
from time import time


GRAY_CODE = {'0':'000','1':'100','2':'110','3':'010','4':'011','5':'111','6':'101','7':'001'}

# definition of the entry operators 0,1,+,- :
hopping_in_entry = {('0', '0'): '0',
                ('1', '1'): '1',
                ('0', '1'): '+',
                ('1', '0'): '-'} ## (from, to)


def binary(num, pre='0b', length=8, spacer=0):
    """ converts int to binary according to given _length_ param 
    """
    return '{0}{{:{1}>{2}}}'.format(pre, spacer, length).format(bin(num)[2:])

def entry_to_pauli(entry):
    """   
    Args:
        entry (_type_): _description_

    Returns:
        _type_: _description_
    """

    if entry == '+':
        return ('X','-iY')
    if entry == '-':
        return ('X','iY')
    if entry == '0':
        return ('I','Z')
    if entry == '1':
        return ('I','-Z')

def convert_tuple_element(pauli_element):
    if pauli_element == 'X':
        return ('X',1)
    if pauli_element == 'I':
        return ('I',1)
    if pauli_element == 'iY':
        return ('Y',1j)
    if pauli_element == 'Z':
        return ('Z',1)
    if pauli_element == '-iY':
        return ('Y',-1j)
    if pauli_element == '-Z':
        return ('Z',-1)


def compute_pauli_and_weight(t1,t2):
    (p1,w1) = t1
    (p2,w2) = t2
    t3 = ((f'{p1}{p2}'),w1*w2)
    return t3

def pauli_tuple_processing(pauli_tuple):
    result = map(convert_tuple_element,pauli_tuple)
    
    (new_pauli_str,pauli_weight) = reduce(compute_pauli_and_weight,list(result))

    Pau = Pauli(new_pauli_str)
    pauli_weight = pauli_weight /(2**len(pauli_tuple))
    return [pauli_weight,Pau]


def entryoperator_to_weightedpauli(operator,num_qubits):
    """_summary_

    Args:
        operator (_type_): a '+-01' operator
        num_qubits (_type_): 

    Returns:
        list: 2^num_qubits pauli strings with coefficients
    """
    # produce all combinations
    pauli_strings_list = [entry_to_pauli(entry) for entry in operator]
    # print(pauli_strings_list)
    # x = list(product(*pauli_strings_list))
    return WeightedPauliOperator(list(map(pauli_tuple_processing,list(product(*pauli_strings_list)))))


class MatrixEncoder:
    def __init__(self,raw_matrix,encoding_type='efficient') -> None:
        """ initiate encoded matrix 

        Args:
            raw_matrix (_type_): _description_
        """
        self.raw_matrix = raw_matrix
        self.raw_matrix_dims = 0 
        self.raw_matrix_sparsity = 0
        self.mapping = 0 # a dictionary that maps each raw matrix |i><j| to entry operator
        self.num_qubits = 0
        if encoding_type == 'efficient':
            self.encoded_matrix_operator = self.encode_matrix()


    def get_raw_matrix(self):
        return self.raw_matrix

    def get_general_info(self):
        general_info = {
            "raw_matrix_dims" : self.raw_matrix_dims,
            "raw_matrix_sparsity": self.raw_matrix_sparsity,
            "num_qubits" : self.num_qubits
        }
        print("\n")
        print("general info:")
        print(general_info)
        print("\n")
        return general_info

    def get_mapping(self):
        return self.mapping

    def get_encoded_qubitOp(self):
        return self.encoded_matrix_operator

    def clean_zero_paulis(self,puali_list):
        THRESHOLD_WEIGHT = 1e-05
        for idx,p in enumerate(puali_list):
            if np.abs(p[0])< THRESHOLD_WEIGHT:
                print(p)
                del puali_list[idx]
        return puali_list

    def generate_gray_code(self):
        n=self.num_qubits
        gray_code = {}
        for i in range(0, 1<<n):
            gray=i^(i>>1)
            gray_code[str(i)]="{0:0{1}b}".format(gray,n)
        self.gray_code = gray_code

    def encode_matrix(self):
        """ encodes a CI Hamiltonian Matrix to a weighted pauli operator with log2(N) qubits. 

        Args:
            M (ndarray): a NxN Configuration Interaction Hamiltonian Matrix that needed to be encoded and later diagonalized.
        """

        print("starting encoding matrix to qubitOp...")
        t1 = time()
        # we need a mapping that for each |i><j| in the matrix encodes it to an entry operator string  '0110-+...'
        mapping = {}
        qubitOp = WeightedPauliOperator(paulis=[])

        M = self.raw_matrix
        (N,N) = M.shape
        self.raw_matrix_dims = (N,N)

        #notice that this limits us to only to matrices that are NxN , N = 2^m
        num_qubits = int(np.log2(N))
        assert 2**num_qubits == N
        self.num_qubits = num_qubits

        num_excluded = 0 # to count the sparsity of the matrixc excluded M(i,j) == 0

        PERFORM_GRAY_ENCODING = False
        if PERFORM_GRAY_ENCODING:
            self.generate_gray_code()
        for i in range(N):
            for j in range(N):
                mapping[(i,j)] = []
                w = M[i,j] # operator weight
                if w == 0:
                    mapping[(i,j)].append(0)
                    num_excluded +=1
                else:
                    entry_op = ''
                    if PERFORM_GRAY_ENCODING:
                        encoded_i = self.gray_code[str(i)]
                        encoded_j = self.gray_code[str(j)]
                    else:
                        encoded_i = binary(i,pre='',length=num_qubits)
                        encoded_j = binary(j,pre='',length=num_qubits)
                    # build entry operator string:
                    # for q0, q in zip(encoded_j, encoded_i):
                    #     entry_op += hopping_in_entry[(q0, q)]
                    
                    entry_op_list = [hopping_in_entry[(q0,q1)] for (q0,q1) in zip(encoded_j,encoded_i)]
                    entry_op = "".join(entry_op_list)
                    mapping[(i,j)].append(entry_op)
                    # newOp = operator2WeightedPauliOperator(operator=entry_op,num_qubits=num_qubits)
                    newOp = entryoperator_to_weightedpauli(entry_op,num_qubits)
                    qubitOp+= w * newOp
                    

        # chop insiginificant paulis by Threshold
        THRESHOLD_PAULI_WEIGHT = 1e-06
        qubitOp.chop(THRESHOLD_PAULI_WEIGHT)
        self.raw_matrix_sparsity = num_excluded
        print("sparsity = ",num_excluded)
        self.mapping = mapping
        print("finished matrix encoding to qubitOp!")
        total_runtime = round(time() - t1,3)
        print("took {0} seconds".format(total_runtime))
        print("number of total paulis = ",len(qubitOp.paulis))
        return qubitOp



    def naive_encoding(self):
        # we need a mapping that for each |i><j| in the matrix encodes it to an entry operator string  '0110-+...'
            qubitOp = WeightedPauliOperator(paulis=[])

            M = self.raw_matrix
            (N,N) = M.shape
            self.raw_matrix_dims = (N,N)
            from qiskit.aqua.operators.primitive_ops import MatrixOp

            p = MatrixOp(M)
            qubitOp = p.to_pauli_op()
            #notice that this limits us to only to matrices that are NxN , N = 2^m
            num_qubits = qubitOp.num_qubits

            self.num_qubits = num_qubits

            return qubitOp
    
def compact_encoding():
    """ 
    taken from  https://arxiv.org/pdf/2111.00627.pdf
    """
    





def test_H2_operator():
    print("\n encoding an H2 minimal CI matrix with HF state and Double Excitation state\n")

    a = -1.8266
    b = 0.1814
    c = 0.1814
    d = -0.2596
    
    # run a simple version of the encoding:
    
    # M = [a,b,c,d]
    # entries = ['0','+','-','1']
    # qubitOp = WeightedPauliOperator(paulis=[])
    # for w,e in zip(M,entries):
    #     qubitOp+= w * operator2WeightedPauliOperator(e, 1, 1)


    M = np.array([[ a,b],[ c,d]])
    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    
    shift = 0.7
    M_lowest_eigenvalue = -1.8
    from qiskit.aqua.algorithms import NumPyEigensolver
    exact_result = NumPyEigensolver(qubitOp).run()
    exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    print("computed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0]))
    print("lowest eigenvalue: ", M_lowest_eigenvalue)


def test_4x4_operator():
    print(" encoding arbitrary CI matrix from water molecule configs")

    # arbitrary CI matrix from water molecule configs
    M = np.array([[-8.28558597e+01 ,2.43717123e-02 ,  0.00000000e+00  ,.00000000e+00],
        [ 2.43717123e-02 , -8.41774586e+01 , -2.82342066e-11 ,-2.49933407e-12],
        [ 0.00000000e+00, -2.82342066e-11, -8.35992573e+01 , 2.15965646e-03],
        [ 0.00000000e+00 ,-2.49933407e-12 , 2.15965646e-03 , -8.28135727e+01]])
    M_lowest_eigenvalue = -84.1779078488533
    shift = 9.2

    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    from qiskit.aqua.algorithms import NumPyEigensolver
    exact_result = NumPyEigensolver(qubitOp).run()
    exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    print("\ncomputed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0]))
    print("lowest eigenvalue: ", M_lowest_eigenvalue)
    print("\n")


def test_naive_encoding():
    a = -1.8266
    b = 0.1814
    c = 0.1814
    d = -0.2596
    M = np.array([[ a,b],[ c,d]])
    M = np.array([[-8.28558597e+01 ,2.43717123e-02 ,  0.00000000e+00  ,.00000000e+00],
        [ 2.43717123e-02 , -8.41774586e+01 , -2.82342066e-11 ,-2.49933407e-12],
        [ 0.00000000e+00, -2.82342066e-11, -8.35992573e+01 , 2.15965646e-03],
        [ 0.00000000e+00 ,-2.49933407e-12 , 2.15965646e-03 , -8.28135727e+01]])
    # from qiskit.quantum_info.operators import Operator, Pauli
    # qubitOp = naive_encoding(M)
    # print(qubitOp)
    from qiskit.aqua.operators.primitive_ops import MatrixOp
    from qiskit.quantum_info.operators.symplectic.sparse_pauli_op import SparsePauliOp
    from qiskit.aqua.algorithms import NumPyEigensolver

    M_op = MatrixOp(M)
    z = M_op.to_pauli_op()

    r = NumPyEigensolver(z).run()
    shift = 9.2
    exact = np.real(r.eigenvalues[0])+ shift
    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    
    r2 = NumPyEigensolver(qubitOp).run()
    exact2 =  np.real(r2.eigenvalues[0])+ shift
    print(z.num_qubits)
    print(z)
    print(exact)
    
    print(qubitOp)
    print(qubitOp.paulis)
    print(exact2)
    
    # s_op = SparsePauliOp.from_operator(M)
    # z = M_op.to_pauli_op()

    # shift = 0.7
    # M_lowest_eigenvalue = -1.8
    # from qiskit.aqua.algorithms import NumPyEigensolver
    # exact_result = NumPyEigensolver(qubitOp).run()
    # exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    # print("computed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0]))
    # print("lowest eigenvalue: ", M_lowest_eigenvalue)


def test_pauli_grouping():
    from qiskit.aqua.operators.legacy import TPBGroupedWeightedPauliOperator
    # M = np.array([[-8.28558597e+01 ,2.43717123e-02 ,  0.00000000e+00  ,.00000000e+00],
    #     [ 2.43717123e-02 , -8.41774586e+01 , -2.82342066e-11 ,-2.49933407e-12],
    #     [ 0.00000000e+00, -2.82342066e-11, -8.35992573e+01 , 2.15965646e-03],
    #     [ 0.00000000e+00 ,-2.49933407e-12 , 2.15965646e-03 , -8.28135727e+01]])
    # M_lowest_eigenvalue = -84.1779078488533
    ci_matrix_filepath = "_O_H10.955_H10.9552105_symmetryc1_/h2o_ci_4.out"
    M = np.loadtxt(ci_matrix_filepath)
    shift = 9.2

    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    TP = TPBGroupedWeightedPauliOperator(qubitOp.paulis,qubitOp.basis)
    print(TP)
    x =TP.sorted_grouping(qubitOp)
    print(x)
    print(x.print_details())
    print("\n")
    # print(qubitOp.paulis)
    from qiskit.aqua.algorithms import NumPyEigensolver
    exact_result = NumPyEigensolver(TP).run()
    exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    print("\ncomputed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0]))
    # print("lowest eigenvalue: ", M_lowest_eigenvalue)
    print("\n")
    # print(len(TP.paulis))
    # print("\n")
    # print(len(qubitOp.paulis))

def test_H2O_operator():
    from qiskit.aqua.operators.legacy import TPBGroupedWeightedPauliOperator
    # M = np.array([[-8.28558597e+01 ,2.43717123e-02 ,  0.00000000e+00  ,.00000000e+00],
    #     [ 2.43717123e-02 , -8.41774586e+01 , -2.82342066e-11 ,-2.49933407e-12],
    #     [ 0.00000000e+00, -2.82342066e-11, -8.35992573e+01 , 2.15965646e-03],
    #     [ 0.00000000e+00 ,-2.49933407e-12 , 2.15965646e-03 , -8.28135727e+01]])
    # M_lowest_eigenvalue = -84.1779078488533
    ci_matrix_filepath = "H2O_tridiagonal/CI_matrices/H2O_tridiagonal_cimat__8.out"
    M = np.loadtxt(ci_matrix_filepath)
    shift = 9.215017810268293

    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    eigenvals, eigenstates = np.linalg.eigh(M)
    diagonalization_energy = np.real(eigenvals[0])
    print("\n")
    print(len(qubitOp.paulis))
    for x in qubitOp.paulis:
        print(x)
    
    print("\n")
    from qiskit.aqua.algorithms import NumPyEigensolver
    exact_result = NumPyEigensolver(qubitOp).run()
    exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    print("\ncomputed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0]))
    print("diagonalization_energy = ",diagonalization_energy)
    # print("lowest eigenvalue: ", M_lowest_eigenvalue)
    print("\n")
    


def test_NH3_operator():

    print(" encoding arbitrary CI matrix from ammonia molecule configs")

    # arbitrary CI matrix from water molecule configs
    ci_matrix_filepath = "__N0.00.00.2851841_H-0.46974900.8136292-0.0950614_H-0.4697490-0.8136292-0.0950614_H0.93949800.0-0.0950614_symmetryc1_/NH3_32.out"
    M = np.loadtxt(ci_matrix_filepath)
    shift = 11.9399594194046816 # Nuclear Repulsion Energy 

    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    from qiskit.aqua.algorithms import NumPyEigensolver
    exact_result = NumPyEigensolver(qubitOp).run()
    exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    fci_value = -55.51941431916908 
    error =  exact_energy - fci_value
    print("\ncomputed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0]))
    print("\ncomputed ground state energy (with repulsion): ",exact_energy)
    print("FCI ground_state energy : ",fci_value)
    print("Error: ",error)
    print("\n")

def test_benzene_op():
    print(" encoding arbitrary CI matrix from ammonia molecule configs")

    # arbitrary CI matrix from water molecule configs
    ci_matrix_filepath = "benzene_nist/benzene_nist_cimat__2.out"
    M = np.loadtxt(ci_matrix_filepath)
    shift = 204.8096752115309 # Nuclear Repulsion Energy 

    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    from qiskit.aqua.algorithms import NumPyEigensolver
    exact_result = NumPyEigensolver(qubitOp).run()
    exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    fci_value = -228.25151419789768
    error =  exact_energy - fci_value
    print("\ncomputed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0])) 
    print("\ncomputed ground state energy (with repulsion): ",exact_energy)
    print("FCI ground_state energy : ",fci_value)
    print("Error: ",error)
    print("\n")


def test_LiH_op():
    print(" encoding arbitrary CI matrix from ammonia molecule configs")

    # arbitrary CI matrix from water molecule configs
    ci_matrix_filepath = "LiH_dist/LiH_1.5/CI_matrices/LiH_1.5_cimat__4.out"
    M = np.loadtxt(ci_matrix_filepath)
    shift = 1.0583544213400002 # Nuclear Repulsion Energy 

    encoded_mat = MatrixEncoder(M)
    encoded_mat.get_general_info()
    qubitOp = encoded_mat.get_encoded_qubitOp()
    from qiskit.aqua.algorithms import NumPyEigensolver
    exact_result = NumPyEigensolver(qubitOp).run()
    exact_energy = np.real(exact_result.eigenvalues)[0] + shift
    fci_value = -7.8823622868109116
    error =  exact_energy - fci_value
    print("\ncomputed lowest eigenvalue (without repulsion): ",np.real(exact_result.eigenvalues[0]))
    print("\ncomputed ground state energy (with repulsion): ",exact_energy)
    print("FCI ground_state energy : ",fci_value)
    print("Error: ",error)
    print("\n")

if __name__ == "__main__":  
    test_H2_operator() # minimal basis endcoding of H2
    # test_4x4_operator() 
    # test_naive_encoding()
    # test_pauli_grouping()
    # test_H2O_operator()
    # test_NH3_operator()
    # test_benzene_op()
    # test_LiH_op()
