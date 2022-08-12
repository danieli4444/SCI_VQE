from matplotlib import pyplot as plt

"""
to understand the minimal basis meaning please see -
https://documents.epfl.ch/users/m/mk/mkilic/public/chapter_3.pdf


or hydrogen (and helium) this means a single s-function. For the first row
in the periodic table it means two s-functions (1s and 2s) and one set of
p-functions (2p x , 2p y and 2p z ). Lithium and beryllium formally only require
two s-functions, but a set of p-functions is usually also added. For the second
row elements, three s-functions (1s, 2s and 3s) and two sets of p-functions
(2p and 3p) are used

"""
"""

addressed molecules:
H2 ,   LiH,  BeH2,       H2O,        NH3,       benzene (C6H6)

numbe of electrons:

2  ,   4  ,  6   ,       10 ,        10,        36

the following minimal basis spin-orbital number with JW encoding without freezing orbitals (or number of qubits):
2*2  , 2*(5 + 1) , 2*(5 + 2) , 2*(5 + 2) , 2*(5 + 3), 2*(5*6 + 6)
4    , 12  , 14        , 14        , 16       , 72


total number of determinants will be (M/2 choose N/2)^2:
4,   , 225 , 1225      , 441        ,3136     , 8.2* 10^19 (and cisdt would be 1.5 *10^7)


the QEE needed qubits log(#determinants) without freezing orbitals:
2    , 8   , 11        , 9          ,12       , 67

CI QEE needed qubits (for achieving chemical accuracy???):
1    , 4   , 5       , 5          , 6        , (cisdt 14)


"""

electrons = [2,4,6,10,10,36]
classic_encoding_qubits = [4, 8, 12, 14, 16, 72]
QEE_qubits = [2    , 6   , 9        , 9          ,12       , 67]
our_qubits = [1, 4, 5, 5, 6, 15]
xlabels = ['H2','LiH','BeH2','H2O','NH3','C6H6(benzene)']

plt.scatter(electrons,classic_encoding_qubits, label="classic_encoding_qubits",color='red')
plt.scatter(electrons,QEE_qubits, label="QEE_qubits",color='blue')
plt.scatter(electrons,our_qubits, label="our_qubits",color='green')
plt.xticks(electrons,xlabels)
plt.legend()
plt.tight_layout()
plt.savefig('q_vs_e.png')
plt.show()



