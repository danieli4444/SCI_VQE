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
H2 ,   LiH       , BeH2,       H2O,        NH3,       C2H4

numbe of electrons:

2  ,   4         ,  6   ,       10 ,        10,        12 + 4 = 16

the following minimal basis spin-orbital number with JW encoding without freezing orbitals (or number of qubits):
2*2  , 2*(5 + 1) , 2*(5 + 2) , 2*(5 + 2) , 2*(5 + 3), 2*(5*2 + 4)
4    , 12        , 14        , 14        , 16       , 28


total number of determinants will be (M/2 choose N/2)^2:
4,   , 225      , 1225      , 441        ,3136           , 9,018,009


the QEE needed qubits log(#determinants) without freezing orbitals - ceil[log2(#dets)]:
2    , 8        , 11        , 9          ,12                     , 24


CI QEE needed qubits (for achieving chemical accuracy???):
1    , 3        , 4       , 5          , 6        , 12


"""

electrons = [2,4,6,10,10,16]
classic_encoding_qubits = [4, 12, 14, 14, 16, 28]
QEE_qubits = [2    , 8   , 11        , 9          ,12       , 24]
our_qubits = [1, 3, 4, 5, 6, 12]
xlabels = ['H2','LiH','BeH2','NH3/H2O','H2O/NH3','C2H4(ethylene)']

# plt.scatter(electrons,classic_encoding_qubits, label="Standard encoding",color='red')
# plt.scatter(electrons,QEE_qubits, label="QEE_qubits",color='blue')
# plt.scatter(electrons,our_qubits, label="Proposed encoding",color='green')
plt.plot(electrons,our_qubits,'-o', label="Proposed encoding",color='green')
plt.plot(electrons,classic_encoding_qubits,'-o', label="Standard encoding",color='red')
plt.xticks(electrons,xlabels)
plt.legend()
plt.tight_layout()
plt.savefig('q_vs_e_withoutQEE.png')
plt.show()


# electrons = [2,4,6,10,16]
# classic_encoding_qubits = [4, 12, 14, 14, 28]
# QEE_qubits = [2    , 8   , 11        , 9             , 24]
# our_qubits = [1, 3, 4, 5, 12]
# xlabels = ['H2','LiH','BeH2','H2O','C2H4(ethylene)']

# plt.plot(electrons,classic_encoding_qubits, label="Standard encoding",color='red')
# # plt.plot(electrons,QEE_qubits, label="QEE_qubits",color='blue')
# plt.plot(electrons,our_qubits, label="Proposed encoding",color='green')
# plt.xticks(electrons,xlabels)
# plt.legend()
# plt.tight_layout()
# plt.savefig('q_vs_e_withoutNH3_withoutQEE.png')
# plt.show()


