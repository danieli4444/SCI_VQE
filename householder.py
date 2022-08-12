## module householder
''' d,c = householder(a).
    Householder similarity transformation of matrix [a] to 
    the tridiagonal form [c\d\c].
    p = computeP(a).
    Computes the acccumulated transformation matrix [p]
    after calling householder(a).
'''    
from numpy import dot,diagonal,outer,identity
from math import sqrt
import numpy as np

def householder(a): 
    n = len(a)
    for k in range(n-2):
        u = a[k+1:n,k]
        uMag = sqrt(dot(u,u))
        if u[0] < 0.0: uMag = -uMag
        u[0] = u[0] + uMag
        h = dot(u,u)/2.0
        v = dot(a[k+1:n,k+1:n],u)/h
        g = dot(u,v)/(2.0*h)
        v = v - g*u
        a[k+1:n,k+1:n] = a[k+1:n,k+1:n] - outer(v,u) \
                         - outer(u,v)
        a[k,k+1] = -uMag
    return diagonal(a),diagonal(a,1)

def computeP(a): 
    n = len(a)
    p = identity(n)*1.0
    for k in range(n-2):
        u = a[k+1:n,k]
        h = dot(u,u)/2.0
        v = dot(p[1:n,k+1:n],u)/h           
        p[1:n,k+1:n] = p[1:n,k+1:n] - outer(v,u)
    return p

def get_tridiagonal_mat(M):
    d,c = householder(M)
    newM = np.diag(d) + np.diag(c, 1) + np.diag(c, -1)
    return newM

if __name__ == '__main__':
    ci_matrix_filepath = "H2O_dev/CI_matrices/H2O_dev_cimat__8.out"
    M = np.loadtxt(ci_matrix_filepath)
    eigenvals, eigenstates = np.linalg.eigh(M)
    print(np.real(eigenvals[0]))
    print('\n')
    d,c = householder(M)
    newM = np.diag(d) + np.diag(c, 1) + np.diag(c, -1)

    eigenvals2, eigenstates = np.linalg.eigh(newM)
    print(np.real(eigenvals2[0]))
