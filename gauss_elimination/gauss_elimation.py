import numpy as np

A = np.array([[1, 0, 0, 0, 0],
           [-2, 3, 0, -1, 0],
           [-1, 0, 5, -2, 0],
           [0, -1, -2, 4, -1],
           [0, 0, 0, -1, 5]], float)

b = np.array([5, 24, -24, 0, 24], float)

tol = 1
n = len(b)
x = np.zeros(n, float)

def Gauss(A, x):
    s = x
    er = 0
    for i in range(0, n-1):
        s[i] = np.abs(A[i, i])
        for j in range(2,n-1):
            if np.abs(A[i, j] > s[i]):
                s[i] = np.abs(A[i, j])
    Elimination(A, b)
    if er != -1:
        Substitue(A, b)
    return x

def Elimination(A, b):

    for k in range(0, n-1):
        if np.abs(A[k, k]) < tol:
            Pivot(A, b, k)

        for i in range(k+1, n):
            factor = A[i, k] / A[k, k]
            for j in range(k+1, n):
                A[i, j] = A[i, j] - factor * A[k, j]
            b[i] = b[i] - factor * b[k]
    return A, b
def Pivot(A, b, k):
    s = x
    p = k
    big = np.abs([A[k, k], s[k]])
    for ii in range(k+1, n):
        dummy = np.abs(np.exp(A[ii, k] / s[ii]))
        if (dummy > big).all:
            big = dummy
            p = ii

    if p != k:
        for jj in range(k, n):
            dummy = A[p, jj]
            A[p, jj] = A[k, jj]
            A[k, jj] = dummy
        dummy = b[p]
        b[p] = b[k]
        b[k] = dummy
        dummy = s[p]
        s[p] = s[k]
        s[k] = dummy
    return s, k

def Substitue(A, b):
    x[n-1] = b[n-1] / A[n-1, n-1]
    for i in range(n-1, -1, -1):
        sum = 0
        for j in range(i+1, n):
            sum = sum + A[i, j] * x[j]
        x[i] = (b[i] - sum) / A[i, i]
    return x


gauss = Gauss(A, x)
print(gauss)