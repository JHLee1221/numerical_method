import numpy as np
from scipy.linalg import solve


def jacobi(A, b, x, n):
    D = np.diag(A)
    R = A - np.diagflat(D)



    for i in range(n):
        x = (b - np.dot(R, x)) / D
    return x


'''___Main___'''

A =  np.array([[1, 0, 0, 0, 0],
           [-2, 3, 0, -1, 0],
           [-1, 0, 5, -2, 0],
           [0, -1, -2, 4, -1],
           [0, 0, 0, -1, 5]], float)
b = np.array([5, 24, -24, 0, 24])
x = [1.0, 1.0, 1.0, 1.0, 1.0]
n = len(b)

x = jacobi(A, b, x, n)
print

solve(A, b)
x = jacobi(A, b, x, n)
print(solve(A, b))