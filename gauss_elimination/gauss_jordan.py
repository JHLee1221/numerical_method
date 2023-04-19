import numpy as np
import scipy as scp
import time

start = time.time()

# setting the tolerance value
tol = 10e-6
np.set_printoptions(precision=4)

A = np.array([[1, 0, 0, 0, 0, 5],
    [-2, 3, 0, -1, 0, 24],
    [-1, 0, 5, -2, 0, -24],
    [0, -1, -2, 4, -1, 0],
    [0, 0, 0, -1, 5, 24]], float)

b = np.array([5, 24, -24, 0, 24])

# printing the equation
def PrintEqs(A, b):
    n = b.size
    for i in range(n):
        for j in range(n):
            print("{0:10.4e} ".format(A[i][j]), end="   ")
        print("|   {0:10.4e}".format(b[i]))

# printing the vector
def PrintVec(b):
    n = b.size
    for j in range(n):
        print("{0:10.4e} ".format(b[j]))
    print("")

# printing the matrix
def PrintMat(A):
    n, m = A.shape
    for i in range(n):
        for j in range(m):
            print("{0:10.4e} ".format(A[i][j]), end="")
        print("")
    print("")

def GaussJordan(A, b, prt=False):
    m, n = A.shape

    # making the augmented matrix
    augMat = np.zeros(shape=(m, n + 1))
    for i in range(0, m):
        for j in range(0, n):
            augMat[i, j] = A[i, j]
    for i in range(0, m):
        augMat[i, n] = b[i]

    # for debugging
    print("Augmented Matrix=")
    PrintMat(augMat)
    m, n = augMat.shape

    # calculating the Gauss-Jordan elimination
    x = np.zeros(m)
    for k in range(0, m):
        factor = augMat[k, k]
        for j in range(0, n):
            augMat[k, j] = augMat[k, j] / factor
        for i in range(0, m):
            if (i != k):
                factor = augMat[i, k]
                for j in range(k, n):
                    augMat[i, j] = augMat[i, j] - factor * augMat[k, j]

    # for debugging
    print(" ")
    print("After calculating Gauss-Jordan elimination:")
    PrintMat(augMat)

    # making the solution vector
    for k in range(0, m):
        x[k] = augMat[k, n - 1]

    return x

matA = A.copy()
vecb = b.copy()

GaussJordan = GaussJordan(A, b, True)

print(" ")
print("GaussJordan=", GaussJordan)
print("")
PrintVec(GaussJordan)
print('Calculation time is', time.time() - start)
