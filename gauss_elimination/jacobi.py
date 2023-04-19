import numpy as np
from scipy.linalg import solve
import time

start = time.time()

A =  np.array([[1, 0, 0, 0, 0],
           [-2, 3, 0, -1, 0],
           [-1, 0, 5, -2, 0],
           [0, -1, -2, 4, -1],
           [0, 0, 0, -1, 5]], float)

b = np.array([5, 24, -24, 0, 24])

tol = 10e-6
np.set_printoptions(precision=6)


def PrintEqs(A, b):
    n = b.size
    for i in range(n):
        for j in range(n):
            print("{0:10.6e} ".format(A[i][j]), end="   ")
        print("|   {0:10.6e}".format(b[i]))


# printing the pretty vector
def PrintVec(b):
    n = b.size
    for j in range(n):
        print("{0:10.6e} ".format(b[j]))
    print("")


# printing the pretty matrix
def PrintMat(A):
    n, m = A.shape
    for i in range(n):
        for j in range(m):
            print("{0:10.6e} ".format(A[i][j]), end="")
        print("")
    print("")

def Jacobi(A, b, tol=1.0e-9, iterNum=500):
    rows, cols = A.shape
    currx = np.zeros((rows), dtype=float)
    prevx = np.zeros((rows), dtype=float)

    invD = np.zeros(A.shape, dtype=float)
    for i in range(0, rows):
        invD[i, i] = 1.0 / (A[i, i])

    L = np.zeros(A.shape, dtype=float)
    for i in range(0, rows):
        for j in range(0, i):
            L[i, j] = -A[i, j]

    U = np.zeros(A.shape, dtype=float)
    for i in range(0, rows):
        for j in range(0, i):
            U[j, i] = -A[j, i]

    count = 0
    while (count < iterNum):
        matFirst = np.matmul(invD, (L + U))
        matSecond = np.matmul(invD, b)
        currx = np.matmul(matFirst, prevx) + matSecond

        err = np.linalg.norm(np.fabs(currx - prevx))
        if err < tol:
            break
        else:
            count = count + 1
            prevx = currx

    return currx, count

matA = A.copy()
vecb = b.copy()
matAA = A.copy()
vecbb = b.copy()

Jacobi, iterNumJacobi = Jacobi(matA, vecb, tol = 1.0e-6, iterNum = 1000)

print("iterNumJacobi=", iterNumJacobi)
print(Jacobi)
print("")
PrintVec(Jacobi)
print("CrossCheck=", np.dot(matAA, Jacobi) - vecbb)
print("")
print('Calculation time is ',time.time() - start)