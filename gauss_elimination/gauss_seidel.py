import numpy as np
import time
start = time.time()

tol = 10e-6
np.set_printoptions(precision=6)

A = np.array([[1, 0, 0, 0, 0],
           [-2, 3, 0, -1, 0],
           [-1, 0, 5, -2, 0],
           [0, -1, -2, 4, -1],
           [0, 0, 0, -1, 5]], float)

b = np.array([5, 24, -24, 0, 24], float)

# printing the equation
def PrintEqs(A, b):
    n = b.size
    for i in range(n):
        for j in range(n):
            print("{0:10.6e} ".format(A[i][j]), end="   ")
        print("|   {0:10.6e}".format(b[i]))


# printing the vector
def PrintVec(b):
    n = b.size
    for j in range(n):
        print("{0:10.6e} ".format(b[j]))
    print("")


# printing the matrix
def PrintMat(A):
    n, m = A.shape
    for i in range(n):
        for j in range(m):
            print("{0:10.6e} ".format(A[i][j]), end="")
        print("")
    print("")


def GaussSeidel(A, b, omega=1.0, tol=1.0e-9, iterNum=500):
    m, n = A.shape
    x = np.zeros(n, dtype=float)
    for i in range(0, n):
        temp = A[i][i]
        for j in range(0, n):
            A[i][j] = A[i][j] / temp
        b[i] = b[i] / temp
    for i in range(0, n):
        sumTemp = b[i]
        for j in range(0, n):
            if (i != j):
                sumTemp = sumTemp - A[i][j] * x[j]
        x[i] = sumTemp

    count = 0
    while (count < iterNum):
        sentinel = 1
        for i in range(0, n):
            old = x[i]
            sumTemp = b[i]
            for j in range(0, n):
                if (i != j):
                    sumTemp = sumTemp - A[i][j] * x[j]
            x[i] = omega * sumTemp + (1.0 - omega) * old
            if ((sentinel == 1) & (x[i] != 0.0)):
                err = np.fabs((x[i] - old) / x[i]) * 100.0
                if (err > tol):
                    sentinel = 0
        count = count + 1
        if (sentinel == 1):
            break

    return x, count

matA = A.copy()
vecb = b.copy()
matAA = A.copy()
vecbb = b.copy()

GaussSeidel, iterNumSeidel = GaussSeidel(A, b, omega = 1.0, tol = 1.0e-6, iterNum = 1000)
print("GaussSeidel=", iterNumSeidel)
print(GaussSeidel)
PrintVec(GaussSeidel)
print("CrossCheck=", np.dot(matAA, GaussSeidel) - vecbb)
print("")
print('Calculation time is ',time.time()- start)
