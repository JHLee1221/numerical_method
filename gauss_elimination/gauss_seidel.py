import numpy as np

def gauss_seidel(A, b, x0, epsilon, max_iterations):
    n = len(A)
    x = x0.copy()

    #Gauss-Seidal Method [By Bottom Science]

    for i in range(max_iterations):
        x_new = np.zeros(n)
        for j in range(n):
            s1 = np.dot(A[j, :j], x_new[:j])
            s2 = np.dot(A[j, j + 1:], x[j + 1:])
            x_new[j] = (b[j] - s1 - s2) / A[j, j]
        if np.allclose(x, x_new, rtol=epsilon):
            return x_new
        x = x_new
    return x

A = np.array([[1, 0, 0, 0, 0],
           [-2, 3, 0, -1, 0],
           [-1, 0, 5, -2, 0],
           [0, -1, -2, 4, -1],
           [0, 0, 0, -1, 5]], float)

b = np.array([5, 24, -24, 0, 24], float)

n = len(b)
x0 = np.zeros(n, float)
eps = 1e-5
max_iter = 100

x = gauss_seidel(A, b, x0, eps, max_iter)
print(x)