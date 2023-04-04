def gauss_jordan(A):
    n = len(A)
    for i in range(n):
        # Find pivot element
        max_row = i
        for j in range(i+1, n):
            if abs(A[j][i]) > abs(A[max_row][i]):
                max_row = j
        A[i], A[max_row] = A[max_row], A[i]

        # Scale pivot row
        pivot = A[i][i]
        if pivot == 0:
            return "No unique solution exists."
        for j in range(i, n+1):
            A[i][j] /= pivot

        # Eliminate other rows
        for j in range(n):
            if j != i:
                factor = A[j][i]
                for k in range(i, n+1):
                    A[j][k] -= factor * A[i][k]

    x = [A[i][n] for i in range(n)]
    return x

A = [[1, 0, 0, 0, 0, 5],
    [-2, 3, 0, -1, 0, 24],
    [-1, 0, 5, -2, 0, -24],
    [0, -1, -2, 4, -1, 0],
    [0, 0, 0, -1, 5, 24]]

x = gauss_jordan(A)

print(x)