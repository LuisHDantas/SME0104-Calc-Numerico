import numpy as np

def lagrange_interp(xi, yi, x):
    xi, yi, x = np.array(xi), np.array(yi), np.array(x)
    if xi.ndim > 1:
        xi, yi, x = xi.T, yi.T, x.T
    n = len(xi)
    L = np.ones((n, len(x)))
    for i in range(n):
        for j in range(n):
            if i != j:
                L[i, :] = L[i, :] * (x - xi[j]) / (xi[i] - xi[j])
    y = yi @ L
    return y

def newton_interp(xi, yi, x):
    xi, yi, x = np.array(xi), np.array(yi), np.array(x)
    if xi.ndim > 1:
        xi, yi, x = xi.T, yi.T, x.T
    n = len(xi)
    ni = len(x)
    N = np.ones((n, ni))
    D = np.zeros((n, n))
    D[:, 0] = yi
    for j in range(n-1):  # table of divided differences
        for i in range(n-j-1):
            D[i, j+1] = (D[i+1, j] - D[i, j]) / (xi[i+j+1] - xi[i])
    for i in range(1, n):  # Newton's form
        N[i, :] = N[i-1, :] * (x - xi[i-1])
    y = D[0, :] @ N
    return y

    