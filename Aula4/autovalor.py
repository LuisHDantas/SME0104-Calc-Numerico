import numpy as np
from scipy.linalg import lu_factor, lu_solve


def clgs(A):
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))
    for j in range(n):
        V = A[:, j].copy()
        for i in range(j):
            R[i, j] = np.dot(Q[:, i], A[:, j])
            V = V - R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(V)
        Q[:, j] = V / R[j, j]
    return Q, R

def mgs(matrix):
    n = matrix.shape[1]
    q1 = np.array(matrix, dtype='float64')
    r1 = np.zeros((n, n))
    for k in range(n):
        a_k = q1[..., k]
        r1[k,k] = np.linalg.norm(a_k)
        a_k /= r1[k, k]
        for i in range(k+1, n):
            a_i = q1[..., i]
            r1[k,i] = np.transpose(a_k) @ a_i
            a_i -= r1[k, i] * a_k
    return q1, r1

def francis(A, tol):
    n = A.shape[0]
    V = np.eye(n)
    erro = np.inf
    
    while erro > tol:
        Q, R = mgs(A)
        A = R @ Q
        V = V @ Q
        erro = np.max(np.abs(np.tril(A, -1)));

    return V, np.diag(A)


# não tinha caso de teste para confirmar se a função está correta
def potencias(A, tol):
    k = 0
    kmax = 1000
    erro = np.inf
    n = A.shape[0]
    y0 = np.zeros(n)
    y0[0] = 1

    while erro > tol and k < kmax:
        x = A @ y0
        y = x / np.linalg.norm(x)
        erro = np.abs(np.abs(y0.T @ y) - 1)
        y0 = y
        k += 1

    lambda_ = y.T @ A @ y
    return lambda_, y, k

# não tinha caso de teste para confirmar se a função está correta
def potencia_inv(A, tol, alpha=0):
    k = 0
    kmax = 1000
    erro = np.inf
    n = A.shape[0]
    y0 = np.zeros(n)
    y0[0] = 1
    B = A - alpha * np.eye(n)
    lu, piv = lu_factor(B)

    while erro > tol and k < kmax:
        x = lu_solve((lu, piv), y0)
        y = x / np.linalg.norm(x)
        erro = np.abs(np.abs(y0.T @ y) - 1)
        y0 = y
        k += 1

    lambda_ = y.T @ A @ y
    return lambda_, y, k
