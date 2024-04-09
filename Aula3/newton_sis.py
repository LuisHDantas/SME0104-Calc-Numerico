import numpy as np

def newton_sis(F, Jac, x, tol=1e-6, kmax=1000):
    for k in range(kmax):
        v = np.linalg.solve(Jac(x), F(x))
        x = x - v
        if np.linalg.norm(v) < tol:
            return x, k
    return x, kmax

def F(x):
    return np.array([x[0]**2 - 2*x[0] - x[1] + 1, x[0]**2 + x[1]**2 - 1])

def J(x):
    return np.array([[2*x[0] - 2, -1], [2*x[0], 2*x[1]]])

# Run the test case
if __name__ == "__main__":
    print(newton_sis(F, J, np.array([1.0, -1.0])))