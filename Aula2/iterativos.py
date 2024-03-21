import numpy as np

def gauss_jacobi(A, b, x0, tol):
    n = A.shape[0]
    D = np.diag(np.diag(A))
    C = np.eye(n) - np.linalg.solve(D,A)
    g = np.linalg.solve(D,b)
    kmax = 10000
    k = 0
    

    while np.linalg.norm(b - A.dot(x0)) > tol and k < kmax:
        k = k + 1
        x0 = C.dot(x0) + g

    if k == kmax:
        print("O método não convergiu")
        return
    
    x = x0
    return (x,k)

def gauss_seidel(A, b, x0, tol):
    L = np.tril(A)
    R = np.triu(A,1)
    C = np.linalg.solve(-L,R)
    g = np.linalg.solve(L,b)
    kmax = 10000
    k = 0

    while np.linalg.norm(b - A.dot(x0)) > tol and k < kmax:
        k = k + 1
        x0 = C.dot(x0) + g

    if k == kmax:
        print("O método não convergiu")
        return
    
    x = x0
    return (x,k)



if __name__ == "__main__":
    print(gauss_jacobi(np.array([[5,-1,1],[2,4,0],[1,1,5]]),np.array([10,12,-1]),np.array([0,0,0]),1e-5))
    print(gauss_seidel(np.array([[5,-1,1],[2,4,0],[1,1,5]]),np.array([10,12,-1]),np.array([0,0,0]),1e-5))