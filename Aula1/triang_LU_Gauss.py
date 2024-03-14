import numpy as np
# L é uma matriz triangular inferior
# b é termo independente
# retorno será solução do sistema
def sub_progressiva(L, b):
    b_len = len(b)
    x = [0] * b_len # inicializa vetor solução com zeros

    for i in range(b_len):
        x[i] = b[i]
        for j in range(i):
            x[i] -= (L[i][j] * x[j])
        x[i] /= L[i][i]

    return x

# U é uma matriz triangular inferior
# y é termo independente
# retorno será solução do sistema
def sub_regressiva(U, y):
    y_len = len(y)
    x = [0] * y_len # inicializa vetor solução com zeros

    for i in range(y_len - 1, -1, -1):
        x[i] = y[i]
        for j in range(i + 1, y_len):
            x[i] -= (U[i][j] * x[j])
        x[i] /= U[i][i]

    return x

# A: matriz quadrada
# retorno: L e U são matrizes triangulares inferior e superior, respectivamente
def decomposicao_LU(A):
    n = len(A)
    L = [[0] * n for i in range(n)]
    for i in range(n):
        L[i][i] = 1

    U = [[0] * n for i in range(n)]

    for k in range(n):
        for j in range(k, n):
            U[k][j] = A[k][j]
            for s in range(k):
                U[k][j] -= (L[k][s] * U[s][j])
        for i in range(k+1, n):
            L[i][k] = A[i][k]
            for s in range(k):
                L[i][k] -= (L[i][s] * U[s][k])
            L[i][k] /= U[k][k]

    return L, U

# A: matriz quadrada
# b: vetor termo independente
def LU(A, b):
    L, U = decomposicao_LU(A)
    y = sub_progressiva(L, b)
    x = sub_regressiva(U, y)
    return x

# A: matriz de coeficientes
# b: vetor termo independente
# x: vetor solução do sistema
def eliminacao_gauss(A,b):
    n = len(A)

    for k in range(n-1):
        for i in range(k+1, n):
            m = -A[i][k]/A[k][k]
            for aux in range(k, n):
                A[i][aux] = A[i][aux] + m*A[k][aux]
            b[i] = b[i] + m*b[k]

    x = sub_regressiva(A,b)

    return x

#A: matriz dos coeficientes
#b: vetor termo independente
#x: vetor solução do sistema
def eliminacao_gauss_pivo(A,b):
    n = len(A)

    for k in range(n-1):
        for aux in range(k, n):
            index = k
            maximum = abs(A[k][k])
            if abs(A[aux][k]) > maximum:
                maximum = abs(A[aux][k])
                index = aux

            A[k], A[index] = A[index], A[k]
            b[k], b[index] = b[index], b[k]

        for i in range(k+1, n):
            m = -A[i][k]/A[k][k]

            for aux in range(k,n):
                A[i][aux] = A[i][aux] + m*A[k][aux]
            
            b[i] = b[i] + m*b[k]

    x = sub_regressiva(A,b)

    return x

#A: matriz não-singular
#L, U: matrizes triangulares inferior e superior, respectivamente
#P: matriz de permutação
def decomposicao_lup(A):
    n = len(A)

    U = np.copy(A)
    L = np.eye(n)
    P = np.copy(L)

    for j in np.arange(n-1):
        # pivô
        pivo = np.argmax(np.abs(U[j:n,j]))
        P[[j,pivo], :] = P[[pivo,j], :]
        U[[j,pivo], j:n] = U[[pivo,j], j:n]
        L[[j,pivo], :j] = L[[pivo,j], :j]

        for i in np.arange(j+1, n):
            L[i,j] = U[i,j]/U[j,j]
            U[i,j+1:n] = U[i,j+1:n] - L[i,j]*U[j,j+1:n]
            U[i,j] = 0

    return L, U, P
    


if __name__ == "__main__":
    A = []
    n = int(input("Enter the size of the matrix: "))
    print("Enter the elements of the matrix:")
    for i in range(n):
        row = list(map(int, input().split()))
        A.append(row)
    
    b = list(map(int, input("Enter the elements of the vector term: ").split()))
    
    L, U, P = decomposicao_lup(A)
    print("Matrix L:")
    for row in L:
        print(row)
    print("Matrix U:")
    for row in U:
        print(row)
    print("Matrix P:")
    for row in P:
        print(row)


