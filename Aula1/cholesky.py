import math 
# A: matriz SPD
# H: triang. inferior, tal que A = H*H'
def cholesky(A):
    n = len(A)
    H = A.copy()

    for i in range(n):
        for j in range(i+1, n):
            H[i][j] = 0

    for k in range(n-1):
        H[k][k] = math.sqrt(H[k][k])
        for aux in range(k+1, n):
            H[aux][k] = H[aux][k]/H[k][k]
        
        for j in range(k+1, n):
            for aux in range(j, n):
                H[aux][j] = H[aux][j] - H[aux][k]*H[j][k]

    H[n-1][n-1] = math.sqrt(H[n-1][n-1])

    return H
            


if __name__ == "__main__":
    matrix = [[6, 15, 55], [15, 55, 225], [55, 225, 979]]
    result = cholesky(matrix)
    for row in result:
        print(row)