import numpy as np

def mmq(x, y, k):
    n = len(x)

    # Create Vandermonde matrix
    X = np.vander(x, k+1)

    # Calculate the least squares solution
    A = np.linalg.inv(X.T @ X) @ X.T @ y

    return A