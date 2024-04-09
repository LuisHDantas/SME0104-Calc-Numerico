def secante(func, x0, x1, tol=1e-6, kmax=100):
    """
    Método da secante para encontrar o zero de uma função.
    """
    f0 = func(x0)

    for k in range(kmax):
        f1 = func(x1)
        ds = f1*(x1 - x0)/(f1 - f0)
        x0 = x1
        x1 -= ds
        if abs(ds) < tol:
            x = x1
            return x, k
        f0 = f1
    
    x = None
    return x, k

if __name__ == "__main__":
    print("Método da secante")
    f = lambda x: x**2 + x - 6
    x0 = 5
    x1 = 5.5
    x, k = secante(f, x0, x1)
    print(f"x = {x} em {k} iterações")