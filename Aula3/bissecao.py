import numpy as np

def bissecao(func, a, b, tol):
    # [a,b] é o intervalo com f(a)*f(b) < 0

    x = (a+b)/2
    erro = np.inf

    while erro > tol:
        if func(a)*func(x) < 0:
            b = x
        else:
            a = x
        
        x0 = x
        x = (a+b)/2
        erro = np.abs(x-x0)

    return (a,b), x

if __name__ == "__main__":
    f = lambda x: (x-2)**2 - np.log(x)
    a = 1
    b = 2
    tol = 1e-5

    raiz = bissecao(f, a, b, tol)
    print(raiz)