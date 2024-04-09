def newton(func, dfunc, x, tol=1e-6, kmax=100):
    for k in range(kmax):
        dx = func(x) / dfunc(x)
        x -= dx
        if abs(dx) < tol:
            return x, k
    x = None
    return x, k


if __name__ == "__main__":
    print(newton(lambda x: x**2 + x - 6, lambda x: 2*x + 1, 5.5))