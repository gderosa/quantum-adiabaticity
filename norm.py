from math import pi, sqrt, sin

from scipy import integrate

def c2(n):
    n = float(n)
    return 4*(pi**2)*(n**2) / ((n+2)**2 * (n-2)**2)

def iterate(N):
    n = 1
    s = 0.0
    while n < N + 0.5:
        s = s + c2(n)
        n = n + 2 # keep odd
        print s

def eigen_i(n, x): # in [0, pi]
    return sqrt(2.0/pi)*sin(n*x)

def eigen_f(n, x): # in [0, 2*pi]
    return sqrt(1.0/pi)*sin((n*x)/2.0)


# iterate(10**4)

for n in range(1, 1000):
    N_i = integrate.quad(
                lambda x: eigen_i(n, x)**2, 0, pi
            )
    N_f = integrate.quad(
                lambda x: eigen_f(n, x)**2, 0, 2*pi
            )
    print N_i, N_f

