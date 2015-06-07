from math import pi

def phi2(n):
    n = float(n)
    return 4*(pi**2)*(n**2) / ((n+2)**2 * (n-2)**2)

def iterate(N):
    n = 1
    s = 0.0
    while n < N + 0.5:
        s = s + phi2(n)
        n = n + 2 # keep odd
        print s

iterate(10)

