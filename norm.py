from math import pi, sqrt, sin

from scipy import integrate

def c2(n):
    n = float(n)
    return 4*(pi**2)*(n**2) / ((n+2)**2 * (n-2)**2)

def eigen_i(n, x): # in [0, pi], 0 elsewhere
    return sqrt(2.0/pi)*sin(n*x)

def eigen_f(n, x): # in [0, 2*pi]
    return sqrt(1.0/pi)*sin((n*x)/2.0)


def check_norms(N):
    for n in range(1, 1000):
        N_i = integrate.quad(
                    lambda x: eigen_i(n, x)**2, 0, pi
                )
        N_f = integrate.quad(
                    lambda x: eigen_f(n, x)**2, 0, 2*pi
                )
        print N_i, N_f

def c(ni, nf): # eigenfunctions are real
    return integrate.quad(
                lambda x: eigen_f(nf, x) * eigen_i(ni, x), # < bra | ket >
                0, pi
            )

def check_convergence(ni, N):
    nf = 1
    s = 0.0
    while nf < N:
        cnf, err = c(ni, nf)
        if 20 * abs(err) > abs(cnf): # discard if more than 5% numeric error
            pass
            # print 'DISCARDED unstable numeric integral: ' + str(cnf) + ' +/- ' + str(err)
        else:
            s = s + cnf**2
            print '+ (' + str(cnf) + ' +/- ' + str(err) + ')**2 = ' + str(s)
        nf = nf + 1

check_convergence(1, 10**4)

