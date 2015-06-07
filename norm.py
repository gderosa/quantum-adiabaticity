from math import pi, sqrt, sin

from scipy import integrate
from skmonaco import mcquad

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
    val, err = 1., 1.
    try:
        c.npoints = c.npoints
    except AttributeError:
        c.npoints = 1e4
    while True:
        print 'npoints=' + str(c.npoints)
        val, err = mcquad(
                lambda x: eigen_f(nf, x[0]) * eigen_i(ni, x[0]), # < bra | ket >
                npoints=c.npoints, xl=[0.], xu=[pi]
            )
        if abs(err / val) < 1e-2:
            return (val, err)
        else:
            c.npoints = c.npoints * 3

def even(x):
    return abs((x / 2.0) - round(x / 2.0)) < 0.25

def check_convergence(ni, N):
    nf = 1
    s = 0.0
    discarded=0
    while nf < N:
        if even(nf) and round(nf) != 2:
            pass
        else:
            cnf, err = c(ni, nf)
            if 20. * abs(err) > abs(cnf): # discard if more than 5% numeric error
                print 'DISCARDED unstable numeric integral: ' + str(cnf) + ' +/- ' + str(err)
                discarded = discarded + 1
            else:
                s = s + cnf**2
                print '+ (' + str(cnf) + ' +/- ' + str(err) + ')**2 = ' + str(s)
            print '  #' + str(nf) + '  discarded=' + str(discarded)
        nf = nf + 1

check_convergence(1, 1000)

