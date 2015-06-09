from cmath import pi, sqrt, sin, exp

from scipy import integrate
from skmonaco import mcquad

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt



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
    return integrate.quadrature(
        lambda x: eigen_f(nf, x) * eigen_i(ni, x), # < bra | ket >
        0., pi,
        rtol=1e-5, maxiter=1000,
        vec_func=False
    )

def c1_analytic(nf):
    if even(nf):
        if abs(nf-2.) < 0.25 :
            return 1./sqrt(2)
        else:
            return 0.0
    else:
        return (4.*sqrt(2)/pi) * ( 1. / ((nf+2.)*(nf-2.)) ) * sin(nf*pi/2.)


def c1_analytic_2(nf):
    return c1_analytic(nf)**2


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

def check_convergence_analytic_1(N):
    nf = 1.
    s = 0.0
    while nf < N:
        add = c1_analytic_2(nf)
        s = s + add
        print '  #' + str(nf) + '  +' + str(add) + ' = ' + str(s)
        nf = nf + 1.

def evolve(x, t, iterations=990):
    s = 0.0
    for nf in range(1, iterations):
        E = nf**2 / 8.
        omega = E
        s = s + c1_analytic(nf) * eigen_f(nf, x) * exp(-1j*omega*t)
    return s

# Returns, X, Y, Z numpy vectors for matplotlib3d; where Z is the
# independent variable and X, Y real and imaginary parts of Psi.
def evolve_v(t):
    Z = np.linspace(0, 2*pi, 960)
    Re = np.ndarray(Z.size)
    Im = np.ndarray(Z.size)
    it = np.nditer(Z, flags=['f_index'])
    while not it.finished:
        i = it.index
        z = it[0]
        Psi = evolve(z, t)
        Re[i], Im[i] = Psi.real, Psi.imag
        it.iternext()
    return Re, Im, Z





mpl.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y, Z = evolve_v(1e-80)
ax.plot(X, Y, Z, label='hello quantum')
ax.legend()

plt.show()





# evolve_all()

# check_convergence(1, 1000)

#check_convergence_analytic_1(10034)

