import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

e = 1.6 * 1e-19
a = 10e-9
m = 0.51* 9.1e-31
b = 3e-9
h = 6.6 * 1e-34


def dirac(V,E):
    alpha, beta = ab(V, E)
    r = 1/a * np.arccos(np.cos(alpha * a) - (beta**2 * b * a) / 2 * np.sin(alpha*a) / (alpha * a))
    return r


def ab(V, E):
    alpha = np.sqrt(8 * np.pi**2 * m * np.abs(E) * e / h**2)
    beta  = np.sqrt(8 * np.pi**2 * m * (V - np.abs(E)) * e / h**2)
    return alpha, beta


def kronig_penney(V, E):
    alpha, beta = ab(V, E)
    r =1/a * np.arccos((np.cos(beta*b) * np.cosh(alpha*(a-b)) ) - (beta**2 - alpha**2) / (2*alpha * beta) * np.sin(beta*b) * np.sinh(alpha*(a-b)))
    return r

V = 0.3
step = 1e-8
E = np.arange(-V, 0, step)

zip_factor = 15
k_kp = kronig_penney(V, E)[::zip_factor]
E = E[::zip_factor]
plt.plot(k_kp, E, 'r')
plt.plot(-k_kp, E, 'r')
plt.grid()
plt.show()

