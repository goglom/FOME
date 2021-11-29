import numpy as np
import matplotlib.pyplot as plt
import warnings
from matplotlib.widgets import Slider

warnings.filterwarnings("ignore")

step = 0.1
e = 1.6 * 1e-19
a = 0.54310 * 1e-9
m = 0.19*9.1 * 1e-31
b = 1/5 * a
h = 6.6 * 1e-34

E = np.arange(0, 50, step)
N = 3

def dirac(V,E):
    alpha, beta = ab(V, E)
    r = 1/a * np.arccos(np.cos(alpha * a) - (beta**2 * b * a) / 2 * np.sin(alpha*a) / (alpha * a))
    return r

def ab(V, E):
    alpha = np.sqrt(8 * np.pi**2 * m * E * e / h**2)
    beta  = np.sqrt(8 * np.pi**2 * m * (E+V) * e / h**2)
    return alpha, beta

def kronig_penney(V, E):
    alpha, beta = ab(V, E)
    r = 1/a * np.arccos((np.cos(beta*b) * np.cos(alpha*(a-b)) ) - (alpha**2 + beta**2) / (2*alpha * beta) * np.sin(beta*b) * np.sin(alpha*(a-b)))
    return r



fig, (ax) = plt.subplots()
fig.subplots_adjust(bottom=0.2, left=0.1)
ax_v= plt.axes([0.10, 0.05, 0.8, 0.03])
slider_v = Slider(ax_v, "V", 0, 30, 10, valstep=0.5)


def update(*args):
    ax.clear()

    V = slider_v.val

    k_kp = kronig_penney(V, E)

    ax.plot(k_kp, E, 'b')
    ax.plot(-k_kp, E, 'b')

    k_d = dirac(V, E)
    ax.plot(k_d, E, 'r')
    ax.plot(-k_d, E, 'r')
    ax.grid()


slider_v.on_changed(update)
update()

plt.show()
