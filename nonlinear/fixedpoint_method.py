import numpy as np
import scipy as scp
import sympy as smp

tol = 10e-9
np.set_printoptions(precision=4)

def Fixpt(g, x0, es, imax):
    xr = x0
    x = []
    e = []
    iterNum = []

    for i in range(imax):
        xrold = xr
        xr = g(xrold)

        if xr != 0:
            ea = np.fabs((xr-xrold)/(xr)) * (100.0)

        if ea < es:
            break

    return iterNum, xr, e, x

def g(x):
    return (1 + (2/x))

iterNum, xr, err, x = Fixpt(g, 1.0, tol, 100)
print(xr)
