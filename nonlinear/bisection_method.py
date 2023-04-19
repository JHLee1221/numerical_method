import numpy as np
import scipy as scp
import sympy as smp

tol = 10e-9
np.set_printoptions(precision=4)

def Bisect(f, xl, xu, es, imax):
    xr = 0.0
    x = []
    e = []
    iterNum = []
    fl = f(xl)

    for i in range(imax):
        xrold = xr

        xr = (xl + xu) / 2
        fr = f(xr)

        if xr != 0:
            ea = np.fabs((xr -xrold) / (xr)) * (100.0)

        iterNum.append(i)
        x.append(xr)
        e.append(ea)

        test = fl*fr

        if test < 0.0:
            xu = xr
        elif test > 0.0:
            xl = xr
            fl = fr
        else:
            ea = 0.0

        if ea < es:
            break
    return iterNum, xr, e, x

def f(x):
    return 2.78966968e-17*x**7 -1.15038952e-13*x**6 + 1.73860468e-10*x**5 -1.15250088e-07*x**4 +2.79743932e-05*x**3 + 1.59246727e-03*x**2 -4.33652768e-01*x + 1.48445801e+01
    #return -2.82930021e-21*x**9 + 1.47564154e-17*x**8 -3.25700158e-14*x**7 + 3.95713155e-11*x**6 -2.88287819e-08*x**5 + 1.28570500e-05*x**4 -3.42816533e-03*x**3 + 5.03661322e-01*x**2 -3.18066570e+01*x + 2.36759241e+02

iterNum, xr, err, x = Bisect(f, 0.0, 1.0, tol, 100)
print(xr)

