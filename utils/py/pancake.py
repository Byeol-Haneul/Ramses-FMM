import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

"""
Parameters for the pancake collapse test
"""
omega_m = 0.3
omega_l = 0.7
n = 1024
del_ini = 0.1
aexp_ini = 0.01

def fy(a):
    """
    Integrand for the integral
    """
    y=omega_m*(1.0/a-1.0) + omega_l*(a*a-1.0) + 1.0
    return 1.0/y**1.5

def d1(a):
    """
    Computes the linear growing mode D1 in a FRW universe.
    """
    y = omega_m * (1.0 / a - 1.0) + omega_l * (a**2 - 1.0) + 1.0
    y12 = np.sqrt(y)
    return y12 / a * sp.integrate.romberg(fy,1e-6,a)

def fpeebl(a):
    """
    Computes the Pebbles growth factor f = d log D1 / d log a.
    """
    y = omega_m * (1.0 / a - 1.0) + omega_l * (a**2 - 1.0) + 1.0
    fact = sp.integrate.romberg(fy,1e-6,a)
    return (omega_l * a**2 - 0.5 * omega_m / a) / y - 1.0 + a * fy(a) / fact

xf = np.linspace(0, 1, n+1)
xp = 0.5 * ( xf[:-1] + xf[1:] )

disp = del_ini * np.sin( 2 * np.pi * xp ) / ( 2 * np.pi )
vfact = fpeebl(aexp_ini) * aexp_ini * np.sqrt( omega_m / aexp_ini + omega_l *aexp_ini**2 )
vp = del_ini * vfact * np.sin( 2 * np.pi * xp ) / ( 2 * np.pi )
mp = np.ones(n) / n
xp = xp + disp - 0.5

data = np.column_stack((xp, vp, mp))
np.savetxt( "ic_part", data, fmt="%.10f", delimiter="  " )





