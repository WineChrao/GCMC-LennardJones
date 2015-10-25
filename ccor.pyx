import numpy as np
cimport numpy as np

__author__ = 'maurizio'

def corp(double r, double rho, double sigma, double epsilon4):
    """
    pressure correction
    """
    cdef double sig3, ri3
    sig3 = sigma**3
    ri3 = sig3/(r*r*r)
    return 4.0*np.pi*epsilon4*(rho**2)*sig3*(2.0*ri3*ri3*ri3/9.0-ri3/3.0)


def coru(double r, double rho, double sigma, double epsilon4):
    """
    energy correction
    """
    cdef double sig3, ri3
    sig3 = sigma**3
    ri3 = sig3/(r*r*r)
    return 2.0*np.pi*epsilon4*(rho*sig3)*(ri3*ri3*ri3/9.0-ri3/3.0)
