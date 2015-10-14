import numpy as np

__author__ = 'maurizio'

def corp(r, rho, sigma, epsilon4):
    """
    pressure correction
    """
    sig3 = sigma**3
    ri3 = sig3/(r*r*r)
    return 4.0*np.pi*epsilon4*(rho**2)*sig3*(2.0*ri3*ri3*ri3/9.0-ri3/3.0)


def coru(r, rho, sigma, epsilon4):
    """
    energy correction
    """
    sig3 = sigma**3
    ri3 = sig3/(r*r*r)
    return 2.0*np.pi*epsilon4*(rho*sig3)*(ri3*ri3*ri3/9.0-ri3/3.0)
