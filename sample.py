import cor
__author__ = 'maurizio'


def sample(icycl, en, vir, npart, box_length, beta, tail_corr, rc, sigma, epsilon4):
    # writes quantities to file
    if npart != 0:
        enp = en/float(npart)
        vol = box_length**3
        rho = npart/vol
        press = rho/beta + vir/(3.0*vol)
        if tail_corr:
            press += cor.corp(rc, rho, sigma, epsilon4)
    else:
        rho = 0.0
        enp = 0.0
        press = 0.0
    # print icycl, float(enp), float(press), float(rho)
    return float(enp), float(press), float(rho)