import numpy as np
import lattice
# import cor
import ccor
# import ener
import adjust
import mcmove
import mcexch
import sample
import atomvis
import matplotlib.pyplot as plt
# import cavity
# import time

import ccavity
import cener

__author__ = 'maurizio'


"""
grand canonical coupled to an ideal gas bath with pressure: p in in_press_list
        chemical potential bath: mu^b = mu^0 + ln(beta*p)/beta
                                      = ln (beta*pid*Lamda^3)/beta
        zz is defined as         zz   = exp(beta*mu^b)/Lambda^3
                                      = beta*pid
        excess chemical pot.     muex = mu^b -mu^id
                                      = mu^0 + ln(beta*pid)/beta - mu^0 - ln(rho)
                                      = ln(zz)/beta - ln <rho>
"""


def main():

    res_en = []
    res_rho = []
    res_press = []
    in_press_list = [0.01, 0.1, 1, 2, 5, 10, 20]
    for press in in_press_list:

        r_en, r_press, r_rho = calc(press)
        res_en.append(r_en)
        res_press.append(r_press)
        res_rho.append(r_rho)
    print res_rho
    print res_press
    plt.plot(res_rho, res_press, 'o')
    plt.xlabel('density')
    plt.ylabel('Pressure')
    plt.show()
    atomvis.Close()


def calc(reservoir_pressure):
    rc = 500  # initial guess for potential cut off distance
    npart = 50  # particles
    temp = 2  # temperature
    rho = 0.6  # initial density of fluid in box
    n_nodes = 1000    # grid nodes for cavity biasing
    equil_cycles = 200
    prod_cycles = 100
    sample_frequency = 2
    displ_attempts = 100
    max_displ = 0.09
    exch_attempts = 10
    epsilon = 1.0
    sigma = 1.0
    mass = 1.0
    i_seed = 68324
    tail_corr = 1
    shift_pot = 0
    use_cavity_bias = True
    box_length = (npart/rho)**(1.0/3.0)
    rc = min(rc, box_length/2.0)
    rc_sq = rc * rc

    beta = 1.0/temp
    zz = beta * reservoir_pressure

    epsilon4 = 4.0 * epsilon
    epsilon48 = 48.0 * epsilon
    sigma_sq = sigma * sigma

    coru = 0.0
    corp = 0.0
    e_cut = 0.0

    part_pos_array = lattice.lattice(box_length, npart)
    press_list = []
    en_list = []
    rho_list = []
    np.random.seed(i_seed)

    if shift_pot:
        e_cut, dummy = cener.ener(rc_sq, rc_sq, sigma_sq, epsilon4, epsilon48, shift_pot, e_cut)

    if tail_corr:
        corp = ccor.corp(rc, rho, sigma, epsilon4)
        coru = ccor.coru(rc, rho, sigma, epsilon4)
    en, vir, rho = cener.totenerg(npart, part_pos_array, box_length, rc, sigma, rc_sq,
                                 sigma_sq, epsilon4, epsilon48, shift_pot, e_cut, tail_corr)
    nmoves = displ_attempts + exch_attempts
    # display particles in real time
    atomvis.Clear()
    atomvis.Init(part_pos_array, box_length)

    cavity_nodes = ccavity.cavity_locations(box_length, n_nodes)

    for ii in range(2):
        """
        ii = 0 equilibration
        ii = 1 production
        """
        if ii == 0:
            ncycles = equil_cycles
        else:
            ncycles = prod_cycles
        attempt = 0
        naccp = 0
        attemptp = 0
        nacc = 0
        atte = 0
        acce = 0
        rhoav = 0.0
        nsampav = 0
        naccp, attemptp, max_displ = adjust.adjust(attempt, nacc, max_displ, box_length/2.0, attemptp, naccp)
        for icycle in range(ncycles):
            for imove in range(nmoves):
                ran = np.random.random()*nmoves
                if ran < displ_attempts and len(part_pos_array) != 0:

                    # attempt to displace a particle
                    en, vir, nacc, attempt, part_pos_array = mcmove.mcmove(npart, part_pos_array, max_displ, beta,
                                                                           en, vir, attempt, nacc, box_length, rc_sq,
                                                                           sigma_sq, epsilon4, epsilon48, shift_pot, e_cut)
                else:
                    av_cav, n_cavs = ccavity.available_cavities(cavity_nodes, part_pos_array, 0.8 * sigma)
                    # attempt to exchange a particle with the reservoir
                    en, vir, nacc, attempt, part_pos_array, npart = mcexch.mcexch(npart, part_pos_array, beta,  en, vir, attempt,
                                                                                  nacc, box_length, rc, rc_sq, sigma, sigma_sq, epsilon4,
                                                                                  epsilon48, shift_pot, tail_corr, e_cut, zz, av_cav, n_cavs,
                                                                                  n_nodes, use_cavity_bias)
                if len(part_pos_array) > 0:
                    atomvis.Update(part_pos_array)

            if ii == 1:
                # sample averages
                if icycle % sample_frequency == 0:
                    samp_enp, samp_press, samp_rho = sample.sample(icycle, en, vir, npart, box_length, beta, tail_corr, rc, sigma, epsilon4)
                    # to determine excess chem. potential
                    nsampav += 1
                    rhoav += npart/(box_length**3)
                    press_list.append(samp_press)
                    en_list.append(samp_enp)
                    rho_list.append(samp_rho)
            if (icycle % ncycles/5) == 0:
                naccp, attemptp, max_displ = adjust.adjust(attempt, nacc, max_displ, box_length/2.0, attemptp, naccp)
            if ncycles != 0:
                en, vir, rho = cener.totenerg(npart, part_pos_array, box_length, rc, sigma, rc_sq,
                                             sigma_sq, epsilon4, epsilon48, shift_pot, e_cut, tail_corr)
                if rhoav != 0 and nsampav != 0:
                    muex = np.log(zz)/beta - np.log(rhoav/nsampav)/beta
                    mu = muex + np.log(rhoav/nsampav)/beta

    print muex, mu
    print attempt, nacc
    print attemptp, naccp
    print npart

    return np.mean(en_list), np.mean(press_list), np.mean(rho_list)


if __name__ == "__main__":
    main()
