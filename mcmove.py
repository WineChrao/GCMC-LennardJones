# import ener
import numpy as np
import cener
import time
__author__ = 'maurizio'


def mcmove(npart, part_pos_array, dr, beta,  en, vir, attempt, nacc, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut):
    attempt += 1
    jb = 0
    # select a particle at random
    o = int(npart*np.random.random())  # +1
    #  calculate energy old configuration
    eno, viro = cener.eneri(npart, part_pos_array, o, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
    # give particle a random displacement and store previous position
    xo = part_pos_array[o, 0]
    yo = part_pos_array[o, 1]
    zo = part_pos_array[o, 2]

    part_pos_array[o, 0] += (np.random.random()-0.5)*dr
    part_pos_array[o, 1] += (np.random.random()-0.5)*dr
    part_pos_array[o, 2] += (np.random.random()-0.5)*dr
    # calculate energy new configuration:
    enn, virn = cener.eneri(npart, part_pos_array, o, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
    # acceptance test
    if np.random.random() < min(1, np.exp(-beta*(enn-eno))):
        # accepted
        nacc += 1
        en = en + (enn-eno)
        vir = vir + (virn-viro)
        # put particle in simulation box
        if part_pos_array[o, 0] < 0:
            part_pos_array[o, 0] += box_length
        if part_pos_array[o, 0] > box_length:
            part_pos_array[o, 0] -= box_length
        if part_pos_array[o, 1] < 0:
            part_pos_array[o, 1] += box_length
        if part_pos_array[o, 1] > box_length:
            part_pos_array[o, 1] -= box_length
        if part_pos_array[o, 2] < 0:
            part_pos_array[o, 2] += box_length
        if part_pos_array[o, 2] > box_length:
            part_pos_array[o, 2] -= box_length
    else:  # if move is rejected return to previous coordinates
        part_pos_array[o, 0] = xo
        part_pos_array[o, 1] = yo
        part_pos_array[o, 2] = zo
    return en, vir, nacc, attempt, part_pos_array
