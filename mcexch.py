import numpy as np
import cor
import ener

__author__ = 'maurizio'


def mcexch(npart, part_pos_array, beta,  en, vir, attempt,
           nacc, box_length, rc, rc2, sigma, sigma2, epsilon4,
           epsilon48, shift_pot, tail_corr, ecut, zz):
    attempt += 1
    vol = box_length**3
    rhoo = npart/vol
    # select to add of delete a particle
    if np.random.random() < 0.5:
        # add a particle at a random position
        xn = np.random.random()*box_length
        yn = np.random.random()*box_length
        zn = np.random.random()*box_length
        o = npart
        jb = 0
        part_pos_array = np.append(part_pos_array, [[xn, yn, zn]], axis=0)
        # determine energy of this particle
        enn, virn = ener.eneri(npart, part_pos_array, o, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
        # tail correction
        if tail_corr:
            rhon = (npart+1)/vol
            enn += ((npart+1)*cor.coru(rc, rhon, sigma, epsilon4)-npart*cor.coru(rc, rhoo, sigma, epsilon4))
        # acceptance test:
        arg = zz*vol*np.exp(-beta*enn)/(npart+1)
        if np.random.random() < arg:
            # accepted
            nacc += 1
            en += enn
            vir += virn
            npart += 1
        else:
            part_pos_array = np.delete(part_pos_array, o, 0)
    else:
        # delete a randomly selected particle
        o = int(npart*np.random.random())
        jb = 0
        eno, viro = ener.eneri(npart, part_pos_array, o, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
        # particle is removed, so new energy
        enn = -eno
        virn = -viro
        # tail correction
        if tail_corr:
            rhon = (npart-1)/vol
            enn += ((npart-1)*cor.coru(rc, rhon, sigma, epsilon4)-npart*cor.coru(rc, rhoo, sigma, epsilon4))
        # acceptance test:
        arg = npart*np.exp(-beta*enn)/(zz*vol)
        if np.random.random() < arg:
            # accepted
            en += enn
            vir += virn
            npart -= 1
            part_pos_array = np.delete(part_pos_array, o, 0)


    return en, vir, nacc, attempt, part_pos_array, npart