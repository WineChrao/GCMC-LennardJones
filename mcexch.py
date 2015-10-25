import numpy as np
#import cor
import ccor
# import ener
import cener

__author__ = 'maurizio'


def mcexch(npart, part_pos_array, beta,  en, vir, attempt,
           nacc, box_length, rc, rc2, sigma, sigma2, epsilon4,
           epsilon48, shift_pot, tail_corr, ecut, zz, av_cav, n_cavs, n_nodes, use_cavity_bias):
    attempt += 1
    vol = box_length**3
    rhoo = npart/vol
    # select to add or delete a particle
    if np.random.random() < 0.5:
        # random insertion
        if len(av_cav) == 0 or not use_cavity_bias:
            # add a particle at a random position
            xn = np.random.random()*box_length
            yn = np.random.random()*box_length
            zn = np.random.random()*box_length
            pc = 1
        else:
            # cavity biased insertion
            if len(av_cav) > 1:
                cav_index = np.random.randint(0, len(av_cav) - 1)
            else:
                cav_index = 0
            xn = av_cav[cav_index, 0]
            yn = av_cav[cav_index, 1]
            zn = av_cav[cav_index, 2]
            pc = float(n_cavs)/float(n_nodes)

        o = npart
        jb = 0
        part_pos_array = np.append(part_pos_array, [[xn, yn, zn]], axis=0)
        # determine energy of this particle
        enn, virn = cener.eneri(npart, part_pos_array, o, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
        # tail correction
        if tail_corr:
            rhon = (npart+1)/vol
            enn += ((npart+1)*ccor.coru(rc, rhon, sigma, epsilon4)-npart*ccor.coru(rc, rhoo, sigma, epsilon4))
        # acceptance test:
        arg = zz*vol*pc*np.exp(-beta*enn)/(npart+1)
        if np.random.random() < min(1, arg):
            # accepted
            nacc += 1
            en += enn
            vir += virn
            npart += 1
        else:
            part_pos_array = np.delete(part_pos_array, o, 0)
    else:
        if np.random.random() < (1 - (float(n_cavs) - 1) / float(n_nodes))**n_nodes:
            pc = 1
        else:
            pc = (float(n_cavs) - 1) / float(n_nodes)
        # delete a randomly selected particle
        o = int(npart*np.random.random())
        jb = 0
        eno, viro = cener.eneri(npart, part_pos_array, o, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
        # particle is removed, so new energy
        enn = -eno
        virn = -viro
        # tail correction
        if tail_corr:
            rhon = (npart-1)/vol
            enn += ((npart-1)*ccor.coru(rc, rhon, sigma, epsilon4)-npart*ccor.coru(rc, rhoo, sigma, epsilon4))
        # acceptance test:
        arg = npart*np.exp(-beta*enn)/(zz*pc*vol)
        if np.random.random() < min(1, arg):
            # accepted
            nacc += 1
            en += enn
            vir += virn
            npart -= 1
            part_pos_array = np.delete(part_pos_array, o, 0)

    return en, vir, nacc, attempt, part_pos_array, npart