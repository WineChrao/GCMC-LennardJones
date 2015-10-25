import numpy as np
cimport numpy as np


FLOATTYPE = np.float
ctypedef np.float_t FLOATTYPE_t
INTTYPE = np.int
ctypedef np.int_t INTTYPE_t

def cavity_locations(double box_length, int n_cavs):
    cdef int n
    cdef double dist
    cdef int cav_index
    cdef double dx, dy, dz
    cdef np.ndarray cavities_pos = np.zeros((n_cavs, 3), dtype=FLOATTYPE)
    n = int(n_cavs**(1.0/3.0))+1
    dist = box_length/float(n)
    cav_index = 0
    dx = 0
    for i in range(n):
        dx += dist
        dy = 0
        for j in range(n):
            dy += dist
            dz = 0
            for k in range(n):
                dz += dist
                if cav_index < n_cavs:
                    cavities_pos[cav_index, 0] = dx
                    cavities_pos[cav_index, 1] = dy
                    cavities_pos[cav_index, 2] = dz
                    cav_index += 1
    return cavities_pos


def available_cavities(np.ndarray[np.float_t, ndim=2] cavities_pos, np.ndarray[np.float_t, ndim=2] part_pos, double cut_off):
    cdef np.ndarray list_to_remove = np.array([], dtype=INTTYPE)
    cdef np.ndarray av_cav = np.array([], dtype=FLOATTYPE)
    cdef int i, j, remove_and_break
    cdef double dx, dy, dz 
    cdef int n_cavs
    for i in range(len(cavities_pos)):
        remove_and_break = 0
        for j in range(len(part_pos)):
            dx = cavities_pos[i, 0] - part_pos[j, 0]
            dy = cavities_pos[i, 1] - part_pos[j, 1]
            dz = cavities_pos[i, 2] - part_pos[j, 2]
            r2 = dx*dx + dy*dy + dz*dz
            if r2 < cut_off*cut_off:
                remove_and_break = 1
                break
        if remove_and_break:
            list_to_remove = np.append(list_to_remove, [i], axis=0)

    av_cav = np.delete(cavities_pos, list_to_remove, 0)
    n_cavs = len(av_cav)
    return av_cav, n_cavs
