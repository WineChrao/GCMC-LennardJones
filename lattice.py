import numpy as np

__author__ = 'maurizio'
"""
place 'npart' particles on a lattice with density 'rho'
"""


def lattice(box_length, npart):
    part_pos = np.zeros((npart, 3), dtype=np.float)
    n = int(npart**(1.0/3.0))+1
    dist = box_length/float(n)
    part_index = 0
    dx = -dist
    for i in range(n):
        dx += dist
        if dx > box_length:
            dx -= box_length
        dy = -dist
        for j in range(n):
            dy += dist
            if dy > box_length:
                dy -= box_length
            dz = -dist
            for k in range(n):
                dz += dist
                if dz > box_length:
                    dz -= box_length
                if part_index < npart:
                    part_pos[part_index, 0] = dx
                    part_pos[part_index, 1] = dy
                    part_pos[part_index, 2] = dz
                    part_index += 1
    return part_pos



