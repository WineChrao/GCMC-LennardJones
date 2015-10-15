import cor
__author__ = 'maurizio'


def ener(r2, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut):
    """
    calculates energy (en) and virial (vir) for given
    distance squared between (r2) two particles
    """
    if r2 <= rc2:
        r2i = sigma2/r2
        r6i = r2i*r2i*r2i
        if shift_pot:
            en = epsilon4*(r6i*r6i-r6i) - ecut
        else:
            en = epsilon4*(r6i*r6i-r6i)

        vir = epsilon48*(r6i*r6i-0.5*r6i)
    else:
        en = 0.0
        vir = 0.0
    return en, vir


def totenerg(npart, part_pos_array, box_length, rc, sigma, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut, rho, tail_corr):
    en = 0.0
    vir = 0.0
    new_rho = rho
    for i in range(npart):
        jb = i + 1
        # calculate energy particle i with particle j=jb,npart
        eni, viri = eneri(npart, part_pos_array, i, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
        en += eni
        vir += viri
    if tail_corr:
        rho = npart/float(box_length**3)
        en += npart * cor.coru(rc, rho, sigma, epsilon4)
    return en, vir, new_rho


def eneri(npart, part_pos_array, i, jb, box_length, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut):
    en = 0.0
    vir = 0.0
    for j in range(jb, npart):
        if j != i:
            dx = part_pos_array[i, 0] - part_pos_array[j, 0]
            dy = part_pos_array[i, 1] - part_pos_array[j, 1]
            dz = part_pos_array[i, 2] - part_pos_array[j, 2]
            if dx > box_length/2.0:
                dx -= box_length
            elif dx < -box_length/2.0:
                dx += box_length
            if dy > box_length/2.0:
                dy -= box_length
            elif dy < -box_length/2.0:
                dy += box_length
            if dz > box_length/2.0:
                dz -= box_length
            elif dz < -box_length/2.0:
                dz += box_length
            r2 = dx*dx + dy*dy + dz*dz
            # calculate energy and virial of pair i,j
            if r2 > 0.0:
                enij, virij = ener(r2, rc2, sigma2, epsilon4, epsilon48, shift_pot, ecut)
                en += enij
                vir += virij
    return en, vir

