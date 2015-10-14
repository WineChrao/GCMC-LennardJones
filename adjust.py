__author__ = 'maurizio'


def adjust(attempt, nacc, dr, half_box_length, attemptp, naccp):
    """
    adjusts maximum displacement such that 50% of the
    moves will be accepted
    """
    if (attempt == 0) or (attemptp > attempt):
        naccp = nacc
        attemptp = attempt
    else:
        frac = float(nacc-naccp)/float(attempt-attemptp)
        dro = dr
        dr = dr*abs(frac/0.50)
        # ---limit the change:
        if dr/dro > 1.5:
            dr = dro*1.5
        if dr/dro < 0.5:
            dr = dro*0.5
        if dr > half_box_length/2.0:
            dr = half_box_length/2.0
        naccp = nacc
        attemptp = attempt
    return naccp, attemptp, dr
