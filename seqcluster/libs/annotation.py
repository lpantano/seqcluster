


def _position_in_feature(pos_a, pos_b):
    """return distance to 3' and 5' end of the feature"""
    strd = "-"
    if pos_a[2] in pos_b[2]:
        strd = "+"
    if pos_a[2] in "+" and pos_b[2] in "+":
        lento5 = pos_a[0] - pos_b[1] + 1
        lento3 = pos_a[1] - pos_b[1] + 1
    if pos_a[2] in "+" and pos_b[2] in "-":
        lento5 = pos_a[1] - pos_b[0] + 1
        lento3 = pos_a[0] - pos_b[1] + 1
    if pos_a[2] in "-" and pos_b[2] in "+":
        lento5 = pos_a[0] - pos_b[1] + 1
        lento3 = pos_a[1] - pos_b[0] + 1
    if pos_a[2] in "-" and pos_b[2] in "-":
        lento3 = pos_a[0] - pos_b[0] + 1
        lento5 = pos_a[1] - pos_b[1] + 1
    else:
        lento5 = pos_a[0] - pos_b[0] + 1
        lento3 = pos_a[1] - pos_b[1] + 1
    return lento5, lento3, strd
