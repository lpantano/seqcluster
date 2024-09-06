try:
    import scipy.stats as stat
except:
    pass

def up_threshold(x, s, p):
    """function to decide if similarity is
    below cutoff"""
    if 1.0 * x/s >= p:
        return True
    elif stat.binomtest(int(x), int(s), p).pvalue > 0.01:
        return True
    return False
