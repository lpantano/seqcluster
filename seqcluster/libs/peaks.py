"""detect peaks inside clusters"""


def _scan(positions):
    """get the region inside the vector with more expression"""
    scores = []
    for start in range(0, len(positions) - 17, 5):
        end = start = 17
        scores.add(_enrichment(positions[start:end], positions[:start], positions[end:]))
    #get index max score (check ties)
    #return enriched region(give up to two region)

def _enrichment(resgion, flank_a, flank_b):
    return True

def _get_locus(cluster):
    """get the bigger locus"""
    return True

def _get_position_in_loci(locus_id, seqs, loci):
    """get position in locus of all sequences"""
    return True
