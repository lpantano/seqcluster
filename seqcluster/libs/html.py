import logging
from read import map_to_precursors


logger = logging.getLogger('html')


def _single_cluster(c, data, args):
    seqs, names = data[0][c]['seqs'].values(), data[0][c]['seqs'].keys()
    loci = data[0][c]['loci']
    logger.info("map all sequences to all loci")
    print "%s" % (loci)
    map_to_precursors(seqs, names, loci, args)
    # map_sequences_w_bowtie(sequences, precursors)
    logger.info("plot sequences on loci")
 
