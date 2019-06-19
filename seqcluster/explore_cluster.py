from __future__ import print_function
#import sys
import os
#from os import listdir
#from os.path import isfile, join
import logging
from seqcluster.libs.read import load_data, get_sequences_from_cluster, map_to_precursors, get_precursors_from_cluster


logger = logging.getLogger('explore')


def explore(args):
    """Create mapping of sequences of two clusters
    """
    logger.info("Reading sequences")
    data = load_data(args.json)
    logger.info("Get sequences from json")
    #get_sequences_from_cluster()
    c1, c2 = args.names.split(",")
    seqs, names = get_sequences_from_cluster(c1, c2, data[0])
    loci = get_precursors_from_cluster(c1, c2, data[0])
    logger.info("Map all sequences to all loci")
    print("%s" % (loci))
    map_to_precursors(seqs, names, loci, os.path.join(args.out, "map.tsv"), args)
    #map_sequences_w_bowtie(sequences, precursors)
    logger.info("Plot sequences on loci")
    #get_matrix_position()
    #plot_sequences()
    logger.info("Done")
