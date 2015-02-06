import logging
import os
# from collections import Counter

import matplotlib
matplotlib.use('Agg')
import pylab
pylab.rcParams['figure.figsize'] = (35.0, 12.0)
import pandas as pd

from read import map_to_precursors
from utils import safe_dirs


logger = logging.getLogger('html')


def make_profile(data, out_dir, args):
    """
    Make html for each cluster
    """
    for c in data[0]:
        logger.debug("creating cluser: {}".format(c))
        safe_dirs(os.path.join(out_dir, c))
        _single_cluster(c, data, os.path.join(out_dir, c, "maps.tsv"), args)


def _convert_to_df(in_file):
    dat = {}
    with open(in_file) as in_handle:
        for line in in_handle:
            cols = line.strip().split(" ")
            counts = int(cols[1].replace("cx", ""))
            dat[4] = counts
    df = pd.data_frame(dat)
    return df


def _single_cluster(c, data, out_file, args):
    # dat_freq = Counter()
    # for s in data[0][c]['freq']:
    #     dat_freq[s.keys()[0]] = round(sum(s.values()[0].values()))

    names = [round(sum(s.values()[0].values())) for s in data[0][c]['freq']]
    seqs = [s.values()[0] for s in data[0][c]['seqs']]

    loci = data[0][c]['loci']

    if loci[0][3] - loci[0][2] > 500:
        logger.info("locus bigger > 500 nt, skipping: %s" % loci)
        return False
    logger.info("map all sequences to all loci %s " % loci)
    map_to_precursors(seqs, names, {loci[0][0]: [loci[0]]}, out_file, args)
    # map_sequences_w_bowtie(sequences, precursors)

    logger.info("plot sequences on loci")
    df = _convert_to_df(out_file)
    plot = graph.plot()
    plot.set_ylabel('% CPU')
