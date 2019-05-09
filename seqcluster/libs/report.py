from __future__ import print_function
import logging
import string
import os

from math import log as mlog2
from collections import Counter, defaultdict

from seqcluster.libs.read import map_to_precursors, precursor_sequence, map_to_precursor_biopython
from seqcluster.libs.utils import safe_dirs
from progressbar import ProgressBar

from seqcluster.libs.utils import file_exists
from seqcluster.function.rnafold import run_rnafold
from seqcluster.html import HTML
from seqcluster import templates

logger = logging.getLogger('html')


def _get_ann(dbs, features):
    """
    Gives format to annotation for html table output
    """
    value = ""
    for db, feature in zip(dbs, features):
        value += db + ":" + feature
    return value


def _parse(profile, size):
    total = Counter()
    for sample in profile:
        for pos in range(len(size)):
            total[pos] += profile[sample][pos]
    return total.values()


def make_profile(data, out_dir, args):
    """
    Make data report for each cluster
    """
    safe_dirs(out_dir)
    main_table = []
    header = ['id', 'ann']
    n = len(data[0])
    bar = ProgressBar(maxval=n).start()
    bar.update(0)
    for itern, c in enumerate(data[0]):
        bar.update(itern)
        logger.debug("creating cluser: {}".format(c))
        safe_dirs(os.path.join(out_dir, c))
        valid, ann, pos_structure = _single_cluster(c, data, os.path.join(out_dir, c, "maps.tsv"), args)
        data[0][c].update({'profile': pos_structure})
        loci = data[0][c]['loci']
        logger.debug("precursor_sequence")
        data[0][c]['precursor'] = {"seq": precursor_sequence(loci[0][0:5], args.ref)}
        logger.debug("parse alignments")
        data[0][c]['precursor']["colors"] = list(_parse(data[0][c]['profile'], data[0][c]['precursor']["seq"]))
        logger.debug("update rnafold")
        data[0][c]['precursor'].update(run_rnafold(data[0][c]['precursor']['seq']))

    return data


def _expand(dat, counts, start, end):
    """
    expand the same counts from start to end
    """
    for pos in range(start, end):
        for s in counts:
            dat[s][pos] += counts[s]
    return dat


def _convert_to_df(in_file, freq, raw_file):
    """
    convert data frame into table with pandas
    """
    dat = defaultdict(Counter)
    if isinstance(in_file, (str, bytes)):
        with open(in_file) as in_handle:
            for line in in_handle:
                cols = line.strip().split("\t")
                counts = freq[cols[3]]
                dat = _expand(dat, counts, int(cols[1]), int(cols[2]))
    else:
        if raw_file:
            out_handle = open(raw_file, "w")
        for name in in_file:
            counts = freq[name]
            if raw_file:
                print("%s\t%s\t%s\t%s\t%s\t%s" % ("chr", in_file[name][0], in_file[name][1], name, sum(counts.values()), "+"), file=out_handle, end="")
            dat = _expand(dat, counts, in_file[name][0], in_file[name][1])

    for s in dat:
        for p in dat[s]:
            dat[s][p] = mlog2(dat[s][p] + 1)
    return dat


def _make(c):
    """
    create html from template, adding figure,
    annotation and sequences counts
    """
    ann = defaultdict(list)

    for pos in c['ann']:
        for db in pos:
            ann[db] += list(pos[db])
    logger.debug(ann)

    valid = [l for l in c['valid']]
    ann_list = [", ".join(list(set(ann[feature]))) for feature in ann if feature in valid]

    return valid, ann_list


def _single_cluster(c, data, out_file, args):
    """
    Map sequences on precursors and create
    expression profile
    """
    valid, ann = 0, 0
    raw_file = None
    freq = defaultdict()
    [freq.update({list(s.keys())[0]: list(s.values())[0]}) for s in data[0][c]['freq']]
    names = [list(s.keys())[0] for s in data[0][c]['seqs']]
    seqs = [list(s.values())[0] for s in data[0][c]['seqs']]
    loci = data[0][c]['loci']

    if loci[0][3] - loci[0][2] > 500:
        logger.info("locus bigger > 500 nt, skipping: %s" % loci)
        return valid, ann, {}
    if not file_exists(out_file):
        if args.razer:
            logger.debug("map with razer all sequences to all loci %s " % loci)
            map_to_precursors(seqs, names, {loci[0][0]: [loci[0][0:5]]}, out_file, args)
        else:
            logger.debug("map with biopython fn all sequences to all loci %s " % loci)
            if args.debug:
                raw_file = out_file
            out_file = map_to_precursor_biopython(seqs, names, loci[0][0:5], args)

    logger.debug("plot sequences on loci")
    df = _convert_to_df(out_file, freq, raw_file)
    
    if df:
        logger.debug("create html")
        valid, ann = _make(data[0][c])
    logger.debug("done single cluster")
    return valid, ann, df
