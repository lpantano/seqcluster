from collections import defaultdict
import pybedtools
import numpy as np
import pandas as pd

import seqcluster.libs.logger as mylog
from seqcluster.libs.classes import sequence
from seqcluster.libs.tool import _normalize_seqs


logger = mylog.getLogger(__name__)


def parse_align_file(file_in):
    """
    Parse sam files with aligned sequences
    """
    loc_id = 1
    bedfile_clusters = ""
    bamfile = pybedtools.BedTool(file_in)
    bed = pybedtools.BedTool.bam_to_bed(bamfile)
    for c, start, end, name, q, strand in bed:
        loc_id += 1
        bedfile_clusters += "%s\t%s\t%s\t%s\t%s\t%s\n" % \
                            (c, start, end, name, loc_id, strand)
    return bedfile_clusters


def parse_ma_file(seq_obj, in_file):
    """
    read seqs.ma file and create dict with
    sequence object
    """
    name = ""
    index = 1
    total = defaultdict(int)
    ratio = list()
    with open(in_file) as handle_in:
        line = handle_in.readline().strip()
        cols = line.split("\t")
        samples = cols[2:]
        for line in handle_in:
            line = line.strip()
            cols = line.split("\t")
            name = int(cols[0].replace("seq_", ""))
            seq = cols[1]
            exp = {}
            for i in range(len(samples)):
                exp[samples[i]] = int(cols[i+2])
                total[samples[i]] += int(cols[i+2])
            ratio.append(np.array(list(exp.values())) / np.mean(list(exp.values())))
            index = index+1
            if name in seq_obj:
                seq_obj[name].set_freq(exp)
                seq_obj[name].set_seq(seq)
            # new_s = sequence(seq, exp, index)
            # seq_l[name] = new_s
    df = pd.DataFrame(ratio)
    df = df[(df.T != 0).all()]
    size_factor = dict(zip(samples, df.median(axis=0)))
    seq_obj = _normalize_seqs(seq_obj, size_factor)
    return seq_obj, total, index


def parse_ma_file_raw(in_file):
    """
    read seqs.ma file and create dict with
    sequence object
    """
    name = ""
    index = 1
    total = defaultdict(int)
    seq_obj = defaultdict(sequence)
    ratio = list()
    with open(in_file) as handle_in:
        line = handle_in.readline().strip()
        cols = line.split("\t")
        samples = cols[2:]
        for line in handle_in:
            line = line.strip()
            cols = line.split("\t")
            name = int(cols[0].replace("seq_", ""))
            seq = cols[1]
            exp = {}
            for i in range(len(samples)):
                exp[samples[i]] = int(cols[i+2])
                total[samples[i]] += int(cols[i+2])
            ratio.append(np.array(list(exp.values())) / np.mean(list(exp.values())))
            index = index+1
            if name not in seq_obj:
                seq_obj[name] = sequence(name)
            seq_obj[name].set_freq(exp)
            seq_obj[name].set_seq(seq)
    df = pd.DataFrame(ratio)
    df = df[(df.T != 0).all()]
    size_factor = dict(zip(samples, df.median(axis=0)))
    seq_obj = _normalize_seqs(seq_obj, size_factor)
    return seq_obj, total, index
