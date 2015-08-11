from collections import defaultdict
import pybedtools

import logger as mylog
from classes import sequence
from tool import _normalize_seqs


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
            index = index+1
            if name in seq_obj:
                seq_obj[name].set_freq(exp)
                seq_obj[name].set_seq(seq)
            # new_s = sequence(seq, exp, index)
            # seq_l[name] = new_s
    seq_obj = _normalize_seqs(seq_obj, total)
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
            index = index+1
            if name not in seq_obj:
                seq_obj[name] = sequence(name)
            seq_obj[name].set_freq(exp)
            seq_obj[name].set_seq(seq)
    seq_obj = _normalize_seqs(seq_obj, total)
    return seq_obj, total, index
