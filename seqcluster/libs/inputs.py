from collections import defaultdict
import pybedtools
from classes import sequence
from tool import _normalize_seqs


def parse_align_file(file_in):
    """parse sam files with aligned sequences"""
    loc_id = 1
    bedfile_clusters = ""
    bamfile = pybedtools.BedTool(file_in)
    bed = pybedtools.BedTool.bam_to_bed(bamfile)
    for c, start, end, name, q, strand in bed:
        loc_id += 1
        bedfile_clusters += "%s\t%s\t%s\t%s\t%s\t%s\n" % \
                            (c, start, end, name, loc_id, strand)
    return bedfile_clusters


def parse_ma_file(in_file):
    """read seqs.ma file and create dict with
    sequence object"""
    name = ""
    seq_l = {}
    index = 1
    total = defaultdict(int)
    with open(in_file) as handle_in:
        line = handle_in.readline().strip()
        cols = line.split("\t")
        samples = cols[2:]
        for line in handle_in:
            line = line.strip()
            cols = line.split("\t")
            exp = {}
            for i in range(len(samples)):
                exp[samples[i]] = int(cols[i+2])
                total[samples[i]] += int(cols[i+2])
            name = cols[0].replace(">", "")
            index = index+1
            seq = cols[1]
            new_s = sequence(seq, exp, index)
            seq_l[name] = new_s
    seq_l = _normalize_seqs(seq_l, total)
    return seq_l


