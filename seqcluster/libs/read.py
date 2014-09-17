"""functions for explore tool"""
import json
from do import find_cmd, run

def load_data(in_file):
    """load json file from seqcluster cluster"""
    with open(in_file) as in_handle:
        return json.load(in_handle)


def get_sequences_from_cluster(c1, c2, data):
    """get all sequences from on cluster"""
    seqs1 = data[c1]['seqs']
    seqs2 = data[c2]['seqs']
    seqs = list(set(seqs1 + seqs2))
    names = []
    for s in seqs:
        if s in seqs1 and s in seqs2:
            names.append("both")
        elif s in seqs1:
            names.append(c1)
        else:
            names.append(c2)
    return seqs, names


def get_precursors_from_cluster(c1, c2, data):
    loci1 = data[c1]['loci']
    loci2 = data[c2]['loci']
    return loci1, loci2


def map_to_precursors(seqs, loci):
    """map sequences to precursors with bowtie"""
    if find_cmd("bowtie"):
        cmd = "bowtie -f {fasta} {index} -S"
        run(cmd.format(**locals()))
    return True


def get_fasta(loci):
    return True


def plot_positions():
	return True
