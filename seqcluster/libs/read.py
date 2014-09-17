"""functions for explore tool"""


def load_data(in_file):
    """load json file from seqcluster cluster"""
    with open(in_file) as in_handle:
        return json.load(in_handle)


def get_sequences_from_cluster(c1, c2, data):
    """get all sequences from on cluster"""
    seqs1 = data[c1]['seqs']
    seqs2 = data[c2]['seqs']
    return seqs1 + seqs2


def get_precursors_from_cluster(c1, c2, data):
    loci1 = data[c1]['loci']
    loci2 = data[c2]['loci']
    return loci1 + loci2


def map_to_precursors(seqs, loci):
	return True


def get_fasta(loci):
    return True


def plot_positions():
	return True
