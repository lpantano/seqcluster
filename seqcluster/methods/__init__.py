"""Functions to access json with classes"""
from collections import defaultdict

from seqcluster.libs.classes import cluster

def read_cluster(data, id=1):
    """Read json cluster and populate as cluster class"""
    cl = cluster(1)

    # seqs = [s.values()[0] for s in data['seqs']]
    names = [s.keys()[0] for s in data['seqs']]
    cl.add_id_member(names, 1)
    freq = defaultdict()
    [freq.update({s.keys()[0]: s.values()[0]}) for s in data['freq']]

