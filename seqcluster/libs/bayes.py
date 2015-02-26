"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""
from collections import defaultdict
from thinkbayes import Pmf
import numpy as np


class _update(Pmf):
    """A map from string bowl ID to probablity."""

    def __init__(self, hypos, loci):
        """Initialize self.
        hypos: sequence of string bowl IDs
        """
        self.loci = loci
        Pmf.__init__(self)
        for hypo in hypos:
            self.Set(hypo, 1)
            self.Normalize()

    def Update(self, data):
        """Updates the PMF with new data.
        data: string cookie type
        """
        for hypo in self.Values():
            like = self.Likelihood(data, hypo)
            self.Mult(hypo, like)
        self.Normalize()

    def Likelihood(self, data, hypo):
        """The likelihood of the data under the hypothesis.
        data: string cookie type
        hypo: string bowl ID
        """
        mix = self.loci[hypo]
        like = mix[data]
        return like


def _transform(seqs):
    seqs_in_c = defaultdict(dict)
    for s, c in seqs:
        seqs_in_c[s].update({c: seqs[(s, c)]})
    return seqs_in_c


def _dict_seq_locus(list_c, loci_obj, seq_obj):
    """
    return dict with sequences = [ cluster1, cluster2 ...]
    """
    seqs = defaultdict(set)
    # n = len(list_c.keys())
    for c in list_c.values():
        for l in c.loci2seq:
            [seqs[s].add(c.id) for s in c.loci2seq[l]]

    common = [s for s in seqs if len(seqs[s]) > 1]
    seqs_in_c = defaultdict(float)
    for c in list_c.values():
        for l in c.loci2seq:
            # total = sum([v for v in loci_obj[l].coverage.values()])
            for s in c.loci2seq[l]:
                if s in common:
                    pos = seq_obj[s].pos[l]
                    # cov = 1.0 * loci_obj[l].coverage[pos] / total
                    cov = 1.0 * loci_obj[l].coverage[pos]
                    if seqs_in_c[(s, c.id)] < cov:
                        seqs_in_c[(s, c.id)] = cov
    seqs_in_c = _transform(seqs_in_c)
    return seqs_in_c


def _bayes(loci):
    data = defaultdict(dict)
    hypos = loci.keys()
    for hypo in hypos:
        data[hypo] = dict(position=loci[hypo])
    pmf = _update(hypos, data)
    pmf.Update('position')
    return pmf


def decide_by_bayes(list_c, s2p):
    # for each cluster get seq ~ loci edges
    loci_obj, seq_obj = s2p
    seqs_in_c = _dict_seq_locus(list_c, loci_obj, seq_obj)
    for s in seqs_in_c:
        total = sum(seqs_in_c[s].values())
        norm_values = 1.0 * np.array(seqs_in_c[s].values()) / total
        prob = _bayes(dict(zip(seqs_in_c[s].keys(), norm_values)))
        for clus in prob:
            list_c[clus].idmembers[s] = prob.loci[clus]["position"]
            # print "%s %s %s" % (s, clus, list_c[clus].idmembers[s])
    return list_c
