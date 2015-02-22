"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""
from collections import defaultdict
from thinkbayes import Pmf



class loci(Pmf):
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


def _normalize(loci, loci_obj):
    """
    normalized number of sequences in that position
    """
    return 0


def _dict_seq_locus(list_c):
    """
    return dict with sequences = [locus1. locus2 ...]
    """
    for c in list_c.values():
        for l in c.loci2seq:
            fake =0 


def decide_by_bayes(list_c, loci_obj):
    # for each cluster get seq ~ loci edges
    loci = _dict_seq_locus(list_c)
    loci = _normalize(loci, loci_obj)
    data = defaultdict(dict)
    hypos = loci.keys()
    for hypo in hypos:
        data[hypo] = dict(position=loci[hypo])
    pmf = loci(hypos, data)
    pmf.Update('position')
    return pmf
    # for hypo, prob in pmf.Items():
    #    print hypo, prob
