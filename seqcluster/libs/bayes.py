"""This file contains code for use with "Think Bayes",
by Allen B. Downey, available from greenteapress.com

Copyright 2012 Allen B. Downey
License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
"""
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


def decide_by_bayes():
    #get common sequences
    #for each
        #get number of loci in common
        #run bayes
        #update proportion of counts
    #returl list of clusters
    hypos = ['Bowl1', 'Bowl2']
    loci = {
    'Bowl1':dict(vanilla=0.75),
    'Bowl2':dict(vanilla=0.5),
    }
    pmf = loci(hypos, loci)
    pmf.Update('vanilla')

    for hypo, prob in pmf.Items():
        print hypo, prob
