from operator import itemgetter

from  seqcluster.libs import pysen
import numpy as np

import seqcluster.libs.logger as mylog
from seqcluster.libs.classes import *
# from seqcluster.function.peakdetect import peakdetect as peakdetect


logger = mylog.getLogger(__name__)


def sort_precursor(c, loci):
    """
    Sort loci according to number of sequences mapped there.
    """
    # Original Py 2.7 code
    #data_loci = map(lambda (x): [x, loci[x].chr, int(loci[x].start), int(loci[x].end), loci[x].strand, len(c.loci2seq[x])], c.loci2seq.keys())
    # 2to3 suggested Py 3 rewrite
    data_loci = [[x, loci[x].chr, int(loci[x].start), int(loci[x].end), loci[x].strand, len(c.loci2seq[x])] for x in list(c.loci2seq.keys())]
    data_loci = sorted(data_loci, key=itemgetter(5), reverse=True)
    return data_loci

def best_precursor(clus, loci):
    """
    Select best precursor asuming size around 100 nt
    """
    data_loci = sort_precursor(clus, loci)
    current_size = data_loci[0][5]
    best = 0
    for item, locus in enumerate(data_loci):
        if locus[3] - locus[2] > 70:
            if locus[5] > current_size * 0.8:
                best = item
                break
    best_loci = data_loci[best]
    del data_loci[best]
    data_loci.insert(0, best_loci)
    return data_loci


def peak_calling(clus_obj):
    """
    Run peak calling inside each cluster
    """
    new_cluster = {}
    for cid in clus_obj.clus:
        cluster = clus_obj.clus[cid]
        cluster.update()
        logger.debug("peak calling for %s" % cid)
        bigger = cluster.locimaxid
        # fn to get longer > 30 < 100 precursor
        # iterate by them mapping reads, detect peak, descart, start_again
        # give quantification of 'matures'
        if bigger in clus_obj.loci:
            s, e = min(clus_obj.loci[bigger].counts.keys()), max(clus_obj.loci[bigger].counts.keys())
            scale = min(s, e)
            logger.debug("bigger %s at %s-%s" % (bigger, s, e))
            dt = np.array([0] * (e - s + 12))
            for pos in clus_obj.loci[bigger].counts:
                ss = int(pos) - scale + 5
                dt[ss] += clus_obj.loci[bigger].counts[pos]
        x = np.array(range(0, len(dt)))
        logger.debug("x %s and y %s" % (x, dt))
        # tab = pd.DataFrame({'x': x, 'y': dt})
        # tab.to_csv( str(cid) + "peaks.csv", mode='w', header=False, index=False)
        if len(x) > 35 + 12:
            peaks = pysen.pysenMMean(x, dt)
            logger.debug(peaks)
        else:
            peaks =  ['short']
        cluster.peaks = peaks
        new_cluster[cid] = cluster
    clus_obj.clus = new_cluster
    return clus_obj
