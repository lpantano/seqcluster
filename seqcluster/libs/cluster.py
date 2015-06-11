import os.path as op
import pickle
from progressbar import ProgressBar

import pysam
import pybedtools
from pypeaks import Data
import numpy as np

from bcbio.utils import file_exists

import logger as mylog
from classes import *
from seqcluster.function.peakdetect import peakdetect as peakdetect
from tool import _get_seqs_from_cluster


logger = mylog.getLogger(__name__)


def detect_complexity(bam_in, genome):
    """
    genome coverage of small RNA
    """
    out_file = bam_in + "_cov.tsv"
    if not genome or file_exists(out_file):
        return None
    fai = genome + ".fai"
    cov = pybedtools.BedTool(bam_in).genome_coverage(g=fai, max=1)
    cov.saveas()
    total = 0
    for region in cov:
        if region[0] == "genome" and int(region[1]) != 0:
            total += float(region[4])
    logger.info("Total genome with sequences: %s " % total)

def clean_bam_file(bam_in, seqs_list):
    """
    Remove from alignment reads with low counts and highly # of hits
    """
    out_file = op.splitext(bam_in)[0] + "_rmlw.bam"
    if file_exists(out_file):
        return out_file
    pysam.index(bam_in)
    bam = pysam.AlignmentFile(bam_in, "rb")
    with pysam.AlignmentFile(out_file, "wb", template=bam) as out_handle:
        for read in bam.fetch():
            seq_name = read.query_name
            try:
                nh = read.get_tag('NH')
            except ValueError:
                nh = 1
                continue
            ratio = seqs_list[seq_name].total() / float(nh)
            if ratio > 0.01:
                out_handle.write(read)
    return out_file

def detect_clusters(c, current_seq, MIN_SEQ):
    """
    Parse the merge file of sequences position to create clusters that will have all
    sequences that shared any position on the genome

    :param c: file from bedtools with merge sequence positions
    :param current_seq: list of sequences
    :param MIN_SEQ: int cutoff to keep the cluster or not. 10 as default

    :return: object with information about:
        * cluster
        * dict with sequences (as keys) and cluster_id (as value)
        * sequences
        * loci

    """
    current_loci = {}
    current_clus = {}
    lindex = 0
    eindex = 0
    previous_id = 0
    for line in c.features():
        c, start, end, name, score, strand, c_id = line
        pos = start if strand == "+" else end
        if c_id != previous_id:
            logger.debug("detect_cluster: %s %s %s" % (c_id, previous_id, name))
            lindex += 1
            eindex += 1
            current_clus[eindex] = cluster(eindex)
            newpos = position(lindex, c, start, end, strand)
            current_loci[lindex] = newpos

        # update locus, sequences in each line
        current_loci[lindex].end = int(end)
        current_loci[lindex].coverage[pos] += 1
        current_clus[eindex].idmembers[name] = 1
        current_clus[eindex].add_id_member([name], lindex)
        current_seq[name].add_pos(lindex, pos)
        previous_id = c_id
    logger.debug("%s clusters read" % eindex)
    # merge cluster with shared sequences  
    cluster_obj, cluster_id = _find_families(current_clus, MIN_SEQ)

    return cluster_info_obj(current_clus, cluster_id, current_loci, current_seq)

def _find_families(clus_obj, min_seqs):
    """
    Mask under same id all clusters that share sequences
    :param clus_obj: cluster object coming from detect_cluster
    :param min_seqs: int cutoff to keep the cluster or not. 10 as default

    :return: updated clus_obj and dict with seq_id: cluster_id
    """
    seen = {}
    c_index = clus_obj.keys()
    with ProgressBar(maxval=len(c_index), redirect_stdout=True) as p:
        for itern, c in enumerate(c_index):
            p.update(itern)
            clus = clus_obj[c]
            if len(clus.idmembers.keys()) < min_seqs:
                del clus_obj[c]
                continue
            logger.debug("reading cluster %s" % c)
            logger.debug("loci2seq  %s" % clus.loci2seq)
            already_in, not_in = _get_seqs_from_cluster(clus.idmembers.keys(), seen)
            logger.debug("seen %s news %s" % (already_in, not_in))
            for s in not_in:
                seen[s] = c
            if len(already_in) > 0:
                logger.debug("seen in %s" % already_in)
                for eindex in already_in:
                    prev_clus = clus_obj[eindex]
                    logger.debug("_find_families: prev %s current %s" % (eindex, clus.id))
                    # add current seqs to seen cluster
                    for s_in_clus in prev_clus.idmembers:
                        seen[s_in_clus] = c
                        clus.idmembers[s_in_clus] = 1
                    # add current locus to seen cluster
                    for loci in prev_clus.loci2seq:
                        logger.debug("adding %s" % loci)
                        # if not loci_old in current_clus[eindex].loci2seq:
                        clus.add_id_member(list(prev_clus.loci2seq[loci]), loci)
                    logger.debug("loci %s" % clus.loci2seq.keys())
                    del clus_obj[eindex]
                clus_obj[c] = clus
                logger.debug("num cluster %s" % len(clus_obj.keys()))
    logger.info("%s clusters merged" % len(clus_obj.keys()))

    return clus_obj, seen

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
        s, e, st = clus_obj.loci[bigger].start, clus_obj.loci[bigger].end, clus_obj.loci[bigger].strand
        scale = min(s, e)
        logger.debug("bigger %s at %s-%s" % (bigger, s, e))
        seqs = cluster.loci2seq[bigger]
        dt = np.array([0] * (e - s + 12))
        for seq in seqs:
            ss = int(clus_obj.seq[seq].pos[bigger]) - scale + 5
            l = clus_obj.seq[seq].len
            se = ss + l
            if st == "-":
                se, ss = ss, ss - l
            dt[ss:se] += clus_obj.seq[seq].total()
        x = np.array(range(0, len(dt)))
        logger.debug("x %s and y %s" % (x, dt))

        if len(x) > 35 + 12:
            profile = Data(x, dt, smoothness=3)
            pickle.dump(profile, open("xy.pickle", 'w'))
            profile.normalize()
            profile.get_peaks(method='slope', lookahead=4, avg_interval=20)
            peaks = list(profile.peaks['peaks'][0])
            # print profile.peaks
            # peaks = peakdetect(dt, x, lookahead=5)[0]
            logger.debug(peaks)
        else:
            peaks = [[s], [e], ['short']]
        cluster.peaks = peaks
        new_cluster[cid] = cluster
    clus_obj.clus = new_cluster
    return clus_obj
