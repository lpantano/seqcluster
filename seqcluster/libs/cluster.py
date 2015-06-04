import os.path as op
import pysam

import logger as mylog
from classes import *
from tool import _get_seqs_from_cluster


logger = mylog.getLogger(__name__)

def clean_bam_file(bam_in, seqs_list):
    """
    Remove from alignment reads with low counts and highly # of hits
    """
    out_file = op.splitext(bam_in)[0] + "_rmlw.bam"
    pysam.index(bam_in)
    bam = pysam.AlignmentFile(bam_in, "rb")
    with pysam.AlignmentFile(out_file, "wb", template=bam) as out_handle:
        for read in bam.fetch():
            seq_name = read.query_name
            nh = read.get_tag('NH')
            ratio = seqs_list[seq_name].total() / float(nh)
            if ratio > 0.1:
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
    for c in c_index:
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


