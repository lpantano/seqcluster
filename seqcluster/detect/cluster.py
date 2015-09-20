import os.path as op
from progressbar import ProgressBar

import pysam
import pybedtools
from  seqcluster.libs import pysen
import numpy as np

from bcbio.utils import file_exists

import seqcluster.libs.logger as mylog
from seqcluster.libs.classes import *
# from seqcluster.function.peakdetect import peakdetect as peakdetect
from seqcluster.detect.metacluster import _get_seqs_from_cluster
from seqcluster.libs.do import run


logger = mylog.getLogger(__name__)


def detect_complexity(bam_in, genome):
    """
    genome coverage of small RNA
    """
    out_file = bam_in + "_cov.tsv"
    if file_exists(out_file):
        return None
    fai = genome + ".fai"
    cov = pybedtools.BedTool(bam_in).genome_coverage(g=fai, max=1)
    cov.saveas(out_file)
    total = 0
    for region in cov:
        if region[0] == "genome" and int(region[1]) != 0:
            total += float(region[4])
    logger.info("Total genome with sequences: %s " % total)

def clean_bam_file(bam_in, mask=None):
    """
    Remove from alignment reads with low counts and highly # of hits
    """
    seq_obj = defaultdict(int)
    if mask:
        mask_file = op.splitext(bam_in)[0] + "_mask.bam"
        if not file_exists(mask_file):
            pybedtools.BedTool(bam_file).intersect(b=mask, v=True).saveas(mask_file)
        bam_in = mask_file
    out_file = op.splitext(bam_in)[0] + "_rmlw.bam"
    # bam.index(bam_in, {'algorithm':{}})
    run("samtools index %s" % bam_in)
    if not file_exists(bam_in + ".bai"):
        raise IOError("Failed to created bam index of %s. Try to do it manually" % bam_in)
    bam_handle = pysam.AlignmentFile(bam_in, "rb")
    with pysam.AlignmentFile(out_file, "wb", template=bam_handle) as out_handle:
        for read in bam_handle.fetch():
            seq_name = int(read.query_name.replace('seq_', ''))
            match_size = [nts for oper, nts in read.cigartuples if oper == 0]
            subs_size = [nts for oper, nts in read.cigartuples if oper == 4]
            if match_size[0] < 17:
                continue
            if subs_size:
                if subs_size[0] > 3:
                    continue
            try:
                nh = read.get_tag('NH')
            except KeyError:
                nh = 1
            seq_obj[seq_name] = sequence(seq_name)
            seq_obj[seq_name].align = nh
            out_handle.write(read)
    return out_file, seq_obj

def detect_clusters(c, current_seq, MIN_SEQ, non_un_gl=False):
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
    # sequence2clusters = [set()] * (max(current_seq.keys()) + 2)
    sequence2clusters = defaultdict(set)
    lindex = 0
    eindex = 0
    previous_id = 0
    for line in c.features():
        c, start, end, name, score, strand, c_id = line
        name = int(name.replace('seq_', ''))
        pos = int(start) if strand == "+" else int(end)
        if name not in current_seq:
            continue
        if c.find('Un_gl') > -1 and non_un_gl:
            continue
        if c_id != previous_id:
            if previous_id > 0:
                if len(current_clus[eindex].idmembers) < MIN_SEQ:
                    for s in current_clus[eindex].idmembers:
                        sequence2clusters[s] = sequence2clusters[s] - set([eindex])
                    del current_clus[eindex]

            logger.debug("detect_cluster: %s %s %s" % (c_id, previous_id, name))
            lindex += 1
            eindex += 1
            current_clus[eindex] = cluster(eindex)
            newpos = position(lindex, c, start, end, strand)
            current_loci[lindex] = newpos

        # update locus, sequences in each line
        current_loci[lindex].end = int(end)
        current_loci[lindex].coverage[pos] += 1
        size = range(pos, pos + current_seq[name].len)
        current_loci[lindex].counts.update(dict(zip(size, [current_seq[name].total()] * current_seq[name].len)))
        current_clus[eindex].idmembers[name] = 1
        current_clus[eindex].add_id_member([name], lindex)
        current_seq[name].add_pos(lindex, pos)
        # current_seq[name].align = 1
        previous_id = c_id
        sequence2clusters[name].add(eindex)
    logger.info("%s Clusters read" % eindex)
    # merge cluster with shared sequences  
    metacluster_obj, cluster_id = _find_metaclusters(current_clus, sequence2clusters, current_seq, MIN_SEQ)

    return cluster_info_obj(current_clus, metacluster_obj, current_loci, current_seq)

def _common(items, seen):
    intersect = map(seen.get, items)
    return filter(None, intersect)

def _update(clusters, idx, hash):
    return hash.update(dict(zip(clusters, [idx] * len(clusters))))

def _find_metaclusters(clus_obj, sequence2clusters, current_seq, min_seqs):
    """
    Mask under same id all clusters that share sequences
    :param clus_obj: cluster object coming from detect_cluster
    :param min_seqs: int cutoff to keep the cluster or not. 10 as default

    :return: updated clus_obj and dict with seq_id: cluster_id
    """
    seen = defaultdict(int)
    metacluster = defaultdict(set)
    c_index = len(sequence2clusters)
    logger.info("Creating meta-clusters based on shared sequences: %s" % c_index)
    meta_idx = 1
    bar = ProgressBar(maxval=c_index)
    bar.start()
    bar.update()
    for itern, name in enumerate(sequence2clusters):
        clusters = sequence2clusters[name]
        if len(clusters) == 0:
            c_index -= 1
            continue
        current_seq[name].align = 1
        meta_idx += 1
        bar.update(itern)
        already_in = _common(clusters, seen)
        _update(clusters, meta_idx, seen)
        metacluster[meta_idx] = metacluster[meta_idx].union(clusters)

        if already_in:
            for seen_metacluster in already_in:
                clusters2merge = metacluster[seen_metacluster]
                metacluster[meta_idx] = metacluster[meta_idx].union(clusters2merge)
                _update(clusters2merge, meta_idx, seen)
                # metacluster[seen_metacluster] = 0
                del metacluster[seen_metacluster]
    logger.info("%s metaclusters from %s sequences" % (len(metacluster), c_index))

    return metacluster, seen

def _find_families_deprecated(clus_obj, min_seqs):
    """
    Mask under same id all clusters that share sequences
    :param clus_obj: cluster object coming from detect_cluster
    :param min_seqs: int cutoff to keep the cluster or not. 10 as default

    :return: updated clus_obj and dict with seq_id: cluster_id
    """
    logger.info("Creating meta-clusters based on shared sequences.")
    seen = defaultdict()
    metacluster = defaultdict(list)
    c_index = clus_obj.keys()
    meta_idx = 0
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
            meta_idx += 1
            metacluster[meta_idx].append(c)
            seen.update(dict(zip(not_in, [meta_idx] * len(not_in))))
            if len(already_in) > 0:
                logger.debug("seen in %s" % already_in)
                for eindex in already_in:
                    for cluster in metacluster[eindex]:
                        metacluster[meta_idx].append(cluster)
                        prev_clus = clus_obj[cluster]
                        logger.debug("_find_families: prev %s current %s" % (eindex, clus.id))
                        # add current seqs to seen cluster
                        seqs_in = prev_clus.idmembers.keys()
                        seen.update(dict(zip(seqs_in, [meta_idx] * len(seqs_in))))
                        # for s_in_clus in prev_clus.idmembers:
                        #    seen[s_in_clus] = meta_idx
                    #    clus.idmembers[s_in_clus] = 1
                    # add current locus to seen cluster
                    # for loci in prev_clus.loci2seq:
                    #    logger.debug("adding %s" % loci)
                        # if not loci_old in current_clus[eindex].loci2seq:
                    #    clus.add_id_member(list(prev_clus.loci2seq[loci]), loci)
                    # logger.debug("loci %s" % clus.loci2seq.keys())
                    del metacluster[eindex]
                # clus_obj[c] = clus

                # logger.debug("num cluster %s" % len(clus_obj.keys()))
    logger.info("%s clusters merged" % len(metacluster))

    return metacluster, seen

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
        if bigger in clus_obj.loci:
            s, e = min(clus_obj.loci[bigger].counts.keys()), max(clus_obj.loci[bigger].counts.keys())
            scale = min(s, e)
            logger.debug("bigger %s at %s-%s" % (bigger, s, e))
            # seqs = cluster.loci2seq[bigger]
            dt = np.array([0] * (e - s + 12))
            for pos in clus_obj.loci[bigger].counts:
                ss = int(pos) - scale + 5
                # l = clus_obj.seq[seq].len
                # if st == "-":
                #    se, ss = ss, ss - l
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
