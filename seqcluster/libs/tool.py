from collections import defaultdict
import operator
import os
import copy
# import time
import math
# import numpy as np
import pybedtools
import logger as mylog
from classes import *
from mystats import up_threshold
from annotation import _position_in_feature, read_gtf_line
from bayes import decide_by_bayes
import parameters

logger = mylog.getLogger(__name__)


def get_distance(p1, p2):
        if p1.strand == "+":
                d = p1.start-p2.start
        else:
                d = (p1.end-(p1.end-p1.start)) - (p2.end-(p2.end-p2.start))
        return d


def get_ini_str(n):
        #get space values
        return "".join(" " for i in range(n))


def readbowtie():
    ##read results from bowtie...no longer needed
    dict = {}
    f = open("temp.map", 'r')
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        numism = line.count(">")
        mism = ""
        if numism >= 1:
            mism = cols[7]
       
        result = "%s:%s %s %s;" % (cols[2],cols[3],cols[1],mism)
        if dict.has_key(cols[0]):
            store = dict[cols[0]]+result
            dict[cols[0]] = store
        else:
            dict[cols[0]] = result
    f.close()
    return dict


def what_is(file):
    with open(file, 'r') as handle_in:
        for line in handle_in:
            if not line.startswith("@"):
                cols =  line.split("\t")
                if (cols[5] == "+" or cols[5] == "-"):
                    return "bed"
                else:
                    return "sam"
                break


def init_numlocidb(beds):
    dict = {}
    for filebed in beds:
        db =  os.path.basename(filebed)
        dict[db] = 0;
    return dict


def calc_complexity(nl):
    ns = len(nl)
    total = 0.0
    for l in nl:
        total += 1.0/l
        #print total
    total /= ns
    return (total)


def calculate_size(vector):
    maxfreq = 0
    zeros = 0
    counts = 0
    total = len(vector.keys())
    for s in vector.keys():
        maxfreq = max(vector[s], maxfreq)
        counts += int(vector[s])
        if vector[s] == 0:
            zeros += 1
    return (counts*((total-zeros)/total))


def show_seq(clus_obj, index):
    """Get the precursor and map sequences to it.
       this way we create a positional map."""
    current = clus_obj.clus
    clus_seqt = clus_obj.seq
    clus_locit = clus_obj.loci

    itern = 0
    for idc in current.keys():
        itern += 1
        timestamp = str(idc)
        seqListTemp = ()
        f = open("/tmp/"+timestamp+".fa","w")
        for idl in current[idc].loci2seq.keys():
            seqListTemp = list(set(seqListTemp).union(current[idc].loci2seq[idl]))
        maxscore = 0
        for s in seqListTemp:
            score = calculate_size(clus_seqt[s].freq)
            maxscore = max(maxscore,score)
            clus_seqt[s].score = score
            seq = clus_seqt[s]
            f.write(">"+s+"\n"+seq.seq+"\n")
        f.close()

        locilen_sorted = sorted(current[idc].locilen.iteritems(), key = operator.itemgetter(1),reverse = True)
        lmax = clus_locit[locilen_sorted[0][0]]
        f = open("/tmp/"+timestamp+".bed","w")
        f.write("%s\t%s\t%s\t.\t.\t%s\n" % (lmax.chr,lmax.start,lmax.end,lmax.strand))
        f.close()
        os.system("bedtools  getfasta -s -fi "+index+" -bed /tmp/"+timestamp+".bed -fo /tmp/"+timestamp+".pre.fa")
        os.system("bowtie2-build /tmp/"+timestamp+".pre.fa  /tmp/"+timestamp+".pre.ind >/dev/null 2>&1")
        os.system("bowtie2 --rdg 7,3 --mp 4 --end-to-end --no-head --no-sq  -D 20 -R 3 -N 0 -i S,1,0.8 -L 3 -f /tmp/"+timestamp+".pre.ind /tmp/"+timestamp+".fa -S /tmp/"+timestamp+".map  >>bowtie.log 2>&1")
        f = open("/tmp/"+timestamp+".map","r")
        seqpos = {}
        minv = 10000000
        for line in f:
            line = line.strip()
            cols = line.split("\t")
            seqpos[cols[0]] = int(cols[3])
            if minv>int(cols[3]):
                minv = int(cols[3])
        f.close()
        seqpos_sorted = sorted(seqpos.iteritems(), key = operator.itemgetter(1),reverse = False)
        showseq = ""
        showseq_plain = ""
        for (s,pos) in seqpos_sorted:
            ratio = (clus_seqt[s].score*1.0/maxscore*100.0)
            realScore = (math.log(ratio,2)*2)
            if realScore<0:
                realScore = 0
            # "score %s max %s  ratio %s real %.0f" % (clus_seqt[s].score,maxscore,ratio,realScore)
            ##calculate the mean expression of the sequence and change size letter
            showseq_plain += "<br>%s<a style = \"font-size:%.0fpx;\"href = javascript:loadSeq(\"%s\")>%s</a>" % ("".join("." for i in range(pos-1)),realScore+10,s,clus_seqt[s].seq)
            #showseq+ = seqviz.addseq(pos-1,clus_seqt[s].len,clus_seqt[s].seq)
        #current[idc].showseq = showseq
        current[idc].showseq_plain = showseq_plain
        os.system("rm /tmp/"+timestamp+"*")
    clus_obj.clus = current
    clus_obj.seq = clus_seqt
    return clus_obj


def anncluster(c, clus_obj, db, type_ann):
    """intersect transcription position with annotation files"""
    id_sa, id_ea, id_id, id_idl, id_sta = 1, 2, 3, 4, 5
    if type_ann == "bed":
        id_sb = 7
        id_eb = 8
        id_stb = 11
        id_tag = 9
    ida = 0
    clus_id = clus_obj.clus
    loci_id = clus_obj.loci
    db = os.path.splitext(db)[0]
    logger.debug("Type:%s\n" % type_ann)
    for cols in c.features():
        if type_ann == "gtf":
            cb, sb, eb, stb, db, tag = read_gtf_line(cols[6:])
        else:
            sb = int(cols[id_sb])
            eb = int(cols[id_eb])
            stb = cols[id_stb]
            tag = cols[id_tag]
        id = int(cols[id_id])
        idl = int(cols[id_idl])
        if (id in clus_id):
            clus = clus_id[id]
            sa = int(cols[id_sa])
            ea = int(cols[id_ea])
            ida += 1
            lento5, lento3, strd = _position_in_feature([sa, ea, cols[id_sta]], [sb, eb, stb])
            if db in loci_id[idl].db_ann:
                ann = annotation(db, tag, strd, lento5, lento3)
                tdb = loci_id[idl].db_ann[db]
                tdb.add_db_ann(ida, ann)
                loci_id[idl].add_db(db, tdb)
            else:
                ann = annotation(db, tag, strd, lento5, lento3)
                tdb = dbannotation(1)
                tdb.add_db_ann(ida, ann)
                loci_id[idl].add_db(db, tdb)
            clus_id[id] = clus
    clus_obj.clus = clus_id
    clus_obj.loci = loci_id
    return clus_obj


def parse_align_file(file_in):
    """parse sam files with aligned sequences"""
    loc_id = 1
    bedfile_clusters = ""
    bamfile = pybedtools.BedTool(file_in)
    bed = pybedtools.BedTool.bam_to_bed(bamfile)
    for c, start, end, name, q, strand in bed:
        loc_id += 1
        bedfile_clusters += "%s\t%s\t%s\t%s\t%s\t%s\n" % \
                    (c, start, end, name, loc_id, strand)
    return bedfile_clusters


def parse_merge_file_deprecated(c, current_seq, MIN_SEQ):
    """
    Parse the merge file of sequences position to create clusters that will have all
    sequences that shared any position on the genome

    :param c: file from bedtools with merge sequence positions
    :param seq_l_in: list of sequences
    :param MIN_SEQ: int cutoff to keep the cluster or not. 10 as default

    :return: object with information about:
        * cluster
        $ dict with sequences (as keys) and cluster_id (as value)
        * sequences
        * loci

    """
    clus_id = {}
    current_loci = {}
    current_clus = {}
    index = 0
    lindex = 0
    eindex = 0
    for line in c.features():
        a = mergealigned(line)
        if len(a.names) >= MIN_SEQ:
            lindex += 1
            already_in, not_in = _get_seqs_from_cluster(a.names, clus_id)
            logger.debug("_merge_aligned: idl %s %s %s %s %s" % (lindex, a.chr, a.start, a.end, a.strand))
            if len(already_in) > 0:
                eindex = already_in[0]
                for toremove in already_in[1:]:
                    prev_clus = current_clus[toremove]
                    for s_in_prev in prev_clus.idmembers:
                        clus_id[s_in_prev] = eindex
                        current_clus[eindex].idmembers[s_in_prev] = 0
                    for loci_old in prev_clus.loci2seq:
                        #if not loci_old in current_clus[eindex].loci2seq:
                        current_clus[eindex].add_id_member(list(prev_clus.loci2seq[loci_old]), loci_old)
                        #if exists loci, merge sequences
                    del current_clus[toremove]
            else:
                index += 1
                eindex = index
            if not eindex in current_clus:
                current_clus[eindex] = cluster(eindex)
            newpos = position(lindex, a.chr, a.start, a.end, a.strand)
            current_loci[lindex] = newpos
            for s in a.names:
                current_seq[s].add_pos(lindex)
                current_clus[eindex].idmembers[s] = 0
                clus_id[s] = eindex
            current_clus[eindex].add_id_member(a.names, lindex)
    return cluster_info_obj(current_clus, clus_id, current_loci, current_seq)


def add_seqs_position_to_loci(fn_bedtools, seqs):
    """return seqL, with exact position in that loci
    :param fn_bedtools: pybedtool object
    :param seqs: object class sequences

    return: dict with locus as key and positions as values"""
    seqs_pos = defaultdict(defaultdict)
    for line in fn_bedtools:
        locus = line[4]
        name = line[9]
        pos = [line[6], line[7], line[8]]
        strand = line[9]
        start = pos[1] if strand == "+" else pos[2]
        if start not in seqs_pos[locus]:
            seqs_pos[locus][start] = []
        seqs_pos[locus][start].append(name)
        #print ("{locus} {name} {pos}").format(**locals())
    return seqs_pos


def _get_seqs_from_cluster(seqs, clus_id):
    """
    Returns the sequences that are already part of the cluster

    :param seqs: list of sequences ids
    :param clus_id: dict of sequences ids that are part of a cluster

    :returns:
        * :code:`already_in`list of cluster id that contained some of the sequences
        * :code:`not_in`list of sequences that don't belong to any cluster yet
    """
    already_in = set()
    not_in = []
    for s in seqs:
        if s in clus_id:
            already_in.add(clus_id[s])
        else:
            not_in.append(s)
    return list(already_in), not_in


def reduceloci(clus_obj,  path):
    """reduce number of loci a cluster has"""
    filtered = {}
    n_cluster = 0
    current = clus_obj.clus
    logger.info("Number of loci: %s" % len(clus_obj.loci.keys()))
    for idc in current:
        logger.debug("_reduceloci: cluster %s" % idc)
        c = copy.deepcopy(current[idc])
        n_loci = len(c.loci2seq)
        if n_loci < 1000:
            filtered, n_cluster = _iter_loci(c, (clus_obj.loci, clus_obj.seq), filtered, n_cluster)
        else:
            n_cluster += 1
            filtered[n_cluster] = _add_complete_cluster(n_cluster, c)
    clus_obj.clus = filtered
    return clus_obj


def _add_complete_cluster(idx, clus1):
    logger.debug("Not resolving cluster %s, too many loci. New id %s" % (clus1.id, idx))
    locilen_sorted = sorted(clus1.locilen.iteritems(), key=operator.itemgetter(1), reverse=True)
    maxidl = locilen_sorted[0][0]
    c = cluster(idx)
    c.add_id_member(clus1.loci2seq[maxidl], maxidl)
    c.id = idx
    c.toomany = len(locilen_sorted)
    return c


def _iter_loci_deprecated(c, filtered, n_cluster, min_seq):
    """Go through all locus and decide if they are part
    of the same TU or not.

    :param idx: int cluster id
    :param filtered: dict with clusters object
    :param n_cluster: int cluster id
    :param min_seq: int min number of sequences inside
    cluster

    :return:
        * filtered: dict of cluster objects
        * n_cluster: int cluster id"""
    n_loci = len(c.loci2seq)
    n_loci_prev = n_loci + 1
    total_seqs = list()
    cicle = 0
    while n_loci < n_loci_prev and n_loci != 0:
        n_loci_prev = n_loci
        cicle += 1
        ma = _calculate_similarity(c)
        if (cicle % 1) == 0:
            logger.debug("_iter_loci:number of cicle: %s with n_loci %s" % (cicle, n_loci))
        locilen_sorted = sorted(c.locilen.iteritems(), key=operator.itemgetter(1), reverse=True)
        maxseq = locilen_sorted[0][1]*1.0
        if maxseq > min_seq:
            logger.debug("_iter_loci:maxseq: %s" % maxseq)
            c, total_seqs, filtered, n_cluster = _solve_loci(c, locilen_sorted, total_seqs, filtered, maxseq, n_cluster)
        else:
            for (idl, lenl) in locilen_sorted:
                logger.debug("_iter_loci:remove locus %s with len %s:" % (idl, lenl))
                c.loci2seq.pop(idl, "None")
                c.locilen.pop(idl, "None")
        n_loci = len(c.loci2seq)
    return filtered, n_cluster


def _iter_loci(c, s2p, filtered, n_cluster):
    """Go through all locus and decide if they are part
    of the same TU or not.

    :param idx: int cluster id
    :param s2p: dict with [loci].coverage[start] = # of sequences there
    :param filtered: dict with clusters object
    :param n_cluster: int cluster id

    :return:
        * filtered: dict of cluster objects
        * n_cluster: int cluster id"""
    n_loci = len(c.loci2seq)
    n_loci_prev = n_loci + 1
    cicle = 0
    # [logger.note("BEFORE %s %s %s" % (c.id, idl, len(c.loci2seq[idl]))) for idl in c.loci2seq]
    internal_cluster = {}
    loci = _convert_to_clusters(c)
    if n_loci == 1:
        n_cluster += 1
        c.id = n_loci
        filtered[n_cluster] = c
    while n_loci < n_loci_prev and n_loci != 1:
        n_loci_prev = n_loci
        cicle += 1
        if (cicle % 1) == 0:
            logger.debug("_iter_loci:number of cicle: %s with n_loci %s" % (cicle, n_loci))
        loci_similarity = _calculate_similarity(loci)
        loci_similarity = sorted(loci_similarity.iteritems(), key=operator.itemgetter(1), reverse=True)
        internal_cluster = _merge_similar(loci, loci_similarity)
        n_loci = len(internal_cluster)
        loci = internal_cluster
        logger.debug("_iter_loci: n_loci %s" % n_loci)
    if n_loci > 1:
        n_internal_cluster = sorted(internal_cluster.keys(), reverse=True)[0]
        internal_cluster = _solve_conflict(internal_cluster, s2p, n_internal_cluster)
    internal_cluster = _clean_cluster(internal_cluster)
    for idc in internal_cluster:
        n_cluster += 1
        logger.debug("_iter_loci: add to filtered %s" % n_cluster)
        filtered[n_cluster] = internal_cluster[idc]
        filtered[n_cluster].id = n_cluster
    logger.debug("_iter_loci: filtered %s" % filtered.keys())
    for new_c in internal_cluster.values():
        [logger.note("%s %s %s %s" % (c.id, new_c.id, idl, len(new_c.loci2seq[idl]))) for idl in new_c.loci2seq]
    return filtered, n_cluster


def _remove_loci(ci, idl):
    for (idl, lenl) in locilen_sorted:
        logger.debug("_remove_loci:remove locus %s with len %s:" % (idl, lenl))
        c.loci2seq.pop(idl, "None")
        c.locilen.pop(idl, "None")


def _convert_to_clusters(c):
    """Return 1 cluster per loci"""
    new_dict = {}
    n_cluster = 0
    logger.debug("_convert_to_cluster: loci %s" % c.loci2seq.keys())
    for idl in c.loci2seq:
        n_cluster += 1
        new_c = cluster(n_cluster)
        #new_c.id_prev = c.id
        new_c.loci2seq[idl] = c.loci2seq[idl]
        new_dict[n_cluster] = new_c
    logger.debug("_convert_to_cluster: new ids %s" % new_dict.keys())
    return new_dict


def _calculate_similarity(c):
    """Get a similarity matrix of % of shared sequence

    :param c: cluster object

    :return ma: similarity matrix
    """
    ma = {}
    for idc in c:
        set1 = _get_seqs(c[idc])
        [ma.update({(idc, idc2): _common(set1, _get_seqs(c[idc2]))}) for idc2 in c if idc != idc2 and (idc2, idc) not in ma]
    logger.debug("_calculate_similarity_ %s" % ma)
    return ma


def _get_seqs(list_idl):
    """get all sequences in a cluster knowing loci"""
    seqs = set()
    for idl in list_idl.loci2seq:
        logger.debug("_get_seqs_: loci %s" % idl)
        [seqs.add(s) for s in list_idl.loci2seq[idl]]
    logger.debug("_get_seqs_: %s" % len(seqs))
    return seqs


def _common(s1, s2):
    """calculate the common % percentage of sequences"""
    c = len(set(s1).intersection(s2))
    t = min(len(s1), len(s2))
    return 1.0 * c / t


def _merge_similar(loci, locilen_sorted):
    """internal function to reduce loci complexity

    :param loci: class cluster
    :param locilen_sorted: list of loci sorted by size

    :return
     c: updated class cluster
    """
    n_cluster = 0
    internal_cluster = {}
    clus_seen = {}
    for pairs, common in locilen_sorted:
        n_cluster += 1
        logger.debug("_merge_similar:try new cluster %s" % n_cluster)
        new_c = cluster(n_cluster)
        logger.debug("_merge_similar:id %s  common %s" % (pairs, common))
        p_seen, p_unseen = [], []
        size = min(len(_get_seqs(loci[pairs[0]])), len(_get_seqs(loci[pairs[1]])))
        if up_threshold(common * size, size, parameters.similar):
            if pairs[0] in clus_seen:
                p_seen.append(pairs[0])
                p_unseen.append(pairs[1])
            if pairs[1] in clus_seen:
                p_seen.append(pairs[1])
                p_unseen.append(pairs[0])
            if len(p_seen) == 0:
                new_c = _merge_cluster(loci[pairs[0]], new_c)
                new_c = _merge_cluster(loci[pairs[1]], new_c)
                [clus_seen.update({p: n_cluster}) for p in pairs]
                internal_cluster[n_cluster] = new_c
            if len(p_seen) == 1:
                idc_seen = clus_seen[p_seen[0]]
                internal_cluster[idc_seen] = _merge_cluster(loci[p_unseen[0]], internal_cluster[idc_seen])
                clus_seen[p_unseen[0]] = idc_seen
        else:
            continue
    internal_cluster.update(_add_unseen(loci, clus_seen, n_cluster))
    logger.debug("_merge_similar: total clus %s" %
                 len(internal_cluster.keys()))
    return internal_cluster


def _merge_cluster(old, new):
    """merge one cluster to another"""
    logger.debug("_merge_cluster: %s to %s" % (old.id, new.id))
    for idl in old.loci2seq:
        logger.debug("_merge_cluster: add idl %s" % idl)
        new.loci2seq[idl] = old.loci2seq[idl]
    return new


def _solve_conflict(list_c, s2p, n_cluster):
    """make sure sequences are counts once.
    Resolve by most-vote or exclussion

    :params list_c: dict of objects cluster
    :param s2p: dict of [loci].coverage = # num of seqs
    :param n_cluster: number of clusters

    return dict: new set of clusters"""
    logger.debug("_solve_conflict: count once")
    if parameters.decision_cluster == "bayes":
        return decide_by_bayes(list_c, s2p)
    loci_similarity = _calculate_similarity(list_c)
    loci_similarity = sorted(loci_similarity.iteritems(), key=operator.itemgetter(1), reverse=True)
    common = sum([score for p, score in loci_similarity])
    while common:
        n_cluster += 1
        logger.debug("_solve_conflict: ma %s" % loci_similarity)
        logger.debug("_solve_conflict: common %s, new %s" % (common, n_cluster))
        pairs = loci_similarity[0][0]
        if parameters.decision_cluster.startswith("most-voted"):
            list_c = _split_cluster_by_most_vote(list_c, pairs)
        else:
            list_c = _split_cluster(list_c, pairs, n_cluster)
        list_c = {k: v for k, v in list_c.iteritems() if len(v.loci2seq) > 0}
        loci_similarity = _calculate_similarity(list_c)
        loci_similarity = sorted(loci_similarity.iteritems(), key=operator.itemgetter(1), reverse=True)
        #logger.note("%s %s" % (pairs, loci_similarity[0][1]))
        common = sum([score for p, score in loci_similarity])
        logger.debug("_solve_conflict: solved clusters %s" % len(list_c.keys()))
    return list_c


def _split_cluster(c , pairs, n):
    """split cluster by exclussion"""
    old = c[p[0]]
    new = c[p[1]]
    new_c = cluster(n)
    common = set(_get_seqs(old)).intersection(_get_seqs(new))
    for idl in old.loci2seq:
        in_common = list(set(common).intersection(old.loci2seq[idl]))
        if len(in_common) > 0:
            logger.debug("_split_cluster: in_common %s with pair 1" % (len(in_common)))
            new_c.loci2seq[idl] = in_common
            old.loci2seq[idl] = list(set(old.loci2seq[idl]) - set(common))
            logger.debug("_split_cluster: len old %s with pair 1" % (len(old.loci2seq)))
    for idl in new.loci2seq:
        in_common = list(set(common).intersection(new.loci2seq[idl]))
        if len(in_common) > 0:
            logger.debug("_split_cluster: in_common %s with pair 2" % (len(in_common)))
            new_c.loci2seq[idl] = in_common
            new.loci2seq[idl] = list(set(new.loci2seq[idl]) - set(common))
            logger.debug("_split_cluster: len old %s with pair 2" % (len(new.loci2seq)))
    old.loci2seq = {k: v for k, v in old.loci2seq.iteritems() if len(v) > 0}
    new.loci2seq = {k: v for k, v in new.loci2seq.iteritems() if len(v) > 0}
    c[n] = new
    c[p[0]] = old
    c[p[1]] = new
    return c


def _split_cluster_by_most_vote(c, p):
    """split cluster by most-vote strategy"""
    old, new = c[p[0]], c[p[1]]
    old_size = _get_seqs(old)
    new_size = _get_seqs(new)
    logger.debug("_most_vote: size of %s %s - %s %s" % (old.id, len(old_size), new.id, len(new_size)))
    if len(old_size) > len(new_size):
        keep, remove = old, new
    else:
        keep, remove = new, old
    common = list(set(old_size).intersection(new_size))
    logger.debug("_most_vote: keep %s remove  %s with common %s" % (keep.id, remove.id, len(common)))
    for idl in remove.loci2seq:
        if len(common) > 0:
            remove.loci2seq[idl] = list(set(remove.loci2seq[idl]) - set(common))
    keep.loci2seq = {k: v for k, v in keep.loci2seq.iteritems() if len(v) > 0}
    remove.loci2seq = {k: v for k, v in remove.loci2seq.iteritems() if len(v) > 0}
    c[keep.id] = keep
    c[remove.id] = remove
    return c


def _add_unseen(loci, clus_seen, n_cluster):
    unseen = {}
    for idc in loci:
        if idc not in clus_seen:
            n_cluster += 1
            loci[idc].id = n_cluster
            unseen[n_cluster] = loci[idc]
            logger.debug("_add_unseen: add  %s as new %s" %
                         (idc, n_cluster))
    return unseen


def _clean_cluster(list_c):
    """Remove cluster with less than 10 sequences and
    loci with size smaller than 60%"""
    list_c = {k: v for k, v in list_c.iteritems() if len(_get_seqs(v)) > parameters.min_seqs}
    logger.debug("_clean_cluster: number of clusters %s " % len(list_c.keys()))
    list_c = {k: _select_loci(v) for k, v in list_c.iteritems()}
    return list_c


def _select_loci(c):
    """Select only loci with most abundant sequences"""
    loci_len = {k: len(v) for k, v in c.loci2seq.iteritems()}
    logger.debug("_select_loci: number of loci %s" % len(c.loci2seq.keys()))
    loci_len_sort = sorted(loci_len.iteritems(), key=operator.itemgetter(1), reverse=True)
    max_size = loci_len_sort[0][1]
    logger.debug("_select_loci: max size %s" % max_size)
    loci_clean = {locus: c.loci2seq[locus] for locus, size in loci_len_sort if size > 0.8 * max_size}
    c.loci2seq = loci_clean
    logger.debug("_select_loci: number of loci %s after cleaning" % len(c.loci2seq.keys()))
    return c


def _calculate_size_enrichment(c):
    """calculate whether there is a size
    enrichment in the cluster

    :param c: cluster object
    """
    seqs = _get_sequences(c)
    return True


def _solve_loci_deprecated(c, locilen_sorted, seen_seqs, filtered, maxseq, n_cluster):
    """internal function to reduce loci complexity

    The function will read the all loci in a cluster of
    sequences and will determine if all loci are part
    of the same transcriptional unit(TU) by most-vote locus
    or by exclusion of common sequence that are the
    minority of two loci.

    :param c: class cluster
    :param locilen_sorted: list of loci sorted by size
    :param seem_seqs: list of seen sequences
    :param filtered: final TU list
    :param maxseq: bigger locus
    "param n_cluster: integer with index of different TU"
    :return
     c: updated class cluster
     seen_seqs: updated list of sequences
     filtered: updated dict of TUs
     n_cluster: updated int with current index of TUs
    """
    first_run = 0
    seen_seqs = list()
    n_cluster += 1
    logger.debug("_solve_loci:new cluster %s" % n_cluster)
    new_c = cluster(n_cluster)
    for idl, lenl in locilen_sorted:
        locus_seqs = c.loci2seq[idl]
        if first_run == 0:
            seen_seqs = locus_seqs
            first_run = 1
            first_idl = idl
        intersect = list(set(seen_seqs).intersection(locus_seqs))
        common = 0
        if intersect:
            common = len(intersect)*1.0/min(len(seen_seqs), len(locus_seqs))
        logger.debug("_sole_loci:id %s idl %s len %s max %s seen %s inter %s common %s " % (c.id, idl, lenl, maxseq, len(seen_seqs), len(intersect), common))
        if common*1.0 >= 0.6:
            if lenl*1.0 >= 0.6*maxseq:
                c, new_c, seen_seqs = _merge_loci_in_cluster(c, new_c, idl, seen_seqs)
            else:
                c, new_c, seen_seqs = _merge_with_first_loci(c, new_c, first_idl, idl, seen_seqs)
        else:
            c = _remove_seqs_from_loci(c, idl, seen_seqs)
    filtered[n_cluster] = new_c
    return c, seen_seqs, filtered, n_cluster


def _merge_loci_in_cluster(c, new_c, idl, current_seqs):
    logger.debug("_merge_loci_in_cluster:join")
    locus_seqs = c.loci2seq[idl]
    common = len(set(locus_seqs).intersection(current_seqs))
    seen = list(set(locus_seqs).union(current_seqs))
    new_c.add_id_member(list(locus_seqs), idl)
    c.loci2seq.pop(idl, "None")
    c.locilen.pop(idl, "None")
    return c, new_c, seen


def _merge_with_first_loci(c, new_c, first_idl, idl, current_seqs):
    logger.debug("_merge_with_first_loci:join first")
    locus_seqs = c.loci2seq[idl]
    common = len(set(locus_seqs).intersection(current_seqs))
    seen = list(set(locus_seqs).union(current_seqs))
    new_c.add_id_member(list(locus_seqs), first_idl)
    c.loci2seq.pop(idl, "None")
    c.locilen.pop(idl, "None")
    return c, new_c, seen


def _remove_seqs_from_loci(c, idl, current_seqs):
    current = c.loci2seq[idl]
    common = len(set(current).intersection(current_seqs))
    seen = list(set(current).intersection(current_seqs))
    unseen = list(set(sorted(current)).difference(sorted(seen)))
    logger.debug("_remove_seqs_from_loci:seen %s unseen %s" % (len(seen), len(unseen)))
    c.locilen[idl] = len(unseen)
    c.loci2seq[idl] = unseen
    logger.debug("_remove_seqs_from_loci:remove; new len %s" % len(unseen))
    if c.locilen[idl] == 0:
        c.loci2seq.pop(idl, "None")
        c.locilen.pop(idl, "None")
    return c


def generate_position_bed(clus_obj):
    ##generate file with positions in bed format
    bedaligned = ""
    clus_id = clus_obj.clus
    for idc in clus_id.keys():
        clus = clus_id[idc]
        for idl in clus.loci2seq.keys():
            pos = clus_obj.loci[idl]
            bedaligned +=  "%s\t%s\t%s\t%s\t%s\t%s\n" % (pos.chr,pos.start,pos.end,idc,idl,pos.strand)
    return bedaligned


def parse_ma_file(in_file):
    """read seqs.ma file and create dict with
    sequence object"""
    name = ""
    seq_l = {}
    index = 1
    total = defaultdict(int)
    with open(in_file) as handle_in:
        line = handle_in.readline().strip()
        cols = line.split("\t")
        samples = cols[2:]
        for line in handle_in:
            line = line.strip()
            cols = line.split("\t")
            exp = {}
            for i in range(len(samples)):
                exp[samples[i]] = int(cols[i+2])
                total[samples[i]] += int(cols[i+2])
            name = cols[0].replace(">", "")
            index = index+1
            seq = cols[1]
            new_s = sequence(seq, exp, index)
            seq_l[name] = new_s
    seq_l = _normalize_seqs(seq_l, total)
    return seq_l


def _normalize_seqs(s, t):
    """Normalize to RPM"""
    for ids in s:
        obj = s[ids]
        [obj.norm_freq.update({sample: 1.0 * obj.freq[sample] / (t[sample]+1) * 1000000}) for sample in obj.norm_freq]
        s[ids] = obj
    return s
