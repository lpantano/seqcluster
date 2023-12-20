from __future__ import print_function

import os
import sys
import os.path as op
from collections import Counter, namedtuple
import pickle
import json
import numpy as np
import pandas as pd

import pybedtools

from seqcluster.libs.utils import file_exists
import seqcluster.libs.logger as mylog
from seqcluster.libs import do
from seqcluster.libs.read import load_data
from seqcluster.libs.mystats import up_threshold
from seqcluster.detect.cluster import detect_clusters, clean_bam_file, peak_calling, detect_complexity
from seqcluster.detect.description import best_precursor
from seqcluster.libs.annotation import anncluster
from seqcluster.libs.inputs import parse_ma_file
from seqcluster.detect.metacluster import reduceloci, _get_seqs
from seqcluster.libs.tool import generate_position_bed
from seqcluster.libs.classes import *
import seqcluster.libs.parameters as param
from seqcluster.db import make_database


logger = mylog.getLogger(__name__)


def cluster(args):
    """
    Creating clusters
    """

    args = _check_args(args)
    read_stats_file = op.join(args.dir_out, "read_stats.tsv")
    if file_exists(read_stats_file):
        os.remove(read_stats_file)

    bam_file, seq_obj = _clean_alignment(args)

    logger.info("Parsing matrix file")
    seqL, y, l = parse_ma_file(seq_obj, args.ffile)
    # y, l = _total_counts(seqL.keys(), seqL)
    logger.info("counts after: %s" % sum(y.values()))
    logger.info("# sequences after: %s" % l)
    dt = pd.DataFrame({'sample': y.keys(), 'counts': y.values()})
    dt['step'] = 'aligned'
    dt.to_csv(read_stats_file, sep="\t", index=False, header=False, mode='a')

    if len(seqL.keys()) < 10:
        logger.error("It seems you have low coverage. Please check your fastq files have enough sequences.")
        raise ValueError("So few sequences.")

    logger.info("Cleaning bam file")
    y, l = _total_counts(list(seqL.keys()), seqL)
    logger.info("counts after: %s" % sum(y.values()))
    logger.info("# sequences after: %s" % l)
    dt = pd.DataFrame({'sample': y.keys(), 'counts': y.values()})
    dt['step'] = 'cleaned'
    dt.to_csv(read_stats_file, sep="\t", index=False, header=False, mode='a')
    clusL = _create_clusters(seqL, bam_file, args)
    # y, l = _total_counts(list(clusL.seq.keys()), clusL.seq, aligned=True)
    y, l = _total_counts(clusL.clus, seqL)
    # import pdb;pdb.set_trace()
    logger.info("counts after: %s" % sum(y.values()))
    logger.info("# sequences after: %s" % l)
    dt = pd.DataFrame({'sample': y.keys(), 'counts': y.values()})
    dt['step'] = 'clusters'
    dt.to_csv(read_stats_file, sep="\t", index=False, header=False, mode='a')

    logger.info("Solving multi-mapping events in the network of clusters")
    clusLred = _cleaning(clusL, args.dir_out)
    y, l = _total_counts(clusLred.clus, seqL)
    # import pdb;pdb.set_trace()
    # y, l = _total_counts(list(clusLred.seq.keys()), clusLred.seq, aligned=True)
    logger.info("counts after: %s" % sum(y.values()))
    logger.info("# sequences after: %s" % l)
    dt = pd.DataFrame({'sample': y.keys(), 'counts': y.values()})
    dt['step'] = 'meta-cluster'
    dt.to_csv(read_stats_file, sep="\t", index=False, header=False, mode='a')
    logger.info("Clusters up to %s" % (len(clusLred.clus.keys())))

    if args.show:
        logger.info("Creating sequences alignment to precursor")
        clusLred = show_seq(clusLred, args.index)

    clusLred = peak_calling(clusLred)

    clusLred = _annotate(args, clusLred)

    logger.info("Creating json and count matrix")
    json_file = _create_json(clusLred, args)

    logger.info("Output file in: %s" % args.dir_out)

    if args.db:
        name = args.db + ".db"
        logger.info("Create database: database/" + name)
        data = load_data(json_file)
        out_dir = op.join(args.dir_out, "database")
        make_database(data, name, out_dir)
    logger.info("Finished")


def _check_args(args):
    """
    check arguments before starting analysis.
    """
    logger.info("Checking parameters and files")
    args.dir_out = args.out
    args.samplename = "pro"
    global decision_cluster
    global similar
    if not os.path.isdir(args.out):
        logger.warning("the output folder doens't exists")
        os.mkdir(args.out)
    if args.bed and args.gtf:
        logger.error("cannot provide -b and -g at the same time")
        raise SyntaxError
    if args.debug:
        logger.info("DEBUG messages will be showed in file.")
    if args.bed:
        args.list_files = args.bed
        args.type_ann = "bed"
    if args.gtf:
        args.list_files = args.gtf
        args.type_ann = "gtf"
    logger.info("Output dir will be: %s" % args.dir_out)
    if not all([file_exists(args.ffile), file_exists(args.afile)]):
        logger.error("I/O error: Seqs.ma or Seqs.bam. ")
        raise IOError("Seqs.ma or/and Seqs.bam doesn't exists.")
    if hasattr(args, 'list_files'):
        beds = args.list_files.split(",")
        for filebed in beds:
            if not file_exists(filebed):
                logger.error("I/O error: {0}".format(filebed))
                raise IOError("%s  annotation files doesn't exist" % filebed)
    param.decision_cluster = args.method
    if args.similar:
        param.similar = float(args.similar)
    if args.min_seqs:
        param.min_seqs = int(args.min_seqs)
    return args


def _total_counts(seqs, seqL, aligned=False):
    """
    Counts total seqs after each step
    """
    total = Counter()
    nseqs = 0
    if isinstance(seqs, list):
        if not aligned:
            nseqs = len([total.update(seqL[s].freq) for s in seqs])
        else:
            nseqs = len([total.update(seqL[s].freq) for s in seqs if seqL[s].align > 0])
    elif isinstance(seqs, dict):
        #[total.update(seqs[s].get_freq(seqL)) for s in seqs]
        # import pdb;pdb.set_trace()
        # !import code; code.interact(local=vars())
        seqs_in=[]
        void=[seqs_in.extend(list(seqs[s].idmembers.keys())) for s in seqs]
        len([total.update(seqL[s].freq) for s in set(seqs_in)])
        nseqs = sum(len(seqs[s].idmembers) for s in seqs)
    return total, nseqs


def _write_size_table(data_freq, data_len, ann_valid, cluster_id):
    dd = Counter()
    for f, l in zip(data_freq, data_len):
        dd[l] += np.mean(list(f.values()))

    table = ""
    for l in sorted(dd):
        table += "%s\t%s\t%s\t%s\n" % (l, dd[l], ann_valid, cluster_id)
    return table

def _get_annotation(c, loci):
    """get annotation of transcriptional units"""
    data_ann_temp = {}
    data_ann = []
    counts = Counter()
    for lid in c.loci2seq:
        # original Py 2.7 code
        #for dbi in loci[lid].db_ann.keys():
        #    data_ann_temp[dbi] = {dbi: map(lambda (x): loci[lid].db_ann[dbi].ann[x].name, loci[lid].db_ann[dbi].ann.keys())}
        # suggestion by 2to3
        for dbi in list(loci[lid].db_ann.keys()):
            data_ann_temp[dbi] = {dbi: [loci[lid].db_ann[dbi].ann[x].name for x in list(loci[lid].db_ann[dbi].ann.keys())]}
            logger.debug("_json_: data_ann_temp %s %s" % (dbi, data_ann_temp[dbi]))
            counts[dbi] += 1
        # original Py 2.7 code
        #data_ann = data_ann + map(lambda (x): data_ann_temp[x], data_ann_temp.keys())
        # suggestion by 2to3
        data_ann = data_ann + [data_ann_temp[x] for x in list(data_ann_temp.keys())]
        logger.debug("_json_: data_ann %s" % data_ann)
    counts = {k: v for k, v in iter(counts.items())}
    total_loci = sum([counts[db] for db in counts])
    valid_ann = [k for k, v in iter(counts.items()) if up_threshold(v, total_loci, 0.7)]
    return data_ann, valid_ann


def _get_counts(list_seqs, seqs_obj, factor):
    scaled = {}
    seq = namedtuple('seq', 'freq norm_freq')
    for s in list_seqs:
        if s not in factor:
            factor[s] = 1
        samples = seqs_obj[s].norm_freq.keys()
        corrected_norm = np.array(list(seqs_obj[s].norm_freq.values())) * factor[s]
        corrected_raw = np.array(list(seqs_obj[s].freq.values())) * factor[s]
        scaled[s] = seq(dict(zip(samples, corrected_raw)), dict(zip(samples, corrected_norm)))
    return scaled


def _sum_by_samples(seqs_freq, samples_order):
    """
    Sum sequences of a metacluster by samples.
    """
    n = len(seqs_freq[list(seqs_freq.keys())[0]].freq.keys())
    y = np.array([0] * n)
    for s in seqs_freq:
        x = seqs_freq[s].freq
        exp = [seqs_freq[s].freq[sam] for sam in samples_order]
        y = list(np.array(exp) + y)
    return y


def _annotate(args, setclus):
    """annotate transcriptional units with
    gtf/bed files provided by -b/g option"""
    logger.info("Creating bed file")
    bedfile = generate_position_bed(setclus)
    a = pybedtools.BedTool(bedfile, from_string=True)
    beds = []
    logger.info("Annotating clusters")
    if hasattr(args, 'list_files'):
        beds = args.list_files.split(",")
        for filebed in beds:
            logger.info("Using %s " % filebed)
            db = os.path.basename(filebed)
            b = pybedtools.BedTool(filebed)
            c = a.intersect(b, wo=True)
            setclus = anncluster(c, setclus, db, args.type_ann, args.feature_id)
    return setclus


def _clean_alignment(args):
    """
    Prepare alignment for cluster detection.
    """
    logger.info("Clean bam file with highly repetitive reads with low counts. sum(counts)/n_hits > 1%")
    bam_file, seq_obj = clean_bam_file(args.afile, args.mask)
    logger.info("Using %s file" % bam_file)
    detect_complexity(bam_file, args.ref, args.out)
    return bam_file, seq_obj


def _create_clusters(seqL, bam_file, args):
    """
    Cluster sequences and
    create metaclusters with multi-mappers.
    """
    clus_obj = []
    cluster_file = op.join(args.out, "cluster.bed")
    if not os.path.exists(op.join(args.out, 'list_obj.pk')):
        if not file_exists(cluster_file):
            logger.info("Parsing aligned file")
            logger.info("Merging sequences")
            bedtools = os.path.join(os.path.dirname(sys.executable), "bedtools")
            bedtools = bedtools if os.path.exists(bedtools) else "bedtools"
            parse_cmd = "awk '{i=i+1;print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"i\"\\t\"$6}'"
            cmd = "{bedtools} bamtobed -i {bam_file} | {parse_cmd} | {bedtools} cluster -s -d 20 -i - > {cluster_file}"
            do.run(cmd.format(**locals()))
        c = pybedtools.BedTool(cluster_file)
        logger.info("Creating clusters")
        clus_obj = detect_clusters(c, seqL, args.min_seqs, args.non_un_gl)
        with open(op.join(args.out, 'list_obj.pk'), 'wb') as output:
            pickle.dump(clus_obj, output, pickle.HIGHEST_PROTOCOL)
    else:
        logger.info("Loading previous clusters")
        with open(op.join(args.out, 'list_obj.pk'), 'rb') as input:
            clus_obj = pickle.load(input)
    # bedfile = pybedtools.BedTool(generate_position_bed(clus_obj), from_string=True)
    # seqs_2_loci = bedfile.intersect(pybedtools.BedTool(aligned_bed, from_string=True), wo=True, s=True)
    # seqs_2_position = add_seqs_position_to_loci(seqs_2_loci, seqL)
    logger.info("%s clusters found" % (len(clus_obj.clusid)))
    return clus_obj


def _cleaning(clusL, path):
    """
    Load saved cluster and jump to next step
    """
    backup = op.join(path, "list_obj_red.pk")
    if not op.exists(backup):
        clus_obj = reduceloci(clusL, path)
        with open(backup, 'wb') as output:
            pickle.dump(clus_obj, output, pickle.HIGHEST_PROTOCOL)
        return clus_obj
    else:
        logger.info("Loading previous reduced clusters")
        with open(backup, 'rb') as in_handle:
            clus_obj = pickle.load(in_handle)
        return clus_obj


def _create_json(clusL, args):
    clus = clusL.clus
    seqs = clusL.seq
    loci = clusL.loci
    data_clus = {}
    out_count = os.path.join(args.dir_out, "counts.tsv")
    out_single_count = os.path.join(args.dir_out, "counts_sequence.tsv")
    out_size = os.path.join(args.dir_out, "size_counts.tsv")
    out_bed = os.path.join(args.dir_out, "positions.bed")
    samples_order = list(seqs[list(seqs.keys())[1]].freq.keys())
    with open(out_count, 'w') as matrix, open(out_size, 'w') as size_matrix, open(out_bed, 'w') as out_bed, open(out_single_count, 'w') as matrix_single:
        matrix.write("id\tnloci\tann\t%s\n" % "\t".join(samples_order))
        matrix_single.write("id\tann\tsequence\t%s\n" % "\t".join(samples_order))
        for cid in clus:
            seqList = []
            c = clus[cid]
            seqList = _get_seqs(c)
            logger.debug("_json_: %s" % seqList)
            data_ann, valid_ann = _get_annotation(c, loci)
            data_loci = best_precursor(c, loci)
            idloci, chrom, s, e, st, size = data_loci[0]
            annotation = valid_ann[0] if valid_ann else "none"

            bed_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, s, e, annotation, cid, st, len(seqList))
            out_bed.write(bed_line)

            # original Py 2.7 code
            #data_seqs = map(lambda (x): {x: seqs[x].seq}, seqList)
            # proposal by 2to3
            data_seqs = [{x: seqs[x].seq} for x in seqList]
            scaled_seqs = _get_counts(seqList, seqs, c.idmembers)
            # original Py 2.7 code
            #data_freq = map(lambda (x): scaled_seqs[x].freq, seqList)
            #data_freq_w_id = map(lambda (x): {x: scaled_seqs[x].norm_freq}, seqList)
            #data_len = map(lambda (x): seqs[x].len, seqList)
            # proposal by 2to3
            data_freq = [scaled_seqs[x].freq for x in seqList]
            data_freq_w_id = [{x: scaled_seqs[x].norm_freq} for x in seqList]
            data_len = [seqs[x].len for x in seqList]

            sum_freq = _sum_by_samples(scaled_seqs, samples_order)

            data_ann_str = [["%s::%s" % (name, ",".join(features)) for name, features in iter(k.items())] for k in data_ann]
            data_valid_str = " ".join(valid_ann)

            for s in seqList:
                f = [seqs[s].freq[so] for so in samples_order]
                if f.count(0) > 0.1 * len(f) and len(f) > 9:
                    continue
                f = map(str, f)
                print("\t".join([str(cid), data_valid_str, seqs[s].seq, "\t".join(f)]), file=matrix_single, end="\n")

            matrix.write("%s\t%s\t%s|%s\t%s\n" % (cid, c.toomany, data_valid_str, ";".join([";".join(d) for d in data_ann_str]), "\t".join(map(str, sum_freq))))
            size_matrix.write(_write_size_table(data_freq, data_len, data_valid_str, cid))

            data_string = {'seqs': data_seqs, 'freq': data_freq_w_id,
                           'loci': data_loci, 'ann': data_ann, 
                           'valid': valid_ann, 'peaks': clus[cid].peaks}
            data_clus[cid] = data_string

    out_file = os.path.join(args.dir_out, "seqcluster.json")
    # import pdb; pdb.set_trace()
    with open(out_file, 'w') as handle_out:
        # https://stackoverflow.com/a/50577730/1772223
        def default(o):
            if isinstance(o, np.int64): return int(o)
            raise TypeError
        handle_out.write(json.dumps([data_clus], default=default, skipkeys=True, indent=2))

    return out_file
