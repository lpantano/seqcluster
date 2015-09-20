import os
import os.path as op
from collections import Counter, namedtuple
from operator import itemgetter
import pickle
import json
import numpy as np
import pandas as pd

import pybedtools

from bcbio.utils import file_exists

import libs.logger as mylog
from libs.read import load_data
from libs.mystats import up_threshold
from detect.cluster import detect_clusters, clean_bam_file, peak_calling, detect_complexity
from libs.annotation import anncluster
from libs.inputs import parse_ma_file, parse_align_file
from detect.metacluster import reduceloci,_get_seqs
from libs.tool import generate_position_bed
from libs.classes import *
import libs.parameters as param
from db import make_database


logger = mylog.getLogger(__name__)


def cluster(args):
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
        logger.error("It seems you have so low coverage. Please check your fastq files have enough sequences.")
        raise ValueError("So few sequences.")

    logger.info("Cleaning bam file")
    y, l = _total_counts(seqL.keys(), seqL)
    logger.info("counts after: %s" % sum(y.values()))
    logger.info("# sequences after: %s" % l)
    dt = pd.DataFrame({'sample': y.keys(), 'counts': y.values()})
    dt['step'] = 'cleaned'
    dt.to_csv(read_stats_file, sep="\t", index=False, header=False, mode='a')

    clusL = _create_clusters(seqL, bam_file, args)
    y, l = _total_counts(clusL.seq.keys(), clusL.seq, aligned=True)
    logger.info("counts after: %s" % sum(y.values()))
    logger.info("# sequences after: %s" % l)
    dt = pd.DataFrame({'sample': y.keys(), 'counts': y.values()})
    dt['step'] = 'clusters'
    dt.to_csv(read_stats_file, sep="\t", index=False, header=False, mode='a')

    logger.info("Solving multi-mapping events in the network of clusters")
    clusLred = _cleaning(clusL, args.dir_out)
    y, l = _total_counts(clusLred.clus, seqL)
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


def _total_counts(seqs, seqL, aligned=False):
    """
    Counts total seqs after each step
    """
    total = Counter()
    if isinstance(seqs, list):
        if not aligned:
            l = len([total.update(seqL[s].freq) for s in seqs])
        else:
            l = len([total.update(seqL[s].freq) for s in seqs if seqL[s].align > 0])
    elif isinstance(seqs, dict):
        [total.update(seqs[s].get_freq(seqL)) for s in seqs]
        l = sum(len(seqs[s].idmembers) for s in seqs)
    return total, l


def _write_size_table(data_freq, data_len, ann_valid, cluster_id):
    dd = Counter()
    for f, l in zip(data_freq, data_len):
        dd[l] += np.mean(f.values())

    table = ""
    for l in sorted(dd):
        table += "%s\t%s\t%s\t%s\n" % (l, dd[l], ann_valid, cluster_id)
    return table


def _create_json(clusL, args):
    clus = clusL.clus
    seqs = clusL.seq
    loci = clusL.loci
    data_clus = {}
    out_count = os.path.join(args.dir_out, "counts.tsv")
    out_size = os.path.join(args.dir_out, "size_counts.tsv")
    out_bed = os.path.join(args.dir_out, "positions.bed")
    samples_order = list(seqs[seqs.keys()[1]].freq.keys())
    with open(out_count, 'w') as matrix, open(out_size, 'w') as size_matrix, open(out_bed, 'w') as out_bed:
        matrix.write("id\tnloci\tann\t%s\n" % "\t".join(samples_order))
        for cid in clus:
            seqList = []
            c = clus[cid]
            data_loci = map(lambda (x): [x, loci[x].chr, int(loci[x].start), int(loci[x].end), loci[x].strand, len(c.loci2seq[x])], c.loci2seq.keys())
            data_loci = sorted(data_loci, key=itemgetter(5), reverse=True)
            seqList = _get_seqs(c)
            logger.debug("_json_: %s" % seqList)
            data_ann, valid_ann = _get_annotation(c, loci)

            idloci, chrom, s, e, st, size = data_loci[0]
            annotation = valid_ann[0] if valid_ann else "none"

            bed_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, s, e, annotation, cid, st, len(seqList))
            out_bed.write(bed_line)


            data_seqs = map(lambda (x): {x: seqs[x].seq}, seqList)
            # data_freq = map(lambda (x): seqs[x].freq, seqList)
            scaled_seqs = _get_counts(seqList, seqs, c.idmembers)
            data_freq = map(lambda (x): scaled_seqs[x].freq, seqList)
            data_freq_w_id = map(lambda (x): {x: scaled_seqs[x].norm_freq}, seqList)
            data_len = map(lambda (x): seqs[x].len, seqList)
            # data_freq_values = map(lambda (x): map(int, scaled_seqs[x].freq.values()), seqList)
            # sum_freq = _sum_by_samples(data_freq_values)
            sum_freq = _sum_by_samples(scaled_seqs, samples_order)

            data_ann_str = [["%s::%s" % (name, ",".join(features)) for name, features in k.iteritems()] for k in data_ann]
            data_valid_str = " ".join(valid_ann)
            matrix.write("%s\t%s\t%s|%s\t%s\n" % (cid, c.toomany, data_valid_str, ";".join([";".join(d) for d in data_ann_str]), "\t".join(map(str, sum_freq))))
            size_matrix.write(_write_size_table(data_freq, data_len, data_valid_str, cid))

            data_string = {'seqs': data_seqs, 'freq': data_freq_w_id,
                    'loci': data_loci, 'ann': data_ann, 'valid': valid_ann, 'peaks': clus[cid].peaks}
            data_clus[cid] = data_string


    out_file = os.path.join(args.dir_out, "seqcluster.json")
    with open(out_file, 'w') as handle_out:
        handle_out.write(json.dumps([data_clus], skipkeys=True, indent=2))

    return out_file


def _get_annotation(c, loci):
    """get annotation of transcriptional units"""
    data_ann_temp = {}
    data_ann = []
    counts = Counter()
    for lid in c.loci2seq:
        for dbi in loci[lid].db_ann.keys():
            data_ann_temp[dbi] = {dbi: map(lambda (x): loci[lid].db_ann[dbi].ann[x].name, loci[lid].db_ann[dbi].ann.keys())}
            logger.debug("_json_: data_ann_temp %s %s" % (dbi, data_ann_temp[dbi]))
            counts[dbi] += 1
        data_ann = data_ann + map(lambda (x): data_ann_temp[x], data_ann_temp.keys())
        logger.debug("_json_: data_ann %s" % data_ann)
    counts = {k: v for k, v in counts.iteritems()}
    total_loci = sum([counts[db] for db in counts])
    valid_ann = [k for k, v in counts.iteritems() if up_threshold(v, total_loci, 0.7)]
    return data_ann, valid_ann


def _get_counts(list_seqs, seqs_obj, factor):
    scaled = {}
    seq = namedtuple('seq', 'freq norm_freq')
    for s in list_seqs:
        if s not in factor:
            factor[s] = 1
        samples = seqs_obj[s].norm_freq.keys()
        corrected_norm = np.array(seqs_obj[s].norm_freq.values()) * factor[s]
        corrected_raw = np.array(seqs_obj[s].freq.values()) * factor[s]
        scaled[s] = seq(dict(zip(samples, corrected_raw)), dict(zip(samples, corrected_norm)))
    return scaled


def _sum_by_samples(seqs_freq, samples_order):
    n = len(seqs_freq[seqs_freq.keys()[0]].freq.keys())
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
            setclus = anncluster(c, setclus, db, args.type_ann)
    return setclus


def _clean_alignment(args):
    """
    Prepare alignment for cluster detection.
    """
    logger.info("Clean bam file with highly repetitive reads with low counts. sum(counts)/n_hits > 1%")
    bam_file, seq_obj = clean_bam_file(args.afile, args.mask)
    logger.info("Using %s file" % bam_file)
    detect_complexity(bam_file, args.ref)
    return bam_file, seq_obj


def _create_clusters(seqL, bam_file, args):
    clus_obj = []
    logger.info("Parsing aligned file")
    cluster_file = op.join(args.out, "cluster.bed")
    if not os.path.exists(args.out + '/list_obj.pk'):
        logger.info("Merging position")
        if not file_exists(cluster_file):
            aligned_bed = parse_align_file(bam_file)
            a = pybedtools.BedTool(aligned_bed, from_string=True)
            # c = a.merge(o="distinct", c="4,5,6", s=True, d=20)
            c = a.cluster(s=True, d=20)
            c.saveas(cluster_file)
        else:
            c = pybedtools.BedTool(cluster_file)
        logger.info("Creating clusters")
        # clus_obj = parse_merge_file(c, seqL, args.MIN_SEQ)
        clus_obj = detect_clusters(c, seqL, args.min_seqs, args.non_un_gl)
        with open(args.out + '/list_obj.pk', 'wb') as output:
            pickle.dump(clus_obj, output, pickle.HIGHEST_PROTOCOL)
    else:
        logger.info("Loading previous clusters")
        with open(args.out + '/list_obj.pk', 'rb') as input:
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


def _check_args(args):
    logger.info("Checking parameters and files")
    args.dir_out = args.out
    args.samplename = "pro"
    global decision_cluster
    global similar
    if not os.path.isdir(args.out):
        logger.warning("the output folder doens't exists")
        os.mkdirs(args.out)
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
        raise IOError("Seqs.ma or Seqs.bam doesn't exists.")
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


#    contentDivC += table.make_div(expchart.getExpDiv(), "expseqs", "css_exp")
#    contentDivC += table.make_div("<pre>" + clus.showseq_plain + "</pre>", "plainseqs", "css_seqs")
#    tmp_cont, allData, exp = _create_sequences_out(seqListTemp, seq)

#def _create_sequences_out(seqListTemp, clus):
#    allData = "var allData = [{"
#    for s in seqListTemp:
#        allData += '"%s":[' % s
#        seq = clus.seq[s]
#        out.write("S %s %s\n" % (seq.seq, seq.len))
#        #showseq+="S %s %s " % (seq.seq,seq.len)
#        colseqs = [seq.seq, seq.len]
#        for sample in samples_list:
#            allData += '{"sample":"%s","reads":"%s","color":"#FF6600"},' % (sample, int(clus.seq[s].freq[sample]))
#            colseqs.append(int(clus.seq[s].freq[sample]))
#            exp += int(clus.seq[s].freq[sample])
#            if (not freq.has_key(sample)):
#                freq[sample] = 0
#            freq[sample] += int(clus.seq[s].freq[sample])
#        allData = allData[:-1] + "],"
#        contentDivS += table.make_line("".join(map(table.make_cell, colseqs)))
#
#    allData = allData[:-1] + "}];"

