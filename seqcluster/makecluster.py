import os
from collections import Counter, namedtuple
from operator import itemgetter
import pybedtools
import pickle
import numpy as np
from bcbio.utils import file_exists
import libs.logger as mylog
from libs.mystats import up_threshold
import json
from libs.cluster import detect_clusters
from libs.tool import parse_ma_file, reduceloci, show_seq, \
    parse_align_file, generate_position_bed, anncluster, _get_seqs, add_seqs_position_to_loci
from libs.classes import *
import libs.parameters as param


logger = mylog.getLogger(__name__)


def cluster(args):
    args = _check_args(args)
    args.MIN_SEQ = 10
    logger.info("Parsing matrix file")
    seqL = parse_ma_file(args.ffile)
    if len(seqL.keys()) < 100:
        logger.error("It seems you have so low coverage. Please check your fastq files have enough sequences.")
        raise ValueError("So few sequences.")
    clusL = _create_clusters(seqL, args)
    logger.info("Solving multi-mapping events in the network of clusters")
    clusLred = reduceloci(clusL, args.dir_out)
    logger.info("Clusters up to %s" % (len(clusLred.clus.keys())))
    if args.show:
        logger.info("Creating sequences alignment to precursor")
        clusLred = show_seq(clusLred, args.index)
    clusLred = _annotate(args, clusLred)
    logger.info("Creating json and count matrix")
    _create_json(clusLred, args)
    logger.info("Output file in: %s" % args.dir_out)
    logger.info("Finished")


def _write_size_table(data_freq, data_len, ann_valid, cluster_id):
    dd = Counter()
    for f, l in zip(data_freq, data_len):
        dd[l] += np.mean(f.values())

    table = ""
    for l in sorted(dd):
        table += "%s %s %s %s\n" % (l, dd[l], ann_valid, cluster_id)
    return table


def _create_json(clusL, args):
    clus = clusL.clus
    seqs = clusL.seq
    loci = clusL.loci
    data_clus = {}
    out_count = os.path.join(args.dir_out, "counts.tsv")
    out_size = os.path.join(args.dir_out, "size_counts.tsv")
    with open(out_count, 'w') as matrix, open(out_size, 'w') as size_matrix:
        matrix.write("id\tann\t%s\n" % "\t".join(list(seqs[seqs.keys()[1]].freq.keys())))
        for cid in clus.keys():
            seqList = []
            c = clus[cid]
            data_loci = map(lambda (x): [x, loci[x].chr, int(loci[x].start), int(loci[x].end), loci[x].strand, len(c.loci2seq[x])], c.loci2seq.keys())
            data_loci = sorted(data_loci, key=itemgetter(5), reverse=True)
            seqList = _get_seqs(c)
            logger.debug("_json_: %s" % seqList)
            data_ann, valid_ann = _get_annotation(c, loci)
            data_seqs = map(lambda (x): {x: seqs[x].seq}, seqList)
            # data_freq = map(lambda (x): seqs[x].freq, seqList)
            scaled_seqs = _get_counts(seqList, seqs, c.idmembers)
            # print data_freq
            data_freq = map(lambda (x): scaled_seqs[x].freq, seqList)
            data_freq_w_id = map(lambda (x): {x: scaled_seqs[x].norm_freq}, seqList)
            data_len = map(lambda (x): seqs[x].len, seqList)
            data_freq_values = map(lambda (x): map(int, scaled_seqs[x].freq.values()), seqList)
            sum_freq = _sum_by_samples(data_freq_values)
            data_ann_str = [["%s::%s" % (name, ",".join(features)) for name, features in k.iteritems()] for k in data_ann]
            data_valid_str = " ".join(valid_ann)
            matrix.write("%s\t%s|%s\t%s\n" % (cid, data_valid_str, ";".join([";".join(d) for d in data_ann_str]), "\t".join(map(str, sum_freq))))
            data_string = {'seqs': data_seqs, 'freq': data_freq_w_id,
                    'loci': data_loci, 'ann': data_ann, 'valid': valid_ann}
            data_clus[cid] = data_string
            size_table = size_matrix.write(_write_size_table(data_freq, data_len, data_valid_str, cid))
    with open(os.path.join(args.dir_out, "seqcluster.json"), 'w') as handle_out:
        handle_out.write(json.dumps([data_clus], skipkeys=True, indent=2))


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


def _sum_by_samples(seqs_freq):
    y = np.array(seqs_freq[0]) * 0
    for x in seqs_freq:
        y = list(np.array(x) + y)
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


def _create_clusters(seqL, args):
    clus_obj = []
    logger.info("Parsing aligned file")
    aligned_bed = parse_align_file(args.afile)
    if not os.path.exists(args.out + '/list_obj.pk'):
        logger.info("Merging position")
        a = pybedtools.BedTool(aligned_bed, from_string=True)
        # c = a.merge(o="distinct", c="4,5,6", s=True, d=20)
        c = a.cluster(s=True, d=20)
        logger.info("Creating clusters")
        # clus_obj = parse_merge_file(c, seqL, args.MIN_SEQ)
        clus_obj = detect_clusters(c, seqL, args.MIN_SEQ)
        with open(args.out + '/list_obj.pk', 'wb') as output:
            pickle.dump(clus_obj, output, pickle.HIGHEST_PROTOCOL)
    else:
        logger.info("Loading previous clusters")
        with open(args.out + '/list_obj.pk', 'rb') as input:
            clus_obj = pickle.load(input)
    bedfile = pybedtools.BedTool(generate_position_bed(clus_obj), from_string=True)
    # seqs_2_loci = bedfile.intersect(pybedtools.BedTool(aligned_bed, from_string=True), wo=True, s=True)
    # seqs_2_position = add_seqs_position_to_loci(seqs_2_loci, seqL)
    logger.info("%s clusters found" % (len(clus_obj.clus.keys())))
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


#def _create_html(args, clus_red, con, log):
#    #db4js={}
#
#    filtered = clus_red.clus
#    clus_locit = clus_red.loci
#    seq = clus_red.seq
#
#    with open(args.ffile) as f:
#        cols = f.readline().strip().split("\t")
#        samples_list = cols[2:]
#
#    # #####################creating html####################
#    con.info("Creating misc for nice html style")
#    os.makedirs(os.path.join(args.dir_out, "html"))
#    os.makedirs(os.path.join(args.dir_out, "html", "clusters"))
#
#    pathscript = os.path.dirname(__file__)
#    shutil.copytree(os.path.join(pathscript, "misc", "js"), os.path.join(args.dir_out, "html"))
#    shutil.copytree(os.path.join(pathscript, "misc", "css"), os.path.join(args.dir_out, "html"))
#    shutil.copytree(os.path.join(pathscript, "misc", "images"), os.path.join(args.dir_out, "html"))
#
#    with open(os.path.join(args.dir_out, "html", "clus.html"), 'w') as chtml:
#        with open(os.path.join(args.dir_out, "html", "annotation.html"), 'w') as dbhtml:
#            ccont = table.make_header("".join(map(table.make_cell_header, ["id", "DB"] + samples_list)))
#            icont = ""
#            con.info("Creating outputs")
#            #####################creating plain text and html files#####################
#            with open(os.path.join(args.dir_out, "clus.json"), 'w') as out:
#                with open(os.path.join(args.dir_out, "ann.tab"), 'w') as outann:
#                    outann.write("\t".join(["id", "ann", "manyloci"] + samples_list))
#                    outann.write("\n")
#                    #outbed=os.path.join(dir_out,"clus.bed")
#                    db4js["None"] = [0, 0, 0]
#                    for id in filtered.keys():
#                        dbsummary = ""
#                        #if (int(id)!=0):
#                        clus = filtered[id]
#                        # "id %s" % (id)
#                        if (len(clus.loci2seq) > 0):
#                            with open(os.path.join(args.dir_out, "html", "clusters", str(id) + ".html"),
#                                      'w') as icluster_html:
#                                ccell = table.make_cell_link("C:%sclusters/%s.html" % (id, id))
#                                ccont += table.make_line(ccell)
#                                tmp_cont, outannD, outD = _create_cluster_out(id, clus, seq, args, con, log)
#                                icluster_html.write(tmp_cont)
#                                out.write(outD)
#                                outann.write(outannD)
#
#                    listbar = []
#                    for db in db4js.keys():
#                        listbar.append({"args": db, "unique": db4js[db][0], "multiple": db4js[db][1],
#                                        "unconsistently": db4js[db][2]})
#                    listkeys = ["unique", "multiple", "unconsistently"]
#                    dbhtml.write(barchart.createhtml(listbar, listkeys))
#
#                    ccont = table.make_table(ccont, args.samplename)
#                    ccont = table.make_jshtml(ccont, args.samplename)
#                    chtml.write(ccont)
#
#
#def _create_clusters_out(id, clus, seq, args, con, log):
#    cluster_id = "C:%s" % id
#    outann = ("%s\t" % cluster_id)
#    ccell = table.make_cell_link(cluster_id, "clusters/%s.html" % id)
#    contentDivC = table.make_div(table.make_a("<b>" + cluster_id + "</b>", cluster_id), "clus", "css_id")
#    out = ("C %s\n" % id)
#    numloci = 0
#    numlocidb = init_numlocidb(beds)
#
#    temp_cont, seqListTemp = _create_loci_out()
#    contentDivC += table.make_div(table.make_table(temp_cont, "loci"), "loci", "css_loci")
#    #contentDivC+=table.make_div(table.make_hs_link("plainseqs"),"linkshowhide","css_seqs_link")
#    contentDivC += table.make_div(expchart.getExpDiv(), "expseqs", "css_exp")
#    contentDivC += table.make_div("<pre>" + clus.showseq_plain + "</pre>", "plainseqs", "css_seqs")
#    #contentDivC+=seqviz.CANVAS
#    tmp_cont, allData, exp = _create_sequences_out(seqListTemp, seq)
#    contentDivC += table.make_div("<pre>" + table.make_table(tmp_cont, "seqs") + "</pre>", "seqs", "css_table")
#
#    #dbsummary = _create_db_out()
#    #ccell+=table.make_cell(dbsummary)
#    icont = table.make_div(contentDivC, cluster_id, "css_clus")
#    icont = table.make_html(icont, expchart.addgraph(allData), args.samplename)# change
#    #idiv=make_div(clus,clus,"cluster")
#    return icont, out, outann
#
#
#def _create_loci_out(clus):
#    seqListTemp = ()
#    contentDivL = ""
#    contentDivL = table.make_header("".join(map(table.make_cell_header, table.HEADER_L)))
#    for idl in clus.loci2seq.keys():
#        # "idl %s" % (idl)
#        seqListTemp = list(set(seqListTemp).union(filtered[id].loci2seq[idl]))
#        tpos = clus_locit[idl]
#        numloci += 1
#        out = ("%sL %s %s %s %s %s\n" % (out, tpos.idl, tpos.chr, tpos.start, tpos.end, tpos.strand))
#        #outtemp.write("%s\t%s\t%s\t%s\t%s\t%s\n" %  (tpos.chr,tpos.start,tpos.end,id,0,tpos.strand))
#        idl = int(tpos.idl)
#        contentA = ""
#        hits = {}
#        for db in clus.dbann.keys():
#            # "DBA %s" % (db)
#            tdb = clus.dbann[db]
#            if not hits.has_key(db):
#                hits[db] = ""
#            if (tdb.ann.has_key(idl)):
#                ann = tdb.ann[idl]
#                numlocidb[db] += 1
#                contentA += "%s(%s %s)-(%s,%s);" % (ann.db, ann.name, ann.strand, ann.to5, ann.to3);
#                out.write("A %s %s %s %s %s\n" % (ann.db, ann.name, ann.strand, ann.to5, ann.to3))
#                hits[db] += ann.name + ","
#                # "%s %s => %s" % (db,ann.name,hits[db])
#                #print result for consistent DB
#        pos = "%s:%s-%s %s" % (tpos.chr, tpos.start, tpos.end, tpos.strand)
#        contentDivL += table.make_line("".join(map(table.make_cell, [pos, contentA])))
#    return contentDivL, seqListTemp
#
#
#def _create_db_out():
#    cons = 0
#    ann = 0
#    tmp_db4js = {}
#    for filebed in args.list_files:
#        db = os.path.basename(filebed)
#        ratio = 1.0 * numlocidb[db] / numloci
#
#        out = ("%sDBstat %s %s %s %s \n" % (out, db, numloci, numlocidb[db], ratio))
#
#        if (numlocidb[db] > 0):
#
#            ann = 1
#            tmp_db4js[db] = 1
#            if (ratio >= 0.33):
#                out.write("DB %s %s %s consistent\n" % (db, numloci, numlocidb[db]))
#                dbsummary += "DB(%s) %s/%s consistent;" % (db, numlocidb[db], numloci)
#                cons += 1
#                outann = ("%s%s:%s;" % (outann, db, hits[db]))
#
#            elif (ratio < 0.33):
#
#                out.write("DB %s %s %s no-consistent\n" % (db, numloci, numlocidb[db]))
#                dbsummary += "DB(%s) %s/%s consistent;" % (db, numlocidb[db], numloci)
#
#    for db in tmp_db4js:
#        if cons == 0:
#            db4js[db][2] += exp
#        if cons > 1:
#            db4js[db][1] += exp
#        if cons == 1:
#            db4js[db][0] += exp
#
#    if (ann == 0):
#        db4js["None"][0] += exp
#
#    return dbsummary
#
#
#def _create_sequences_out(seqListTemp, clus):
#    seq_header = "".join(map(table.make_cell_header, table.HEADER_SEQ + samples_list))
#    contentDivS = table.make_header(seq_header)
#    exp = 0
#    freq = {}
#    showseq = ""
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
#
#    colseqs = []
#    outann = "%s\t%s\t" % (outann, clus.toomany)
#    for sample in samples_list:
#        colseqs.append(freq[sample])
#        outann.write("%s\t" % freq[sample])
#
#    outann = ("%s\n" % outann)
#    ccell += "".join(map(table.make_cell, colseqs))
#    ccont += table.make_line(ccell)
#
#    return contentDivS, allData, exp


