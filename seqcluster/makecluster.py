import sys
import os
import operator
import pybedtools
import pickle
import shutil
import numpy as np
import libs.logger as mylog
import json
from libs.tool import parse_ma_file, reduceloci, show_seq, what_is, \
    parse_merge_file, parse_align_file, generate_position_bed, anncluster
from libs.classes import *


logger = mylog.getLogger(__name__)


def cluster(args):
    args = _check_args(args)
    args.MIN_SEQ = 10
    logger.info("Parsing matrix file")
    seqL = parse_ma_file(args.ffile)
    clusL = _create_clusters(seqL, args)
    logger.info("Solving multi-mapping events in the network of clusters")
    clusLred = reduceloci(clusL, args.MIN_SEQ, args.dir_out)
    logger.info("Clusters up to %s" % (len(clusLred.clus.keys())))
    if args.show:
        logger.info("Creating sequences alignment to precursor")
        clusLred = show_seq(clusLred, args.index)
    clusLred = _annotate(args, clusLred)
    logger.info("Creating json and count matrix")
    _create_json(clusLred, args)
    logger.info("Output file in: %s" % args.dir_out)
    logger.info("Finished")

def _create_json(clusL, args):
    clus = clusL.clus
    seqs = clusL.seq
    loci = clusL.loci
    data_clus = {}
    with open(os.path.join(args.dir_out, "counts.tsv"),'w') as matrix:
        matrix.write( "id\tann\t%s\n" % "\t".join(list(seqs[seqs.keys()[1]].freq.keys())))
        for cid in clus.keys():
            seqList = []
            c = clus[cid]
            data_loci = map(lambda (x): [loci[x].chr, loci[x].start, loci[x].end], c.loci2seq.keys())
            data_ann_temp = {}
            data_ann= []
            for lid in c.loci2seq:
                loci[lid].chr
                seqList = list(set(seqList).union(c.loci2seq[lid]))
                for dbi in loci[lid].db_ann.keys():
                    data_ann_temp[dbi] = {dbi: map(lambda (x): loci[lid].db_ann[dbi].ann[x].name, loci[lid].db_ann[dbi].ann.keys())}
                data_ann = data_ann + map(lambda (x): data_ann_temp[x], data_ann_temp.keys())
            data_seqs = map(lambda (x): seqs[x].seq, seqList)
            data_freq = map(lambda (x): seqs[x].freq, seqList)
            data_freq_values = map(lambda (x): map(int, seqs[x].freq.values()), seqList)
            sum_freq = _sum_by_samples(data_freq_values)
            data_ann_str = [k.keys() for k in data_ann]
            matrix.write("%s\t%s\t%s\n" % (cid, ";".join([ ";".join(d) for d in data_ann_str]), "\t".join(map(str, sum_freq))))
            data_string = {'seqs': data_seqs, 'freq': data_freq,
                'loci': data_loci, 'ann': data_ann}
            data_clus[cid] = data_string
    with open(os.path.join(args.dir_out, "seqcluster.json"), 'w') as handle_out:
        handle_out.write(json.dumps([data_clus], skipkeys=True, indent=2))


def _sum_by_samples(seqs_freq):
    y = np.array(seqs_freq[0]) * 0
    for x in seqs_freq:
        y = list(np.array(x) + y) 
    return y


def _annotate(args, setclus):
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
            #db4js[db] = [0, 0, 0]
            b = pybedtools.BedTool(filebed)
            c = a.intersect(b, wo=True)
            setclus = anncluster(c, setclus, db, args.type_ann)
    return setclus


def _create_clusters(seqL, args):
    clus_obj = []
    if not os.path.exists(args.out + '/list_obj.pk'):
        #merge positions to create clusters: everything connected by
        #positions on the genome.
        #assumption: minimun number of common loci
        logger.info("Parsing aligned file")
        bed_obj = parse_align_file(args.afile, args.format)
        logger.info("Merging position")
        a = pybedtools.BedTool(bed_obj, from_string=True)
        c = a.merge(o="distinct", c="4,5,6", s=True, d=20)
        logger.info("Creating clusters")
        clus_obj = parse_merge_file(c, seqL, args.MIN_SEQ)
        with open(args.out + '/list_obj.pk', 'wb') as output:
            pickle.dump(clus_obj, output, pickle.HIGHEST_PROTOCOL)
    else:
        logger.info("Loading previous clusters")
        with open(args.out + '/list_obj.pk', 'rb') as input:
            clus_obj = pickle.load(input)
    logger.info("%s clusters found" % (len(clus_obj.clus.keys())))
    return clus_obj


def _check_args(args):
    logger.info("Checking parameters and files")
    args.dir_out = args.out
    args.samplename = "pro"
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
    if hasattr(args, 'list_files'):
        try:
            f = open(args.ffile, 'r')
            f.close()
            f = open(args.afile, 'r')
            f.close()
            beds = args.list_files.split(",")
            for filebed in beds:
                f = open(filebed, 'r')
                f.close()
        except IOError as e:
            logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
    args.format = what_is(args.afile)
    logger.info("Aligned file is in: %s" % args.format)
    if not args.format:
        logging.error("Format of aligned reads not in sam or bed")
        raise "Format of aligned reads not in sam or bed"

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


