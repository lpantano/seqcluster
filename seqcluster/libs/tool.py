from collections import OrderedDict
import operator
import os
import copy
from sam2bed import makeBED
#from . import table,barchart,seqviz,expchart
import time
import math
import pysam
import logging
from classes import *

logger = logging.getLogger(__name__)


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


def show_seq(clus_obj,index):
    ##create vizualization of sequence along precursor for the results
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


def anncluster(c,clus_obj,db,type_ann):
    ##intersect positions of clusters with annotation files
    if type_ann == "gtf":
        id_sa = 1
        id_ea = 2
        id_sb = 9
        id_eb = 10
        id_id = 3
        id_idl = 4
        id_sta = 5
        id_stb = 12
        id_tag = 8
    elif type_ann == "bed":
        id_sa = 1
        id_ea = 2
        id_sb = 7
        id_eb = 8
        id_id = 3
        id_idl = 4
        id_sta = 5
        id_stb = 11
        id_tag = 9

    ida  =  0
    clus_id = clus_obj.clus
    loci_id = clus_obj.loci
    for cols in c.features():
        id = int(cols[id_id])
        idl = int(cols[id_idl])
        if (clus_id.has_key(id)):
            clus = clus_id[id]
            strd = "-"
            sa = int(cols[id_sa])
            ea = int(cols[id_ea])
            sb = int(cols[id_sb])
            eb = int(cols[id_eb])
            if cols[id_sta] in cols[id_stb]:
                strd = "+"
            if cols[id_sta] in "+" and cols[id_stb] in "+":
                lento5 = sa-sb+1
                lento3 = ea-eb+1
            if cols[id_sta] in "+" and cols[id_stb] in "-":
                lento5 = ea-sb+1
                lento3 = sa-eb+1
            if cols[id_sta] in "-" and cols[id_stb] in "+":
                lento5 = sa-eb+1
                lento3 = ea-sb+1
            if cols[id_sta] in "-" and cols[id_stb] in "-":
                lento3 = sa-sb+1
                lento5 = ea-eb+1
            else:
                lento5 = sa-sb+1
                lento3 = ea-eb+1
            #"EXISTS clus %s with db DBA %s" % (cols[9],db)
            ida +=  1
            if db in loci_id[idl].db_ann:
                # "NEW clus %s with db DBA %s" % (cols[9],db)
                ann = annotation(db,cols[id_tag],strd,lento5,lento3)
                tdb = loci_id[idl].db_ann[db]
                tdb.add_db_ann(ida,ann)
                loci_id[idl].add_db(db,tdb)
                #clus.adddb(db,tdb)
            else:
                # "UPDATE clus %s with db DBA %s" % (cols[9],db)
                ann = annotation(db,cols[id_tag],strd,lento5,lento3)
                tdb = dbannotation(1)
                tdb.add_db_ann(ida,ann)
                loci_id[idl].add_db(db,tdb)
                #clus.adddb(db,tdb)

            clus_id[id] = clus
    clus_obj.clus = clus_id
    clus_obj.loci = loci_id
    return clus_obj


def parse_align_file(file_in, format):
    #parse sam files with aligned sequences
    loc_id = 1
    #clus_obj = {}
    bedfile_clusters = ""
    if format == "sam":
        samfile = pysam.Samfile( file_in, "r" )
        for a in samfile.fetch():
            loc_id += 1
            a = makeBED(a)
            if a:
                bedfile_clusters += "%s\t%s\t%s\t%s\t%s\t%s\n" % \
                    (samfile.getrname(int(a.chr)), a.start, a.end, a.name, loc_id, a.strand)
        samfile.close()
    elif format == "bed":
        f = open(file_in, 'r')
        for line in f:
            loc_id += 1
            #line = processSAM(line)
            a = bedaligned(line)
            #print "%s\t%s\t%s\t%s\t%s\t%s\n" % (a.chr,a.start,a.end,a.name,loc_id,a.strand)
            bedfile_clusters += "%s\t%s\t%s\t%s\t%s\t%s\n" % \
                (a.chr, a.start, a.end, a.name, loc_id, a.strand)
        f.close()
    return bedfile_clusters


def parse_merge_file(c, current_seq, MIN_SEQ):
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
        #only keep locus with 10 or more sequences
        if len(a.names) >= MIN_SEQ:
            lindex += 1
            already_in, not_in = _get_seqs_from_cluster(a.names, clus_id)
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


def reduceloci(clus_obj,min_seq,path):
    ##reduce number of loci a cluster has
    filtered = {}
    idcNew = 0
    current = clus_obj.clus
    clus_idt = clus_obj.clusid
    clus_locit = clus_obj.loci
    clus_seqt = clus_obj.seq
    #saved = list()
    #nodes = ""
    #nodesInfo = ""
    for idc in current.keys():
        logger.debug("Cluster analized: %s " % idc)
        clus1 = copy.deepcopy(current[idc])
        nElements = len(clus1.loci2seq)
        logger.debug("Number of elements: %s " % nElements)
        currentElements = nElements+1
        removeSeqs = list()
        cicle = 0
        seqfound = 0
        if nElements < 1000:
            while nElements<currentElements and nElements != 0:
                cicle += 1
                if (cicle % 10) == 0:
                    logger.debug("Number of cicle: %s with nElements %s" % (cicle,nElements))
                ##get loci with more number of sequences
                locilen_sorted = sorted(clus1.locilen.iteritems(), key=operator.itemgetter(1), reverse=True)
                ##get gigest loci
                maxseq = locilen_sorted[0][1]*1.0
                if maxseq > min_seq:
                    clus1, removeSeqs, filtered, idcNew = _solve_loci(clus1, locilen_sorted, removeSeqs, filtered, maxseq, idcNew)
                else:
                    ##if number of sequences < 10, just remove it
                    for (idl, lenl) in locilen_sorted:
                        clus1.loci2seq.pop(idl, "None")
                        clus1.locilen.pop(idl, "None")
                nElements = len(clus1.loci2seq)
        else:
            idcNew += 1 ##creating new cluster id
            logger.debug("Not resolving cluster %s, too many loci. New id %s" % (idc, idcNew))
            locilen_sorted = sorted(clus1.locilen.iteritems(), key=operator.itemgetter(1), reverse=True)
            maxidl = locilen_sorted[0][0]
            #filtered[idcNew] = copy.deepcopy(clus1)
            filtered[idcNew] = cluster(idcNew)
            filtered[idcNew].add_id_member(clus1.loci2seq[maxidl], maxidl)
            filtered[idcNew].id = idcNew
            filtered[idcNew].toomany = nElements
    clus_obj.clus = filtered
    #nodesFiles = open(path+"/nodes.txt",'w')
    #nodesFiles.write(nodes)
    #nodesFiles.close()
    #nodesFiles = open(path+"/nodesInfo.txt",'w')
    #nodesFiles.write(nodesInfo)
    #nodesFiles.close()
    return clus_obj


def _solve_loci(clus1, locilen_sorted, removeSeqs, filtered, maxseq, idcNew):
    first_run = 0
    seqListTemp = list()
    for (idl,lenl) in locilen_sorted:
        print clus1.loci2seq.keys()
        #nodesInfo += "o%s\tblue\n" % (idl)
        tempList = clus1.loci2seq[idl]
        ##this should be remove to merge more clusters
        #tempList = list(set(sorted(tempList)).difference(sorted(removeSeqs)))
        if first_run == 0:
            seqListTemp = tempList
        intersect = list(set(seqListTemp).intersection(tempList))
        common = 0
        if intersect:
            ##this could be change to min to merge more clusters
            common = len(intersect)*1.0/min(len(seqListTemp), len(tempList))
        ##check if number of common sequences is 60% or greater than maxseq
        print "%s %s %s %s %s " % (clus1.id, idl, lenl, maxseq, common)
        if (first_run == 0 and lenl*1.0>0) or (first_run !=0 and lenl*1.0 >= 0.6*maxseq and common*1.0 >= 0.6):
            if first_run == 0:
                idcNew += 1 ##creating new cluster id
                #nodesInfo += "n%s\tyellow\n" % (idcNew)
            if not idcNew in filtered:
                filtered[idcNew] = cluster(idcNew)
            clus1, removeSeqs = _merge_loci_in_cluster(clus1, idl, removeSeqs, filtered[idcNew])
            #nodes += "o%s\tn%s\t-1\t%s\n" % (idl,idcNew, len(removeSeqs))
        else:
            clus1 = _remove_seqs_from_loci(clus1, idl, removeSeqs)
        first_run = 1
    return clus1, removeSeqs, filtered, idcNew

def _merge_loci_in_cluster(clus1, idl, current_seqs, new_clus):
    ##remove sequences from cluster
    tempList = clus1.loci2seq[idl]
    removeSeqs = list(set(tempList).union(current_seqs))
    ##adding sequences to new cluster
    new_clus.add_id_member(list(tempList), idl)
    ##remove locus from cluster
    clus1.loci2seq.pop(idl, "None")
    clus1.locilen.pop(idl, "None")
    return clus1, removeSeqs


def _remove_seqs_from_loci(clus1, idl, current_seqs):
    ##remove all sequences in the other cluster
    tempList = clus1.loci2seq[idl]
    removeSeqs = list(set(tempList).union(current_seqs))
    newList = list(set(sorted(tempList)).difference(sorted(removeSeqs)))
    ##update sequences in this locus
    #seqsGiven = lenl-len(newList)
    #if len(newList) != lenl and len(newList)>0:
        #nodes += "o%s\tn%s\t%s\t%s\n" % (idl,idcNew,seqsGiven,len(removeSeqs))
    clus1.locilen[idl] = len(newList)
    clus1.loci2seq[idl] = newList
    #if no more sequence in cluster, remove it
    if clus1.locilen[idl] == 0:
        #nodes += "o%s\tn%s\t%s\t%s\n" % (idl,idcNew,lenl,len(removeSeqs))
        clus1.loci2seq.pop(idl,"None")
        clus1.locilen.pop(idl,"None")
    return clus1

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


def parse_ma_file(file):
    f  =  open(file, 'r')
    name = ""
    seq_l = {}
    index = 1
    line = f.readline()
    line = line.strip()
    cols = line.split("\t")
    samples = cols[2:]
    for line in f:
        line = line.strip()
        cols = line.split("\t")
        exp = {}
        for i in range(len(samples)):
            exp[samples[i]] = cols[i+2] 
        name = cols[0].replace(">","")
        index = index+1
        seq = cols[1]      
        new_s = sequence(seq,exp,index)
        seq_l[name] = new_s      
    f.close()
    return seq_l
