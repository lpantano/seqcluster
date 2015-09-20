import operator
import os

import math

import logger as mylog

logger = mylog.getLogger(__name__)


REMOVED = 0


def get_distance(p1, p2):
        if p1.strand == "+":
                d = p1.start-p2.start
        else:
                d = (p1.end-(p1.end-p1.start)) - (p2.end-(p2.end-p2.start))
        return d


def get_ini_str(n):
        # get space values
        return "".join(" " for i in range(n))


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


def _normalize_seqs(s, t):
    """Normalize to RPM"""
    for ids in s:
        obj = s[ids]
        [obj.norm_freq.update({sample: 1.0 * obj.freq[sample] / (t[sample]+1) * 1000000}) for sample in obj.norm_freq]
        s[ids] = obj
    return s
