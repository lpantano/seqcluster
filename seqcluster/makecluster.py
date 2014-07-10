import sys
import os
import re
import operator
#import pybedtools
import logging
import pickle

from libs.tool import *


def  cluster(args,con,log):

    args = _check_args(args,con,log)

    MIN_SEQ=10
    db4js={}
       

    # ############read input files##################################
    con.info("Parsing matrix file")
    seq_l=parse_ma_file(args.ffile)
    # ###############################################################
    clus_obj = _create_clusters(args,con,log)
    # #####################reduce loci when possible#############################
    con.info("Solving multi-mapping events in the network of clusters")
    clus_red=reduceloci(clus_obj,MIN_SEQ,dir_out,log)
    con.info("Clusters up to %s" % (len(clus_red.clus.keys())))
    # ###########################################################################
    # #####################create sequences overview ############################
    if args.show:
        con.info("Creating sequences alignment to precursor")
        clus_red=show_seq(clus_red,args.index)
    # ###########################################################################
    # #####################overlap with features#################################
    clus_red=_annotate(args,clus_red,con,log)
    # #############################################################
    _create_html(args,con,log)
    con.info("Finished")

def _create_html(args,clus_red,con,log):
    filtered=clus_red.clus
    clus_locit=clus_red.loci

    f=open(args.ffile) 
    line=f.readline()
    line=line.strip()
    cols=line.split("\t")
    samples_list=cols[2:]
    f.close

    # #####################creating html####################
    os.makedirs(os.path.join(args.dir_out,"html"))
    os.makedirs(os.path.join(args.dir_out,"html","clusters"))

    pathscript=os.path.abspath(__file__).replace("seqcluster.py","")
    cmdcp="rsync -a -u --exclude='- *.' "+pathscript+"js "+pathscript+"css "+pathscript+"images "+dir_out+"/html/."
    #print bcolors.OKBLUE+cmdcp+bcolors.ENDC 
    os.system(cmdcp)

    chtml = open(os.path.join(args.dir_out,"html","clus.html"), 'w')
    dbhtml = open(os.path.join(args.dir_out,"html","annotation.html"), 'w')
    ccont=table.make_header("".join(map(table.make_cell_header,["id","DB"]+samples_list)))
    icont=""

    con.info("Creating outputs")
    #####################creating plain text and html files#####################
    out = open(os.path.join(args.dir_out,"clus.parse.txt"), 'w')
    outann = open(os.path.join(args.dir_out,"ann.tab"), 'w')
    outann.write("\t".join(["id","ann","manyloci"]+samples_list))
    outann.write("\n")
    outbed=os.path.join(dir_out,"clus.bed")

    db4js["None"]=[0,0,0]
    for id in filtered.keys():
        dbsummary=""
        #if (int(id)!=0):
        clus=filtered[id]
        # "id %s" % (id)
        if (len(clus.loci2seq)>0):
            icluster_html = open(os.path.join(dir_out,"html","clusters",str(id)+".html"), 'w')
            cluster_id="C:%s" % id
            outann.write("%s\t" % cluster_id)
            ccell=table.make_cell_link(cluster_id,"clusters/%s.html" % id)
            contentDivC=table.make_div(table.make_a("<b>"+cluster_id+"</b>",cluster_id),"clus","css_id")
            out.write("C %s\n" % id)
            numloci=0
            numlocidb=init_numlocidb(beds)
            
            contentDivL=table.make_header("".join(map(table.make_cell_header,table.HEADER_L)))
            seqListTemp=()

            for idl in clus.loci2seq.keys():
                # "idl %s" % (idl)
                seqListTemp=list(set(seqListTemp).union(filtered[id].loci2seq[idl]))
                tpos=clus_locit[idl]
                numloci+=1
                out.write("L %s %s %s %s %s\n" % (tpos.idl,tpos.chr,tpos.start,tpos.end,tpos.strand))
                #outtemp.write("%s\t%s\t%s\t%s\t%s\t%s\n" %  (tpos.chr,tpos.start,tpos.end,id,0,tpos.strand))
                idl=int(tpos.idl)
                contentA=""
                hits={}
                for db in clus.dbann.keys():
                    # "DBA %s" % (db)
                    tdb=clus.dbann[db]
                    if not hits.has_key(db):
                        hits[db]=""
                    if (tdb.ann.has_key(idl)):
                        ann=tdb.ann[idl]
                        numlocidb[db]+=1
                        contentA+="%s(%s %s)-(%s,%s);" % (ann.db,ann.name,ann.strand,ann.to5,ann.to3);
                        out.write("A %s %s %s %s %s\n" % (ann.db,ann.name,ann.strand,ann.to5,ann.to3))
                        hits[db]+=ann.name+","
                        # "%s %s => %s" % (db,ann.name,hits[db])
                        #print result for consistent DB
                pos="%s:%s-%s %s" % (tpos.chr,tpos.start,tpos.end,tpos.strand)
                contentDivL+=table.make_line("".join(map(table.make_cell,[pos,contentA])))

            contentDivC+=table.make_div(table.make_table(contentDivL,"loci"),"loci","css_loci")
            #contentDivC+=table.make_div(table.make_hs_link("plainseqs"),"linkshowhide","css_seqs_link")
            contentDivC+=table.make_div(expchart.getExpDiv(),"expseqs","css_exp")
            contentDivC+=table.make_div("<pre>"+clus.showseq_plain+"</pre>","plainseqs","css_seqs")
            
            #contentDivC+=seqviz.CANVAS
            ##print seq_header
            exp=0
            seq_header="".join(map(table.make_cell_header,table.HEADER_SEQ+samples_list)) 
            contentDivS=table.make_header(seq_header)
            freq={}
            showseq=""
            allData="var allData = [{"
            for s in seqListTemp:
                allData+='"%s":[' % s
                seq=clus_red.seq[s]
                out.write("S %s %s\n" % (seq.seq,seq.len))
                #showseq+="S %s %s " % (seq.seq,seq.len)
                colseqs=[seq.seq,seq.len]
                for sample in samples_list:
                    allData+='{"sample":"%s","reads":"%s","color":"#FF6600"},' % (sample,int(clus_red.seq[s].freq[sample]))
                    colseqs.append(int(clus_red.seq[s].freq[sample]))
                    exp+=int(clus_red.seq[s].freq[sample])
                    if (not freq.has_key(sample)):
                        freq[sample]=0
                    freq[sample]+=int(clus_red.seq[s].freq[sample])
                allData=allData[:-1]+"],"
                contentDivS+=table.make_line("".join(map(table.make_cell,colseqs)))
            
            allData=allData[:-1]+"}];"

            contentDivC+=table.make_div("<pre>"+table.make_table(contentDivS,"seqs")+"</pre>","seqs","css_table")
            cons=0
            ann=0
            tmp_db4js={}
            for filebed in beds:
                db= os.path.basename(filebed)
                ratio=1.0*numlocidb[db]/numloci

                out.write("DBstat %s %s %s %s \n" % (db,numloci,numlocidb[db],ratio))
                
                if (numlocidb[db]>0):
                    
                    ann=1
                    tmp_db4js[db]=1
                    if(ratio>=0.33):
                        out.write("DB %s %s %s consistent\n" % (db,numloci,numlocidb[db]))
                        dbsummary+="DB(%s) %s/%s consistent;" % (db,numlocidb[db],numloci)
                        cons+=1
                        outann.write("%s:%s;" % (db,hits[db]))
                        
                    elif (ratio<0.33):
                       
                        out.write("DB %s %s %s no-consistent\n" % (db,numloci,numlocidb[db]))
                        dbsummary+="DB(%s) %s/%s consistent;" % (db,numlocidb[db],numloci)

            for db in tmp_db4js:
                if cons==0: 
                    db4js[db][2]+=exp
                if cons>1: 
                    db4js[db][1]+=exp
                if cons==1: 
                    db4js[db][0]+=exp

            if (ann==0): 
                db4js["None"][0]+=exp


            ccell+=table.make_cell(dbsummary)
            colseqs=[]
            outann.write("\t%s\t" % clus.toomany)
            for sample in samples_list:
                colseqs.append(freq[sample])
                outann.write("%s\t" % freq[sample])

            outann.write("\n")
            ccell+="".join(map(table.make_cell,colseqs))
            ccont+=table.make_line(ccell)
            icont=table.make_div(contentDivC,cluster_id,"css_clus")
            #idiv=make_div(clus,clus,"cluster")
        icont=table.make_html(icont,expchart.addgraph(allData),samplename)         
        icluster_html.write(icont)
        icluster_html.close()

    listbar=[]
    for db in db4js.keys():
        listbar.append({"args":db,"unique":db4js[db][0],"multiple":db4js[db][1],"unconsistently":db4js[db][2]})
    listkeys=["unique","multiple","unconsistently"]

    dbhtml.write(barchart.createhtml(listbar,listkeys))
    dbhtml.close()
    ccont=table.make_table(ccont,samplename)
    ccont=table.make_jshtml(ccont,samplename)         
    chtml.write(ccont)
    chtml.close()
    out.close()
    outann.close()

def _annotate(args,setclus,con,log):
    con.info("Creating bed file")
    bedfile=generate_position_bed(setclus)
    a = pybedtools.BedTool(bedfile,from_string=True)
    beds=[]
    con.info("Annotating clusters")
    if args.list_files:
        beds=list_files.split(",")
        for filebed in beds:
            db= os.path.basename(filebed)
            db4js[db]=[0,0,0]
            b = pybedtools.BedTool(filebed)
            c=a.intersect(b, wo=True)
            setclus=anncluster(c,setclus,db,args.type_ann) 
    return setclus
def _create_clusters(args,con,log):
    clus_obj = []
    if not os.path.exists(args.out+'/list_obj.pk'):
        #merge positions to create clusters: everything connected by 
        #positions on the genome. 
        #assumption: minimun number of common loci
        con.info("Parsing aligned file")
        bed_obj=parse_align_file(args.afile,format)
        a = pybedtools.BedTool(bed_obj,from_string=True)
        c = a.merge(nms=True,d=20,s=True,scores="collapse")
        con.info("Creating clusters")
        clus_obj = parse_merge_file(c,seq_l,MIN_SEQ)
        with open(args.out+'/list_obj.pk', 'wb') as output:
            pickle.dump(clus_obj, output, pickle.HIGHEST_PROTOCOL)
    else:
        con.info("Loading previous clusters")
        with open(args.out+'/list_obj.pk', 'rb') as input:
            clus_obj = pickle.load(input)

    con.info("%s clusters found" % (len(clus_obj.clus.keys())))
    return clus_obj

def _check_args(args,con,log):
    ########################################################
    args.dir_out = args.out
    samplename="pro"

    if not os.path.isdir(args.out):
        print bcolors.FAIL + "the output folder doens't exists" + bcolors.ENDC
        sys.exit(1)
    if args.bed and args.gtf:
        print bcolors.FAIL + "cannot provide -b and -g at the same time" +  bcolors.ENDC
        sys.exit(1)


    #####################define variables####################    

    if args.debug:
        con.info("DEBUG messages will be showed in file.")

    if args.bed:
        args.list_files=args.bed
        args.type_ann="bed"
    if args.gtf:
        args.list_files=args.gtf
        args.type_ann="gtf"

    con.info("Output dir will be: %s" % dir_out) 

    # #try to open all files to avoid future I/O errors
    try:
        f=open(args.ffile,'r')
        f.close()
        f=open(args.afile,'r')
        f.close()
        beds=list_files.split(",")
        for filebed in beds:
            f=open(filebed,'r')
            f.close()
    except IOError as e:
        con.error("I/O error({0}): {1}".format(e.errno, e.strerror))
 
    #####################read aligned sequences#####################
    format = what_is(options.afile)
    if not format:
        logging.error("Format of aligned reads not in sam or bed")
        raise("Format of aligned reads not in sam or bed")

    return args
