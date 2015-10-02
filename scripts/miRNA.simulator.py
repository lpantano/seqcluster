from optparse import OptionParser
import sys
import os
import re
import random
import numpy
from libs.objectsSim import *

def createMirna(string):
        cols = string.split("[")
        idx = 0
        obj=mirna()
        for c in cols[1:]:
                idx += 1
                c=c.replace("]","")
                p=c.split(":")
                chrom=p[0]
                p=p[1].split("-")
                start=int(p[0])
                end=int(p[1])
                if idx == 1:
                        obj.addp5(chrom,start,end)
                else:
                        obj.addp3(chrom,start,end)
        return obj


usagetxt = "usage: %prog  -f precurso.fa -m miRNA.str -n 10 -s hsa"
parser = OptionParser(usage=usagetxt, version="%prog 1.0")
parser.add_option("-f", "--fa", dest="fasta",
                  help="", metavar="FILE")
parser.add_option("-m", "--ma", dest="strfile",
                                help="", metavar="FILE")
parser.add_option("-n", "--num", dest="numsim",
                                help="")
parser.add_option("-s", "--species", dest="species",
                                help="", metavar="FILE")
parser.add_option("-e", "--exp", dest="exp", action="store_true",
                                help="give expression", default=False)

if len(sys.argv)<4:
    parser.print_help()
    sys.exit(1)

(options, args) = parser.parse_args()
species=options.species

pos=open(options.strfile,'r')
listmirna={}
for line in pos:
        if line.find("[")>0:
                line = line.strip()
                name = line.split(" ")
                name = name[0].replace(">","")
                listmirna[name]=createMirna(line)
pos.close()

data = {}
fas = open(options.fasta, 'r')
name = ""
seq = ""
nt = ['A', 'T', 'G', 'C']
for line in fas:
        line = line.strip()
        if (line.find(">")>=0):
                if (len(name)!=0):
                        mir=listmirna[name].p5
                        if mir!=0 and name.find(species)>=0:
                                for rand in range(int(options.numsim)):
                                        randS=random.randint(mir.s-3,mir.s+3)+1
                                        randE=random.randint(mir.e-3,mir.e+3)+1
                                        if randS<1:
                                                randS=1
                                        if randE>mir.e:
                                                randE=mir.e-1
                                        randSeq=seq[randS-1:randE]
                                        randName=name+"_"+mir.id+"_"+str(randS)+":"+str(randE)
                                        randName=randName+"_"+str(randS-mir.s)+":"+str(randE-mir.e)
                                        isMut=random.randint(0,3)
                                        mut_tag = "_mut:null"
                                        if isMut==3:
                                                ntMut=random.randint(0,3)
                                                posMut=random.randint(0,len(randSeq)-1)
                                                tempList=list(randSeq)
                                                if tempList[posMut] == nt[ntMut]:
                                                    ntMut -= 1
                                                    if ntMut < 0 :
                                                        ntMut +=2
                                                tempList[posMut]=nt[ntMut]
                                                randSeq=''.join(tempList)
                                                mut_tag="_mut:"+str(posMut+1)+nt[ntMut]
                                        randName+=mut_tag
                                        isAdd=random.randint(0,2)
                                        add_tag = "_add:null"
                                        if isAdd==2:
                                                posAdd=random.randint(1,3)
                                                randAdd=""
                                                add_tag="_add:"
                                                for numadd in range(posAdd):
                                                        ntAdd=random.randint(0,1)
                                                        randSeq+=nt[ntAdd]
                                                        add_tag+=nt[ntAdd]
                                        randName+=add_tag
                                        exp = ""
                                        if not data.has_key(randSeq):
                                            if options.exp:
                                                trial =random.randint(1, 100)
                                                p = random.randint(1, 50) / 50.0
                                                exp = "_x%s" % numpy.random.negative_binomial(trial, p, 1)[0]
                                            print ">%s%s" % (randName, exp)
                                            print randSeq
                                            data[randSeq]=1
                name = line.split(" ")
                name = name[0].replace(">","")
                seq=""
        else:
                seq+=line.replace("U","T")
