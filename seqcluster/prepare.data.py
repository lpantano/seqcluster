import sys
import os
from os import listdir
from os.path import isfile, join
import re
from optparse import OptionParser


#define class object to store seq info
class seq_obj:
    def __init__(self,idx,seq):
    	self.idx=idx
    	self.seq=seq
        self.group={}
    def addExp(self,gr,exp):
    	self.group[gr]=exp



try:
	f = open(options.dir, 'r')
	out = open(options.out+"/seqs.fa", 'w')
	maout = open(options.out+"/seqs.ma", 'w')
except IOError as e:
    print "I/O error({0}): {1}".format(e.errno, e.strerror)
    sys.exit(1)


name=""


#read file with list of fasta files
seq_l={}
list_s=[]
idx=1
for line1 in f:
	line1=line1.strip()
	cols=line1.split("\t")
	fasta = open(cols[0],'r')
	list_s.append(cols[1])
	print line1
	for line in fasta:
	    if (">" in line):
	    	idx+=1
	    	line=line.strip()
	        counts = re.search("x([0-9]+)", line)
	        name=line.replace(">","")
	    else:
	    	line=line.strip()
	        seq=line
	        #add seq to repository       
	        if not seq_l.has_key(seq):
	        	seq_l[seq]=seq_obj(idx,seq)
	        seq_l[seq].addExp(cols[1],counts.group(1))

f.close()

#create matrix of counts for each seq in each sample
#create reduced fasta files, with only unique sequences
maout.write("id\tseq")
for g in list_s:
	maout.write("\t%s" % g)
for s in seq_l.keys():
	maout.write("\n>seq_%s\t%s" % (seq_l[s].idx,seq_l[s].seq))
	for g in list_s:
		if seq_l[s].group.has_key(g):
			maout.write("\t%s" % seq_l[s].group[g])
		else:
			maout.write("\t0")
	out.write(">seq_%s\n%s\n" % (seq_l[s].idx,seq_l[s].seq))
out.close()
maout.close()