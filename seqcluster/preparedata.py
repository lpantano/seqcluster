import sys
import os
from os import listdir
from os.path import isfile, join
import re
import logging


#define class object to store seq info
class seq_obj:
    def __init__(self,idx,seq):
    	self.idx=idx
    	self.seq=seq
        self.group={}
    def addExp(self,gr,exp):
    	self.group[gr]=exp


def prepare(args, con, log):
	"""
	Read all seq.fa files and create a matrix and unique fasta files
	
	:param args: options parsed from command line
	"""
	try:
		f = open(args.dir, 'r')
		seq_out = open(os.path.join(args.out,"seqs.fa"), 'w')
		ma_out = open(os.path.join(args.out,"seqs.ma"), 'w')
	except IOError as e:
	    con.error("I/O error({0}): {1}".format(e.errno, e.strerror))

	name=""
	#read file with list of fasta files
	con.info("reading sequeces")
	seq_l, list_s = _read_fasta_files(f)
	#create matrix of counts for each seq in each sample
	#create reduced fasta files, with only unique sequences
	con.info("creating matrix with unique sequences")
	_create_matrix_uniq_seq(list_s,seq_l,ma_out,seq_out)
	con.info("Finish preprocessing. Get an SAM file of seqs.fa and run seqcluster cluster.")

def _read_fasta_files(f):
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
		        if not seq_l.has_key(seq):
		        	seq_l[seq]=seq_obj(idx,seq)
		        seq_l[seq].addExp(cols[1],counts.group(1))

	f.close()
	return (seq_l,list_s)



def _create_matrix_uniq_seq(list_s,seq_l,maout,out):
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
