#import sys
import os
#from os import listdir
#from os.path import isfile, join
import re
import logging
from libs.classes import sequence_unique


logger = logging.getLogger('seqbuster')


def prepare(args):
    """
	Read all seq.fa files and create a matrix and unique fasta files.
	The information is 

	:param args: options parsed from command line
	:param con: logging messages going to console
	:param log: logging messages going to console and file

	:returns: files - matrix and fasta files that should be used with
	 and aligner (as bowtie) and run `seqcluster cluster`
	"""
    try:
        f = open(args.dir, 'r')
        seq_out = open(os.path.join(args.out, "seqs.fa"), 'w')
        ma_out = open(os.path.join(args.out, "seqs.ma"), 'w')
    except IOError as e:
        logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
    #name = ""
    logger.info("reading sequeces")
    seq_l, sample_l = _read_fasta_files(f)
    logger.info("creating matrix with unique sequences")
    _create_matrix_uniq_seq(sample_l, seq_l, ma_out, seq_out)
    logger.info("Finish preprocessing. Get an SAM file of seqs.fa and run seqcluster cluster.")


def _read_fasta_files(f):
    """ read fasta files of each sample and generate a seq_obj
	with the information of each unique sequence in each sample

	:param f: file containing the path for each fasta file and
	 the name of the sample. Two column format with `tab` as field
	 separator

	:returns: * :code:`seq_l`: is a list of seq_obj objects, containing
	            the information of each sequence
	          * :code:`sample_l`: is a list with the name of the samples
	            (column two of the config file)
	"""
    seq_l = {}
    sample_l = []
    idx = 1
    for line1 in f:
        line1 = line1.strip()
        cols = line1.split("\t")
        fasta = open(cols[0], 'r')
        sample_l.append(cols[1])
        #print line1
        for line in fasta:
            if line.startswith(">"):
                idx += 1
                line = line.strip()
                counts = re.search("x([0-9]+)", line)
                #name = line.replace(">", "")
            else:
                seq = line.strip()
                if not seq in seq_l:
                    seq_l[seq] = sequence_unique(idx, seq)
                seq_l[seq].add_exp(cols[1], counts.group(1))

    f.close()
    return seq_l, sample_l


def _create_matrix_uniq_seq(sample_l, seq_l, maout, out):
    """ create matrix counts for each different sequence in all the fasta files

	:param sample_l: :code:`list_s` is the output of :code:`_read_fasta_files`
	:param seq_l: :code:`seq_s` is the output of :code:`_read_fasta_files`
	:param maout: is a file handler to write the matrix count information
	:param out: is a file handle to write the fasta file with unique sequences

	:returns: Null
	"""
    maout.write("id\tseq")
    for g in sample_l:
        maout.write("\t%s" % g)
    for s in seq_l.keys():
        maout.write("\n>seq_%s\t%s" % (seq_l[s].idx, seq_l[s].seq))
        for g in sample_l:
            if seq_l[s].group.has_key(g):
                maout.write("\t%s" % seq_l[s].group[g])
            else:
                maout.write("\t0")
        out.write(">seq_%s\n%s\n" % (seq_l[s].idx, seq_l[s].seq))
    out.close()
    maout.close()
