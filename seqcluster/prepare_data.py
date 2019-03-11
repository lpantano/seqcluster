from __future__ import print_function
#import sys
import os
import os.path as op
import traceback
#from os import listdir
#from os.path import isfile, join
import re
import logging
from seqcluster.libs.classes import sequence_unique
from seqcluster.libs.classes import quality
from seqcluster.libs.fastq import is_fastq, open_fastq


logger = logging.getLogger('prepare')


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
        f = open(args.config, 'r')
        seq_out = open(op.join(args.out, "seqs.fastq"), 'w')
        ma_out = open(op.join(args.out, "seqs.ma"), 'w')
    except IOError as e:
        traceback.print_exc()
        raise IOError("Can not create output files: %s, %s or read %s" % (op.join(args.out, "seqs.ma"), op.join(args.out, "seqs.fastq"), args.config))
    logger.info("Reading sequeces")
    seq_l, sample_l = _read_fastq_files(f, args)
    logger.info("Creating matrix with unique sequences")
    logger.info("Filtering: min counts %s, min size %s, max size %s, min shared %s" % (args.minc, args.minl, args.maxl, args.min_shared))
    _create_matrix_uniq_seq(sample_l, seq_l, ma_out, seq_out, args.min_shared)
    logger.info("Finish preprocessing. Get a sorted BAM file of seqs.fa and run seqcluster cluster.")


def _read_fasta_files(f, args):
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
        with open(cols[0], 'r') as fasta:
            sample_l.append(cols[1])
            for line in fasta:
                if line.startswith(">"):
                    idx += 1
                    counts = int(re.search("x([0-9]+)", line.strip()).group(1))
                else:
                    seq = line.strip()
                    seq = seq[0:int(args.maxl)] if len(seq) > int(args.maxl) else seq
                    if counts > int(args.minc) and len(seq) > int(args.minl):
                        if seq not in seq_l:
                            seq_l[seq] = sequence_unique(idx, seq)
                        seq_l[seq].add_exp(cols[1], counts)
    return seq_l, sample_l


def _read_fastq_files(f, args):
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
    p = re.compile("^[ATCGNU]+$")
    with open(op.join(args.out, "stats_prepare.tsv"), 'w') as out_handle:
        for line1 in f:
            line1 = line1.strip()
            cols = line1.split("\t")
            # if not is_fastq(cols[0]):
            #    raise ValueError("file is not fastq: %s" % cols[0])
            with open_fastq(cols[0]) as handle:
                sample_l.append(cols[1])
                total = added = 0
                line = handle.readline()
                while line:
                    if line.startswith("@") or line.startswith(">"):
                        seq = handle.readline().strip()
                        if not p.match(seq):
                            continue
                        idx += 1
                        total += 1
                        keep = {}
                        counts = int(re.search("x([0-9]+)", line.strip()).group(1))
                        if is_fastq(cols[0]):
                            handle.readline().strip()
                            qual = handle.readline().strip()
                        else:
                            qual = "I" * len(seq)
                        qual = qual[0:int(args.maxl)] if len(qual) > int(args.maxl) else qual
                        seq = seq[0:int(args.maxl)] if len(seq) > int(args.maxl) else seq
                        if counts > int(args.minc) and len(seq) > int(args.minl):
                            added += 1
                            if seq in keep:
                                keep[seq].update(qual)
                            else:
                                keep[seq] = quality(qual)
                            if seq not in seq_l:
                                seq_l[seq] = sequence_unique(idx, seq)
                            seq_l[seq].add_exp(cols[1], counts)
                            seq_l[seq].quality = keep[seq].get()
                    line=handle.readline()
                print("total\t%s\t%s" % (idx, cols[1]), file=out_handle, end="")
                print("added\t%s\t%s" % (len(seq_l), cols[1]), file=out_handle, end="")
                logger.info("%s: Total read %s ; Total added %s" % (cols[1], idx, len(seq_l)))
    return seq_l, sample_l


def _create_matrix_uniq_seq(sample_l, seq_l, maout, out, min_shared):
    """ create matrix counts for each different sequence in all the fasta files

    :param sample_l: :code:`list_s` is the output of :code:`_read_fasta_files`
    :param seq_l: :code:`seq_s` is the output of :code:`_read_fasta_files`
    :param maout: is a file handler to write the matrix count information
    :param out: is a file handle to write the fasta file with unique sequences

    :returns: Null
    """
    skip = 0
    if int(min_shared) > len(sample_l):
        min_shared = len(sample_l)
    maout.write("id\tseq")
    for g in sample_l:
        maout.write("\t%s" % g)
    for s in seq_l.keys():
        seen = sum([1 for g in seq_l[s].group if seq_l[s].group[g] > 0])
        if seen < int(min_shared):
            skip += 1
            continue
        maout.write("\nseq_%s\t%s" % (seq_l[s].idx, seq_l[s].seq))
        for g in sample_l:
            if g in seq_l[s].group:
                maout.write("\t%s" % seq_l[s].group[g])
            else:
                maout.write("\t0")
        qual = "".join(seq_l[s].quality)
        out.write("@seq_%s\n%s\n+\n%s\n" % (seq_l[s].idx, seq_l[s].seq, qual))
    out.close()
    maout.close()
    logger.info("Total skipped due to --min-shared parameter (%s) : %s" % (min_shared, skip))
