from __future__ import print_function
import os
from collections import Counter, defaultdict
from seqcluster.libs.classes import quality, umi
from itertools import product
import gzip
import re
import logging


logger = logging.getLogger('seqbuster')

def collapse(in_file):
    """collapse identical sequences and keep Q"""
    keep = Counter()
    with open_fastq(in_file) as handle:
        line = handle.readline()
        while line:
            if line.startswith("@"):
                if line.find("UMI") > -1:
                    logger.info("Find UMI tags in read names, collapsing by UMI.")
                    return collapse_umi(in_file)
                seq = handle.readline().strip()
                handle.readline()
                qual = handle.readline().strip()
                if seq in keep:
                    keep[seq].update(qual)
                else:
                    keep[seq] = quality(qual)
            if line.startswith(">"):
                seq = handle.readline().strip()
                if seq not in keep:
                    keep[seq] = quality("".join(["I"] * len(seq)))
                else:
                    keep[seq].update("".join(["I"] * len(seq)))
            line = handle.readline()
    logger.info("Sequences loaded: %s" % len(keep))
    return keep


def collapse_umi(in_file):
    """collapse reads using UMI tags"""
    keep = defaultdict(dict)
    with open_fastq(in_file) as handle:
        line = handle.readline();
        while line:
            if line.startswith("@"):
                m = re.search('UMI_([ATGC]*)', line.strip())
                umis = m.group(0)
                seq = handle.readline().strip()
                handle.readline()
                qual = handle.readline().strip()
                if (umis, seq) in keep:
                    keep[(umis, seq)][1].update(qual)
                    keep[(umis, seq)][0].update(seq)
                else:
                    keep[(umis, seq)] = [umi(seq), quality(qual)]
            line = handle.readline();
    logger.info("Sequences loaded: %s" % len(keep))
    return keep


def open_fastq(in_file):
    """ open a fastq file, using gzip if it is gzipped
    from bcbio package
    """
    _, ext = os.path.splitext(in_file)
    if ext == ".gz":
        return gzip.open(in_file, 'rt')
    if ext in [".fastq", ".fq", ".fasta", ".fa"]:
        return open(in_file, 'r')
    return ValueError("File needs to be fastq|fasta|fq|fa [.gz]")


def is_fastq(in_file):
    """copy from bcbio package"""
    fastq_ends = [".txt", ".fq", ".fastq"]
    zip_ends = [".gzip", ".gz"]
    base, first_ext = os.path.splitext(in_file)
    second_ext = os.path.splitext(base)[1]
    if first_ext in fastq_ends:
        return True
    elif (second_ext, first_ext) in product(fastq_ends, zip_ends):
        return True
    else:
        return False


def splitext_plus(f):
    """Split on file extensions, allowing for zipped extensions.
    copy from bcbio
    """
    base, ext = os.path.splitext(f)
    if ext in [".gz", ".bz2", ".zip"]:
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base, ext


def write_output(out_file, seqs, minimum=1, size=15):
    idx =0
    logger.info("Writing %s sequences to %s" % (len(seqs.keys()), out_file))
    with open(out_file, 'w') as handle:
        for s in seqs:
            idx += 1
            if isinstance(seqs[s], list):
                seq = seqs[s][0].get()
                qual = "".join(seqs[s][1].get())
                counts = seqs[s][0].times
            else:
                seq = s
                qual = "".join(seqs[s].get())
                counts = seqs[s].times
            if int(counts) >= minimum and len(seq) > size:
                handle.write(("@seq_{idx}_x{counts}\n{seq}\n+\n{qual}\n").format(**locals()))
    return out_file
