# Re-aligner small RNA sequence from SAM/BAM file (miRBase annotation)
import os.path as op

import pysam
from bcbio import bam

import seqcluster.libs.logger as mylog
from realign import *

logger = mylog.getLogger(__name__)

def _get_pos(string):
    name = string.split(":")[0][1:]
    pos = string.split(":")[1][:-1].split("-")
    return name, pos

def _read_mature(matures, sps):
    mature = {}
    with open(matures) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                name = line.strip().replace(">", " ").split()
                mir5p = _get_pos(name[2])
                if len(name) > 3:
                    mir3p = _get_pos(name[3])
    return mature



def _read_precursor(precursor, sps):
    """
    read precurso file for that species
    """
    hairpin = {}
    with open(precursor) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                name = line.strip().replace(">", " ").split()[0]
            else:
                hairpin[name] = line.strip()
    return hairpin

def _download_mirbase(precursor, reference, version="22"):
    """
    Download files from mirbase
    """
    if not precursors or not reference:
        print "Working with version %s" % version

def _annotate(read, mirbase_ref):
    """
    Using SAM/BAM coordinates, mismatches and realign to annotate isomiRs
    """

def _realign(seq, precursor, start):
    """
    The actual fn that will realign the sequence
    """
    # print seq
    # print precursor[start:]

def _realign_hits(read):
    """
    realign read to the list precursors
    """
    for genome in read.precursors:
        new_hit = _realign(read.sequence, genome)

def _sort_by_name(bam_fn):
    """
    sort bam file by name sequence
    """

def _read_bam(bam_fn, precursors):
    """
    read bam file and perform realignment of hits
    """
    handle = bam.open_samfile(bam_fn)
    reads = defaultdict(realign)
    for line in handle:
        chrom = handle.getrname(line.reference_id)
        # print "%s %s %s %s" % (line.query_name, line.reference_start, line.query_sequence, chrom)
        if line.query_name not in reads:
            reads[line.query_name].sequence = line.query_sequence
        _realign(reads[line.query_name].sequence, precursors[chrom], line.reference_start)


def miraligner(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    config = {"algorithm": {"num_cores":1}}
    precursors = _read_precursor(args.hairpin, args.sps)
    matures = _read_mature(args.mirna, args.sps)
    for bam_fn in args.files:
        logger.info("Reading %s" % bam_fn)
        bam_fn = bam.sam_to_bam(bam_fn, config)
        bam_sort_by_n = op.splitext(bam_fn)[0] + "_sort"
        pysam.sort("-n", bam_fn, bam_sort_by_n)
        reads = _read_bam(bam_sort_by_n + ".bam", precursors)
