"""functions for explore tool"""
from __future__ import print_function
import json, itertools
import tempfile, os, contextlib, shutil
from sam2bed import makeBED
import logging
from do import find_cmd, run
import pysam

logger = logging.getLogger('read')


@contextlib.contextmanager
def make_temp_directory():
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


def load_data(in_file):
    """load json file from seqcluster cluster"""
    with open(in_file) as in_handle:
        return json.load(in_handle)


def get_sequences_from_cluster(c1, c2, data):
    """get all sequences from on cluster"""
    seqs1 = data[c1]['seqs']
    seqs2 = data[c2]['seqs']
    seqs = list(set(seqs1 + seqs2))
    names = []
    for s in seqs:
        if s in seqs1 and s in seqs2:
            names.append("both")
        elif s in seqs1:
            names.append(c1)
        else:
            names.append(c2)
    return seqs, names


def get_precursors_from_cluster(c1, c2, data):
    loci1 = data[c1]['loci']
    loci2 = data[c2]['loci']
    return {c1:loci1, c2:loci2}


def map_to_precursors(seqs, names, loci, args):
    """map sequences to precursors with bowtie"""
    with make_temp_directory() as temp:
        pre_fasta = os.path.join(temp, "pre.fa")
        seqs_fasta = os.path.join(temp, "seqs.fa")
        out_sam = os.path.join(temp, "out.fa")
        pre_fasta = get_loci_fasta(loci, pre_fasta, args.ref)
        seqs_fasta = get_seqs_fasta(seqs, names, seqs_fasta)
        if find_cmd("bowtie-build"):
            cmd = "bowtie-build -f {pre_fasta} {temp}/pre"
            run(cmd.format(**locals()))
            cmd = "bowtie {temp}/pre -f {seqs_fasta} -S > {out_sam}"
            run(cmd.format(**locals()))
            run("head {0}".format(out_sam))
            read_alignment(out_sam, loci, seqs)
    return True


def get_seqs_fasta(seqs, names, out_fa):
    """get fasta from sequences"""
    with open(out_fa, 'w') as fa_handle:
        for s, n in itertools.izip(seqs, names):
            print(">c{1}-{0}\n{0}".format(s, n), file=fa_handle)
    return out_fa


def get_loci_fasta(loci, out_fa, ref):
    """get fasta from precursor"""
    if not find_cmd("bedtools"):
        raise ValueError
    with make_temp_directory() as temp:
        bed_file = os.path.join(temp, "file.bed")
        with open(bed_file, 'w') as bed_handle:
            for nc, loci in loci.iteritems():
                for l in loci:
                    logger.info("get_fasta: loci %s" % l)
                    c, s, e = l
                    print("{0}\t{1}\t{2}\t{3}".format(c, s, e, nc), file=bed_handle)
        cmd = "bedtools getfasta -fi {ref} -bed {bed_file} -fo {out_fa}"
        run(cmd.format(**locals()))
    return out_fa


def read_alignment(out_sam, loci, seqs):
    """read which seqs map to which loci and
    return a tab separated file"""
    samfile = pysam.Samfile(out_sam, "r")
    for a in samfile.fetch():
        a = makeBED(a)


def plot_positions():
	return True
