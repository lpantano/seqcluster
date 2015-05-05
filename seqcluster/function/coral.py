"""prepare data for CoRaL"""
import os.path as op
import pybedtools

from bcbio import utils
from bcbio.distributed.transaction import tx_tmpdir

from seqcluster.libs.do import run

min_trimmed_read_len = 14
max_trimmed_read_len = 50
seg_threshold = 2
seg_maxgap = 44
seg_minrun = min_trimmed_read_len
antisense_min_reads = 0

def prepare_bam(bam_in, precursors):
    """
    Clean BAM file to keep only position inside the bigger cluster
    """
    # use pybedtools to keep valid positions
    # intersect option with -b bigger_cluster_loci
    a = pybedtools.BedTool(bam_in)
    b = pybedtools.BedTool(precursors)
    c = a.intersect(b, u=True)
    out_file = utils.splitext_plus(op.basename(bam_in))[0] + "_clean.bam"
    c.saveas(out_file)
    return op.abspath(out_file)


def detect_regions(bam_in, out_dir, prefix):
    """
    Detect regions using first CoRaL module
    """
    bam2bigwig_cmd = ("bam_to_bigwig.sh {bam_in} {out_dir}/loci {prefix}")
    segment_bigwig_cmd = ("segment_bigwig_into_loci.sh "
                          "{out_dir}/loci.pos.bigwig "
                          "{out_dir}/locid.neg.bigwig "
                          "{seg_thr} {seg_maxgap} "
                          "{seg_minrun} {out_dir}/loci.bed")
    # with tx_tmpdir() as temp_dir:
    with utils.chdir(out_dir):
        run(bam2bigwig_cmd.format(**locals()), "run bam2bigwig")
        run(segment_bigwig_cmd.format(seg_thr=seg_threshold, seg_minrun=seg_minrun, seg_maxgap=seg_maxgap, **locals()), "run segment_bigwig")


def prepare_ann_file(args):
    """
    Create custom ann_file for Coral
    """


def download_hsa_file(args):
    """
    In case of human, download from server
    """
