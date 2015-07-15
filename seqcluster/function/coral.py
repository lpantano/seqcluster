"""prepare data for CoRaL"""
import os.path as op
from collections import Counter
import pybedtools

from bcbio import utils
from bcbio.distributed.transaction import tx_tmpdir

from seqcluster.libs.do import run, find_cmd
from seqcluster.function.rnafold import calculate_structure

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


def _select_anno(annotation):
    """
    Select annotation from multiple elements
    """
    options = annotation.split(",")
    return options[0]

def _reorder_columns(bed_file):
    """
    Reorder columns to be compatible with CoRaL
    """
    new_bed = utils.splitext_plus(bed_file)[0] + '_order.bed'
    with open(bed_file) as in_handle:
        with open(new_bed, 'w') as out_handle:
            for line in in_handle:
                cols = line.strip().split("\t")
                cols[3] = _select_anno(cols[3]) + "_" + cols[4]
                cols[4] = "0"
                print >>out_handle, "\t".join(cols)
    return new_bed


def _fix_score_column(cov_file):
    """
    Move counts to score columns in bed file
    """
    new_cov = utils.splitext_plus(cov_file)[0] + '_fix.cov'
    with open(cov_file) as in_handle:
        with open(new_cov, 'w') as out_handle:
            for line in in_handle:
                cols = line.strip().split("\t")
                cols[4] = cols[6]
                print >>out_handle, "\t".join(cols[0:6])
    return new_cov


def detect_regions(bam_in, bed_file, out_dir, prefix):
    """
    Detect regions using first CoRaL module
    """
    bed_file = _reorder_columns(bed_file)
    counts_reads_cmd = ("coverageBed -s -counts -b {bam_in} "
                        "-a {bed_file} | sort -k4,4 "
                        "> {out_dir}/loci.cov")
    # with tx_tmpdir() as temp_dir:
    with utils.chdir(out_dir):
        run(counts_reads_cmd.format(min_trimmed_read_len=min_trimmed_read_len, max_trimmed_read_len=max_trimmed_read_len, **locals()), "Run counts_reads")
        loci_file = _fix_score_column(op.join(out_dir, "loci.cov"))
        return loci_file


def _order_antisense_column(cov_file, min_reads):
    """
    Move counts to score columns in bed file
    """
    new_cov = op.join(op.dirname(cov_file), 'feat_antisense.txt')
    with open(cov_file) as in_handle:
        with open(new_cov, 'w') as out_handle:
            print >>out_handle, "name\tantisense"
            for line in in_handle:
                cols = line.strip().split("\t")
                cols[6] = 0 if cols[6] < min_reads else cols[6]
                print >>out_handle, "%s\t%s" % (cols[3], cols[6])
    return new_cov


def _reads_per_position(bam_in, loci_file, out_dir):
    """
    Create input for compute entropy
    """
    data = Counter()
    a = pybedtools.BedTool(bam_in)
    b = pybedtools.BedTool(loci_file)
    c = a.intersect(b, s=True, bed=True, wo=True)
    for line in c:
        end = int(line[1]) + 1 + int(line[2]) if line[5] == "+" else int(line[1]) + 1
        start = int(line[1]) + 1 if line[5] == "+" else int(line[1]) + 1 + int(line[2])
        side5 = "%s\t5p\t%s" % (line[15], start)
        side3 = "%s\t3p\t%s" % (line[15], end)
        data[side5] += 1
        data[side3] += 1

    counts_reads = op.join(out_dir, 'locus_readpos.counts')
    with open(counts_reads, 'w') as out_handle:
        for k in data:
            print >>out_handle, k

    return counts_reads


def create_features(bam_in, loci_file, reference, out_dir):
    """
    Use feature extraction module from CoRaL
    """
    lenvec_plus = op.join(out_dir, 'genomic_lenvec.plus')
    lenvec_minus = op.join(out_dir, 'genomic_lenvec.minus')
    compute_genomic_cmd = ("compute_genomic_lenvectors "
                           "{bam_in} {lenvec_plus} "
                           "{lenvec_minus} "
                           "{min_len} "
                           "{max_len} ")
    index_genomic_cmd = ("index_genomic_lenvectors "
                         "{lenvec} ")
    genomic_lenvec = op.join(out_dir, 'genomic_lenvec')
    feat_len_file = op.join(out_dir, 'feat_lengths.txt')
    compute_locus_cmd = ("compute_locus_lenvectors "
                         "{loci_file} "
                         "{genomic_lenvec} "
                         "{min_len} "
                         "{max_len} "
                         "> {feat_len_file}")
    cov_S_file = op.join(out_dir, 'loci.cov_anti')
    coverage_anti_cmd = ("coverageBed -S -counts -b "
                         "{bam_in} -a {loci_file} "
                         "> {cov_S_file}")
    feat_posentropy = op.join(out_dir, 'feat_posentropy.txt')
    entropy_cmd = ("compute_locus_entropy.rb "
                   "{counts_reads} "
                   "> {feat_posentropy}")
    with utils.chdir(out_dir):
        run(compute_genomic_cmd.format(min_len=min_trimmed_read_len, max_len=max_trimmed_read_len, **locals()), "Run compute_genomic")
        run(index_genomic_cmd.format(lenvec=lenvec_plus), "Run index in plus")
        run(index_genomic_cmd.format(lenvec=lenvec_minus), "Run index in minus")
        run(compute_locus_cmd.format(min_len=min_trimmed_read_len, max_len=max_trimmed_read_len, **locals()), "Run compute locus")
        run(coverage_anti_cmd.format(**locals()), "Run coverage antisense")
        feat_antisense = _order_antisense_column(cov_S_file, min_trimmed_read_len)

        counts_reads = _reads_per_position(bam_in, loci_file, out_dir)
        run(entropy_cmd.format(**locals()), "Run entropy")

        rnafold = calculate_structure(loci_file, reference)


def prepare_ann_file(args):
    """
    Create custom ann_file for Coral
    """


def download_hsa_file(args):
    """
    In case of human, download from server
    """
