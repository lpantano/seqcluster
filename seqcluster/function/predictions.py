"""Implementation of open source tools to predict small RNA functions"""

from os import path as op
import os
import shutil

from bcbio.distributed import transaction
from bcbio.utils import chdir, safe_makedir

from seqcluster.libs import utils, logger as mylog
# import logger as mylog
from seqcluster.libs.read import get_loci_fasta, make_temp_directory
from seqcluster.libs.do import run
import coral

logger = mylog.getLogger(__name__)


def run_coral(clus_obj, out_dir, args):
    """
    Run some CoRaL modules to predict small RNA function
    """
    if not args.bed:
        raise ValueError("This module needs the bed file output from cluster subcmd.")
    workdir = op.abspath(op.join(args.out, 'coral'))
    safe_makedir(workdir)
    bam_in = op.abspath(args.bam)
    bed_in = op.abspath(args.bed)
    reference = op.abspath(args.ref)
    with chdir(workdir):
        bam_clean = coral.prepare_bam(bam_in, bed_in)
        out_dir = op.join(workdir, "regions")
        safe_makedir(out_dir)
        prefix = "seqcluster"
        loci_file = coral.detect_regions(bam_clean, bed_in, out_dir, prefix)
        coral.create_features(bam_clean, loci_file, reference, out_dir)


def is_tRNA(clus_obj, out_dir, args):
    """
    Iterates through cluster precursors to predict sRNA types
    """
    ref = os.path.abspath(args.reference)
    utils.safe_dirs(out_dir)
    for nc in clus_obj[0]:
        c = clus_obj[0][nc]
        loci = c['loci']
        out_fa = "cluster_" + nc
        if loci[0][3] - loci[0][2] < 500:
            with make_temp_directory() as tmpdir:
                os.chdir(tmpdir)
                get_loci_fasta({loci[0][0]: [loci[0][0:5]]}, out_fa, ref)
                summary_file, str_file = _run_tRNA_scan(out_fa)
                if "predictions" not in c:
                    c['predictions'] = {}
                c['predictions']['tRNA'] = _read_tRNA_scan(summary_file)
                score = _read_tRNA_scan(summary_file)
                logger.debug(score)
                shutil.move(summary_file, op.join(out_dir, summary_file))
                shutil.move(str_file, op.join(out_dir, str_file))
        else:
            c['errors'].add("precursor too long")
        clus_obj[0][nc] = c

    return clus_obj


def _read_tRNA_scan(summary_file):
    """
    Parse output from tRNA_Scan
    """
    score = 0
    if os.path.getsize(summary_file) == 0:
        return 0
    with open(summary_file) as in_handle:
        # header = in_handle.next().strip().split()
        for line in in_handle:
            if not line.startswith("--"):
                pre = line.strip().split()
                score = pre[-1]
    return score


def _run_tRNA_scan(fasta_file):
    """
    Run tRNA-scan-SE to predict tRNA
    """
    out_file = fasta_file + "_trnascan"
    se_file = fasta_file + "_second_str"
    cmd = "tRNAscan-SE -q -o {out_file} -f {se_file} {fasta_file}"
    run(cmd.format(**locals()))
    return out_file, se_file
