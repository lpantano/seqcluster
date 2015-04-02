"""Implementation of open source tools to predict small RNA functions"""
from os import path as op
from bcbio.distributed import transaction
import os
import shutil

import utils
from read import get_loci_fasta
from do import run


def make_predictions(clus_obj, out_dir, args):
    """
    Iterates through cluster precursors to predict sRNA types
    """
    ref = args.reference
    utils.safe_dirs(out_dir)
    for c in clus_obj[0]:
        loci = c['loci']
        out_fa = "cluster_" + c.id
        if loci[0][3] - loci[0][2] > 500:
            with transaction.tx_tmpdir as tmpdir:
                os.chdir(tmpdir)
                get_loci_fasta({loci[0][0]: [loci[0][0:5]]}, out_fa, ref)
                summary_file, str_file = _run_tRNA_scan(out_fa)
                shutil.move(summary_file, op.join(out_dir, summary_file))
                shutil.move(str_file, op.join(out_dir, str_file))


def _run_tRNA_scan(fasta_file):
    """
    Run tRNA-scan-SE to predict tRNA
    """
    out_file = fasta_file + "_trnascan"
    se_file = fasta_file + "_second_str"
    cmd = "tRNAscan-SE {fasta_file} -o {out_file} -f {se_file}"
    run(cmd.format(**locals()))
    return out_file, se_file
