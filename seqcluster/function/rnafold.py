"""
Wrap RNAfold command
"""

import subprocess
import pybedtools

from bcbio.utils import splitext_plus


def run_rnafold(seqs):
    out = 0
    cmd = ("echo {seqs} | RNAfold").format(**locals())
    print cmd
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    for line in iter(process.stdout.readline, ''):
        if line.find(" ") > -1:
            out = float(line.split(" ")[1][1:7])
    return out

def calculate_structure(loci_file):
    structure_file = splitext_plus(loci_file, "-fold.tsv")
    loci_bed = pybedtools.BedTool(loci_file)
    # getfasta with reference fasta
    # for line: get fasta -> run_rnafold -> add 
