"""
Wrap RNAfold command
"""
import os
import subprocess
import pybedtools


def run_rnafold(seqs):
    out = 0
    cmd = ("echo {seqs} | RNAfold").format(**locals())
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    for line in iter(process.stdout.readline, ''):
        if line.find(" ") > -1:
            out = float(line.split(" ")[1][1:7])
    # print out
    if not out:
        return ""
    return out

def calculate_structure(loci_file):
    structure_file = os.path.splitext(loci_file)[0] + "-fold.tsv"
    loci_bed = pybedtools.BedTool(loci_file)
    # getfasta with reference fasta
    # for line: get fasta -> run_rnafold -> add value to out_file
    return structure_file
