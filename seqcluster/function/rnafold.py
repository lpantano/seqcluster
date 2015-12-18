"""
Wrap RNAfold command
"""
import os
import subprocess
import pybedtools

from seqcluster.libs.utils import safe_run
from seqcluster.libs import do

def run_rnafold(seqs):
    out = 0
    cmd = ("echo {seqs} | RNAfold").format(**locals())
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    assert do.find_cmd("RNAfold"), "Install RNAfold first."
    for line in iter(process.stdout.readline, ''):
        if line.find(" ") > -1:
            out = float(line.split(" (")[1][0:6].strip())
    return out

def calculate_structure(loci_file, genome):
    """
    Get fasta sequence for each loci and calculate structure
    """
    structure_file = os.path.splitext(loci_file)[0] + "-fold.tsv"
    with safe_run(structure_file) as out_tx:
        with open(out_tx, 'w') as out_handle:
            loci_bed = pybedtools.BedTool(loci_file)
            print >>out_handle, "id\tfold"
            for region in loci_bed:
                seq = pybedtools.BedTool(str(region), from_string=True).sequence(fi=genome)
                seq = open(seq.seqfn).read().split("\n")[1]
                g = run_rnafold(seq)
                print >>out_handle, "%s\t%s" % (region[3], g)
    return structure_file
