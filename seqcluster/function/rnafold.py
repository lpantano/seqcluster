"""
Wrap RNAfold command
"""

import subprocess


def run_rnafold(seqs):
    out = 0
    cmd = ("echo {seqs} | RNAfold").format(**locals())
    print cmd
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    for line in iter(process.stdout.readline, ''):
        if line.find(" ") > -1:
            out = float(line.split(" ")[1][1:7])
    return out
