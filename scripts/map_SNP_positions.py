from __future__ import print_function
from argparse import ArgumentParser
import os
import gzip

import pybedtools

def _open_file(in_file):
    """From bcbio code"""
    _, ext = os.path.splitext(in_file)
    if ext == ".gz":
        return gzip.open(in_file, 'rb')
    if ext in [".fastq", ".fq"]:
        return open(in_file, 'r')
    # default to just opening it
    return open(in_file, "r")

def _lift_positions(line):
    mi_s, mi_e = line[11], line[12]
    snp_p = [e for e, f in enumerate(line) if f.startswith("rs")][0] - 1
    if line[14] == "-":
        rel_p = int(mi_e) - int(line[snp_p]) + 1
    else:
        rel_p = int(line[snp_p]) - int(mi_s) + 1
    return rel_p

def _complement(nt, strand):
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return nt if strand == "+" else pairs[nt]

def _get_mirna_name(attrs):
    items = attrs.split(";")
    return [i.split("=")[1] for i in items if i.startswith("Name")][0]

def _create_header(mirna, snp, out):
    names = set()
    h = []
    last_h = ""
    with _open_file(snp) as in_handle:
        for line in in_handle:
            if line.startswith("##") and not line.startswith("##contig"):
                h.append(line.strip())
            elif line.startswith("#CHROM"):
                last_h = line.strip()
            else:
                break
    with open(mirna) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                cols = line.strip().split("\t")
                if cols[2] == "miRNA":
                    names.add(_get_mirna_name(cols[8]))
        for name in names:
            h.append("##contig=<ID=%s>" % name)
    h.append(last_h)
    return "\n".join(h)

def select_snps(mirna, snp, out):
    """
    Use bedtools to intersect coordinates
    """
    with open(out, 'w') as out_handle:
        print(_create_header(mirna, snp, out), file=out_handle, end="")
        snp_in_mirna = pybedtools.BedTool(snp).intersect(pybedtools.BedTool(mirna), wo=True)
        for single in snp_in_mirna:
            if single[10] == "miRNA" and len(single[3]) + len(single[4]) == 2:
                line = []
                rel_p = _lift_positions(single)
                line.append(_get_mirna_name(single[16]))
                line.append(str(rel_p))
                line.append(single[2])
                line.append(_complement(single[3], single[14]))
                line.append(_complement(single[4], single[14]))
                line.append(single[5])
                line.append(single[6])
                line.append(single[7])
                print("\t".join(line), file=out_handle, end="")
    return out

if __name__ == "__main__":
    parser = ArgumentParser(description="task related to allele methylation specific")
    parser.add_argument("--gtf", help="gtf file with miRNA annotation", required=1)
    parser.add_argument("--vcf", help="vcf file with SNP variants", required=1)
    parser.add_argument("--out", help="out put file.", required=1)
    args = parser.parse_args()

    mirna_snp = select_snps(args.gtf, args.vcf, args.out)
    print("%s mapped to %s: %s" % (args.vcf, args.gtf, mirna_snp))
