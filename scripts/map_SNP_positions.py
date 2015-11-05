from argparse import ArgumentParser

import pybedtools

def _lift_positions(line):
    print line
    mi_s, mi_e = line[3], line[4]
    snp_p = [e for e, f in enumerate(line) if f.startswith("rs")][0] - 1
    if line[6] == "-":
        rel_p = int(mi_e) - int(line[snp_p]) + 1
        print rel_p

def select_snps(mirna, snp):
    """
    Use bedtools to intersect coordinates
    """
    snp_in_mirna = pybedtools.BedTool(mirna).intersect(pybedtools.BedTool(snp), wo=True)
    for single in snp_in_mirna:
        if single[2] == "miRNA":
            _lift_positions(single)

if __name__ == "__main__":
    parser = ArgumentParser(description="task related to allele methylation specific")
    parser.add_argument("--gtf", help="gtf file with miRNA annotation", required=1)
    parser.add_argument("--vcf", help="vcf file with SNP variants", required=1)
    args = parser.parse_args()

    mirna_snp = select_snps(args.gtf, args.vcf)
    print "%s mapped to %s: %s" % (args.vcf, args.gtf, mirna_snp)
