from argparse import ArgumentParser

import pybedtools

def _lift_positions(line):
    mi_s, mi_e = line[3], line[4]
    snp_p = [e for e, f in enumerate(line) if f.startswith("rs")][0] - 1
    if line[6] == "-":
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

def select_snps(mirna, snp):
    """
    Use bedtools to intersect coordinates
    """
    snp_in_mirna = pybedtools.BedTool(mirna).intersect(pybedtools.BedTool(snp), wo=True)
    for single in snp_in_mirna:
        if single[2] == "miRNA":
            rel_p = _lift_positions(single)
            single[10] = str(rel_p)
            single[12] = _complement(single[12], single[6])
            single[13] = _complement(single[13], single[6])
            single[9] = _get_mirna_name(single[8])
            print single[9:]

if __name__ == "__main__":
    parser = ArgumentParser(description="task related to allele methylation specific")
    parser.add_argument("--gtf", help="gtf file with miRNA annotation", required=1)
    parser.add_argument("--vcf", help="vcf file with SNP variants", required=1)
    args = parser.parse_args()

    mirna_snp = select_snps(args.gtf, args.vcf)
    print "%s mapped to %s: %s" % (args.vcf, args.gtf, mirna_snp)
