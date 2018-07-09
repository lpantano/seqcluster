from __future__ import print_function

import sys

import seqcluster.libs.logger as mylog

STDOUT = sys.stdout
logger = mylog.getLogger(__name__)

def _parse_mut(mut):
    """
    Parse mutation field to get position and nts.
    """
    multiplier = 1
    if mut.startswith("-"):
        mut = mut[1:]
        multiplier = -1
    nt = mut.strip('0123456789')
    pos = int(mut[:-2]) * multiplier
    return nt, pos

def _get_reference_position(isomir):
    """
    Liftover from isomir to reference mature
    """
    mut = isomir.split(":")[1]
    if mut == "0":
        return mut
    nt, pos = _parse_mut(mut)
    trim5 = isomir.split(":")[-2]
    off = -1 * len(trim5)
    if trim5.islower():
        off = len(trim5)
    if trim5 == "NA" or trim5 == "0":
        off = 0
    # print(isomir)
    # print([mut, pos, off, nt])
    return "%s%s" % (pos + off, nt)

def _get_pct(isomirs, mirna):
    """
    Get pct of variants respect to the reference
    using reads and different sequences
    """
    pass_pos = []
    for isomir in isomirs.iterrows():
        mir = isomir[1]["chrom"]
        mut = isomir[1]["sv"]
        mut_counts = isomir[1]["counts"]
        total = mirna.loc[mir, "counts"] * 1.0 - mut_counts
        mut_diff = isomir[1]["diff"]
        ratio = mut_counts / total
        if mut_counts > 10 and ratio  > 0.4 and mut != "0" and mut_diff > 1:
            isomir[1]["ratio"] = ratio
            pass_pos.append(isomir[1])
    return pass_pos

def _genotype(data):
    """Simple decision about genotype."""
    if  data['ratio'] > 0.9:
        return "1/1"
    return "1/0"

def _print_header(data):
    """
    Create vcf header to make
    a valid vcf.
    """
    print("##fileformat=VCFv4.2", file=STDOUT, end="")
    print("##source=seqbuster2.3", file=STDOUT, end="")
    print("##reference=mirbase", file=STDOUT, end="")
    for pos in data:
        print("##contig=<ID=%s>" % pos["chrom"], file=STDOUT, end="")
    print('##INFO=<ID=ID,Number=1,Type=String,Description="miRNA name">', file=STDOUT, end="")
    print('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">', file=STDOUT, end="")
    print('##FORMAT=<ID=NR,Number=A,Type=Integer,Description="Total reads supporting the variant">', file=STDOUT, end="")
    print('##FORMAT=<ID=NS,Number=A,Type=Float,Description="Total number of different sequences supporting the variant">', file=STDOUT, end="")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP001", file=STDOUT, end="")

def print_vcf(data):
    """Print vcf line following rules."""
    id_name = "."
    qual = "."
    chrom = data['chrom']
    pos = data['pre_pos']
    nt_ref = data['nt'][1]
    nt_snp = data['nt'][0]
    flt = "PASS"
    info = "ID=%s" % data['mature']
    frmt = "GT:NR:NS"
    gntp = "%s:%s:%s" % (_genotype(data), data["counts"], data["diff"])
    print("\t".join(map(str, [chrom, pos, id_name, nt_ref, nt_snp, qual, flt, info, frmt, gntp])), file=STDOUT, end="")

def _make_header():
    """
    Make vcf header for SNPs in miRs
    """

def liftover(pass_pos, matures):
    """Make position at precursor scale"""
    fixed_pos = []
    _print_header(pass_pos)
    for pos in pass_pos:
        mir = pos["mature"]
        db_pos = matures[pos["chrom"]]
        mut = _parse_mut(pos["sv"])
        print([db_pos[mir], mut, pos["sv"]])
        pos['pre_pos'] = db_pos[mir][0] + mut[1] - 1
        pos['nt'] = list(mut[0])
        fixed_pos.append(pos)
        print_vcf(pos)
    return fixed_pos

def create_vcf(isomirs, matures, gtf, vcf_file=None):
    """
    Create vcf file of changes for all samples.
    PASS will be ones with > 3 isomiRs supporting the position
         and > 30% of reads, otherwise LOW
    """
    global STDOUT
    isomirs['sv'] = [_get_reference_position(m) for m in isomirs["isomir"]]
    mirna = isomirs.groupby(['chrom']).sum()
    sv = isomirs.groupby(['chrom', 'mature', 'sv'], as_index=False).sum()
    sv["diff"] = isomirs.groupby(['chrom', 'mature', 'sv'], as_index=False).size().reset_index().loc[:,0]
    pass_pos = _get_pct(sv, mirna)
    if vcf_file:
        with open(vcf_file, 'w') as out_handle:
            STDOUT = out_handle
            pass_pos = liftover(pass_pos, matures)

    if gtf:
        vcf_genome_file = vcf_file.replace(".vcf", "_genome.vcf")
        with open(vcf_genome_file, 'w') as out_handle:
            STDOUT = out_handle
            pass_pos = liftover_to_genome(pass_pos, gtf)

def liftover_to_genome(pass_pos, gtf):
    """Liftover from precursor to genome"""

    fixed_pos = []
    for pos in pass_pos:
        if pos["chrom"] not in gtf:
            continue
        db_pos = gtf[pos["chrom"]][0]
        mut = _parse_mut(pos["sv"])
        print([db_pos, pos])
        if db_pos[3] == "+":
            pos['pre_pos'] = db_pos[1] + pos["pre_pos"] + 1
        else:
            pos['pre_pos'] = db_pos[2] - (pos["pre_pos"] - 1)
        pos['chrom'] = db_pos[0]
        pos['nt'] = list(mut[0])
        fixed_pos.append(pos)

    _print_header(fixed_pos)
    for pos in fixed_pos:
        print_vcf(pos)
