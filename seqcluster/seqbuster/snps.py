def _parse_mut(mut):
    multiplier = 1
    if mut.startswith("-"):
        mut = mut[1:]
        multiplier = -1
    nt = mut.strip('0123456789')
    pos = int(mut[:len(nt)-1]) * multiplier
    return nt, pos

def _get_reference_position(isomir):
    mut = isomir.split(":")[1]
    if mut == "0":
        return mut
    nt, pos = _parse_mut(mut)
    trim5 = isomir.split(":")[-2]
    off = -1 * len(trim5)
    if trim5.isupper():
        off = len(trim5)
    return "%s%s" % (pos + off, nt)

def _get_pct(isomirs, mirna):
    pass_pos = []
    for isomir in isomirs.iterrows():
        mir = isomir[1]["chrom"]
        mut = isomir[1]["sv"]
        total = mirna.loc[mir, "counts"] * 1.0
        mut_counts = isomir[1]["counts"]
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

def print_vcf(data):
    """Print vcf line following rules."""
    chrom = data['chrom']
    pos = data['pre_pos']
    nt_ref = data['nt'][1]
    nt_snp = data['nt'][0]
    flt = "PASS"
    info = "id=%s" % data['mature']
    frmt = "%s:%s:%s" % (data["counts"], data["diff"], _genotype(data))
    print "\t".join(map(str, [chrom, pos, nt_ref, nt_snp, flt, info, frmt]))

def liftover(pass_pos, matures):
    """Make position at precursor"""
    fixed_pos = []
    for pos in pass_pos:
        mir = pos["mature"]
        db_pos = matures[pos["chrom"]]
        mut = _parse_mut(pos["sv"])
        pos['pre_pos'] =  db_pos[mir][0] + mut[1]
        pos['nt'] = list(mut[0])
        fixed_pos.append(pos)
        print_vcf(pos)
    return fixed_pos

def create_vcf(isomirs, matures):
    """
    Create vcf file of changes for all samples.
    PASS will be ones with > 3 isomiRs supporting the position
         and > 30% of reads, otherwise LOW
    """
    isomirs['sv'] = [_get_reference_position(m) for m in isomirs["isomir"]]
    mirna = isomirs.groupby(['chrom']).sum()
    sv = isomirs.groupby(['chrom', 'mature', 'sv'], as_index=False).sum()
    sv["diff"] = isomirs.groupby(['chrom', 'mature', 'sv'], as_index=False).size().reset_index().loc[:,0]
    # print isomirs
    # print sv
    pass_pos = _get_pct(sv, mirna)
    pass_pos = liftover(pass_pos, matures)
    # annotate with bedtools? or maybe better another pypackage
