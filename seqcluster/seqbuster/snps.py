
def _get_reference_position(isomir):
    mut = isomir.split(":")[1]
    if mut == "0":
        return mut
    nt = mut.strip('0123456789')
    pos = int(mut[:len(nt)])
    trim5 = isomir.split(":")[-2]
    off = -1 * len(trim5)
    if trim5.isupper():
        off = len(trim5)
    return "%s%s" % (pos + off, nt)

def _get_pct(isomirs, mirna):
    for isomir in isomirs.iterrows():
        mir = isomir[1]["chrom"]
        total = mirna.loc[mir, "counts"] * 1.0
        mut_counts = isomir[1]["counts"]
        if mut_counts > 10 and mut_counts / total > 0.4:
            print isomir

def create_vcf(isomirs):
    """
    Create vcf file of changes for all samples.
    PASS will be ones with > 3 isomiRs supporting the position
         and > 30% of reads, otherwise LOW
    """
    isomirs['sv'] = [_get_reference_position(m) for m in isomirs["isomir"]]
    mirna = isomirs.groupby(['chrom']).sum()
    sv = isomirs.groupby(['chrom', 'sv'], as_index=False).sum()
    sv["diff"] = isomirs.groupby(['chrom', 'sv'], as_index=False).size().reset_index().loc[:,0]
    print isomirs
    print sv
    _get_pct(sv, mirna)
    # fn to get counts supporting SNPs from the total
    # fn to get sequences supporting the SNPs
    # annotate with bedtools? or maybe better another pypackage

