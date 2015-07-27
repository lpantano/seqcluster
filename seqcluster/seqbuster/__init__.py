# Re-aligner small RNA sequence from SAM/BAM file (miRBase annotation)
import os.path as op
import pandas as pd
import pysam
from bcbio import bam
from bcbio.provenance import do
import seqcluster.libs.logger as mylog
from realign import *

logger = mylog.getLogger(__name__)

def _download_mirbase(args, version="CURRENT"):
    """
    Download files from mirbase
    """
    if not args.hairpin or not args.mirna:
        logger.info("Working with version %s" % version)
        hairpin_fn = op.join(op.abspath(args.out), "hairpin.fa.gz")
        mirna_fn = op.join(op.abspath(args.out), "miRNA.str.gz")
        if not file_exists(hairpin_fn):
            cmd_h = "wget ftp://mirbase.org/pub/mirbase/%s/hairpin.fa.gz -O %s &&  gunzip -f !$" % (version, hairpin_fn)
            do.run(cmd_h, "download hairpin")
        if not file_exists(mirna_fn):
            cmd_m = "wget ftp://mirbase.org/pub/mirbase/%s/miRNA.str.gz -O %s && gunzip -f !$" % (version, mirna_fn)
            do.run(cmd_m, "download mirna")
    else:
        return args.hairpin, args.mirna

def _get_pos(string):
    name = string.split(":")[0][1:]
    pos = string.split(":")[1][:-1].split("-")
    return name, pos

def _read_mature(matures, sps):
    mature = defaultdict(dict)
    with open(matures) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                name = line.strip().replace(">", " ").split()
                mir5p = _get_pos(name[2])
                if len(name) > 3:
                    mir3p = _get_pos(name[3])
                mature[name[0]] = {mir5p[0]: map(int, mir5p[1]),
                                   mir3p[0]: map(int, mir3p[1])}
    return mature

def _read_precursor(precursor, sps):
    """
    read precurso file for that species
    """
    hairpin = {}
    with open(precursor) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                name = line.strip().replace(">", " ").split()[0]
            else:
                hairpin[name] = line.strip()
    return hairpin

def _coord(sequence, start, mirna, precursor, iso):
    """
    Define t5 and t3 isomirs
    """
    dif = abs(mirna[0] - start)
    if start < mirna[0]:
        iso.t5 = "I" + sequence[:dif]
    elif start > mirna[0]:
        iso.t5 = "D" + precursor[mirna[0] - 1:mirna[0] - 1 + dif]

    end = start + (len(sequence) - len(iso.add)) - 1
    logger.debug("%s %s %s " % (start, len(sequence), end))
    dif = abs(mirna[1] - end)
    if dif > 3:
        return None
    if end > mirna[1]:
        iso.t3 = "I" + sequence[-dif:]
    elif end < mirna[1]:
        iso.t3 = "D" + precursor[mirna[1] - dif:mirna[1]]

def _annotate(reads, mirbase_ref, precursors):
    """
    Using SAM/BAM coordinates, mismatches and realign to annotate isomiRs
    """
    for r in reads:
        for p in reads[r].precursors:
            start = reads[r].precursors[p].start + 1  # convert to 1base
            for mature in mirbase_ref[p]:
                mi = mirbase_ref[p][mature]
                if start < mi[0] + 4 and start > mi[0] - 4:
                    logger.debug(("{r} {start} {start} {mi} {mature_s}").format(s=reads[r].sequence, mature_s=precursors[p][mi[0]-1:mi[1]], **locals()))
                    _coord(reads[r].sequence, start, mi, precursors[p], reads[r].precursors[p])
                    reads[r].precursors[p].mirna = mature
                    break
    return reads

def _realign(seq, precursor, start):
    """
    The actual fn that will realign the sequence
    """
    error = set()
    pattern_addition = [[1, 1, 0], [1, 0, 1], [0, 1, 0], [0, 1, 1], [0, 0, 1], [1, 1, 1]]
    for pos in range(0, len(seq)):
        if seq[pos] != precursor[(start + pos)]:
            error.add(pos)

    subs, add = [], []
    for e in error:
        if e < len(seq) - 3:
            subs.append([e, precursor[start + e]])

    pattern, error_add = [], []
    for e in range(len(seq) - 3, len(seq)):
        if e in error:
            pattern.append(1)
            error_add.append(e)
        else:
            pattern.append(0)
    for p in pattern_addition:
        if pattern == p:
            add = seq[error_add[0]:]
            break
    if not add and error_add:
        for e in error_add:
            subs.append([e, precursor[start + e]])

    return subs, add

def _clean_hits(reads):
    """
    Select only best matches
    """
    new_reads = defaultdict(realign)
    for r in reads:
        world = {}
        sc = 0
        for p in reads[r].precursors:
            world[p] = reads[r].precursors[p].get_score(len(reads[r].sequence))
            if sc < world[p]:
                sc = world[p]
        new_reads[r] = reads[r]
        for p in world:
            if sc != world[p]:
                new_reads[r].remove_precursor(p)

    return new_reads

def _sort_by_name(bam_fn):
    """
    sort bam file by name sequence
    """

def _read_bam(bam_fn, precursors):
    """
    read bam file and perform realignment of hits
    """
    handle = bam.open_samfile(bam_fn)
    reads = defaultdict(realign)
    for line in handle:
        chrom = handle.getrname(line.reference_id)
        # print "%s %s %s %s" % (line.query_name, line.reference_start, line.query_sequence, chrom)
        if line.query_name not in reads:
            reads[line.query_name].sequence = line.query_sequence
        iso = isomir()
        iso.align = line
        iso.start = line.reference_start
        iso.subs, iso.add = _realign(reads[line.query_name].sequence, precursors[chrom], line.reference_start)
        reads[line.query_name].set_precursor(chrom, iso)

    reads = _clean_hits(reads)
    return reads

def _get_freq(name):
    """
    Check if name read contains counts (_xNumber)
    """
    try:
        counts = name.split("_x")[1]
    except:
        return "NA"
    return counts

def _tab_output(reads, out_file, sample):
    seen = set()
    lines = []
    with open(out_file, 'w') as out_handle:
        print >>out_handle, "name\tfreq\tchrom\tsubs\tadd\tt5\tt3"
        for r, read in reads.iteritems():
            hits = set()
            [hits.add(mature.mirna) for mature in read.precursors.values() if mature.mirna]
            hits = len(hits)
            for p, iso in read.precursors.iteritems():
                if (r, iso.mirna) not in seen:
                    seen.add((r, iso.mirna))
                    chrom = iso.mirna
                    if not chrom:
                        chrom = p
                    count = _get_freq(r)
                    res = ("{r}\t{count}\t{chrom}\t{format}\t{hits}")
                    annotation = "%s:%s" % (chrom, iso.format(":"))
                    lines.append([annotation, chrom, count, sample, hits])
                    print >>out_handle, res.format(format=iso.format(), **locals())
    dt = pd.DataFrame(lines)
    dt.columns = ["isomir", "chrom", "counts", "sample", "hits"]
    return out_file, dt

def _merge(dts):
    """
    merge multiple samples in one matrix
    """
    df = None
    for dt in dts:
        if not df:
            df = dt
        else:
            df.join(dt)

    ma = df.pivot(index='isomir', columns='sample', values='counts')
    ma_mirna = ma
    ma_mirna['mirna'] = [m.split(":")[0] for m in ma.index.values]
    ma_mirna = ma_mirna.groupby(['mirna']).sum()

    return ma

def miraligner(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    config = {"algorithm": {"num_cores": 1}}
    hairpin, mirna = _download_mirbase(args)
    precursors = _read_precursor(args.hairpin, args.sps)
    matures = _read_mature(args.mirna, args.sps)
    out_dts = []
    for bam_fn in args.files:
        logger.info("Reading %s" % bam_fn)
        sample = op.splitext(op.basename(bam_fn))[0]
        out_file = op.join(args.out, sample + ".mirna")
        bam_fn = bam.sam_to_bam(bam_fn, config)
        bam_sort_by_n = op.splitext(bam_fn)[0] + "_sort"
        pysam.sort("-n", bam_fn, bam_sort_by_n)
        reads = _read_bam(bam_sort_by_n + ".bam", precursors)
        _annotate(reads, matures, precursors)
        out_file, dt = _tab_output(reads, out_file, sample)
        out_dts.append(dt)

    ma = _merge(out_dts)
    out_ma = op.join(args.out, "counts.tsv")
    ma.to_csv(out_ma, sep="\t")
    # _summarize(out_dts)
