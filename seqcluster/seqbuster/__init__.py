# Re-aligner small RNA sequence from SAM/BAM file (miRBase annotation)

from __future__ import print_function
import os.path as op
import re
import shutil
import pandas as pd
import pysam
import argparse

from seqcluster.libs import do
from seqcluster.libs.utils import file_exists
import seqcluster.libs.logger as mylog
from seqcluster.install import _get_miraligner
from seqcluster.seqbuster.snps import create_vcf
from seqcluster.collapse import collapse_fastq
from seqcluster.seqbuster.realign import *
from mirtop.gff import reader

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


def _make_unique(name, idx):
    """Make name unique in case only counts there"""
    p = re.compile(".[aA-zZ]+_x[0-9]+")
    if p.match(name):
        tags = name[1:].split("_x")
        return ">%s_%s_x%s" % (tags[0], idx, tags[1])
    return name.replace("@", ">")


def _filter_seqs(fn):
    """Convert names of sequences to unique ids"""
    out_file = op.splitext(fn)[0] + "_unique.fa"
    idx = 0
    if not file_exists(out_file):
        with open(out_file, 'w') as out_handle:
            with open(fn) as in_handle:
                line = in_handle.readline()
                while line:
                    if line.startswith("@") or line.startswith(">"):
                        fixed_name = _make_unique(line.strip(), idx)
                        seq = in_handle.readline().strip()
                        counts = _get_freq(fixed_name)
                        if len(seq) < 26 and (counts > 1 or counts == 0):
                            idx += 1
                            print(fixed_name, file=out_handle, end="\n")
                            print(seq, file=out_handle, end="\n")
                        if line.startswith("@"):
                            in_handle.readline()
                            in_handle.readline()
                    line = in_handle.readline()
    return out_file


def _convert_to_fasta(fn):
    out_file = op.splitext(fn)[0] + ".fa"
    with open(out_file, 'w') as out_handle:
        with open(fn) as in_handle:
            line = in_handle.readline()
            while line:
                if line.startswith("@"):
                    seq = in_handle.readline()
                    _ = in_handle.readline()
                    qual = in_handle.readline()
                elif line.startswith(">"):
                    seq = in_handle.readline()
                count = 2
                if line.find("_x"):
                    count = int(line.strip().split("_x")[1])
                if count > 1:
                    print(">%s" % line.strip()[1:], file=out_handle, end="")
                    print(seq.strip(), file=out_handle, end="")
                line = in_handle.readline()
    return out_file


def _get_pos(string):
    name = string.split(":")[0][1:]
    pos = string.split(":")[1][:-1].split("-")
    return name, map(int, pos)


def _read_mature(matures, sps):
    mature = defaultdict(dict)
    with open(matures) as in_handle:
        for line in in_handle:
            if line.startswith(">") and line.find(sps) > -1:
                name = line.strip().replace(">", " ").split()
                mir5p = _get_pos(name[2])
                mature[name[0]] = {mir5p[0]: mir5p[1]}
                if len(name) > 3:
                    mir3p = _get_pos(name[3])
                    mature[name[0]].update({mir3p[0]: mir3p[1]})
    return mature


def _read_precursor(precursor, sps):
    """
    Load precursor file for that species
    """
    hairpin = defaultdict(str)
    name = None
    with open(precursor) as in_handle:
        for line in in_handle:
            if line.startswith(">"):
                if hairpin[name]:
                    hairpin[name] = hairpin[name] + "NNNNNNNNNNNN"
                name = line.strip().replace(">", " ").split()[0]
            else:
                hairpin[name] += line.strip()
        hairpin[name] = hairpin[name] + "NNNNNNNNNNNN"
    return hairpin


def _read_gtf(gtf):
    """
    Load GTF file with precursor positions on genome
    """
    if not gtf:
        return gtf
    db = defaultdict(list)
    with open(gtf) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            name = [n.split("=")[1] for n in cols[-1].split(";") if n.startswith("Name")]
            chrom, start, end, strand = cols[0], cols[3], cols[4], cols[6]
            if cols[2] == "miRNA_primary_transcript":
                db[name[0]].append([chrom, int(start), int(end), strand])
    return db


def _coord(sequence, start, mirna, precursor, iso):
    """
    Define t5 and t3 isomirs
    """
    dif = abs(mirna[0] - start)
    if start < mirna[0]:
        iso.t5 = sequence[:dif].upper()
    elif start > mirna[0]:
        iso.t5 = precursor[mirna[0] - 1:mirna[0] - 1 + dif].lower()
    elif start == mirna[0]:
        iso.t5 = "NA"
    if dif > 4:
        logger.debug("start > 3 %s %s %s %s %s" % (start, len(sequence), dif, mirna, iso.format()))
        return None

    end = start + (len(sequence) - len(iso.add)) - 1
    dif = abs(mirna[1] - end)
    if iso.add:
        sequence = sequence[:-len(iso.add)]
    # if dif > 3:
    #    return None
    if end > mirna[1]:
        iso.t3 = sequence[-dif:].upper()
    elif end < mirna[1]:
        iso.t3 = precursor[mirna[1] - dif:mirna[1]].lower()
    elif end == mirna[1]:
        iso.t3 = "NA"
    if dif > 4:
        logger.debug("end > 3 %s %s %s %s %s" % (len(sequence), end, dif, mirna, iso.format()))
        return None
    logger.debug("%s %s %s %s %s %s" % (start, len(sequence), end, dif, mirna, iso.format()))
    return True


def _annotate(reads, mirbase_ref, precursors):
    """
    Using SAM/BAM coordinates, mismatches and realign to annotate isomiRs
    """
    for r in reads:
        for p in reads[r].precursors:
            start = reads[r].precursors[p].start + 1  # convert to 1base
            end = start + len(reads[r].sequence)
            for mature in mirbase_ref[p]:
                mi = mirbase_ref[p][mature]
                is_iso = _coord(reads[r].sequence, start, mi, precursors[p], reads[r].precursors[p])
                logger.debug(("{r} {p} {start} {is_iso} {mature} {mi} {mature_s}").format(s=reads[r].sequence, mature_s=precursors[p][mi[0]-1:mi[1]], **locals()))
                if is_iso:
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
            subs.append([e, seq[e], precursor[start + e]])

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
            subs.append([e, seq[e], precursor[start + e]])

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
            logger.debug("score %s %s %s" % (r, p, world[p]))
            if sc != world[p]:
                logger.debug("remove %s %s %s" % (r, p, world[p]))
                new_reads[r].remove_precursor(p)

    return new_reads


def _sort_by_name(bam_fn):
    """
    sort bam file by name sequence
    """


def _sam_to_bam(bam_fn):
    if bam_fn.endswith("bam"):
        bam_out = "%s.bam" % os.path.splitext(bam_fn)[0]
        cmd = "samtools view -Sbh {bam_fn} -o {bam_out}"
        do.run(cmd)
        return bam_out
    return bam_fn


def _read_bam(bam_fn, precursors):
    """
    read bam file and perform realignment of hits
    """
    mode = "r" if bam_fn.endswith("sam") else "rb"
    handle = pysam.Samfile(bam_fn, mode)
    reads = defaultdict(realign)
    for line in handle:
        chrom = handle.getrname(line.reference_id)
        # print("%s %s %s %s" % (line.query_name, line.reference_start, line.query_sequence, chrom))
        query_name = line.query_name
        if query_name not in reads:
            reads[query_name].sequence = line.query_sequence
        iso = isomir()
        iso.align = line
        iso.start = line.reference_start
        iso.subs, iso.add = _realign(reads[query_name].sequence, precursors[chrom], line.reference_start)
        reads[query_name].set_precursor(chrom, iso)

    reads = _clean_hits(reads)
    return reads


def _collapse_fastq(in_fn):
    """
    collapse reads into unique sequences
    """
    args = argparse.Namespace()
    args.fastq = in_fn
    args.minimum = 1
    args.out = op.dirname(in_fn)
    return collapse_fastq(args)


def _read_pyMatch(fn, precursors):
    """
    read pyMatch file and perform realignment of hits
    """
    with open(fn) as handle:
        reads = defaultdict(realign)
        for line in handle:
            query_name, seq, chrom, reference_start, end, mism, add = line.split()
            reference_start = int(reference_start)
            # chrom = handle.getrname(cols[1])
            # print("%s %s %s %s" % (line.query_name, line.reference_start, line.query_sequence, chrom))
            if query_name not in reads:
                reads[query_name].sequence = seq
            iso = isomir()
            iso.align = line
            iso.start = reference_start
            iso.subs, iso.add = _realign(reads[query_name].sequence, precursors[chrom], reference_start)
            logger.debug("%s %s %s %s %s" % (query_name, reference_start, chrom, iso.subs, iso.add))
            if len(iso.subs) > 1:
                continue
            reads[query_name].set_precursor(chrom, iso)

        reads = _clean_hits(reads)
    return reads


def _parse_mut(subs):
    """
    Parse mutation tag from miraligner output
    """
    if subs!="0":
        subs = [[subs.replace(subs[-2:], ""),subs[-2], subs[-1]]]
    return subs


def _read_miraligner(fn):
    """Read ouput of miraligner and create compatible output."""
    reads = defaultdict(realign)
    with open(fn) as in_handle:
        in_handle.readline()
        for line in in_handle:
            cols = line.strip().split("\t")
            iso = isomir()
            query_name, seq = cols[1], cols[0]
            chrom, reference_start = cols[-2], cols[3]
            iso.mirna = cols[3]
            subs, add, iso.t5, iso.t3 = cols[6:10]
            if query_name not in reads:
                reads[query_name].sequence = seq
            iso.align = line
            iso.start = reference_start
            iso.subs, iso.add = _parse_mut(subs), add
            logger.debug("%s %s %s %s %s" % (query_name, reference_start, chrom, iso.subs, iso.add))
            reads[query_name].set_precursor(chrom, iso)
    return reads


def _cmd_miraligner(fn, out_file, species, hairpin, out):
    """
    Run miraligner for miRNA annotation
    """
    tool = _get_miraligner()
    path_db = op.dirname(op.abspath(hairpin))
    cmd = "{tool} -freq -i {fn} -o {out_file} -s {species} -db {path_db} -sub 1 -trim 3 -add 3"
    if not file_exists(out_file):
        logger.info("Running miraligner with %s" % fn)
        do.run(cmd.format(**locals()), "miraligner with %s" % fn)
        shutil.move(out_file + ".mirna", out_file)
    return out_file


def _mirtop(out_files, hairpin, gff3, species, out):
    """
    Convert miraligner to mirtop format
    """
    args = argparse.Namespace()
    args.hairpin = hairpin
    args.sps = species
    args.gtf = gff3
    args.add_extra = True
    args.files = out_files
    args.format = "seqbuster"
    args.out_format = "gff"
    args.out = out
    reader(args)


def _get_freq(name):
    """
    Check if name read contains counts (_xNumber)
    """
    try:
        counts = int(name.split("_x")[1])
    except:
        return 0
    return counts


def _tab_output(reads, out_file, sample):
    seen = set()
    lines = []
    lines_pre = []
    seen_ann = {}
    dt = None
    with open(out_file, 'w') as out_handle:
        print("name\tseq\tfreq\tchrom\tstart\tend\tsubs\tadd\tt5\tt3\ts5\ts3\tDB\tprecursor\thits", file=out_handle, end="")
        for (r, read) in reads.items():
            hits = set()
            [hits.add(mature.mirna) for mature in read.precursors.values() if mature.mirna]
            hits = len(hits)
            for (p, iso) in read.precursors.items():
                if len(iso.subs) > 3 or not iso.mirna:
                    continue
                if (r, iso.mirna) not in seen:
                    seen.add((r, iso.mirna))
                    chrom = iso.mirna
                    if not chrom:
                        chrom = p
                    count = _get_freq(r)
                    seq = reads[r].sequence
                    if iso.get_score(len(seq)) < 1:
                        continue
                    if iso.subs:
                        iso.subs = [] if "N" in iso.subs[0] else iso.subs
                    annotation = "%s:%s" % (chrom, iso.format(":"))
                    res = ("{seq}\t{r}\t{count}\t{chrom}\tNA\tNA\t{format}\tNA\tNA\tmiRNA\t{p}\t{hits}").format(format=iso.format().replace("NA", "0"), **locals())
                    if annotation in seen_ann and seq.find("N") < 0 and seen_ann[annotation].split("\t")[0].find("N") < 0:
                        raise ValueError("Same isomir %s from different sequence: \n%s and \n%s" % (annotation, res, seen_ann[annotation]))
                    seen_ann[annotation] = res
                    lines.append([annotation, chrom, count, sample, hits])
                    lines_pre.append([annotation, chrom, p, count, sample, hits])
                    print(res, file=out_handle, end="")

    if lines:
        dt = pd.DataFrame(lines)
        dt.columns = ["isomir", "chrom", "counts", "sample", "hits"]
        dt = dt[dt['hits']>0]
        dt = dt.loc[:, "isomir":"sample"]
        dt = dt.groupby(['isomir', 'chrom', 'sample'], as_index=False).sum()
        dt.to_csv(out_file + "_summary")
        dt_pre = pd.DataFrame(lines_pre)
        dt_pre.columns = ["isomir", "mature", "chrom", "counts", "sample", "hits"]
        dt_pre = dt_pre[dt_pre['hits']==1]
        dt_pre = dt_pre.loc[:, "isomir":"sample"]
        dt_pre = dt_pre.groupby(['isomir', 'chrom', 'mature', 'sample'], as_index=False).sum()
        return out_file, dt, dt_pre
    return None


def _merge(dts):
    """
    merge multiple samples in one matrix
    """
    df = pd.concat(dts)

    ma = df.pivot(index='isomir', columns='sample', values='counts')
    ma_mirna = ma
    ma = ma.fillna(0)
    ma_mirna['mirna'] = [m.split(":")[0] for m in ma.index.values]
    ma_mirna = ma_mirna.groupby(['mirna']).sum()
    ma_mirna = ma_mirna.fillna(0)
    return ma, ma_mirna


def _create_counts(out_dts, out_dir):
    """Summarize results into single files."""
    ma, ma_mirna = _merge(out_dts)
    out_ma = op.join(out_dir, "counts.tsv")
    out_ma_mirna = op.join(out_dir, "counts_mirna.tsv")
    ma.to_csv(out_ma, sep="\t")
    ma_mirna.to_csv(out_ma_mirna, sep="\t")
    return out_ma_mirna, out_ma


def miraligner(args):
    """
    Realign BAM hits to miRBAse to get better accuracy and annotation
    """
    hairpin, mirna = _download_mirbase(args)
    precursors = _read_precursor(args.hairpin, args.sps)
    matures = _read_mature(args.mirna, args.sps)
    gtf = _read_gtf(args.gtf)
    out_dts = []
    out_files = []
    for bam_fn in args.files:
        sample = op.splitext(op.basename(bam_fn))[0]
        logger.info("Reading %s" % bam_fn)
        if bam_fn.endswith("bam") or bam_fn.endswith("sam"):
            bam_fn = _sam_to_bam(bam_fn)
            bam_sort_by_n = op.splitext(bam_fn)[0] + "_sort"
            pysam.sort("-n", bam_fn, bam_sort_by_n)
            reads = _read_bam(bam_sort_by_n + ".bam", precursors)
        elif bam_fn.endswith("fasta") or bam_fn.endswith("fa") or \
                bam_fn.endswith("fastq"):
            if args.collapse:
                bam_fn = _collapse_fastq(bam_fn)
            out_file = op.join(args.out, sample + ".premirna")
            bam_fn = _filter_seqs(bam_fn)
            if args.miraligner:
                _cmd_miraligner(bam_fn, out_file, args.sps, args.hairpin, args.out)
                reads = _read_miraligner(out_file)
                out_files.append(out_file)
        else:
            raise ValueError("Format not recognized.")

        if args.miraligner:
            _mirtop(out_files, args.hairpin, args.gtf, args.sps, args.out)

        if not args.miraligner:
            reads = _annotate(reads, matures, precursors)

        out_file = op.join(args.out, sample + ".mirna")
        out_file, dt, dt_pre = _tab_output(reads, out_file, sample)
        try:
            vcf_file = op.join(args.out, sample + ".vcf")
            if not file_exists(vcf_file):
                # if True:
                create_vcf(dt_pre, matures, gtf, vcf_file)
            try:
                import vcf
                vcf.Reader(filename=vcf_file)
            except Exception as e:
                logger.warning(e.__doc__)
                logger.warning(e)
        except Exception as e:
            # traceback.print_exc()
            logger.warning(e.__doc__)
            logger.warning(e)
        if isinstance(dt, pd.DataFrame):
            out_dts.append(dt)

    if out_dts:
        _create_counts(out_dts, args.out)
    else:
        print("No files analyzed!")
