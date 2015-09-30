"""functions for explore tool"""
from __future__ import print_function
from collections import defaultdict
import json, itertools
import tempfile, os, contextlib, shutil

from sam2bed import makeBED
import logging
from do import find_cmd, run

import pysam
import pybedtools
try:
    from seqcluster.align import pyMatch
except:
    pass

logger = logging.getLogger('read')


@contextlib.contextmanager
def make_temp_directory(remove=True):
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)

def load_data(in_file):
    """load json file from seqcluster cluster"""
    with open(in_file) as in_handle:
        return json.load(in_handle)

def write_data(data, out_file):
    """write json file from seqcluster cluster"""
    with open(out_file, 'w') as handle_out:
        handle_out.write(json.dumps([data], skipkeys=True, indent=2))

def get_sequences_from_cluster(c1, c2, data):
    """get all sequences from on cluster"""
    seqs1 = data[c1]['seqs']
    seqs2 = data[c2]['seqs']
    seqs = list(set(seqs1 + seqs2))
    names = []
    for s in seqs:
        if s in seqs1 and s in seqs2:
            names.append("both")
        elif s in seqs1:
            names.append(c1)
        else:
            names.append(c2)
    return seqs, names


def get_precursors_from_cluster(c1, c2, data):
    loci1 = data[c1]['loci']
    loci2 = data[c2]['loci']
    return {c1:loci1, c2:loci2}


def map_to_precursors(seqs, names, loci, out_file, args):
    """map sequences to precursors with razers3"""
    with make_temp_directory() as temp:
        pre_fasta = os.path.join(temp, "pre.fa")
        seqs_fasta = os.path.join(temp, "seqs.fa")
        out_sam = os.path.join(temp, "out.sam")
        pre_fasta = get_loci_fasta(loci, pre_fasta, args.ref)
        out_precursor_file = out_file.replace("tsv", "fa")
        seqs_fasta = get_seqs_fasta(seqs, names, seqs_fasta)

        # print(open(pre_fasta).read().split("\n")[1])
        if find_cmd("razers3"):
            cmd = "razers3 -dr 2 -i 80 -rr 90 -f -o {out_sam} {temp}/pre.fa  {seqs_fasta}"
            run(cmd.format(**locals()))
            out_file = read_alignment(out_sam, loci, seqs, out_file)
            shutil.copy(pre_fasta, out_precursor_file)
    return out_file


def map_to_precursors_on_fly(seqs, names, loci, args):
    """map sequences to precursors with franpr algorithm to avoid writting on disk"""
    region = "%s\t%s\t%s\t.\t.\t%s" % (loci[1], loci[2], loci[3], loci[4])
    precursor = pybedtools.BedTool(str(region), from_string=True).sequence(fi=args.ref, s=True)
    precursor = open(precursor.seqfn).read().split("\n")[1]
    dat = dict()
    for s, n in itertools.izip(seqs, names):
        res = pyMatch.Match(precursor, str(s), 1, 3)
        if res > -1:
            dat[n] = [res, res + len(s)]
    return dat


def deprecated_map_to_precursors(seqs, names, loci, out_file, args):
    """map sequences to precursors with bowtie"""
    with make_temp_directory() as temp:
        pre_fasta = os.path.join(temp, "pre.fa")
        seqs_fasta = os.path.join(temp, "seqs.fa")
        out_sam = os.path.join(temp, "out.sam")
        pre_fasta = get_loci_fasta(loci, pre_fasta, args.ref)
        out_precursor_file = out_file.replace("tsv", "fa")
        seqs_fasta = get_seqs_fasta(seqs, names, seqs_fasta)
        if find_cmd("bowtie2-build"):
            cmd = "bowtie2-build -f {pre_fasta} {temp}/pre"
            run(cmd.format(**locals()))
            cmd = "bowtie2 -a --rdg 7,3 --mp 4 --end-to-end -D 20 -R 3 -N 0 -i S,1,0.8 -L 3 -f -x  {temp}/pre -U {seqs_fasta} -S {out_sam}"
            run(cmd.format(**locals()))
            out_file = read_alignment(out_sam, loci, seqs, out_file)
            shutil.copy(pre_fasta, out_precursor_file)
    return out_file


def get_seqs_fasta(seqs, names, out_fa):
    """get fasta from sequences"""
    with open(out_fa, 'w') as fa_handle:
        for s, n in itertools.izip(seqs, names):
            print(">cx{1}-{0}\n{0}".format(s, n), file=fa_handle)
    return out_fa


def get_fasta(bed_file, ref, out_fa):
    """Run bedtools to get fasta from bed file"""
    cmd = "bedtools getfasta -s -fi {ref} -bed {bed_file} -fo {out_fa}"
    run(cmd.format(**locals()))


def get_loci_fasta(loci, out_fa, ref):
    """get fasta from precursor"""
    if not find_cmd("bedtools"):
        raise ValueError("Not bedtools installed")
    with make_temp_directory() as temp:
        bed_file = os.path.join(temp, "file.bed")
        for nc, loci in loci.iteritems():
            for l in loci:
                with open(bed_file, 'w') as bed_handle:
                    logger.debug("get_fasta: loci %s" % l)
                    nc, c, s, e, st = l
                    print("{0}\t{1}\t{2}\t{3}\t{3}\t{4}".format(c, s, e, nc, st), file=bed_handle)
                get_fasta(bed_file, ref, out_fa)
    return out_fa


def read_alignment(out_sam, loci, seqs, out_file):
    """read which seqs map to which loci and
    return a tab separated file"""
    hits = defaultdict(list)
    with open(out_file, "w") as out_handle:
        samfile = pysam.Samfile(out_sam, "r")
        for a in samfile.fetch():
            if not a.is_unmapped:
                nm = int([t[1] for t in a.tags if t[0] == "NM"][0])
                a = makeBED(a)
                if not a:
                    continue
                ref, locus = get_loci(samfile.getrname(int(a.chr)), loci)
                hits[a.name].append((nm, "%s %s %s %s %s %s" % (a.name, a.name.split("-")[0], locus, ref, a.start, a.end)))
        for hit in hits.values():
            nm = hit[0][0]
            for l in hit:
                if nm == l[0]:
                    print(l[1], file=out_handle)
    return out_file


def get_loci(name, loci):
    for nc in loci:
        for l in loci[nc]:
            lname = "{0}:{1}-{2}({3})".format(l[1], l[2], l[3], l[4])
            if name == lname:
                return nc, l[0]


def plot_positions():
	return True
