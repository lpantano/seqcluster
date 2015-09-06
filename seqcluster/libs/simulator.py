"""simulate cluster over the genome"""
import random
from read import get_fasta


def simulate(args):
    """Main function that manage simulatin of small RNAs"""
    if args.fasta:
        name = None
        seq = ""
        reads = dict()
        with open(args.fasta) as in_handle:
            for line in in_handle:
                if line.startswith(">"):
                    if name:
                        reads.update(_generate_reads(seq))
                    seq = ""
                    name = line[1:-1]
                else:
                    seq += line.strip()

        reads.update(_generate_reads(seq))
    _write_reads(reads, args.out)


def _generate_reads(seq):
    """Main function that create reads from precursors"""
    reads = dict()
    if len(seq) < 130 and len(seq) > 70:
        reads.update(_mature(seq[:40], 0))
        reads.update(_mature(seq[-40:], len(seq) - 40))
        reads.update(_noise(seq))
        reads.update(_noise(seq, 20))
    return reads


def _mature(subseq, absolute, size=33, total=5000):
    """Create mature sequences around start/end"""
    reads = dict()
    probs = [0.1, 0.2, 0.4, 0.2, 0.1]
    end = 5 + size
    error = [-2, -1, 0, 1, 2]
    for error5 in error:
        for error3 in error:
            s = 5 - error5
            e = end - error3
            seen = subseq[s:e]
            counts = int(probs[error5 + 2] * probs[error3 + 2] * total) + 1
            name = "seq_%s_%s_x%s" % (s, e, counts)
            reads[name] = (seen, counts)
    return reads


def _noise(seq, size=33, total=1000):
    """Create mature sequences around start/end"""
    reads = dict()
    seen = 0
    while seen < total:
        s = random.randint(0, len(seq) - size)
        e = s + size + random.randint(-5,5)
        p = random.uniform(0, 0.1)
        counts = int(p * total) + 1
        seen += counts
        name = "seq_%s_%s_x%s" % (s, e, counts)
        reads[name] = (seq[s:e], counts)
    return reads


def _write_reads(reads, prefix):
    """
    Write fasta file, ma file and real position
    """
    out_ma = prefix + ".ma"
    out_fasta = prefix + ".fasta"
    out_real = prefix + ".txt"
    with open(out_ma, 'w') as ma_handle:
        print >>ma_handle, "id\tseq\tsample"
        with open(out_fasta, 'w') as fa_handle:
            with open(out_real, 'w') as read_handle:
                for idx, r in enumerate(reads):
                    print >>ma_handle, "seq_%s\t%s\t%s" % (idx, reads[r][0], reads[r][1])
                    print >>fa_handle, ">seq_%s\n%s" % (idx, reads[r][0])
                    print >>read_handle, "%s\t%s\t%s\t%s" % (idx, r, reads[r][0], reads[r][1])


def _get_precursor(bed_file, reference, out_fa):
    """
    get sequence precursor from position
    """
    get_fasta(bed_file, reference, out_fa)
    return 0


def _get_spot(precursor):
    """
    get spot that will be enriched
    """
    return 0


def _get_type(pob):
    """
    randomly decide if is small rna or degradation
    """
    return 0


def _random_sequences(precursor, start= None, end=None):
    """
    randomly get sequences around some nucleotides.
    It could be enriched in some positions
    """
    return 0
