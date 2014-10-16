import os
import pysam
import logging
from libs.tool import parse_ma_file
from libs.sam2bed import makeBED
from collections import defaultdict, Counter

logger = logging.getLogger('stats')


def stats(args):
    """Create stats from the analysis
    """
    logger.info("Reading sequeces")
    data = parse_ma_file(args.ma)
    logger.info("Get sequences from sam")
    is_align = _read_sam(args.sam)
    res = _summarise_sam(data, is_align)
    _write_suma(res, os.path.join(args.out, "stats_align.dat"))
    logger.info("Done")


def _read_sam(sam):
    is_align = set()
    with pysam.Samfile(sam, "r") as samfile:
        for a in samfile.fetch():
            a = makeBED(a)
            if a:
                is_align.add(a.name)
    return is_align


def _summarise_sam(counts, is_align):
    suma = defaultdict(Counter)
    for s, o in counts.iteritems():
        l = len(o.seq)
        if s in is_align:
            suma[l].update(o.freq)
    return suma

def _write_suma(d, fn):
    with open(fn, 'w') as handle:
        for l, c in d.iteritems():
            for e, f in c.iteritems():
                handle.write("%s %s %s\n" % (l, e, f))
