import os
from libs.fastq import collapse, splitext_plus
import logging


logger = logging.getLogger('seqbuster')


def collapse_fastq(args):
    """collapse fasq files after adapter trimming
    """
    idx = 0
    try:
        seqs = collapse(args.fastq)
        out_file = splitext_plus(os.path.basename(args.fastq))[0] + "_trimmed.fastq"
    except IOError as e:
        logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
        raise "Can not read file"
    logger.info("writing output")
    with open(os.path.join(args.out, out_file), 'w') as handle:
        for seq in seqs:
            idx += 1
            qual = "".join(seqs[seq].get())
            counts = seqs[seq].times
            handle.write(("@seq_{idx}_x{counts}\n{seq}\n+\n{qual}\n").format(**locals()))
