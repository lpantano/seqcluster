import os
from libs.fastq import collapse, splitext_plus, write_output
import logging


logger = logging.getLogger('seqbuster')


def collapse_fastq(args):
    """collapse fasq files after adapter trimming
    """
    try:
        seqs = collapse(args.fastq)
        out_file = splitext_plus(os.path.basename(args.fastq))[0] + "_trimmed.fastq"
    except IOError as e:
        logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
        raise "Can not read file"
    logger.info("writing output")
    out_file = os.path.join(args.out, out_file)
    write_output(out_file, seqs)
