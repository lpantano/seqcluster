#import sys
import os
#from os import listdir
#from os.path import isfile, join
import re
import logging
from libs.classes import sequence_unique
from libs.tools import parse_ma_file

logger = logging.getLogger('explore')


def explore(args):
    """Create mapping of sequences of two clusters
    """
    logger.info("reading sequeces")
    seq_l = parse_ma_file(args.f)
    logger.info("get sequences from json")
    logger.info("map all sequences to all loci")
    logger.info("plot sequences on loci")
    logger.info("Done")
