from libs.logger import initialize_logger
import sys
from libs.parse import parse_cl
from preparedata import prepare
from makecluster import cluster


def main(**kwargs):
    kwargs = parse_cl(sys.argv[1:])
    initialize_logger(kwargs['args'].out)
    if "prepare" in kwargs:
		print "run prepare"
		prepare(kwargs["args"])
    if "cluster" in kwargs:
		print "run cluster"
		cluster(kwargs["args"])


