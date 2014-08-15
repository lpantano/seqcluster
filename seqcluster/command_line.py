from libs.logger import initialize_logger
import sys
from libs.parse import parse_cl
from preparedata import prepare
from makecluster import cluster
import libs.logger as mylog
import time


def main(**kwargs):
    kwargs = parse_cl(sys.argv[1:])
    initialize_logger(kwargs['args'].out)
    logger = mylog.getLogger(__name__)
    start = time.time()
    if "prepare" in kwargs:
		logger.info("Run prepare")
		prepare(kwargs["args"])
    elif "cluster" in kwargs:
		logger.info("Run cluster")
		cluster(kwargs["args"])
    logger.info('It took %.3f minutes' % ((time.time()-start)/60))

