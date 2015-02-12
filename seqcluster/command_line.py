from libs.logger import initialize_logger
import sys
from libs.parse import parse_cl
from preparedata import prepare
from makecluster import cluster
from create_report import report
from explore_cluster import explore
from collapse import collapse_fastq
from stats import stats
import libs.logger as mylog
import time


def main(**kwargs):
    kwargs = parse_cl(sys.argv[1:])
    initialize_logger(kwargs['args'].out, kwargs['args'].debug, kwargs['args'].print_debug)
    logger = mylog.getLogger(__name__)
    start = time.time()
    if "prepare" in kwargs:
        logger.info("Run prepare")
        prepare(kwargs["args"])
    elif "cluster" in kwargs:
        logger.info("Run cluster")
        cluster(kwargs["args"])
    elif "report" in kwargs:
        logger.info("Run report")
        report(kwargs["args"])
    elif "explore" in kwargs:
        logger.info("Run explore")
        explore(kwargs["args"])
    elif "stats" in kwargs:
        logger.info("Run stats")
        stats(kwargs["args"])
    elif "collapse" in kwargs:
        logger.info("Run collapse")
        collapse_fastq(kwargs["args"])
    logger.info('It took %.3f minutes' % ((time.time()-start)/60))
