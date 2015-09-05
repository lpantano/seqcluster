from libs.logger import initialize_logger
import sys
from libs.parse import parse_cl
from prepare_data import prepare
from make_clusters import cluster
from create_report import report
from make_predictions import predictions
from explore_cluster import explore
from collapse import collapse_fastq
from seqbuster import miraligner
from libs.simulator import simulate
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
    elif "predict" in kwargs:
        logger.info("Run predictions")
        predictions(kwargs["args"])
    elif "seqbuster" in kwargs:
        logger.info("Run seqbuster")
        miraligner(kwargs["args"])
    elif "explore" in kwargs:
        logger.info("Run explore")
        explore(kwargs["args"])
    elif "stats" in kwargs:
        logger.info("Run stats")
        stats(kwargs["args"])
    elif "collapse" in kwargs:
        logger.info("Run collapse")
        collapse_fastq(kwargs["args"])
    elif "simulator" in kwargs:
        logger.info("Run simulator")
        simulate(kwargs["args"])
    logger.info('It took %.3f minutes' % ((time.time()-start)/60))
