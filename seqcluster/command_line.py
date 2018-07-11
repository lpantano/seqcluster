from seqcluster.libs.logger import initialize_logger
import sys
from seqcluster.libs.parse import parse_cl
from seqcluster.prepare_data import prepare
from seqcluster.make_clusters import cluster
from seqcluster.create_report import report
from seqcluster.make_predictions import predictions
from seqcluster.explore_cluster import explore
from seqcluster.collapse import collapse_fastq
from seqcluster.seqbuster import miraligner
from seqcluster.libs.simulator import simulate
from seqcluster.function.target import targets_enrichment
from seqcluster.stats import stats
import seqcluster.libs.logger as mylog
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
    elif "target" in kwargs:
        logger.info("Run target annotation")
        targets_enrichment(kwargs["args"])
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
