from unittest import TestCase
from collections import namedtuple
import os
import inspect
import seqcluster.libs.parameters
import seqcluster.libs.logger as mylog

from seqcluster.libs.read import load_data
from seqcluster.libs.predictions import make_predictions
from seqcluster.libs.logger import initialize_logger


class TestPredict(TestCase):
    def test_predict(self):
        mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/", "")
        os.chdir(os.path.join(mod_dir, "data/examples"))
        arg = namedtuple('args', 'debug print_debug MIN_SEQ json ref')
        args = arg(True, True, 1, "seqcluster.json", "../genomes/genome.fa")
        initialize_logger(".", args.debug, args.print_debug)
        logger = mylog.getLogger(__name__)
        logger.info("reading data")
        load_data(args.json)
        # make_predictions(c, out_dir, args)
        # self.assertTrue(os.path.exists("pre/seqs.fa"))
