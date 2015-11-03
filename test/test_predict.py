from unittest import TestCase
from collections import namedtuple
import os
import inspect
import seqcluster.libs.parameters
import seqcluster.libs.logger as mylog

from seqcluster.libs.do import find_cmd
from seqcluster.libs.read import load_data
from seqcluster.function.predictions import is_tRNA
from seqcluster.libs.logger import initialize_logger


class TestPredict(TestCase):
    def test_predict(self):
        if find_cmd("tRNAscan-SE"):
            mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/", "")
            os.chdir(os.path.join(mod_dir, "data/examples"))
            out_dir = os.path.join(mod_dir, "data/examples/predictions")
            arg = namedtuple('args', 'debug print_debug MIN_SEQ json reference')
            args = arg(True, True, 1, "seqcluster.json", "../genomes/genome.fa")
            initialize_logger(".", args.debug, args.print_debug)
            logger = mylog.getLogger(__name__)
            logger.info(args)
            logger.info("Reading data")
            data = load_data(args.json)
            logger.info("Start prediction")
            is_tRNA(data, out_dir, args)
        # self.assertTrue()
