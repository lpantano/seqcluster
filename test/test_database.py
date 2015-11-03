from unittest import TestCase
from collections import namedtuple
import os
import inspect
import seqcluster.libs.parameters
import seqcluster.libs.logger as mylog

from seqcluster.libs.do import find_cmd
from seqcluster.libs.read import load_data
from seqcluster.db import make_database
from seqcluster.libs.logger import initialize_logger


class TestDatabase(TestCase):
    def test_databse(self):
        if find_cmd("sqlite3"):
            mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/", "")
            os.chdir(os.path.join(mod_dir, "data/examples"))
            arg = namedtuple('args', 'debug print_debug json ')
            args = arg(True, True, "seqcluster.json")
            initialize_logger(".", args.debug, args.print_debug)
            logger = mylog.getLogger(__name__)
            logger.info(args)
            logger.info("Reading data")
            data = load_data(args.json)
            logger.info("Create databse")
            make_database(data)
        # self.assertTrue()
