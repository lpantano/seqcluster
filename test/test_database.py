from unittest import TestCase
from collections import namedtuple
import os
import inspect
import contextlib
import shutil
import seqcluster.libs.parameters
import seqcluster.libs.logger as mylog

from seqcluster.libs.do import find_cmd
from seqcluster.libs.read import load_data
from seqcluster.db import make_database
from seqcluster.libs.logger import initialize_logger
from nose.plugins.attrib import attr

@contextlib.contextmanager
def make_workdir():
    remove_old_dir = True
    dirname = os.path.join(os.path.dirname(__file__), "test_automated_output")
    if remove_old_dir:
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield dirname
    finally:
        os.chdir(orig_dir)


class TestDatabase(TestCase):
    @attr(database=True)
    def test_database(self):
        if find_cmd("sqlite3"):
            with make_workdir() as workdir:
                arg = namedtuple('args', 'debug print_debug json ')
                args = arg(True, True, "../../data/examples/seqcluster.json")
                initialize_logger(".", args.debug, args.print_debug)
                logger = mylog.getLogger(__name__)
                logger.info(args)
                logger.info("Reading data")
                data = load_data(args.json)
                logger.info("Create database")
                make_database(data)
        # self.assertTrue()
