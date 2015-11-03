from unittest import TestCase
from collections import namedtuple
import os
import os.path as op
import inspect
import pybedtools
import seqcluster.libs.parameters
import seqcluster.libs.logger as mylog

from seqcluster.detect.metacluster import reduceloci
from seqcluster.detect.cluster import detect_clusters, peak_calling
from seqcluster.libs.logger import initialize_logger
from seqcluster.libs.inputs import parse_ma_file
from seqcluster.make_clusters import _create_json
from seqcluster.align import pyMatch

class TestCluster(TestCase):
    def test_cluster(self):
        mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/", "")
        os.chdir(os.path.join(mod_dir, "data/examples"))
        if not op.exists("res_cluster"):
            os.mkdir("res_cluster")
        arg = namedtuple('args', 'debug print_debug MIN_SEQ dir_out')
        args = arg(True, True, 1, "res_cluster")
        # seqL = parse_ma_file("seqs_set.ma")
        # c = pybedtools.BedTool("2_clusters_2_seqs_shared")
        # initialize_logger(".", args.debug, args.print_debug)
        # logger = mylog.getLogger(__name__)
        # logger.info("Start reduceloci test")
        # clus_obj = detect_clusters(c, seqL, args.MIN_SEQ)
        # clus_red = reduceloci(clus_obj, ".")
        # logger.info("Peak calling")
        # clus_red = peak_calling(clus_red)
        # logger.info("Write output")
        # _create_json(clus_red, args)
        # self.assertTrue()

class TestAlign(TestCase):
    def test_align(self):
        pyMatch.Match("AAGAGTAAACAGCCTTCTCCCAGCTTTCTTACTTTCCACAGCTGAGAGTGTAGGATGTTTACA", "TGTAAACATCCTACACTCTTAGCT", 1, 3)
