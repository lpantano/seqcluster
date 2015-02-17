from unittest import TestCase
from collections import namedtuple
import os
import inspect
import pybedtools
import seqcluster.libs.parameters
import seqcluster.libs.logger as mylog

from seqcluster.libs.tool import reduceloci
from seqcluster.libs.cluster import detect_clusters
from seqcluster.libs.logger import initialize_logger
from seqcluster.libs.tool import parse_ma_file


class TestCluster(TestCase):
    def test_cluster(self):
        mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/", "")
        os.chdir(os.path.join(mod_dir, "data/examples"))
        arg = namedtuple('args', 'debug print_debug MIN_SEQ')
        args = arg(True, True, 1)
        seqL = parse_ma_file("seqs_set.ma")
        c = pybedtools.BedTool("2_clusters_2_seqs_shared")
        initialize_logger(".", args.debug, args.print_debug)
        logger = mylog.getLogger(__name__)
        clus_obj = detect_clusters(c, seqL, args.MIN_SEQ)
        clusLred = reduceloci(clus_obj, ".")
        # self.assertTrue(os.path.exists("pre/seqs.fa"))
