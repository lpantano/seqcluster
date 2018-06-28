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
from seqcluster.libs.read import _align


class TestAlign(TestCase):
    def test_align(self):
        _align("TGTAAACATCCTACACTCTTAGCT", "AAGAGTAAACAGCCTTCTCCCAGCTTTCTTACTTTCCACAGCTGAGAGTGTAGGATGTTTACA")
