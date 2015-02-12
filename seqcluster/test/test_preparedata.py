from unittest import TestCase
from collections import namedtuple
import os
import inspect
from seqcluster.preparedata import _read_fastq_files, _create_matrix_uniq_seq
import seqcluster


class TestPreparedata(TestCase):
    def test_preparedata(self):
        mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/","")
        os.chdir(os.path.join(mod_dir, "data/test_collapse"))
        if os.path.exists("seqs.ma"):
            os.remove("seqs.ma")
        if os.path.exists("seqs.fa"):
            os.remove("seqs.fa")
        arg = namedtuple('args', 'minl maxl minc')
        args = arg(15, 40, 1)
        seq_l, list_s = _read_fastq_files(open("config"), args)
        ma_out = open("seqs.ma", 'w')
        seq_out = open("seqs.fa", 'w')
        _create_matrix_uniq_seq(list_s, seq_l, ma_out, seq_out)
        self.assertTrue(os.path.exists("seqs.ma"))
        self.assertTrue(os.path.exists("seqs.fa"))
