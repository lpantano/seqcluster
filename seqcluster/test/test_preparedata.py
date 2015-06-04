from unittest import TestCase
from collections import namedtuple
import os
import shutil
import inspect
from seqcluster.prepare_data import _read_fastq_files, _create_matrix_uniq_seq
import seqcluster


class TestPreparedata(TestCase):
    def test_preparedata(self):
        mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/","")
        os.chdir(os.path.join(mod_dir, "data/test_collapse"))
        if os.path.exists("prepare"):
            shutil.rmtree("prepare")
        os.mkdir("prepare")
        arg = namedtuple('args', 'minl maxl minc out')
        args = arg(15, 40, 1, "prepare")
        seq_l, list_s = _read_fastq_files(open("config"), args)
        ma_out = open("prepare/seqs.ma", 'w')
        seq_out = open("prepare/seqs.fa", 'w')
        _create_matrix_uniq_seq(list_s, seq_l, ma_out, seq_out, 1)
        self.assertTrue(os.path.exists("prepare/seqs.ma"))
        self.assertTrue(os.path.exists("prepare/seqs.fa"))
