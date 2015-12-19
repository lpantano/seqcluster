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
        os.chdir(os.path.join(mod_dir, "data/examples/collapse"))
        out_dir = "test_out_prepare"
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)
        arg = namedtuple('args', 'minl maxl minc out')
        args = arg(15, 40, 1, out_dir)
        seq_l, list_s = _read_fastq_files(open("config"), args)
        ma_out = open(os.path.join(out_dir, "seqs.ma"), 'w')
        seq_out = open(os.path.join(out_dir, "seqs.fa"), 'w')
        _create_matrix_uniq_seq(list_s, seq_l, ma_out, seq_out, 1)
        os.chdir(out_dir)
        self.assertTrue(os.path.exists("seqs.ma"))
        self.assertTrue(os.path.exists("seqs.fa"))
