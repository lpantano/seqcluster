from unittest import TestCase
from collections import namedtuple
import os
import shutil
import inspect
from seqcluster.prepare_data import _read_fastq_files, _create_matrix_uniq_seq
import seqcluster
from nose.plugins.attrib import attr


class TestPreparedata(TestCase):
    @attr(collapse=True)
    def test_preparedata(self):
        out_dir = "test/test_out_prepare"
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)
        arg = namedtuple('args', 'minl maxl minc out')
        args = arg(15, 40, 1, out_dir)
        seq_l, list_s = _read_fastq_files(open("data/examples/collapse/config"), args)
        ma_out = open(os.path.join(out_dir, "seqs.ma"), 'w')
        seq_out = open(os.path.join(out_dir, "seqs.fa"), 'w')
        _create_matrix_uniq_seq(list_s, seq_l, ma_out, seq_out, 1)
        self.assertTrue(os.path.exists(os.path.join(out_dir, "seqs.ma")))
        self.assertTrue(os.path.exists(os.path.join(out_dir, "seqs.fa")))
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)

    @attr(umis=True)
    def test_umis(self):
        from seqcluster.libs.fastq import collapse, write_output
        umis = collapse(os.path.abspath("data/examples/umis/sample.fastq"))
        if len(umis.keys()) != 2:
            raise ValueError("umis didn't detect two unique sequences")
        out_dir = "test/test_automated_output"
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)
        os.mkdir(out_dir)
        write_output(os.path.join(out_dir, "umis.fastq"), umis)
