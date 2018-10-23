"""This directory is setup with configurations to run the main functional test.

Inspired in bcbio-nextgen code
"""
from __future__ import print_function
import os
import subprocess
import unittest
import shutil
import contextlib
import collections
import functools

from nose import SkipTest
from nose.plugins.attrib import attr
import yaml


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

def expected_failure(test):
    """Small decorator to mark tests as expected failure.
    Useful for tests that are work-in-progress.
    """
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner

class AutomatedAnalysisTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "examples")

    def _install_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        #       self._download_to_dir(url, dirname)

    def _download_to_dir(self, url, dirname):
        print(dirname)
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        shutil.move(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))

    @attr(complete=True)
    @attr(cluster=True)
    def test_srnaseq_cluster(self):
        """Run cluster analysis
        """
        with make_workdir() as workdir:
            cl = ["seqcluster",
                  "cluster",
                  "-m", "../../data/examples/cluster/seqs_set.ma",
                  "-a", "../../data/examples/cluster/seqs_map.bam",
                  "--gtf", "../../data/examples/cluster/ann_reduced.gtf",
                  "-r", "../../data/genomes/genome.fa",
                  "-o", "test_out_res"]
            print(" ".join(cl))
            subprocess.check_call(cl)
            cl = ["seqcluster",
                  "report",
                  "-j", "test_out_res/seqcluster.json",
                  "-r", "../../data/genomes/genome.fa",
                  "-o", "test_out_report"]
            print(" ".join(cl))
            subprocess.check_call(cl)

    @attr(complete=True)
    @attr(miraligner=True)
    def test_srnaseq_miraligner(self):
        """Run miraligner analysis
        """
        with make_workdir() as workdir:
            cl = ["seqcluster",
                  "seqbuster",
                  "--sps", "hsa",
                  "--hairpin", "../../data/examples/miraligner/hairpin.fa",
                  "--mirna", "../../data/examples/miraligner/miRNA.str",
                  "--gtf", "../../data/examples/miraligner/hsa.gff3",
                  "-o", "test_out_mirs_fasta",
                  "--miraligner",
                  "../../data/examples/miraligner/sim_isomir.fa"]
            print(" ".join(cl))
            # subprocess.check_call(cl)
