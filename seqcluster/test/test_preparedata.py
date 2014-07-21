from unittest import TestCase
import os 
import inspect
from seqcluster.preparedata import _read_fasta_files,_create_matrix_uniq_seq
import seqcluster

class TestPreparedata(TestCase):
    def test_preparedata(self):
    	mod_dir = os.path.dirname(inspect.getfile(seqcluster)).replace("seqcluster/","")
    	os.chdir(os.path.join(mod_dir,"data"))
    	os.remove("pre/seqs.ma")
    	os.remove("pre/seqs.fa")
        seq_l, list_s = _read_fasta_files(f=open("config"))
        ma_out = open("pre/seqs.ma",'w')
        seq_out = open("pre/seqs.fa",'w')
        _create_matrix_uniq_seq(list_s,seq_l,ma_out,seq_out)
        self.assertTrue(os.path.exists("pre/seqs.ma"))
        self.assertTrue(os.path.exists("pre/seqs.fa"))
