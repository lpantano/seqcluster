"""simulate cluster over the genome"""
from read import get_fasta

def _get_precursor(bed_file, reference, out_fa):
    """
    get sequence precursor from position
    """
    get_fasta(bed_file, reference, out_fa)
    return 0


def _get_spot(precursor):
    """
    get spot that will be enriched
    """
    return 0


def _get_type(pob):
    """
    randomly decide if is small rna or degradation
    """
    return 0


def _random_sequences(precursor, start= None, end=None):
    """
    randomly get sequences around some nucleotides.
    It could be enriched in some positions
    """
    return 0
