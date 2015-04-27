"""prepare data for CoRaL"""


def prepare_bam_file(args):
    """
    Clean BAM file to keep only position inside the bigger cluster
    """
    # use pybedtools to keep valid positions
    # intersect option with -b bigger_cluster_loci


def prepare_ann_file(args):
    """
    Create custom ann_file for Coral
    """


def download_hsa_file(args):
    """
    In case of human, download from server
    """
