"""utils from http://www.github.com/chapmanb/bcbio-nextgen.git"""
import os
from bcbio import utils


def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False
