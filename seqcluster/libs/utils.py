"""utils from http://www.github.com/chapmanb/bcbio-nextgen.git"""
import os
from bcbio import utils


def safe_dirs(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    return dirs


def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False
