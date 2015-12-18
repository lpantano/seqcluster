"""utils from http://www.github.com/chapmanb/bcbio-nextgen.git"""
import os
import shutil

import contextlib

@contextlib.contextmanager
def chdir(p):
    cur_dir =  os.getcwd()
    os.chdir(p)
    try:
        yield
    finally:
        os.chdir(cur_dir)

def safe_dirs(dirs):
    if not os.path.exists(dirs):
        os.makedirs(dirs)
    return dirs

def safe_remove(fn):
    if os.path.exists(fn):
        if os.path.isfile(fn):
            os.remove(fn)
        elif os.path.isdir(fn):
            shutil.rmtree(fn)

@contextlib.contextmanager
def safe_run(fn):
    safe_remove(fn)
    try:
        yield
    except:
        safe_remove(fn)
        raise

def file_exists(fname):
    """Check if a file exists and is non-empty.
    """
    try:
        return fname and os.path.exists(fname) and os.path.getsize(fname) > 0
    except OSError:
        return False
