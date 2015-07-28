"""
Some commands to install common databases like mirbase and some human/mouse annotation
"""
import os.path as op
import os
import sys
import subprocess
import contextlib

def _mkdir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


@contextlib.contextmanager
def chdir(new_dir):
    """
    stolen from bcbio.
    Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    _mkdir(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)


def main(args):
    db = set(args)
    if "mirbase" in db:
        for fn in ["hairpin.fa.gz", "miRNA.str.gz"]:
            out_file = op.join("mirbase", fn)
            _mkdir("mirbase")
            url = "ftp://mirbase.org/pub/mirbase/CURRENT/%s" % fn
            cmd = ["wget", "-O", out_file, "--no-check-certificate", url]
            subprocess.check_call(cmd)
            subprocess.check_call(["gunzip", "-f", out_file])
    if "hg19" in db:
        with chdir("hg19"):
            subprocess.check_call(["wget", "--no-check-certificate", "-p", "https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/hsapiens.sh", "-O", "hg19.sh", ])
            subprocess.check_call(["bash", "hg19.sh"])
    if "mm10" in db:
        with chdir("mm10"):
            subprocess.check_call(["wget", "--no-check-certificate", "-p", "https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/mmusculus.sh", "-O", "mm10.sh", ])
            subprocess.check_call(["bash", "mm10.sh"])


if __name__ == "__main__":
    main(sys.argv[1:])
