"""
Some commands to install common databases like mirbase and some human/mouse annotation
"""
from __future__ import print_function
import os.path as op
import os
import sys
from argparse import ArgumentParser
import subprocess
import contextlib
import yaml

from seqcluster.libs import do, utils

try:
    import bcbio
except:
    print ("Probably this will fail, you need bcbio-nextgen for many installation functions.")
    pass

REMOTES = {
            "requirements": "https://raw.github.com/lpantano/seqcluster/master/requirements.txt",
            "gitrepo": "https://github.com/lpantano/seqcluster.git",
            "gitrepo-bcbio": "https://github.com/chapmanb/bcbio-nextgen.git"
           }

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

def _get_miraligner():
    opts = "-Xms750m -Xmx4g"
    try:
        tool = "miraligner"
        ret = os.system(tool)
        if ret != 0:
            raise SystemExit("%s not installed." % tool)
    except SystemExit:
        tool = None
        print("miraligner not found. I'll try to download it.")
        pass
    if not tool:
        if not utils.file_exists(op.abspath("miraligner.jar")):
            url = "https://raw.githubusercontent.com/lpantano/seqbuster/miraligner/modules/miraligner/miraligner.jar"
            cmd = ["wget", "-O miraligner.jar", "--no-check-certificate", url]
            do.run(" ".join(cmd), "Download miraligner.")
        tool = "java -jar {opts} %s" % op.abspath("miraligner.jar")
    else:
        tool = "%s {opts}" % tool
    return tool.format(**locals())

def _get_flavor():
    """
    Download flavor from github
    """
    target = op.join("seqcluster", "flavor")
    url = "https://github.com/lpantano/seqcluster.git"
    if not os.path.exists(target):
    #   shutil.rmtree("seqcluster")
        subprocess.check_call(["git", "clone","-b", "flavor", "--single-branch", url])
    return op.abspath(target)

def _install(path, args):
    """
    small helper for installation in case outside bcbio
    """
    try:
       from bcbio import install as bcb
    except:
        raise ImportError("It needs bcbio to do the quick installation.")

    path_flavor = _get_flavor()
    s = {"fabricrc_overrides": {"system_install": path,
                                "local_install": os.path.join(path, "local_install"),
                                "use_sudo": "false",
                                "edition": "minimal"}}
    s = {"flavor": path_flavor,
         # "target": "[brew, conda]",
         "vm_provider": "novm",
         "hostname": "localhost",
         "fabricrc_overrides": {"edition": "minimal",
                                "use_sudo": "false",
                                "keep_isolated": "true",
                                "conda_cmd": bcb._get_conda_bin(),
                                "distribution": "__auto__",
                                "dist_name": "__auto__"}}


    s["actions"] = ["install_biolinux"]
    s["fabricrc_overrides"]["system_install"] = path
    s["fabricrc_overrides"]["local_install"] = os.path.join(path, "local_install")
    cbl = bcb.get_cloudbiolinux(bcb.REMOTES)
    sys.path.insert(0, cbl["dir"])
    cbl_deploy = __import__("cloudbio.deploy", fromlist=["deploy"])
    cbl_deploy.deploy(s)

def _install_data(data_dir, path_flavor, args):
    """Upgrade required genome data files in place.
    """
    try:
       from bcbio import install as bcb
    except:
        raise ImportError("It needs bcbio to do the quick installation.")

    bio_data = op.join(path_flavor, "../biodata.yaml")
    s = {"flavor": path_flavor,
         # "target": "[brew, conda]",
         "vm_provider": "novm",
         "hostname": "localhost",
         "fabricrc_overrides": {"edition": "minimal",
                                "use_sudo": "false",
                                "keep_isolated": "true",
                                "conda_cmd": bcb._get_conda_bin(),
                                "distribution": "__auto__",
                                "dist_name": "__auto__"}}
    s["actions"] = ["setup_biodata"]
    s["fabricrc_overrides"]["data_files"] = data_dir
    s["fabricrc_overrides"]["galaxy_home"] = os.path.join(data_dir, "galaxy")
    cbl = bcb.get_cloudbiolinux(bcb.REMOTES)
    s["genomes"] = _get_biodata(bio_data, args)
    sys.path.insert(0, cbl["dir"])
    cbl_deploy = __import__("cloudbio.deploy", fromlist=["deploy"])
    cbl_deploy.deploy(s)

def _get_biodata(base_file, args):
    with open(base_file) as in_handle:
        config = yaml.load(in_handle)
    config["install_liftover"] = False
    config["genome_indexes"] = args.aligners
    config["genomes"] = [g for g in config["genomes"] if g["dbkey"] in args.genomes]
    return config

def _install_mirbase():
    for fn in ["hairpin.fa.gz", "miRNA.str.gz"]:
        out_file = op.join("mirbase", fn)
        _mkdir("mirbase")
        url = "ftp://mirbase.org/pub/mirbase/CURRENT/%s" % fn
        cmd = ["wget", "-O", out_file, "--no-check-certificate", url]
        subprocess.check_call(cmd)
        subprocess.check_call(["gunzip", "-f", out_file])
    return "mirbase/hairpin.fa", "mirbase/miRNA.str"

def _upgrade():
    conda_dir = os.path.join(os.path.dirname(os.path.realpath(sys.executable)))
    try:
        import bcbio.install as install
        install._set_pip_ssl(conda_dir)
    except ImportError:
        print("bcbio was not found, this may error.")
        pass
    pip_bin = os.path.join(conda_dir, "pip")
    cmd = [pip_bin, "install", "--upgrade", "--no-deps",
           "git+%s#egg=seqcluster" % REMOTES["gitrepo"]]
    print(" ".join(cmd))
    subprocess.check_call(cmd)
    #subprocess.check_call([pip_bin, "install", "--upgrade", "--no-deps",
    #                       "git+%s#egg=bcbio-nextgen" % REMOTES["gitrepo-bcbio"]])

def actions(args):
    if args.upgrade:
        _upgrade()
    if False:
        db = set(args.data) if isinstance(args.data, list) else [args.data]
        if "mirbase" in db:
            _install_mirbase()
        if "hg19" in db:
            with chdir("hg19"):
                subprocess.check_call(["wget", "--no-check-certificate", "-p", "https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/hg19.sh", "-O", "hg19.sh", ])
                subprocess.check_call(["bash", "hg19.sh"])
        if "mm10" in db:
            with chdir("mm10"):
                subprocess.check_call(["wget", "--no-check-certificate", "-p", "https://raw.githubusercontent.com/lpantano/seqcluster/master/scripts/mm10.sh", "-O", "mm10.sh", ])
                subprocess.check_call(["bash", "mm10.sh"])
    if args.tools:
        _install(op.abspath(args.tools), args)
    if args.data:
        path_flavor = _get_flavor()
        _install_data(args.data, path_flavor, args)

def main(**kwargs):
    parser = ArgumentParser(description="small RNA analysis installer")
    parser.add_argument("--tools", help="install tools")
    parser.add_argument("--data", help="path install data")
    parser.add_argument("--upgrade", action="store_true", help="upgrade seqcluster", default=[])
    parser.add_argument("--genomes", default=[], action="append")
    parser.add_argument("--aligners", default=[], action="append")
    args = parser.parse_args()
    actions(args)
