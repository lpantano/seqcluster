import argparse
import sys


def parse_cl(in_args):
    print in_args
    sub_cmds = {"prepare": add_subparser_prepare,
                "cluster": add_subparser_cluster,
                "seqbuster": add_subparser_mirbuster,
                "report": add_subparser_report,
                "predict": add_subparser_predict,
                "explore": add_subparser_explore,
                "collapse": add_subparser_collapse,
                "simulator": add_subparser_simulator,
                "stats": add_subparser_stats}
    parser = argparse.ArgumentParser(description="small RNA analysis")
    sub_cmd = None
    if len(in_args) > 0 and in_args[0] in sub_cmds:
        subparsers = parser.add_subparsers(help="seqcluster supplemental commands")
        sub_cmds[in_args[0]](subparsers)
        sub_cmd = in_args[0]
    else:
        print "use %s" % sub_cmds.keys()
        sys.exit(0)
    args = parser.parse_args()

    assert sub_cmd is not None
    kwargs = {"args": args, sub_cmd: True}
    return kwargs


def _add_debug_option(parser):
    parser.add_argument("-d", "--debug", action="store_true",
                        dest="debug", help="max verbosity mode", default=False)
    parser.add_argument("-vd", "--print_debug", action="store_true",
                        help="print debug messageson terminal", default=False)
    return parser


def add_subparser_report(subparsers):
    parser = subparsers.add_parser("report", help="report data")
    parser.add_argument("-j", "--json", dest="json", required=1,
            help="json file from seqcluster")
    parser.add_argument("-o", "--out", dest="out", required=1,
            help="dir of output files")
    parser.add_argument("-r", "--reference", dest="ref", required=1,
            help="reference fasta file with index"),
    parser.add_argument("--razer", action="store_true",
            help="map sequences with razer3")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_mirbuster(subparsers):
    parser = subparsers.add_parser("seqbuster", help="realign miRNA BAM file")
    parser.add_argument("files", nargs="*", help="Bam files.")
    parser.add_argument("-o", "--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("--sps", required=1,
                        help="species")
    parser.add_argument("--hairpin", help="hairpin.fa")
    parser.add_argument("--mirna", help="miRNA.str")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_explore(subparsers):
    parser = subparsers.add_parser("explore", help="explore data")
    parser.add_argument("-j", "--json", dest="json", required=1,
            help="json file from seqcluster")
    parser.add_argument("-n", "--names", dest="names", required=1,
            help="comma-separeted id clusters"),
    parser.add_argument("-r", "--reference", dest="ref", required=1,
            help="reference fasta file with index"),
    parser.add_argument("-o", "--out", dest="out", required=1,
            help="dir of output files")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_prepare(subparsers):
    parser = subparsers.add_parser("prepare", help="prepare data")
    parser.add_argument("-c", "--conf", dest="config", required=1,
            help="file with config file:1st column:path_to_fasta_file ; 2nd column:name")
    parser.add_argument("-o", "--out", dest="out", required=1,
            help="output dir")
    parser.add_argument("-l", "--minl", dest="minl", required=0,
            help="minimum length", default=18)
    parser.add_argument("-u", "--maxl", dest="maxl", required=0,
            help="maximum length", default=35)
    parser.add_argument("-e", "--minc", dest="minc", required=0,
            help="minimum counts", default=10)
    parser.add_argument("--min-shared", dest="min_shared", required=0,
            help="minimum shamples with same sequences", default=2)
    parser = _add_debug_option(parser)
    return parser


def add_subparser_cluster(subparsers):
    parser = subparsers.add_parser("cluster", help="cluster data")
    parser.add_argument("-a", "--afile", dest="afile", required=1,
                      help="aligned file in bam format")
    parser.add_argument("-m", "--ma", dest="ffile", required=1,
                      help="matrix file with sequences and counts for each sample")
    parser.add_argument("-g", "--gtf",
                       dest="gtf", help="annotate with gtf_file. It will use the 3rd column as the tag to annotate" +
                       "\nchr1    source  intergenic      1       11503   .       +       .       ")
    parser.add_argument("-b", "--bed",
                       dest="bed", help="annotate with bed_file. It will use the 4rd column as the tag to annotate" +
                       "\nchr1    157783  157886  snRNA   0       -")
    parser.add_argument("-o", "--out",
                       dest="out", help="output dir", required=1)
    parser.add_argument("-ref", required=1,
                       dest="ref", help="reference fasta")
    parser.add_argument("--mask",
                        help="bed file with regions to mask")
    parser = _add_debug_option(parser)
    parser.add_argument("-s", "--show", action="store_true",
                       dest="show", help="no show sequences", default=False)
    parser.add_argument("--non-un-gl", action="store_true",
                        help="remove Un_gl chromosomes", default=False)
    parser.add_argument("--method", choices=["most-voted", "split", "bayes"],
                       dest="method", help="most-voted, split, bayes", default='most-voted')
    parser.add_argument("--similar",
                       dest="similar", help="threshold to consider two clusters identicals", default=0.8)
    parser.add_argument("--min_seqs",
                       dest="min_seqs", help="threshold to consider a cluster as valid", default=10)
    parser.add_argument("--db",
                        help="prefix for sqlite3 database with results to use htmlViz plugin (in dev).")
    return parser


def add_subparser_stats(subparsers):
    parser = subparsers.add_parser("stats", help="stats data")
    parser.add_argument("-j", "--json", dest="json", required=0,
            help="json file from seqcluster")
    parser.add_argument("-m", "--ma", dest="ma", required=0,
            help="seqs.ma from prepare"),
    parser.add_argument("-a", "--sam", dest="sam", required=0,
            help="aligned file")
    parser.add_argument("-o", "--out",
                       dest="out", help="output dir", required=1)
    parser = _add_debug_option(parser)
    return parser


def add_subparser_collapse(subparsers):
    parser = subparsers.add_parser("collapse", help="collapse data")
    parser.add_argument("-f", "--fastq", dest="fastq", required=1,
                        help="fastq file"),
    parser.add_argument("-o", "--out",
                        dest="out", help="output file", required=1)
    parser = _add_debug_option(parser)
    return parser

    parser = _add_debug_option(parser)
    return parser


def add_subparser_predict(subparsers):
    parser = subparsers.add_parser("predict", help="predict smallRNA types")
    parser.add_argument("-j", "--json", dest="json", required=1,
            help="json file from seqcluster")
    parser.add_argument("--bed", help="BED ouput from cluster to clean BAM file")
    parser.add_argument("--bam", help="BAM file used in cluster subcmd.")
    parser.add_argument("-o", "--out", dest="out", required=1,
            help="dir of output files")
    parser.add_argument("--reference", required=1,
            help="reference fasta file with index")
    parser.add_argument("--coral", action='store_true',
            help="Run CoRaL pipeline")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_explore(subparsers):
    parser = subparsers.add_parser("explore", help="explore data")
    parser.add_argument("-j", "--json", dest="json", required=1,
            help="json file from seqcluster")
    parser.add_argument("-n", "--names", dest="names", required=1,
            help="comma-separeted id clusters"),
    parser.add_argument("-r", "--reference", dest="ref", required=1,
            help="reference fasta file with index"),
    parser.add_argument("-o", "--out", dest="out", required=1,
            help="dir of output files")
    parser = _add_debug_option(parser)
    return parser


def add_subparser_prepare(subparsers):
    parser = subparsers.add_parser("prepare", help="prepare data")
    parser.add_argument("-c", "--conf", dest="config", required=1,
            help="file with config file:1st column:path_to_fasta_file ; 2nd column:name")
    parser.add_argument("-o", "--out", dest="out", required=1,
            help="output dir")
    parser.add_argument("-l", "--minl", dest="minl", required=0,
            help="minimum length", default=18)
    parser.add_argument("-u", "--maxl", dest="maxl", required=0,
            help="maximum length", default=35)
    parser.add_argument("-e", "--minc", dest="minc", required=0,
            help="minimum counts", default=10)
    parser.add_argument("--min-shared", dest="min_shared", required=0,
            help="minimum shamples with same sequences", default=2)
    parser = _add_debug_option(parser)
    return parser


def add_subparser_cluster(subparsers):
    parser = subparsers.add_parser("cluster", help="cluster data")
    parser.add_argument("-a", "--afile", dest="afile", required=1,
                      help="aligned file in bam format")
    parser.add_argument("-m", "--ma", dest="ffile", required=1,
                      help="matrix file with sequences and counts for each sample")
    parser.add_argument("-g", "--gtf",
                       dest="gtf", help="annotate with gtf_file. It will use the 3rd column as the tag to annotate" +
                       "\nchr1    source  intergenic      1       11503   .       +       .       ")
    parser.add_argument("-b", "--bed",
                       dest="bed", help="annotate with bed_file. It will use the 4rd column as the tag to annotate" +
                       "\nchr1    157783  157886  snRNA   0       -")
    parser.add_argument("-o", "--out",
                       dest="out", help="output dir", required=1)
    parser.add_argument("-ref",
                       dest="ref", help="reference fasta")
    parser.add_argument("--mask",
                        help="bed file with regions to mask")
    parser = _add_debug_option(parser)
    parser.add_argument("-s", "--show", action="store_true",
                       dest="show", help="no show sequences", default=False)
    parser.add_argument("--non-un-gl", action="store_true",
                        help="remove Un_gl chromosomes", default=False)
    parser.add_argument("--method", choices=["most-voted", "split", "bayes"],
                       dest="method", help="most-voted, split, bayes", default='most-voted')
    parser.add_argument("--similar",
                       dest="similar", help="threshold to consider two clusters identicals", default=0.8)
    parser.add_argument("--min_seqs",
                       dest="min_seqs", help="threshold to consider a cluster as valid", default=10)
    parser.add_argument("--db",
                        help="prefix for sqlite3 database with results to use htmlViz plugin (in dev).")
    return parser


def add_subparser_stats(subparsers):
    parser = subparsers.add_parser("stats", help="stats data")
    parser.add_argument("-j", "--json", dest="json", required=0,
            help="json file from seqcluster")
    parser.add_argument("-m", "--ma", dest="ma", required=0,
            help="seqs.ma from prepare"),
    parser.add_argument("-a", "--sam", dest="sam", required=0,
            help="aligned file")
    parser.add_argument("-o", "--out",
                       dest="out", help="output dir", required=1)
    parser = _add_debug_option(parser)
    return parser


def add_subparser_collapse(subparsers):
    parser = subparsers.add_parser("collapse", help="collapse data")
    parser.add_argument("-f", "--fastq", dest="fastq", required=1,
                        help="fastq file"),
    parser.add_argument("-o", "--out",
                        dest="out", help="output file", required=1)
    parser = _add_debug_option(parser)
    return parser


def add_subparser_simulator(subparsers):
    parser = subparsers.add_parser("simulator", help="simulate small RNA  from bed file")
    parser.add_argument("--bed",
                        help="bed file with position of precursors <=200 nt")
    parser.add_argument("--fasta", help = "fasta with precursors.")
    parser.add_argument("--out", dest="out", required=1,
                        help="dir of output files")
    parser.add_argument("-r", "--reference", dest="ref",
                        help="reference fasta file with index"),
    parser = _add_debug_option(parser)
    return parser
