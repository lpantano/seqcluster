
def parse_cl(in_args):
	usagetxt = "usage: %prog  -c file -o output"
	parser = argparse.ArgumentParser(usage=usagetxt, version="%prog 0.1")


	

	args = parser.parse_args()


	assert sub_cmd is not None
    kwargs = {"args": args,
                  sub_cmd: True}
    return kwargs

def add_subparser_prepare(subparsers):
	parser = subparsers.add_parser("prepare", help="prepare data")
	parser.add_option("-c", "--conf", dest="dir",
                  help="file with fasta format paths:1st column:path 2nd column:name", metavar="FILE")
	parser.add_option("-o", "--output", dest="out",
                  help="dir of output files", metavar="FILE")
	return parser

def add_subparser_cluster(subparsers):
	parser = subparsers.add_parser("prepare", help="prepare data")
	parser.add_option("-a", "--afile", dest="afile",
	                  help="aligned file in bed/sam format", metavar="FILE")
	parser.add_option("-m", "--ma", dest="ffile",
	                  help="matrix file with sequences and counts for each sample", metavar="FILE")
	parser.add_option("-g", "--gtf",
	                   dest="gtf", help="annotate with gtf_file. It will use the 3rd column as the tag to annotate" +
	                   "\nchr1    source  intergenic      1       11503   .       +       .       ",metavar="GTF")
	parser.add_option("-b", "--bed",
	                   dest="bed", help="annotate with bed_file. It will use the 4rd column as the tag to annotate" +
	                   "\nchr1    157783  157886  snRNA   0       -",metavar="BED")
	parser.add_option("-o", "--out",
	                   dest="out", help="output dir",metavar="FILE")
	parser.add_option("-i", "--index",
	                   dest="index", help="reference fasta",metavar="FILE")
	parser.add_option("-d", "--debug", action="store_true",
	                   dest="debug", help="max verbosity mode",default=False)
	parser.add_option("-s", "--show", action="store_false",
	                   dest="show", help="no show sequences",default=True)
	return parser

