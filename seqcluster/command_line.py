import logging
import  sys
from libs.parse import parse_cl
from preparedata import prepare
from makecluster import cluster
import argparse


def main(**kwargs):
	# set up logging to file - see previous section for more details
	logging.basicConfig(level=logging.DEBUG,
	    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
	    datefmt='%m-%d %H:%M',
	    filename='log',
	    filemode='w')
	# define a Handler which writes INFO messages or higher to the sys.stderr
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	# set a format which is simpler for console use
	formatter = logging.Formatter('%(asctime)s %(name)-12s: %(levelname)-8s %(message)s')
	# tell the handler to use this format
	console.setFormatter(formatter)
	# add the handler to the root logger
	logging.getLogger('').addHandler(console)
	# Now, define a couple of other loggers which might represent areas in your
	# application:

	con = logging.getLogger('console')
	log = logging.getLogger('file')

	if "prepare" in kwargs:
		print "run prepare"
		print kwargs["args"]
		#run.prepare(kwargs["args"],con,log)
	if "cluster" in kwargs:
		print "run cluster"
		#run.cluster(kwargs["args"],con,log)


if __name__ == "__main__":
	kwargs = parse_cl(sys.argv[1:])
	main(**kwargs)
	

