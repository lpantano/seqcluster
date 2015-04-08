from optparse import OptionParser

usagetxt = "usage: %prog  -r res_seqcluster -c clusters_id -o outputdir"
parser = OptionParser(usage=usagetxt, version="%prog 0.1")
parser.add_option("-c", "--clusterid", dest="clus",
                  help="", metavar="FILE")
parser.add_option("-r", "--results", dest="res",
                  help="", metavar="FILE")
parser.add_option("-o", "--output", dest="out",
                  help="", metavar="FILE")
(options, args) = parser.parse_args()

clus = {}
c = open(options.clus, 'r')
for line in c:
    line = line.strip()
    clus[line] = 0
c.close()

isin = 0
f = open(options.res, 'r')
for line in f:
    line = line.strip()
    if not "#" in line:
        if "C:" in line:
            if isin==1:
                out.close()
                outfa.close()
                isin = 0
            if clus.has_key(line):
                out = open(options.out+"/"+line, 'w')
                outfa = open(options.out+"/"+line+".fa", 'w')
                out.write(line+"\n")
                isin = 1
        else:
            if isin == 1:
                out.write(line+"\n")
            if "seq_" in line:
                cols = line.split(" ")
                outfa.write(">"+cols[0]+"\n"+cols[1]+"\n")
