"""
sam2bed.py

Created by Aaron Quinlan on 2009-08-27.
Copyright (c) 2009 Aaron R. Quinlan. All rights reserved.
"""

from classes import bedaligned


def processSAM(line):
    """
    Load a SAM file and convert each line to BED format.
    """
    samLine = splitLine(line.strip())
    return makeBED(samLine)


def makeBED(samFields):
    samFlag = int(samFields.flag)
    # Only create a BED entry if the read was aligned
    if (not (samFlag & 0x0004)) and samFields.pos:
            name = samFields.qname
            seq = name.split("-")[1]
            chrom = samFields.rname
            start = str(int(samFields.pos))
            end = str(int(samFields.pos) + len(seq) - 1)
            strand = getStrand(samFlag)
            # Write out the BED entry
            return bedaligned("%s\t%s\t%s\t%s\t.\t%s\n" %  (chrom, start, end, name, strand))
    else:
            return False


def splitLine(line, delim="\t"):
    splitline = line.split(delim)
    return splitline


def getStrand(samFlag):
    strand = "+"
    if (samFlag & (0x10)):  # minus strand if true.
            strand = "-"
    return strand
