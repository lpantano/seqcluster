"""
Main classes used in seqcluster
"""

import copy
import numpy as np
from operator import add
from collections import Counter, defaultdict


class sequence_unique:
    """
    Object to store the sequence information like: **counts**, **sequence**, **id**
    """
    def __init__(self, idx, seq):
        self.idx = idx
        self.seq = seq
        self.group = {}
        self.quality = ""
        self.total = 0
    def add_exp(self,gr,exp):
        """Function to add the counts for each sample

        :param gr: name of the sample
        :param exp: counts of sample **gr**

        :returns: dict with key,values equally to name,counts.
        """
        self.group[gr] = exp
        self.total = sum(self.group.values())


class quality:

    def __init__(self, q):
        self.qual = [ord(value) for value in q]
        self.times = 1

    def update(self, q, counts = 1):
        now = self.qual
        self.qual = map(add, now, [ord(value) for value in q])
        self.times += counts

    def get(self):
        average = map(lambda x: int(round(x/self.times)), self.qual)
        return [str(unichr(char)) for char in average]


class cluster_info_obj:
    """
    Object containing information about clusters(:code:`clus_obj`), 
    positions(:code:`positions`) and sequences(:code:`sequences`)
    """
    def __init__(self, clus_obj, clus_id, loci_obj, seq_obj):
        self.clus = clus_obj
        self.clusid = clus_id
        self.loci = loci_obj
        self.seq = seq_obj


class sequence:
    """
    Object with information about sequences, counts, size, position, id and score
    """
    def __init__(self, seq_id, seq=None, freq=None):
        # self.seq = seq
        # self.freq = copy.deepcopy(freq)
        # self.norm_freq = copy.deepcopy(freq)
        self.pos = {}
        self.id = seq_id
        self.align = 0
        self.score = 0
        self.factor = {}

    def set_seq(self, seq):
        self.seq = seq
        self.len = len(seq)

    def set_freq(self, freq):
        self.freq = copy.deepcopy(freq)
        self.norm_freq = copy.deepcopy(freq)

    def add_pos(self, pos_id, pos):
        self.pos[pos_id] = pos

    def total(self):
        return sum(self.freq.values())


class position:
    """
    Object with information about position: chr,start,end,strand
    as well, with annotation information throuhg :code:`dbannotation` object 
    """
    def __init__(self, idl, chr, start, end, strand):
        self.idl = idl
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.coverage = Counter()
        self.counts = Counter()
        self.db_ann = {}

    def list(self):
        return map(str, [self.chr, self. start, self.end, self.idl, self.strand])

    def add_db(self, db, ndb):
        self.db_ann[db] = ndb


class annotation:
    """
    Object with information about annotation: database, name of the feature, strand,
    distance to 5' end, distance to 3' end
    """
    def __init__(self,db,name,strand,to5,to3):
        self.db = db
        self.name = name
        self.strand = strand
        self.to5 = to5
        self.to3 = to3


class dbannotation:
    """
    Object with information about annotation: containg one dict that 
    store all features for each database type
    """
    def __init__(self,na):
        self.ann = {}
    def add_db_ann(self,ida,ndba):
        self.ann[ida] = ndba


class cluster:
    """
    Object with cluster information. This is the main object.
    """
    def __init__(self, id):
        self.id = id
        self.idmembers = defaultdict(int)
        self.locimax = None
        self.locimaxid = None
        self.locilen = {}
        self.loci2seq = {}
        self.ref = 0
        self.score = 0
        self.peaks = []
        self.showseq = ""
        self.showseq_plain = ""
        self.toomany = 0
        self.predictions = {}
        self.errors = []
        self.freq = []

    def normalize(self, seq, factor):
        return dict(zip(seq.freq.keys(), list(np.array(seq.freq.values()) * factor)))

    def set_freq(self, seqL):
        total = Counter()
        [total.update(self.normalize(seqL[s], f)) for s, f in self.idmembers.iteritems()]
        self.freq = total
        return total

    def get_freq(self, seqL, force=False):
        self.update()
        if self.freq and not force:
            return self.freq
        else:
            return self.set_freq(seqL)

    def set_ref(self, r):
        self.ref = r

    def update(self, id=None):
        if id:
            self.id = id
        # self.idmembers = defaultdict(int)
        seen = set()
        self.locimax = 0
        for idl in self.loci2seq:
            l = len(self.loci2seq[idl])
            # self.idmembers.update(dict(zip(self.loci2seq[idl], [1] * l)))
            seen = seen.union(set(self.loci2seq[idl]))
            if l > self.locimax:
                self.locimax = l
                self.locimaxid = idl
        remove = set(self.idmembers.keys()) - seen
        add = seen - set(self.idmembers.keys())
        self.idmembers.update(dict(zip(add, [1] * len(add))))
        map(self.idmembers.__delitem__, remove)

    def add_id_member(self, ids, idl):
        for s in ids:
            self.idmembers[s] = 1
            if idl not in self.loci2seq:
                self.loci2seq[idl] = []
            self.loci2seq[idl].append(s)
        self.loci2seq[idl] = list(set(self.loci2seq[idl]))
        lenid = len(self.loci2seq[idl])
        self.locilen[idl] = lenid
        if lenid > self.locimax:
            self.locimax = lenid
            self.locimaxid = idl


class bcolors:
    HEADER  =  '\033[95m'
    OKBLUE  =  '\033[94m'
    OKGREEN  =  '\033[92m'
    WARNING  =  '\033[93m'
    FAIL  =  '\033[91m'
    ENDC  =  '\033[0m'


class bedaligned:
    """
    Object that has the bed format attributes
    """
    def __init__(self,l):
        l = l.strip()
        cols = l.split("\t")
        self.chr = cols[0]
        self.start = cols[1]
        self.end = cols[2]
        self.name = cols[3]
        self.att = cols[4]
        self.strand = cols[5]


class mergealigned:
    """
    Object that has bed format after merge sequence positions
    """
    def __init__(self,l):
            self.chr = l[0]
            self.strand = l.strand
            self.start = l.start
            self.end = l.end
            self.names = list(l.name.split(","))
            self.loci = l.score.split(",")
