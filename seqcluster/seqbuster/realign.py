from collections import defaultdict

class realign:

    def __init__(self):
        self.sequence = ""
        self.precursors = defaultdict(isomir)
        self.score = []
        self.best_hits = [] # maybe sam object?

    def set_precursor(self, precursor, isomir):
        self.precursors[precursor] = isomir

    def remove_precursor(self, precursor):
        del self.precursors[precursor]

class isomir:

    def __init__(self):
        self.t5 = 0
        self.t3 = 0
        self.add = []
        self.subs = []
        self.align = None
        self.start = 0
        self.mirna = None

    def format(self):
        return "subs %s add %s %s %s" % (self.subs, self.add,
                                         self.t5, self.t3)

    def get_score(self, sc):
        for a in self.add:
            if a in ['A', 'T']:
                sc -= 0.25
            else:
                sc -= 0.75
        for e in self.subs:
            sc -= 1
        return sc
