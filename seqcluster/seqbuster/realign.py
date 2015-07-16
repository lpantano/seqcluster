from collections import defaultdict

class realign:

    def __init__(self):
        self.sequence = ""
        self.precursor = defaultdict(isomir)
        self.score = []
        self.best_hits = [] # maybe sam object?

    def set_precursor(self, precursor, isomir):
        self.precursor[precursor] = isomir

class isomir:

    def __init__(self):
        self.t5 = ""
        self.t3 = ""
        self.add = ""
        self.subs = ""
