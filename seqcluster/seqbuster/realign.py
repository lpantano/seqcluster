class realign:

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.precursor = []
        self.score = []
        self.best_hits = [] # maybe sam object?
