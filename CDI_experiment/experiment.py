import material as m

class Experiment(m.Material):
    def __init__(self, name, start_date, numReps):
        m.Material.__init__(self, name)
        self.start_date = start_date
        self.numReps = numReps
