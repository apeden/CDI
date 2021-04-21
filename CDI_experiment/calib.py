import material as m

class Calib(m.Material):
    def __init__(self, name, conc):
        m.Material.__init__(self, name)
        self.conc = conc
    def getConc(self):
        return self.conc
