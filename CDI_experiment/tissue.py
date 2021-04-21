import material as m

class Tissue(m.Material):
    def __init__(self, name, dis_state, organ):
        m.Material.__init__(self, name)
        self.dis_state = dis_state
        self.organ = organ
        self.PKtreat = 0.0
    def getDis_state(self):
        return self.dis_state
    def setPKtreat(self, val):
        self.PKtreat = val
    def getPKtreat(self):
        return self.PKtreat
    def getOrgan (self):
        return self.organ
