import material as m
import pandas as pd
import numpy as np
import math

class Plate(m.Material):
    vol_per_well = 200
    def __init__(self, name):
        m.Material.__init__(self, name)
        self.comp = pd.DataFrame(
            {
                "ROW": pd.Series(list("ABCDEFGH"*12)),
                "COLUMN": pd.Series([math.ceil(x/8) for x in range(1,97)]),
                "WELL": pd.Series(list(range(1,97)), dtype="int32"),
            }
        )
        self.plate_summary = {}
        self.triblock = 1 ## CDI plate filled three wells at a time
    def addAnalyte(self, compon):
        tissue, conc = compon
        point = 1
        point += (self.triblock-1)*4 
        self.comp.loc[point: (point+5), "Sample"] = tissue.getName()
        self.comp.loc[point: (point+5), "PK"] = conc
        self.comp.loc[point: (point+2), "Denatured"] = False
        self.comp.loc[point+3: (point+5), "Denatured"] = True
        try:
            self.plate_summary[tissue.getName()].append(conc)
        except KeyError:
            self.plate_summary[tissue.getName()] = [conc]
        self.triblock += 2
    def addCalibrant(self, compon):
        point = 1
        if self.triblock%2 != 0:
            point += (self.triblock-1)*4
        else:
            point += (self.triblock-2)*4 + 3
        self.comp.loc[point: (point+2), "Sample"] = compon.getName()
        self.comp.loc[point: (point+2), "ug/well recPrP"] = compon.getConc()
        self.comp.loc[point: (point+2), "Denatured"] = True
        self.triblock += 1
    def addData(self, file, start = (9,4)):
        row, column = start
        inData = pd.read_csv(file)
        self.comp["Data"] = inData.iloc[row:,column:].reset_index(drop = True)
        self.comp["Data"] = self.comp["Data"].astype(float)
    def calibrate(self, calibrant_name):
        calibrant_df = self.comp[self.comp["Sample"].isin([calibrant_name])]   
        polynomial_coeff = np.polyfit(calibrant_df["ug/well recPrP"],
                                      calibrant_df["Data"] ,1)
        self.comp["Slope"] = polynomial_coeff[0]
        self.comp["ug PrP per g brain"] = (self.comp["Data"]/polynomial_coeff[0])*400
    def get_capacity(self):
        return 25 - self.triblock
    def get_buff_vol(self):
        return math.ceil((self.triblock*200*3*1.1)/1000)
    def getComps(self): ##get plate components
        return self.comp
    def get_ab_dil(self, factor):
        printOut = (str(round(self.get_buff_vol() * (1000/factor),1))+ " ul Ab in " + \
str(self.get_buff_vol()) + "mls")
        return printOut
    def get_summary(self):
        return self.plate_summary
    def __str__(self):
        printOut = "\n\nPlate " + str(self.name) + "\n=====\n" \
        + "Vol of buffers: "+ str(self.get_buff_vol()) + "\n" \
        + "biotin 3F4 and Eu-strep dils: "\
        + "\n "+ self.get_ab_dil(5000) + "\n "\
        + self.comp.__str__()
        return printOut
