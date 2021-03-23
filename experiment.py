import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.set_option('display.max_rows', None)

FFI = {"44/2019(MV)":"FC",
       "82/2017(MM)":"CNS",
       "38/2014(MM)":"CNS",
       "89/2011(MV)":"OC"}
fCJD = {"30/2002":"CNS",
       "124/2001":"FC",
       "119/2011":"FC"}
sFI = {"05/2000":"FC"}

tissues = []
PK_concs = [0.0, 1.0, 2.5, 5.0, 10.0, 50.0] 

class Material(object):
    def __init__(self, name):
        self.name = name
    def getName(self):
        return self.name

class Plate(Material):
    vol_per_well = 200
    def __init__(self, name):
        Material.__init__(self, name)
        self.comp = pd.DataFrame(
            {
                "ROW": pd.Series(list("ABCDEFGH"*12)),
                "COLUMN": pd.Series([(x//8)+1 for x in range(1,97)]),
                "WELL": pd.Series(list(range(1,97)), dtype="int32"),
                "ug per well recPrP": 0.0,
            }
        )
        self.triblock = 1
    def addAnalyte(self, compon):
        tissue, conc = compon
        point = 1
        point += (self.triblock-1)*4
        self.comp.loc[point: (point+2), "Denatured"] = False
        self.comp.loc[point: (point+2), "Sample"] = tissue.getName()
        self.comp.loc[point: (point+2), "PK"] = conc
        self.comp.loc[point+3: (point+5), "Denatured"] = True
        self.comp.loc[point+3: (point+5), "Sample"] = tissue.getName()
        self.comp.loc[point+3: (point+5), "PK"] = conc
        self.triblock += 2
    def addCalibrant(self, compon):
        point = 1
        if self.triblock%2 != 0:
            point += (self.triblock-1)*4
        else:
            point += (self.triblock-2)*4 + 3
        self.comp.loc[point: (point+2), "Sample"] = compon.getName()
        self.comp.loc[point: (point+2), "ug per well recPrP"] = compon.getConc()
        self.triblock += 1
    def get_capacity(self):
        return 25 - self.triblock
    def get_buff_vol(self):
        return math.ceil((self.triblock*200*3*1.1)/1000)
    def getComps(self): ##get plate components
        return self.comp
    def get_ab_dil(self, factor):
        printOut = (str(self.get_buff_vol() * (1000/factor))+ " ul Ab in " + \
str(self.get_buff_vol()) + "mls")
        return printOut
    def __str__(self):
        printOut = "\nPlate " + str(self.name) + "\n=====\n"
        printOut += "Vol of buffers: "+ str(self.get_buff_vol()) + "\n"
        printOut += "biotin 3F4 and Eu-strep dils: "
        printOut += "\n "+ self.get_ab_dil(5000) + "\n "
        printOut += self.comp.__str__()
        return printOut

class Tissue(Material):
    def __init__(self, name, dis_state, organ):
        Material.__init__(self, name)
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

class Calib(Material):
    def __init__(self, name, conc):
        Material.__init__(self, name)
        self.conc = conc
    def getConc(self):
        return self.conc

class Experiment(Material):
    def __init__(self, name, start_date, numReps):
        Material.__init__(self, name)
        self.start_date = start_date
        self.numReps = numReps

class CDI(Experiment):
    def __init__(self, name, start_date, numReps, tissues,
                 PK_concs, calib_name, calib_concs):
        Experiment.__init__(self, name, start_date, numReps)
        self.tissues = tissues
        self.PK_concs = PK_concs
        self.calib_name = calib_name
        self.calib_concs = calib_concs
        self.tiss_vol = self.numReps * len(PK_concs) * 60
        self.samples , self.plates, self.calibs = [],[],[]
    def print_PKprep(self):
        printOut = "To prepare PK dils, for adding 1.5ul per 60ul brain homog:\n"
        for c in self.PK_concs:
            PK_vol = 50
            printOut += \
str((1-(c/50))*PK_vol)+"ul water plus " + \
str((c/50)*PK_vol) + "ul 2mg/ml PK\n"
        print (printOut)
    def get_tiss_vol(self):
        return self.tiss_vol
    def getPlates(self):
        return self.plates
    def get_tissues(self):
        return self.tissues
    def set_samples(self):
        for r in range(self.numReps):
            for t in self.tissues:
                for c in self.PK_concs:
                    self.samples.append((t,c))
    def set_calibs(self):
        for i in self.calib_concs:
            c = Calib(self.calib_name, i)
            self.calibs.append(c)
    def fillPlatesCalib(self):
        self.set_samples()
        self.set_calibs()
        plate_num = 1
        p = Plate(plate_num)
        for t,c in self.samples:
            if p.get_capacity() > 5:
                p.addAnalyte((t,c))
            else:
                for c in self.calibs:
                    p.addCalibrant(c)
                self.plates.append(p)
                plate_num += 1
                p = Plate(plate_num)
##                p.addAnalyte((t,c))
        self.plates.append(p)
    def __str__(self):
        printOut = "Plates:\n"
        printOut += "Tissue vols required: "+str(self.get_tiss_vol()) +"\n"
        for p in self.plates:
            printOut += p.__str__()
        return printOut
    
def get_cal_dils(subname, subconc, cal_ug_vals, stock_dilFac, 
                 numCurves = 6, sample_vol = 25):
    cal_vol = sample_vol * (numCurves + 1)
    cal_concs, cal_vals, buf_vals = [],[],[]
    dil_stock_conc = stock_dilFac*subconc    
    for u in cal_ug_vals:
        cal_concs.append(u/sample_vol)
    for conc in cal_concs:
        cal_vals.append(round(cal_vol*(conc/dil_stock_conc),1))
    for val in cal_vals:
        buf_vals.append(cal_vol-val)
    print(subname + "\n========")
    print("numCurves: ",str(numCurves))
    print("Stock dil Fac: " + str(stock_dilFac))
    for i in range(len(cal_ug_vals)):    
        printOut = str(cal_ug_vals[i])+ ": "
        printOut += "Vol 1 in 10 stock dil: " + str(cal_vals[i])
        printOut += " Vol of buf: " + str(buf_vals[i])
        print(printOut)
    return cal_vals, buf_vals

for key, value in FFI.items():
    tissues.append(Tissue(key, "FFI", value))
for key, value in fCJD.items():
    tissues.append(Tissue(key, "fCJD", value))
for key, value in sFI.items():
    tissues.append(Tissue(key, "sFI", value))

subconc = {"431":0.2515, "432":0.2253333}
cal_ug_vals = [0.25, 0.50, 0.75, 1.00]

####get_cal_dils("old HuRePrP", 0.691, cal_ug_vals, 0.1, numCurves = 2)
##get_cal_dils("431", 0.2515, concs, 1, numCurves = 3)
##get_cal_dils("432", 0.2253333, concs, 1, numCurves = 3)

##calib_test = []
##reps = 2
##subs = "431" , "432"
##concs = [0.0, 0.25, 0.50, 0.75, 1.00]
##for r in range(reps):
##    for s in subs:
##        for c in concs:
##            calib_test.append(Calib(s,c))    
##e = CDI("1", "12-3-21", 2, tissues, PK_concs)
##e.fillPlates(calib_test)
##p = e.getPlates()[0]
##p_df = p.getComps()
##print(p_df.iloc[0:45,:])

e =CDI(1, "", 2, tissues, PK_concs, "431", cal_ug_vals)
e.fillPlatesCalib()
print(e)
