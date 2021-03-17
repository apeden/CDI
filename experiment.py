import math

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
        self.comp = {"11BD":"HaRecPrP 1 Denat",
                     "11EG":"HaRecPrP 0.75 Denat",
                     "12BD":"HaRecPrP 0.5 Denat",
                     "12EG":"HaRecPrP 0.25 Denat"}
        self.next_col = 1
    def addAnalyte(self, compon):
        self.comp[str(self.next_col)+"BD"] = compon + " Nat"
        self.comp[str(self.next_col)+"EG"] = compon + " Denat"
        self.next_col += 1
    def has_space(self):
        return (self.next_col < 11)
    def get_buff_vol(self):
        return math.ceil((len(self.comp)*200*3*1.1)/1000)
    def getComps(self):
        return self.comp
    def print_mar1_dil(self):
        pass
    def get_ab_dil(self, factor):
        printOut = (str(self.get_buff_vol() * (1000/factor))+ " ul Ab in " + \
str(self.get_buff_vol()) + "mls")
        return printOut
    def __str__(self):
        printout = str(self.name) + "\n=====\n"
        printout += "Vol of buffers: "+ str(self.get_buff_vol()) + "\n"
        for key, value in self.comp.items():
            printout += key + ": "+ value + "\n"
        return printout

class Tissue(Material):
    def __init__(self, name, dis_state, organ):
        Material.__init__(self, name)
        self.dis_state = dis_state
        self.organ = organ
    def getDis_state(self):
        return self.dis_state
    def getOrgan (self):
        return self.organ

class Experiment(Material):
    def __init__(self, name, start_date, numReps):
        Material.__init__(self, name)
        self.start_date = start_date
        self.numReps = numReps

class CDI(Experiment):
    def __init__(self, name, start_date, numReps, tissues, PK_concs):
        Experiment.__init__(self, name, start_date, numReps)
        self.tissues = tissues
        self.PK_concs = PK_concs
        self.tiss_vol = self.numReps * len(PK_concs) * 60  
        self.plates = []
    def print_PKprep(self):
        printOut = "To prepare PK dils, for adding 1.5ul per 60ul brain homog:\n"
        for c in PK_concs:
            PK_vol = 50
            printOut += \
str((1-(c/50))*PK_vol)+"ul water plus " + \
str((c/50)*PK_vol) + "ul 2mg/ml PK\n"
        print (printOut)
    def get_tiss_vol(self):
        return self.tiss_vol
    def get_tissues(self):
        return self.tissues
    def fillPlates(self):
        plate_num = 1
        p = Plate(plate_num)
        for r in range(self.numReps):
            for t in self.tissues:
                for c in self.PK_concs:
                    if p.has_space():
                        p.addAnalyte(t.getName()+ ": " + str(c))
                    else:
                        self.plates.append(p)
                        plate_num += 1
                        p = Plate(plate_num)
                        p.addAnalyte(t.getName()+ ": " + str(c))
    def print_sandwichAb_dils(self):
        printOut = "Mar1, biotin 3F4 and Eu-strep dils\n====\n"
        for p in self.plates:          
            printOut += "\nPlate  " + str(p.getName())
            printOut += "\n "+ p.get_ab_dil(5000)
        print(printOut)
    def __str__(self):
        printOut  = "Plates:\n"
        printOut += "Tissue vols required: "+str(self.get_tiss_vol()) +"\n"
        for p in self.plates:
            printOut += p.__str__()
        return printOut
    
for key, value in FFI.items():
    tissues.append(Tissue(key, "FFI", value))
for key, value in fCJD.items():
    tissues.append(Tissue(key, "fCJD", value))
for key, value in sFI.items():
    tissues.append(Tissue(key, "sFI", value))


e = CDI("1", "12-3-21", 2, tissues, PK_concs)
e.fillPlates()
e.print_PKprep()
e.print_sandwichAb_dils()
print(e)


subconc = {"431":0.2515, "432":0.2253333}
cal_ug_vals = [0.25, 0.50, 0.75, 1.00]

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

##get_cal_dils("old HuRePrP", 0.691, cal_ug_vals, 0.1, numCurves = 2)
##get_cal_dils("431", 0.2515, cal_ug_vals, 1, numCurves = 3)
##get_cal_dils("432", 0.2253333, cal_ug_vals, 1, numCurves = 3)




