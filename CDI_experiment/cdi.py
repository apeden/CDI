import experiment as e
import plate as p
import calib as cal
import pandas as pd

class CDI(e.Experiment):
    def __init__(self, name, start_date, numReps, tissues,
                 PK_concs, calib_name, calib_concs):
        e.Experiment.__init__(self, name, start_date, numReps)
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
        tissue_table = pd.DataFrame(
            {
                "NAME": pd.Series([t.getName() for t in self.tissues]),
                "DIAGNOSIS": pd.Series([t.getDis_state() for t in self.tissues]),
                "TISSUE": pd.Series([t.getOrgan() for t in self.tissues]),
                "PK CONCS":pd.Series([self.PK_concs.__str__() for t in self.tissues]),
                "REPEATS": pd.Series([self.numReps for t in self.tissues])
            }
        )
        return tissue_table
    def set_samples(self):
        for r in range(self.numReps):
            for t in self.tissues:
                for c in self.PK_concs:
                    self.samples.append((t,c))
    def set_calibs(self):
        for i in self.calib_concs:
            c = cal.Calib(self.calib_name, i)
            self.calibs.append(c)
    def fillPlatesCalib(self):
        self.set_samples()
        self.set_calibs()
        plate_num = 1
        pl = p.Plate(plate_num)
        for t,conc in self.samples:
            if pl.get_capacity() > 5:
                pl.addAnalyte((t,conc))
            else:
                for c in self.calibs:
                    pl.addCalibrant(c)
                self.plates.append(pl)
                plate_num += 1
                pl = p.Plate(plate_num)
                pl.addAnalyte((t,conc))
        for c in self.calibs:
            pl.addCalibrant(c)
        self.plates.append(pl)
    def get_cal_dils(self, subname, subconc, cal_ug_vals, stock_dilFac, 
                     numCurves = 2, sample_vol = 25):
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
    def __str__(self):
        printOut = "\n\n===============\n\n" \
        + "Experiment design\n" \
        + "Tissue vols required: "+str(self.get_tiss_vol()) +"\n\n" \
        + self.get_tissues().__str__()\
        +"\nPlate Summary\n===========\n"
        for pl in self.plates:
            printOut += "Plate number "  + str(pl.getName())+"\n"
            for k, v in pl.get_summary().items():
                printOut += k + ": " + v.__str__() + "\n"
        for pl in self.plates:
            printOut += "\nPlates in detail\n===========\n" + pl.__str__()      
        return printOut
