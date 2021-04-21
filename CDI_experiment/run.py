import pandas as pd
import plate 
import tissue as t
import cdi as c
import calib
import pk_data as pPK
import plotBar as pB

pd.set_option('display.max_rows',  None)

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
    tissues.append(t.Tissue(key, "FFI", value))
for key, value in fCJD.items():
    tissues.append(t.Tissue(key, "fCJD", value))
for key, value in sFI.items():
    tissues.append(t.Tissue(key, "sFI", value))

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
##            calib_test.append(calib.Calib(s,c))    
##
##cal_p = pl.Plate("caib plate")
##for te in calib_test:
##    cal_p.addCalibrant(te)
##print(cal_p)

data_files = ("CDI_21_002_column",
              "CDI_21_003_column",
              "CDI_21_004_column",
              "CDI_21_005_column",
              "CDI_21_006_column",
              "CDI_21_007_column",
              "CDI_21_008_column",
              "CDI_21_009_column",
              "CDI_21_010_column",
              "CDI_21_011_column")

e = c.CDI(1, "March 2021", 2, tissues, PK_concs, "431", cal_ug_vals)
e.fillPlatesCalib()
plates = e.getPlates()
assert len(plates) == len(data_files)
for i in range(len(plates)):
    plates[i].addData("data_to_analyse/"+data_files[i]+".csv")
    plates[i].calibrate("431")
#first_plate = e.getPlates()[0]
#first_plate.addData("data_to_analyse/CDI_21_002_column.csv")
#first_plate.calibrate("431")
#print (first_plate)
frames = []
for p in plates:
    frames.append(p.getComps())
df = pd.concat(frames)

chart = pB.PlotBar(df, PK_concs, tissues)
chart.assemblePlot()
