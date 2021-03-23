import csv
import numpy as np
import pandas as pd

calib_concs = [0.0, 0.25, 0.5, 0.75, 1.0]
calib_vals = np.array([[4,3,4],[1,2,1]])
#np.append(calib_vals, [[3,2,3],[2,4,4]]
print(calib_vals)
class Plate_data_csv(object):
    def __init__(self, file):
        self.file = open(file,'r')
        self.contents = csv.reader(self.file)
        self.intMatrix = []
    def getContents(self):
        return self.contents
    def getFile(self):
        return self.file
    def setMatrix(self):
        for r in self.contents:
            self.intMatrix.append(r)
    def getMatrix(self):
        return self.intMatrix
    def __str__(self):
        toPrint = ""
        for r in self.contents:
            toPrint += r.__str__() + "\n"
        return toPrint

class Plate_data_pandas(Plate_data_csv):
    def __init__(self, file):
        Plate_data_csv.__init__(self, file)
        self.contents = pd.read_csv(file)
    def __str__(self):
        return self.contents.__str__()
        

d = pd.read_csv("CDI_21_001_column.csv")
print(d)
s = d.iloc[9:,4]
s = s.reset_index(drop = True)
s = pd.to_numeric(s)
print(s[:8])
p_df["CDI vals"]= s
print(p_df.iloc[:20,:])

def calib_plot(df , sub):
    df_sub = df[df["Sample"]==sub]
    polynomial_coeff=np.polyfit(df_sub["ug per well recPrP"], df_sub["CDI vals"],1)
    print(polynomial_coeff)
    plt.scatter(x=df_sub["ug per well recPrP"], y=df_sub["CDI vals"])
    plt.xlabel("ug per well recPrP")
    plt.ylabel("CDI values")
    plt.title("HuRecPrP " + sub + "; Slope: "+ str(polynomial_coeff.round()[0]))
    xnew=np.linspace(0,1,100)
    ynew=np.poly1d(polynomial_coeff)
    plt.plot(xnew,ynew(xnew),'r')
    plt.show()

calib_plot(p_df , "432")














p = Plate_data_pandas('CDI_21_001.csv')
print(p)




##
##print("Program to demonstrate csv.reader() function:")
##print("\n")
##c = open('CDI_21_001.csv','r')
##print("The csv file is opened for reading:")
##print("\n")
##o = csv.reader(c)
##print("The contents of the above file is as follows:")
##for r in o:
##    print (r)
##c.close()
