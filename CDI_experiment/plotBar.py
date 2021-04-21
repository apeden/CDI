import matplotlib.pyplot as plt
import pk_data as pPK
import numpy as np
import math

class PlotBar(object):
    def __init__(self, df, PK_concs, tissues):
        self.df = df
        self.PK_concs = PK_concs
        self.tissues = tissues
        self.bar = None ##will be a PKdata object
    def assembleData(self, tissue):
        self.bar = pPK.PKdata(self.df, self.PK_concs, tissue.getName())
        self.bar.set_df_nat()
        self.bar.set_df_denat()
        self.bar.set_DN()
        self.bar.set_plot_data()
    def assemblePlot(self):
        fig, axs = plt.subplots(2,4)
        ##axs is an array of axes
        plt.subplots_adjust(left=0.1,
                            bottom=0.1, 
                            right=0.9, 
                            top=0.9, 
                            wspace=0.3, 
                            hspace=0.3)
        for i in range(len(self.tissues)):
            self.assembleData(self.tissues[i])
            subplot = axs[i//4,i%4] 
            subplot.set_ylim(bottom=0.)
            N = len(self.bar.get_conc_list())
            ind = np.arange(N)    # the x locations for the groups
            width = 0.35       # the width of the bars: can also be len(x) sequence
            subplot.bar(ind,
                        self.bar.get_mean_vals(),
                        width,
                        yerr = self.bar.get_DNstd())
            subplot.set_ylabel(chr(956)+"g PrPSc /gram Brain")
            subplot.set_xlabel("PK "+chr(956)+"g/ml")
            subplot.set_title(self.tissues[i].getDis_state()+", "+self.tissues[i].getName())
            subplot.set_xticks(ind)
            subplot.set_xticklabels(self.bar.get_conc_list())
            heighest_val = max(self.bar.get_mean_vals()) + max(self.bar.get_DNstd())
            gradation = 10
            y_scale = 100 + gradation
            if heighest_val > 100:
                y_scale = (round(heighest_val, -int(math.log(heighest_val,10))))
                gradation = y_scale/10
                y_scale += gradation
            subplot.set_yticks(np.arange(0, y_scale, gradation))
        plt.show()
