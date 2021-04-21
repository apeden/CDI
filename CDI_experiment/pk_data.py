import matplotlib.pyplot as plt
import numpy as np
import math

class PKdata (object):
    def __init__(self, df, conc_list, tissue):
        self.df = df
        self.df_denat = None
        self.df_nat = None
        self.conc_list = conc_list
        self.tissue = tissue
        self.mean_vals = []
        self.DNstd =  []
        self.n = []
    def get_conc_list(self):
        return self.conc_list
    def set_df_nat(self):
        self.df_nat = self.df[self.df["Denatured"]== False]
        self.df_nat = self.df_nat[self.df_nat["Sample"] == self.tissue]
        self.df_nat.reset_index(inplace = True)
    def set_df_denat(self):
        self.df_denat = self.df[self.df["Denatured"]== True]
        self.df_denat = self.df_denat[self.df_denat["Sample"] == self.tissue]
        self.df_denat.reset_index(inplace = True)
    def set_DN(self):
        self.df_denat["D-N"] = self.df_denat["Data"]- self.df_nat["Data"]
    def get_df_nat(self):
        return self.df_nat
    def get_df_denat(self):
        return self.df_denat
    def set_plot_data(self):
        mean, std, n = 0.0, 0.0, 0
        for conc in self.conc_list:
            DNvals = self.df_denat[self.df_denat["PK"] == conc]["D-N"]
            mean, std, n = DNvals.mean(), DNvals.std(), DNvals.size
            self.mean_vals.append(mean)
            self.DNstd.append(std)
            self.n.append(n)
    def get_mean_vals(self):
        return self.mean_vals
    def get_DNstd(self):
        return self.DNstd
    def get_n(self):
        return self.n

