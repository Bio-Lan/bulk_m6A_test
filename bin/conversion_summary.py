#!/usr/bin/env python

import argparse
import os
import sys
from collections import defaultdict
import pandas as pd

import utils

class Conversion_summary:
    """
    filter conversion site,get m6A loci
    get DRACH motif and gene m6A site counts
    """
    def __init__(self,args):
        # input
        self.args = args
        self.sample = args.sample
        well_loci = args.well_loci.split(",")

        # set
        self.well_inf = {}
        for i in well_loci:
            well = os.path.basename(i).split(".")[0]
            self.well_inf[well] = i
        
        self.num_dict = {}
        self.check_dict = {}

        # output
        os.makedirs(args.sample, exist_ok=True)
        self.loci_suffix = ".m6A_loci.csv"
        self.DRACH_suffix = ".loci_motif.csv"
        self.gene_suffix = ".loci_in_gene.csv"

    def run(self):
        for well in self.well_inf.keys():
            loci_file = self.well_inf[well]
            df = pd.read_csv(loci_file,header=0,sep=",")
            snp_num = df.shape[0]
            if snp_num == 0:
                continue
            
            m6A_df = self.m6A_confirm(well,df)
            m6A_num = m6A_df.shape[0]
            if m6A_num > 0:
                self.DRACH_gene_counts(well,m6A_df)
            
            self.num_dict[well] = {"conversion_site":snp_num,"m6A_site":m6A_num}
            
        loci_df = pd.DataFrame.from_dict(self.num_dict,orient='index')
        outfile = f"{self.sample}.conversion.csv"
        loci_df.to_csv(outfile,index=True)

        if self.args.check:
            check_df = pd.DataFrame.from_dict(self.check_dict,orient='index')
            outfile = f"{self.sample}.check_base.csv"
            check_df.to_csv(outfile,index=True)

    def m6A_confirm(self,well,df):
        df["posratio"] = df["convs"] / df["covers"]
        df = df[df["convs"] >= self.args.min_m6A_depth]
        df = df[df["covers"] >= self.args.min_coverage]
        df = df[df["posratio"] >= self.args.min_m6A_threshold]
        df = df[df["posratio"] <= self.args.max_m6A_threshold]
        
        if self.args.check:
            df = self.m6A_check(well,df)

        outfile = f"{self.sample}/{well}{self.loci_suffix}"
        df.to_csv(outfile, index=False)
        return df
    
    def m6A_check(self,well,df):
        check_dict = {"1":0,"2":0,"3":0,"filtered":0}
        df["check"] = 0
        for index,row in df.iterrows():
            DRACH = row["DRACH"]
            if "A" in DRACH[0:3] and DRACH[3] == 'C':
                df.loc[index,"check"] = 1
                if DRACH[2] == "A":
                    check_dict["1"] += 1
                elif DRACH[1] == "A":
                    check_dict["2"] += 1
                elif DRACH[0] == "A":
                    check_dict["3"] += 1                
            else:
                check_dict["filtered"] += 1

        df = df[ df["check"]==1 ]
        df = df.drop('check', axis=1)
        self.check_dict[well] = check_dict
        return df

    def DRACH_gene_counts(self,well,df):
        DRACH_dict = defaultdict(int)
        gene_dict = defaultdict(int)
        for _,row in df.iterrows():
            DRACH_dict[row["DRACH"]] += 1
            gene_dict[row['gene_name']] += 1

        DRACH_df = pd.DataFrame(list(DRACH_dict.items()),columns=['DRACH','counts'])
        DRACH_df = DRACH_df.sort_values(by='counts',ascending=False)
        outfile = f"{self.sample}/{well}{self.DRACH_suffix}"
        DRACH_df.to_csv(outfile, index=False)

        gene_df = pd.DataFrame(list(gene_dict.items()),columns=['gene','counts'])
        gene_df = gene_df.sort_values(by='counts',ascending=False)
        outfile = f"{self.sample}/{well}{self.gene_suffix}"
        gene_df.to_csv(outfile, index=False)

if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",required=True)
    parser.add_argument('--well_loci',required=True)
    parser.add_argument("--min_m6A_depth",default=2, type=int,
                        help='Minimum depth to call a m6A site')
    parser.add_argument("--min_coverage",default=10, type=int,
                        help='Minimum coverage  reads of a m6A site')
    parser.add_argument("--min_m6A_threshold", type=float, default=0.1,
                        help='m6A threshold filter, greater than min_m6A_threshold will be recognized as m6A site')
    parser.add_argument("--max_m6A_threshold", type=float, default=0.95,
                        help='m6A threshold filter, less than min_m6A_threshold will be recognized as m6A')
    parser.add_argument('--check',action='store_true')

    args = parser.parse_args()
    
    runner = Conversion_summary(args)
    runner.run()