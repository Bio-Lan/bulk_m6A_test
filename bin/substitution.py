#!/bin/env python

import argparse
import os
import re

import pandas as pd
import pysam
import utils

sorted_header=['A_to_C','A_to_G','A_to_T','C_to_A','C_to_G','C_to_T','G_to_A','G_to_C','G_to_T','T_to_A','T_to_C','T_to_G']

class Substitution:
    """
    ## Features
    - Computes the overall conversion rates in reads and plots a barplot.

    ## Output
    - `{sample}.substitution.txt` Tab-separated table of the overall conversion rates.
    """

    def __init__(self, args):
        self.args = args
        # input
        self.sample = args.sample
        bams = args.conv_bam.split(",")
        
        conv_sample = pd.read_csv(args.conv_csv, index_col=0)
        well_list = conv_sample.index.to_list()

        self.well_bam = {}
        for well in well_list:
            temp = [j for j in bams if well in j]
            if len(temp) > 1:
                sys.exit(f"ERROR:input file wrong")
            self.well_bam[well] = temp[0]
        
        self.well_dict = {}
        # output files
        self.outstat = f"{self.sample}.substitution.csv"

    def run(self):
        for i in self.well_bam.keys():
            for_base, rev_base, is_forward, is_reverse = self.get_sub_tag(self.well_bam[i])
            xout = self.sub_stat(for_base, rev_base, is_forward, is_reverse)
            self.well_dict[i] = xout
        
        df = pd.DataFrame.from_dict(self.well_dict, orient='index', dtype=float)
        df = df[sorted_header]
        df = df.sort_values(by="C_to_T", ascending=False)
        df.to_csv(self.outstat)

    def get_sub_tag(self,bam):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, "rb")
        pysam.set_verbosity(save)
        is_reverse = {
            "cA": 0,"gA": 0,"tA": 0,
            "aC": 0,"gC": 0,"tC": 0,
            "aG": 0,"cG": 0,"tG": 0,
            "aT": 0,"cT": 0,"gT": 0}
        is_forward = {
            "cA": 0,"gA": 0,"tA": 0,
            "aC": 0,"gC": 0,"tC": 0,
            "aG": 0,"cG": 0,"tG": 0,
            "aT": 0,"cT": 0,"gT": 0}
        for_base = {"a": 0, "c": 0, "g": 0, "t": 0}
        rev_base = {"a": 0, "c": 0, "g": 0, "t": 0}
        snp_tags = ["", "cA", "gA", "tA", "aC", "gC", "tC", "aG", "cG", "tG", "aT", "cT", "gT"]
        ref_tags = ["", "a", "c", "g", "t"]
        for read in bamfile.fetch(until_eof=True):
            snpmatch = re.match(
                r"cA(\d+);gA(\d+);tA(\d+);aC(\d+);gC(\d+);tC(\d+);aG(\d+);cG(\d+);tG(\d+);aT(\d+);cT(\d+);gT(\d+);",
                read.get_tag("SC"),
                re.M)
            totmatch = re.match(r"a(\d+);c(\d+);g(\d+);t(\d+)", read.get_tag("TC"), re.M)
            if snpmatch and totmatch:
                if read.is_reverse:
                    for j in range(1, len(ref_tags)):
                        rev_base[ref_tags[j]] += int(totmatch.group(j))
                    for i in range(1, len(snp_tags)):
                        is_reverse[snp_tags[i]] += int(snpmatch.group(i))
                else:
                    for j in range(1, len(ref_tags)):
                        for_base[ref_tags[j]] += int(totmatch.group(j))
                    for i in range(1, len(snp_tags)):
                        is_forward[snp_tags[i]] += int(snpmatch.group(i))
        bamfile.close()

        return for_base, rev_base, is_forward, is_reverse

    def sub_stat(self, for_base, rev_base, is_forward, is_reverse):
        convertdict = {
            "a": ["aC", "aG", "aT"],
            "c": ["cA", "cG", "cT"],
            "g": ["gA", "gC", "gT"],
            "t": ["tA", "tC", "tG"],
        }
        subdict = {
            "a": "t","t": "a","c": "g","g": "c",
            "aC": "tG","aG": "tC","aT": "tA",
            "cA": "gT","cG": "gC","cT": "gA",
            "gA": "cT","gC": "cG","gT": "cA",
            "tA": "aT","tC": "aG","tG": "aC",
        }
        outdict = {
            "aC": "A_to_C",
            "aG": "A_to_G",
            "aT": "A_to_T",
            "cA": "C_to_A",
            "cG": "C_to_G",
            "cT": "C_to_T",
            "gA": "G_to_A",
            "gC": "G_to_C",
            "gT": "G_to_T",
            "tA": "T_to_A",
            "tC": "T_to_C",
            "tG": "T_to_G",
        }

        outw = {}
        for x in ["a", "c", "g", "t"]:
            fbase = for_base[x]
            rbase = rev_base[subdict[x]]
            for y in convertdict[x]:
                fcov = is_forward[y] * 100 / float(fbase) if float(fbase) > 0 else 0
                rcov = is_reverse[subdict[y]] * 100 / float(rbase) if float(rbase) > 0 else 0
                outw[outdict[y]] = "%.3f" % (fcov + rcov)
        
        return outw

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--sample", help="Sample name.", required=True)
    parser.add_argument("--conv_csv", required=True)
    parser.add_argument("--conv_bam", help="Required. bam file from step conversion.", required=True)
    args = parser.parse_args()
    
    runner = Substitution(args)
    runner.run()
