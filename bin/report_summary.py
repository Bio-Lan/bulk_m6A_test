#!/usr/bin/env python

import argparse
from collections import defaultdict

import os
import pandas as pd

import utils
from __init__ import ASSAY

class Report_Summary:
    def __init__(self, args):
        self.sample = args.sample
        self.args = args
        self.stats = {}
    
    def run(self):
        self.parse_read_stats()
        self.parse_summary()
        utils.write_multiqc(self.stats, self.args.sample, ASSAY, "reads.stats")

        # get well BC
        df_raw = pd.read_csv(self.args.sample_raw,header=0,sep=",",index_col= 0)
        df_filter = pd.read_csv(self.args.sample_filter,header=0,sep=",",index_col= 0)
        df = df_filter if df_filter.shape[0] > 0 else df_raw
        self.well_list = df.index.to_list()

        self.parse_conversion_motif()
        self.parse_substitution_rate()
        self.parse_quant_result(df)

    def parse_read_stats(self):
        dtypes = defaultdict(lambda: "int")
        dtypes["CB"] = "object"
        df = pd.read_csv(
            self.args.read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
        )  # skip first line cb not pass whitelist
        df = df.loc[
            :,
            [
                "cbMatch",
                "cbPerfect",
                "genomeU",
                "genomeM",
                "exonic",
                "intronic",
                "exonicAS",
                "intronicAS",
                "countedU",
                "nUMIunique",
                "nGenesUnique",
            ],
        ]
        s = df.sum()
        # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
        valid = int(s["cbMatch"])
        perfect = int(s["cbPerfect"])
        corrected = valid - perfect
        genome_uniq = int(s["genomeU"])
        genome_multi = int(s["genomeM"])
        mapped = genome_uniq + genome_multi
        exonic = int(s["exonic"])
        intronic = int(s["intronic"])
        antisense = int(s["exonicAS"] + s["intronicAS"])
        intergenic = mapped - exonic - intronic - antisense
        counted_uniq = int(s["countedU"])
        data_dict = {
            "Corrected Barcodes": corrected / valid,
            "Reads Mapped To Unique Loci": genome_uniq / valid,
            "Reads Mapped To Multiple Loci": genome_multi / valid,
            "Reads Mapped Uniquely To Transcriptome": counted_uniq / valid,
            "Mapped Reads Assigned To Exonic Regions": exonic / mapped,
            "Mapped Reads Assigned To Intronic Regions": intronic / mapped,
            "Mapped Reads Assigned To Intergenic Regions": intergenic / mapped,
            "Mapped Reads Assigned Antisense To Gene": antisense / mapped,
        }
        for k in data_dict:
            data_dict[k] = utils.get_frac(data_dict[k])
        self.stats.update(data_dict)        

    def parse_summary(self):
        data = utils.csv2dict(self.args.summary)
        origin_new = {
            "Number of Reads": "Raw Reads",
            "Reads With Valid Barcodes": "Valid Reads",
            "Q30 Bases in CB+UMI": "Q30 Bases in CB+UMI",
            "Q30 Bases in RNA read": "Q30 Bases in RNA read"
        }
        parsed_data = {}
        for origin, new in origin_new.items():
            parsed_data[new] = data[origin]
        frac_names = {"Valid Reads","Q30 Bases in CB+UMI","Q30 Bases in RNA read"}
        for k in frac_names:
            parsed_data[k] = utils.get_frac(parsed_data[k])
        for k in set(origin_new.values()) - frac_names:
            parsed_data[k] = int(parsed_data[k])
        self.stats.update(parsed_data)
    
    def parse_quant_result(self,df):
        json_dict = df.to_dict(orient="index")
        utils.write_multiqc(json_dict, self.sample, ASSAY, "quant.well_inf")
        
        stats = df.describe()
        stats.columns = ['UMI','Genes','labeled_UMI','labeled_gene']

        df['UMI_rate'] = df['labeled_UMI']/df['UMI']
        df['gene_rate'] = df['labeled_gene']/df['gene']
        labeled_rate = round( df['UMI_rate'].median() * 100,2)
        labeled_rate2 = round( df['UMI_rate'].mean() * 100,2)

        data_dict = {}
        for item in ['UMI', 'Genes']:
            temp = f'Median {item} across Well'
            data_dict[temp] = int(stats.loc['50%',item])
            temp = f'Mean {item} across Well'
            data_dict[temp] = int(stats.loc['mean',item])
        
        temp = 'Median labeled rate across wells'
        data_dict[temp] = labeled_rate
        temp = 'Mean labeled rate across wells'
        data_dict[temp] = labeled_rate2
        utils.write_multiqc(data_dict, self.sample, ASSAY, "quant.stats")
        
        df = df.loc[:,['UMI_rate','gene_rate']]
        df.columns = ['UMI','Gene']
        label_dict = df.to_dict(orient="list")
        utils.write_multiqc(label_dict, self.sample, ASSAY, "quant.labeled_rate")

    def parse_conversion_motif(self):
        DRACH_dict = {}
        motif_file = self.args.conv_motif.split(",")
        for file in motif_file:
            well = os.path.basename(file).split(".")[0]
            if well in self.well_list:
                df = pd.read_csv(file,header=0,sep=",")
                top_DRACH = df.head(self.args.DRACH_display)
                top_DRACH = top_DRACH.set_index('DRACH')['counts'].to_dict()
                DRACH_dict[well] = top_DRACH
        utils.write_multiqc(DRACH_dict, self.sample, ASSAY, "conversion.motif")
    
    def parse_substitution_rate(self):
        df = pd.read_csv(self.args.sub_stats,header = 0 ,sep=",",index_col=0)
        df = df.loc[self.well_list,]
        df = df.sort_values(by="C_to_T",ascending=False)
        
        box_dict = df.to_dict(orient='list')
        bar_dict = df.to_dict(orient="index")
        utils.write_multiqc(box_dict, self.sample, ASSAY, "substitution.boxplot")
        utils.write_multiqc(bar_dict, self.args.sample, ASSAY, "substitution.barplot")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Result summary")
    parser.add_argument("--sample", help="sample name")
    parser.add_argument("--read_stats", help="cellReadsStats file")
    parser.add_argument("--summary", help="summary file")
    parser.add_argument("--conv_motif", help="motif for each well")
    parser.add_argument("--sub_stats", help="substitution csv file")
    parser.add_argument("--sample_raw", help="qutan sample raw csv file")
    parser.add_argument("--sample_filter", help="qutan sample filter csv file")
    parser.add_argument("--DRACH_display", default=3, type=int, help="how many DRACH display")
    args = parser.parse_args()

    Report_Summary(args).run()