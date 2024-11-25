#!/bin/env python

import argparse
import gzip
import os
import sys
from multiprocessing import Pool

import pandas as pd
import pysam
import scipy.io
import scipy.sparse

import utils
from __init__ import ASSAY, BARCODE_FILE_NAME, FEATURE_FILE_NAME, MATRIX_FILE_NAME

BULK_M6A_MATRIX_DIR_SUFFIX = ["raw", "labeled", "unlabeled"]

class Quant:
    """
    Features
    - Quantify total RNA, unlabeled and labeled RNA.
    Output
    - `{sample}_{RNA}` The expression matrix of total/labeled/unlabeled in Matrix Market Exchange Formats.
    """
    def __init__(self, args):
        self.args = args
        #input
        self.sample = args.sample
        self.matrix_dir = args.matrix_dir

        #set
        barcodes_file = os.path.join(self.matrix_dir,BULK_M6A_MATRIX_DIR_SUFFIX[0], BARCODE_FILE_NAME)
        features_file = os.path.join(self.matrix_dir,BULK_M6A_MATRIX_DIR_SUFFIX[0], FEATURE_FILE_NAME)
        self.features = pd.read_csv(features_file, sep="\t", header=None, index_col=0)
        self.barcodes = utils.read_one_col(barcodes_file)
        self.totaldf = pd.DataFrame()

        self.well_dict = {}
        bam_file = args.conv_bam.split(",")
        loci_file = args.conv_loci.split(",")
        df = pd.read_csv(args.conv_sample, index_col=0)
        self.well_list = df.index.to_list()
        for well in self.well_list:
            temp1 = [x for x in bam_file if well in x]
            temp2 = [y for y in loci_file if well in y]
            if len(temp1) > 1 or len(temp2) > 1:
                sys.exit(f"ERROR:input file wrong")
            self.well_dict[well] = {"bam":temp1[0], "m6A":temp2[0]}

        #output
        self.dir_labeled = f"{self.matrix_dir}/{BULK_M6A_MATRIX_DIR_SUFFIX[1]}"
        self.dir_unlabeled = f"{self.matrix_dir}/{BULK_M6A_MATRIX_DIR_SUFFIX[2]}"
        self.detail_txt = f"{self.sample}.labeled_detail.txt"
        self.rawcsv = f'{self.sample}_raw.csv'
        self.fltcsv = f'{self.sample}_filtered.csv'
    
    def run(self):
        # get backgroud m6a and snp
        self.get_m6A_sie()

        well_dfs = self.run_quant()
        for i in well_dfs:
            self.totaldf = pd.concat([self.totaldf,i], axis=0)
        self.totaldf.to_csv(self.detail_txt, sep="\t", index=False)
        
        self.get_labeled_unlabeled()
        self.write_csv_file()
    
    def get_m6A_sie(self):
        m6A_snp_dict = {}
        if self.args.m6A_file or self.args.snp_file:
            m6A_snp_dict = self.additional_file()
        for well in self.well_dict.keys():
            m6A_site = self.get_loci_list(self.well_dict[well]["m6A"])
            if self.args.well_file:
                if well in m6A_snp_dict.keys():
                    add_site = m6A_snp_dict[well]
            else:
                add_site = m6A_snp_dict
                if "add_m6A" in add_site.keys():
                    m6A_site += add_site["add_m6A"]
                if "add_snp" in add_site.keys():
                    m6A_site = [loci for loci in m6A_site if loci not in add_site["add_snp"]]
            self.well_dict[well]["m6A"] = m6A_site

    def additional_file(self):
        m6A_snp = {}
        if self.args.well_file:
            if self.args.m6A_file:
                m6A_file = self.args.m6A_file
                self.file_check(m6A_file)
                m6A_snp = self.get_each_well_file(m6A_file,"add_m6A",m6A_snp)
            if self.args.snp_file:
                snp_file = self.args.snp_file
                self.file_check(snp_file)
                m6A_snp = self.get_each_well_file(snp_file,"add_snp",m6A_snp)
        else:
            if self.args.m6A_file:
                self.file_check(self.args.m6A_file)
                add_m6A = self.get_loci_list(self.args.m6A_file)
                m6A_snp["add_m6A"] = add_m6A
            if self.args.snp_file:
                self.file_check(self.args.snp_file)
                add_snp = self.get_loci_list(self.args.snp_file)
                m6A_snp["add_snp"] = add_snp
        return m6A_snp
    
    def get_each_well_file(self,file,label,dict1):
        with open(file) as f:
            for i in f:
                temp = i.strip().split(',')
                self.well_file_check(temp,m6A_file)
                site_list = self.get_loci_list(temp[1])
                if temp[0] not in dict1.keys():
                    dict1[ temp[0] ] = {}
                    dict1[ temp[0] ][label] = site_list
                else:
                    if label in dict1[ temp[0] ].keys():
                        sys.exit(f"[ERROR] same well [{temp[0]}] in {file}.")
        return dict1

    def file_check(self,file):
        if not os.path.exists(file):
            sys.exit(f"[Error] {file}  not exists, please check the file.\n")
        if not file.endswith('.csv'):
            sys.exit(f"[Error] {file} format error. Only csv format is allowed.")
    
    def well_file_check(self,well,file_name):
        if len(well) == 1:
            sys.exit(f"[Error] {well[0]} does not have matched file. If not, please modify it in {file_name} !!\n")
        if well[0] not in self.well_list:
            sys.exit(f"[Error] {well[0]} is not in well list, please check the well BC in {file_name} !!\nWell BC list can be found in {self.sample}.conversion.csv")
        self.file_check(well[1])

    def get_loci_list(self,file):
        loci = pd.read_csv(file,header=0,sep=",")
        if 'chrom' in loci.columns and 'loci' in loci.columns:
            loci['chr_loci'] = loci['chrom'].astype(str)+ '_' + loci['loci'].astype(str)
            loci_list = list(loci['chr_loci'])
            return loci_list
        else:
            sys.exit(f"[Error]: Column name in {file} must have 'chrom' and 'loci'.")
    
    def run_quant(self):
        in_bam_list, m6A_list = [], []
        for well in self.well_list:
            in_bam_list.append(f"{self.well_dict[well]['bam']}")
            m6A_list.append(f"{self.well_dict[well]['m6A']}")        
        
        mincpu = min(len(self.well_list),self.args.thread)
        with Pool(mincpu) as pool:
            results = pool.starmap(Quant.quant,zip(in_bam_list,m6A_list)) 
        return results
    
    @staticmethod
    def quant(bam,m6A):
        save = pysam.set_verbosity(0)
        bamfile = pysam.AlignmentFile(bam, "rb")
        pysam.set_verbosity(save)
        countdict = {}

        for read in bamfile.fetch(until_eof=True):
            cb = read.get_tag("CB")
            chro = read.reference_name
            ub = read.get_tag('UB')
            gene = read.get_tag('GX')
            
            if read.get_tag("ST") == "+":
                stag = read.get_tag("CL")
            else:
                stag = read.get_tag("GL")
            if stag == '-':
                m6a_num = 0
            else:
                m6a_pos = []
                for si in range(0, len(stag)):
                    pos = str(chro) + '_' + str(stag[si])
                    if pos in m6A:
                        m6a_pos.append(int(stag[si]))
                m6a_num = len(m6a_pos)
            readid = ":".join([cb, ub, gene])
            if readid not in countdict:
                countdict[readid] = m6a_num
            else:
                if m6a_num > countdict[readid]:
                    countdict[readid] = m6a_num
        bamfile.close()
        ct_df = pd.DataFrame.from_dict(countdict, orient="index", columns=["CT"])
        ct_df = ct_df.reset_index()
        ct_df[['Barcode', 'UMI',"geneID"]] = ct_df['index'].str.split(':', expand=True)
        ct_df = ct_df[['Barcode', 'UMI',"geneID","CT"]]
        
        return ct_df

    def get_labeled_unlabeled(self):
        self.labeled = self.totaldf[self.totaldf["CT"] > 0]
        unlabeled = self.totaldf[self.totaldf["CT"] == 0]
        # matrix
        self.write_sparse_matrix(self.labeled.drop("CT", axis=1), self.dir_labeled)
        self.write_sparse_matrix(unlabeled.drop("CT", axis=1), self.dir_unlabeled)

    def write_sparse_matrix(self,df,matrix_dir):
        count_matrix = self.dataframe_to_matrix(df, features=self.features.index, barcodes=self.barcodes)
        self.to_matrix_dir(count_matrix, matrix_dir)

    def dataframe_to_matrix(self, df, features, barcodes, barcode_column="Barcode", feature_column="geneID", value="UMI"):
        if df.shape[0] > 0:
            series_grouped = df.groupby([barcode_column, feature_column], observed=True).size()
            series_grouped.name = value
            df_grouped = pd.DataFrame(series_grouped)
        else:
            empty_matrix = scipy.sparse.coo_matrix((len(features), len(barcodes)))
            return empty_matrix

        feature_index_dict = {}
        for index, gene_id in enumerate(features):
            feature_index_dict[gene_id] = index
        barcode_index_dict = {}
        for index, barcode in enumerate(barcodes):
            barcode_index_dict[barcode] = index

        # use all barcodes
        barcode_codes = [barcode_index_dict[barcode] for barcode in df_grouped.index.get_level_values(level=0)]
        # use all gene_id from features even if it is not in df
        gene_id_codes = [feature_index_dict[gene_id] for gene_id in df_grouped.index.get_level_values(level=1)]
        mtx = scipy.sparse.coo_matrix(
            (df_grouped[value], (gene_id_codes, barcode_codes)), shape=(len(features), len(barcodes))
        )

        return mtx

    def to_matrix_dir(self, count_matrix, matrix_dir):
        self.check_mkdir(dir_name=matrix_dir)
        self.features.to_csv(f"{matrix_dir}/{FEATURE_FILE_NAME}", sep="\t", header=False)
        pd.Series(self.barcodes).to_csv(f"{matrix_dir}/{BARCODE_FILE_NAME}", index=False, sep="\t", header=False)
        matrix_path = f"{matrix_dir}/{MATRIX_FILE_NAME}"
        with gzip.open(matrix_path, "wb") as f:
            scipy.io.mmwrite(f, count_matrix)

    def check_mkdir(self, dir_name):
        """if dir_name is not exist, make one"""
        if not os.path.exists(dir_name):
            os.system(f"mkdir -p {dir_name}")
    
    def write_csv_file(self):
        df_sum = self.totaldf.groupby('Barcode').agg({
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df_labeled_sum = self.labeled.groupby('Barcode').agg({
            'UMI': 'count',
            'geneID': 'nunique'
        })
        df = pd.concat([df_sum,df_labeled_sum], axis=1)
        df.columns=['UMI','gene','labeled_UMI','labeled_gene']
        df = df.fillna(0)
        df = df.astype(int)
        df.to_csv(self.rawcsv)

        df_filter = df[ (df['UMI']>=self.args.umi_cutoff) & (df['gene']>=self.args.gene_cutoff) ]
        df_filter.to_csv(self.fltcsv)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",required=True)
    parser.add_argument("--matrix_dir",required=True)
    parser.add_argument('--conv_sample',required=True)
    parser.add_argument('--conv_bam',required=True)
    parser.add_argument('--conv_loci',required=True)
    parser.add_argument('--umi_cutoff', default=500, type=int,
        help='If the UMI number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--gene_cutoff', default=0, type=int,
        help='If the gene number exceeds the threshold, it is considered a valid well and reported.'
    )
    parser.add_argument('--thread', type=int)   
    parser.add_argument('--m6A_file', type=str,
                        help="""m6A position file.If this option is set, it is valid for all wells""")
    parser.add_argument('--snp_file', type=str,
                        help="""backgroud snp file.If this option is set, it is valid for all wells""")
    parser.add_argument('--well_file',action='store_true')
    args = parser.parse_args()

    runner = Quant(args)
    runner.run()
