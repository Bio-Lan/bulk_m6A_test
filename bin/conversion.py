#!/usr/bin/env python
import argparse
import os

import pysam
import pandas as pd

import utils
import filter_gtf
import func_conversion

BASE_DICT = {
    "A":"T","T":"A","C":"G","G":"C",
    "a":"t","t":"a","c":"g","g":"c"
}

class Conversion:
    """
    Get conversion for each read and add tags:
        - get a list of sites of the specified conversion type
        - get statistics for all conversion types   
    """
    def __init__(self,args):
        #input
        self.args = args
        gp = filter_gtf.GtfParser(self.args.GTF)
        self.strand = gp.get_id_strand()

        #output
        os.makedirs(args.sample, exist_ok=True)
        wellBC = os.path.splitext(os.path.basename(args.wellBAM))[0]
        self.output_prefix = args.sample + "/"+ wellBC
        self.outbam = self.output_prefix + '.PosTag.bam'
        self.outloci = self.output_prefix + ".conversion_loci.csv"
    
    def run(self):
        loci_dict = func_conversion.addTags(self.strand, self.outbam, self.args)

        # bam index
        utils.index_bam(self.outbam)
        if len(loci_dict) == 0:
            loci_df = pd.DataFrame(columns=['chrom', 'loci', 'read_strand', 'gene_strand', 'type', 'convs','covers','DRACH', 'gene_name'])
        else:
            loci_df = self.loci_count(loci_dict)
        loci_df.to_csv(self.outloci,index=False)

    def loci_count(self,loci_dict):
        save = pysam.set_verbosity(0)
        bam = pysam.AlignmentFile(self.outbam, 'rb')
        pysam.set_verbosity(save)

        ref_fa = pysam.FastaFile(self.args.FASTA)
        for key in loci_dict.keys():
            loci_dict[key]['gene_name'] = ';'.join([str(i) for i in loci_dict[key]['gene_name']])

            chrom = str(loci_dict[key]["chrom"])
            loci = loci_dict[key]["loci"]
            # counts read with loci
            loci_dict[key]["covers"] = bam.count(chrom,loci,loci+1)
           
            # DRACH sequence
            ref_seq = ref_fa.fetch(chrom,loci-3,loci+2)
            if loci_dict[key]["read_strand"] == "+":
                DRACH = ref_seq
            else:
                DRACH = self.seq_reverse_complement(ref_seq)

            loci_dict[key]["DRACH"] = DRACH

        bam.close()
        ref_fa.close()
        df = pd.DataFrame.from_dict(loci_dict, orient='index')
        return df
    
    def seq_reverse_complement(self,seq):
        reversed_complement = ''.join(BASE_DICT.get(base, base) for base in reversed(seq))
        return reversed_complement


if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--wellBAM", required=True,
                        help='STARsolo split BAM(sortedByCoord),CB in tag must only one ,must have "MD" tag')
    parser.add_argument("--GTF", required=True)
    parser.add_argument("--FASTA",required=True,
                        help='Fasta file path')
    parser.add_argument("--conversion_type", type=str,default="CT",
                        help='conversion type, CT for bulk_m6A', required=False)
    parser.add_argument("--basequalilty", default=20, type=int,
                        help='min base quality of the read sequence', required=False)
    parser.add_argument("--MAPQ",default=0, type=int,
                        help='read MAPQ threshold')

    # add version
    parser.add_argument("--version", action="version", version="1.0")

    args = parser.parse_args()

    runner = Conversion(args)
    runner.run()







