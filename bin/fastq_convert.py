#!/usr/bin/env python

import argparse
import sys
import os
import gzip
import subprocess

import pysam
import utils

BASE_DICT = {
    "A":"T","T":"A","C":"G","G":"C",
    "a":"t","t":"a","c":"g","g":"c"
}

logger = utils.get_logger(__name__)

class Fastq_Convert:
    def __init__(self, args):
        self.args = args
        #output
        os.makedirs(f"{args.sample}", exist_ok=True)
        self.out_fq1 = f"{args.sample}/{args.sample}_bulk_m6A_R1.fastq.gz"
        # reverse complement fastq
        self.out_fq2RC = f"{args.sample}/{args.sample}_bulk_m6A_RC_R2.fastq.gz"
        
        self.fq1_list = args.fq1.split(",")
        self.fq2_list = args.fq2.split(",")

        self.fq1_number = len(self.fq1_list)
        fq2_number = len(self.fq2_list)
        if self.fq1_number != fq2_number:
            sys.exit("fastq1 and fastq2 do not have same file number!")
    
    def run(self):
        logger.info(f"reverse complement fastq")
        fq2out = gzip.open(self.out_fq2RC, 'wb')
        if self.fq1_number == 1:
            cmd = f"mv {self.fq1_list[0]} {self.out_fq1} 2>&1"
            self.cmd_run(cmd)
            logger.info(f"convert {self.fq2_list[0]}")
            with pysam.FastxFile(self.fq2_list[0], persist=False) as fq2in:
                for entry2 in fq2in:
                    header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                    seq22 = self.seq_reverse_complement(seq2)
                    qual22 = qual2[::-1]
                    fq2out.write(f"@{header2}\n{seq22}\n+\n{qual22}\n".encode())
            fq2out.close()
        else:
            fq1out = gzip.open(self.out_fq1, 'wb')
            for i in range(self.fq1_number):
                logger.info(f"convert {self.fq1_list[i]}")
                logger.info(f"convert {self.fq2_list[i]}")
                with pysam.FastxFile(self.fq1_list[i], persist=False) as fq1in:
                    for entry1 in fq1in:
                        header1, seq1, qual1 = entry1.name, entry1.sequence, entry1.quality
                        fq1out.write(f"@{header1}\n{seq1}\n+\n{qual1}\n".encode())
                with pysam.FastxFile(self.fq2_list[i], persist=False) as fq2in:
                    for entry2 in fq2in:
                        header2, seq2, qual2 = entry2.name, entry2.sequence, entry2.quality
                        seq22 = self.seq_reverse_complement(seq2)
                        qual22 = qual2[::-1]
                        fq2out.write(f"@{header2}\n{seq22}\n+\n{qual22}\n".encode())
            fq1out.close()
            fq2out.close()
        logger.info(f"reverse complement fastq finish")

    def cmd_run(self,cmd):
        logger.info(cmd)
        subprocess.check_call(cmd, shell=True)
    
    def seq_reverse_complement(self,seq):
        reversed_complement = ''.join(BASE_DICT.get(base, base) for base in reversed(seq))
        return reversed_complement

if __name__ == "__main__":
    """
    Split fastq based on information provided by the user.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--fq1', required=True)
    parser.add_argument('--fq2', required=True)
    args = parser.parse_args()
    
    runner = Fastq_Convert(args)
    runner.run()