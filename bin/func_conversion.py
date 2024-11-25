###### step Func
import re
import sys
from collections import defaultdict

import pysam
import pandas as pd

rdict = {'A':'T','T':'A','G':'C','C':'G',
         'a':'t','t':'a','g':'c','c':'g'}
sub_types = ['cA','gA','tA','aC','gC','tC','aG','cG','tG','aT','cT','gT','aN','cN','gN','tN']
def getConvTypes(convtype):
    """
    get forward and reverse conversion type
    return type: ref_base + read_base
    """
    ftypes = convtype[0].lower()+convtype[1].upper()
    rtypes = rdict[convtype[0]].lower()+rdict[convtype[1]].upper()
    return ftypes, rtypes

def MD_split(MD):
    """
    split MD tag into MD items
    example:
        MD:"75M0CD3^CD0CD5M"
        return:['75', 'M', '0', 'CD', '3', '^CD', '0', 'CD', '5', 'M']
    """
    temp = [i for i in MD]
    tag = []
    i = 0
    while i < len(temp):
        if temp[i] == '^':
            deleted = '^' + re.match(r'[A-Z]+', MD[i+1:]).group()
            tag.append(deleted)
            i += len(deleted)
        elif temp[i].isalpha():
            diff = re.match(r'[A-Z]+', MD[i:]).group()
            tag.append(diff)
            i += len(diff)
        else:
            num = re.match(r'\d+',MD[i:]).group()
            tag.append(num)
            i += len(num)
    return tag

def diff_base(read):
    """
    parse MD tag,return different base in read(1-based)
    example: 
        MD:'2A5A5A5A11A77'
        return: {3: 'A', 9: 'A', 15: 'A', 21: 'A', 33: 'A'}
    """
    MD_list = MD_split(read.get_tag("MD"))
    base_dict = {}
    base_pos = 0
    for i in MD_list:
        if i.isdigit():
            base_pos += int(i)
        elif '^' in i:
            base_pos += 0
        else:
            for j in i:
                base_dict[base_pos+1] = j
                base_pos += 1
    return base_dict

def createTag(data1,separator1,separator2=''):
    if isinstance(data1, dict):
        return separator1.join([key + separator2 + str(data1[key]) for key in data1.keys()])
    if isinstance(data1, list):
        return separator1.join([str(i) for i in data1])

def ConvInRead(read,ftypes,rtypes,qual,MAPQ):
    # set
    total_content = {'a': 0, 'c': 0, 'g': 0, 't': 0}
    specific_conversions = {}
    for i in sub_types:
        specific_conversions[i] = 0
    conv_loci = defaultdict(list)

    # ref base summary
    refseq = read.get_reference_sequence().lower()
    for base in total_content.keys():
        total_content[base] += refseq.count(base)
    
    base_dict = diff_base(read) # 1-based
    if len(base_dict) == 0 or read.mapping_quality < MAPQ:
        SC_tag = createTag(specific_conversions,';')
        TC_tag = createTag(total_content,';')
        conv_loci[ftypes] = "-"
        conv_loci[rtypes] = "-"
        return SC_tag, TC_tag, conv_loci
    else:
        match_pair = read.get_aligned_pairs(with_seq=True,matches_only=True) # 0-based
        read_seq = read.query_sequence
        read_qual = read.query_qualities
        read_name = read.query_name
        for i in base_dict.keys():
            pair = match_pair[i-1]
            seq_pos, ref_pos, ref_base = pair[0], pair[1], pair[2]
            if base_dict[i].lower() != ref_base:
                sys.exit(f"Wrong in this read, please check: {read_name}")
            
            read_conv = base_dict[i].lower()+read_seq[seq_pos].upper()
            if read_conv in specific_conversions.keys() and read_qual[seq_pos] > qual:
                specific_conversions[read_conv] += 1
                if read_conv == ftypes or read_conv == rtypes:
                    conv_loci[read_conv].append(ref_pos)
        SC_tag = createTag(specific_conversions,';')
        TC_tag = createTag(total_content,';')

        if len(conv_loci[ftypes]) == 0:
            conv_loci[ftypes] = '-' 
        if len(conv_loci[rtypes]) == 0:
            conv_loci[rtypes] = '-'

        return SC_tag, TC_tag, conv_loci

def addTags(strand_dict,outbam,args):
    #set
    inbam = args.wellBAM
    conv_type = args.conversion_type
    qual = args.basequalilty
    MAPQ = args.MAPQ

    conv_loci={}
    ftypes, rtypes = getConvTypes(conv_type)

    save = pysam.set_verbosity(0)
    bamfile = pysam.AlignmentFile(inbam, 'rb')
    header = bamfile.header
    mod_bamfile = pysam.AlignmentFile(outbam, mode='wb', header=header,check_sq=False)
    pysam.set_verbosity(save)

    for read in bamfile.fetch(until_eof=True):
        gene_ID = read.get_tag('GX')
        if (not gene_ID) or gene_ID == '-':
            continue
            
        tags = ConvInRead(read,ftypes,rtypes,qual,MAPQ)
        """
        example return:
            tags[0] SC_tag: str 'cA:0gA:0tA:0aC:0gC:0tC:0aG:0cG:0tG:0aT:0cT:0gT:0aN:0cN:0gN:0tN:0'
            tags[1] TC_tag: str 'a:0c:0g:0t:0'
            tags[2] conv_loc: dict {'cT':[int],'gA':[int]} #0-based
        """
        read_strand = '-' if read.is_reverse else '+'
        gene_strand = strand_dict[gene_ID]
        read_chrom = read.reference_name
        read.set_tag('SC', tags[0], 'Z') #special conversion
        read.set_tag('TC', tags[1], 'Z') # read ref sequence base count
        read.set_tag('ST', gene_strand) # gene strand
        read.set_tag('CL', tags[2][ftypes])
        read.set_tag('GL', tags[2][rtypes])

        # conversion loci in read - tag CL
        read_conv = ftypes if gene_strand== "+" else rtypes # Determine the conversion type based on the direction of gene
        read_loci = tags[2][read_conv]
        
        if read_loci != "-":
            for i in read_loci:
                key1 = f'{read_chrom}_{i}_{gene_strand}'
                if key1 not in conv_loci.keys():
                    loci_dict = {
                        "chrom":read_chrom,
                        "loci":i,
                        "read_strand":read_strand,
                        "gene_strand":gene_strand,
                        "type":read_conv,
                        "convs":1,
                        "covers":0,
                        "DRACH":0,
                        "gene_name":[read.get_tag('GN')]
                    }
                    conv_loci[key1] = loci_dict
                else:
                    conv_loci[key1]["convs"] += 1
                    if read.get_tag('GN') not in conv_loci[key1]["gene_name"]:
                        conv_loci[key1]["gene_name"].append(read.get_tag('GN'))
        mod_bamfile.write(read)
    bamfile.close()
    mod_bamfile.close()
    return conv_loci
