## xQTL_pGWASsigsnpExtract.py ##
## extract all snps in xQTL regions from pGWAS tassel output file, determine the gene location ##
## Zhang Fei <zhangfei-123@foxmail.com> ##
## 2018-05-20 ##

import sys
from itertools import islice

## pGWAS tassel file, could be all SNPs or sigsnps, depends on 
fl_pGWAS = sys.argv[1]

## xQTL_summary file, could be xQTL_summary or xQTL_final_summary file
fl_summary = sys.argv[2]

## flanking region 
flanking = int(sys.argv[3])

with open(fl_summary, 'r') as fh_summary:
    dic_xQTLsummary = {}
    for line in islice(fh_summary, 1, None):
        line = line.strip('\n').split('\t')
        xQTL, xQTL_chr, xQTL_location = line[2], int(line[3]), int(line[12])
        ## {xQTL: [xQTL_chr, xQTL_location, xQTL_start, xQTL_end]}
        if dic_xQTLsummary.get(xQTL) is None:
            dic_xQTLsummary[xQTL] = [xQTL_chr, xQTL_location, xQTL_location - flanking, xQTL_location + flanking]

with open(fl_pGWAS, 'r') as fh_pGWAS:
    all_snps = []
    for line in islice(fh_pGWAS, 1, None):
        line = line.strip('\n').split('\t')
        all_snps.append(line)

selected_snps = []

for key in dic_xQTLsummary.keys():
    value = dic_xQTLsummary[key]
    xQTL_chr, xQTL_left, xQTL_right = value[0], value[2], value[3]
    for entry in all_snps:
        if entry[3] != '':
            snp_chr, snp_location = int(entry[2]), int(entry[3])
            if snp_chr == xQTL_chr and snp_location <= xQTL_right and snp_location >= xQTL_left:
                selected_snps.append(entry)

print selected_snps


