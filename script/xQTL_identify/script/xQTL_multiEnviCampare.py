## xQTL_multiEnvirCompar.py ##
## Zhang Fei <zhangfei-123@foxmail.com> ##
## 2018-05-13 ##

## Given xQTL-trait pair, compare leadp and effect in normal and drought condition 
## the tassel file should be names as: mlm_GC_drought_00012_sorted.txt

import sys, re, os, math
from itertools import islice
#import numpy as np

path_drought = sys.argv[1]
path_normal = sys.argv[2]
## xQT_summary file come from xQTL_group.py
fl_xQTL_summary = sys.argv[3]
fl_out = sys.argv[4]

def dicCpdResult_build(path_cpd, prefix):
    dic_cpd = {}
    fls_cpd = os.listdir(path_cpd)
    fls_cpd = sorted(filter(lambda x: re.search("_sorted", x), fls_cpd))
    for fl in fls_cpd:
        cpd = prefix + fl.split('_')[3][0: 3]
        if dic_cpd.get(cpd) is None:
            dic_cpd[cpd] = [path_cpd + fl]
        else:
            dic_cpd[cpd].append(path_cpd + fl)
    return dic_cpd

def snpExtract(dic_cpd, trait, chr_xQTL, start_xQTL, end_xQTL):
    fls = dic_cpd[trait]
    snp_list = []
    for fl in fls:
        with open(fl, 'r') as fh:
            for line in islice(fh, 2, None):
                line = line.strip('\n').split('\t')
                snp, snp_chr, snp_location = line[1], int(line[2]), int(line[3])
                if snp_chr == int(chr_xQTL) and snp_location >= int(start_xQTL) and snp_location <= int(end_xQTL):
                    p = float(line[6])
                    entry = [snp, p]
                    snp_list.append(entry)
    snp_list.sort(key = lambda x: x[1])
    snp_min = snp_list[0]
    return snp_min

dic_cpd_drought = dicCpdResult_build(path_drought, 'GC_')
dic_cpd_normal = dicCpdResult_build(path_normal, 'GC_')

with open(fl_xQTL_summary, 'r') as fh_xQTL_summary:
    result = [['trait', 'xQTL', 'Location', 'leadsnp_normal', 'leadp_normal', 'leadsnp_drought', 'leadp_drought']]
    for line in islice(fh_xQTL_summary, 1, None):
        line = line.strip('\n').split('\t')
        trait, chr_xQTL, start_xQTL, end_xQTL, leadp, xQTL, leadsnp = line[0], int(line[3]), int(line[4]), int(line[5]), float(line[9]), line[2], line[8]
        xQTL_location = str(chr_xQTL) + ':' + str(start_xQTL) + '-' + str(end_xQTL)
        # things in line: [trait, condition, xQTL, chr, start, end, length, snp_number, leadsnp, leadp, gene_list, allsnp, lead_location]
        if line[1] == 'normal':
            snp_normal, leadp_normal = leadsnp, -math.log(leadp, 10)
            condition = 'drought'
            list_leadsnp_drought = snpExtract(dic_cpd_drought, trait, chr_xQTL, start_xQTL, end_xQTL)
            snp_drought, leadp_drought = list_leadsnp_drought[0], -math.log(list_leadsnp_drought[1], 10)
            entry = [trait, xQTL, xQTL_location, snp_normal, leadp_normal, snp_drought, leadp_drought]
            entry = [str(i) for i in entry]
            print '\t'.join(entry)
            result.append(entry)
        elif line[1] == 'drought':
            snp_drought, leadp_drought = leadsnp, -math.log(leadp, 10)
            condition = 'normal'
            list_leadsnp_normal = snpExtract(dic_cpd_normal, trait, chr_xQTL, start_xQTL, end_xQTL)
            snp_normal, leadp_normal = list_leadsnp_normal[0], -math.log(list_leadsnp_normal[1], 10)
            entry = [trait, xQTL, xQTL_location, snp_normal, leadp_normal, snp_drought, leadp_drought]
            entry = [str(i) for i in entry]
            print '\t'.join(entry)
            result.append(entry)

with open(fl_out, 'w') as fh_out:
    result = ['\t'.join(i) for i in result]
    fh_out.write('\n'.join(result))

