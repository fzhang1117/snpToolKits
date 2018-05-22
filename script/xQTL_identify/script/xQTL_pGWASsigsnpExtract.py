## xQTL_pGWASsigsnpExtract.py ##
## extract all snps in xQTL regions from pGWAS tassel output file, determine the gene location ##
## Zhang Fei <zhangfei-123@foxmail.com> ##
## 2018-05-20 ##

import sys, re
from itertools import islice
from operator import itemgetter, attrgetter

## pGWAS tassel file, could be all SNPs or sigsnps, depends on 
fl_pGWAS = sys.argv[1]

## xQTL_summary file, could be xQTL_summary or xQTL_final_summary file
fl_summary = sys.argv[2]

## flanking region 
flanking = int(sys.argv[3])

## pGWAS_threshold
pGWAS_threshold = float(sys.argv[4])

## gff3 gene only
fl_gff3 = sys.argv[5]

## prefix of output
prefix = sys.argv[6]
fl_generesult = prefix + '_pGWAS_result.txt'

## output file 
#fl_out = sys.argv[5]

with open(fl_gff3, 'r') as fh_gff3:
    dic_gff = {}
    for line in fh_gff3:
        line = line.strip('\n').split('\t')
        chr_gene, start_gene, end_gene = line[0][3:], int(line[3]), int(line[4])
        gene = re.split('=|;', line[8])[1]
        if dic_gff.get(gene) is None:
            dic_gff[gene] = [chr_gene, start_gene, end_gene]

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
        if line[0] == 'Trait':
            pass
        else:
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

filtered_snps = []

for snp in selected_snps:
    p = float(snp[6])
    if p <= pGWAS_threshold:
        filtered_snps.append(snp)

filtered_snps = sorted(filtered_snps, key = itemgetter(0, 2, 3))
gene_mapping = [['triat', 'snp', 'chr', 'location', 'p', 'line', 'MarkerR2', 'Gene']]
#print dic_gff

for line in filtered_snps:
    snp_chr, snp_location = line[2], int(line[3])
    for key in dic_gff.keys():
        value = dic_gff[key]
        if snp_chr == value[0] and snp_location <= value[2] and snp_location >= value[1]:
            result = [line[0], line[1], line[2], line[3], line[6], line[13], line[14]] + [key]
            gene_mapping.append(result)
            break


#print filtered_snps
print 'sigsnps less than', pGWAS_threshold, len(filtered_snps)

with open(fl_generesult, 'w') as fh_generesult:
    gene_mapping = ['\t'.join(i) for i in gene_mapping]
    fh_generesult.write('\n'.join(gene_mapping))
