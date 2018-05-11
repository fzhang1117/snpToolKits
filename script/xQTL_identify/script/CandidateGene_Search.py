## CandidateGene_Search.py ##
## zhangfei <zhangfei-123@foxmail.com> ##
## 2018-05-10 ##

import sys, re
from itertools import islice

## fl_xQTL_final: the xQTL_final info after single snp and LD filter
## fl_anno: GFF3 annotation file

## output:
## trait, condition, xQTL, chr, start, end, lead_snp, leadp, gene

fl_xQTL_final = sys.argv[1]
fl_anno = sys.argv[2]
flanking_length = int(sys.argv[3])
fl_out = sys.argv[4]

## build gene dictionary
## {chrom:{gene:[gene, chrom, gene_start, gene_end]}}
with open(fl_anno, 'r') as fh_anno:
	dic_anno = {}
	for line in fh_anno:
		line = line.strip('\n').split('\t')
		chrom, gene_start, gene_end = line[0][3:], int(line[3]), int(line[4])
		gene = re.split(';|=', line[8])[1]
		if dic_anno.get(chrom) == None:
			dic_anno[chrom] = {}
			dic_anno[chrom][gene] = [gene, chrom, gene_start, gene_end]
		else:
			dic_anno[chrom][gene] = [gene, chrom, gene_start, gene_end]

with open(fl_out, 'a') as fh_out:
	with open(fl_xQTL_final, 'r') as fh_xQTL_final:
#		allgenes_flanking = []
		for line in islice(fh_xQTL_final, 1, None):
			line = line.strip('\n').split('\t')
			trait, condition, xQTL, xQTL_chr, lead_snp, xQTL_start, xQTL_end, leadp, location = line[0], line[1], line[2], line[3], line[8], int(line[4]), int(line[5]), line[9], int(line[12])
			if dic_anno.get(xQTL_chr) is not None:
				dic_xQTL_chr = dic_anno[xQTL_chr]
				bound_left = location - flanking_length
				bound_right = location + flanking_length
				for gene in dic_xQTL_chr.keys():
					gene_start, gene_end = dic_xQTL_chr[gene][2], dic_xQTL_chr[gene][3]
					if gene_start >= bound_left and gene_end <= bound_right:
						entry = [trait, condition, xQTL, xQTL_chr, xQTL_start, xQTL_end, lead_snp, leadp, gene]
#						allgenes_flanking.append(entry)
						output = [str(i) for i in output]
						fh_out.writelines('\t'.join(output))
						fh_out.writelines('\n')

