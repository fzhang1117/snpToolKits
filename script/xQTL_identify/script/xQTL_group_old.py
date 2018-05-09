## xQTL_group.py ##
## Zhang Fei <zhangfei-123@foxmail.com> ##
## 2018-04-11 ##
## update: 2018-04-16
## merge xQTL together before set snp by trait

"""
The input file is the significant tassel output stat file,
You can get the file by the line below:
	awk -F '\t' '$7 < (threshold) {print $0}' *sorted.txt(all the tassel output file) > <your output file>

this script are used to group group snps as QTLs

"""

import sys
from operator import itemgetter, attrgetter
import time

## The significant tassel file, don't need to be sorted, and have no title line
fl_tassel = sys.argv[1]
## The snps which intrevl less than the wing_length will be merged together 
wing_length = int(sys.argv[2])
## Output Prefix, could add path but no '/' in the end
index_export = sys.argv[3]
## annotation, gff file
fl_annotation = sys.argv[4]

with open(fl_annotation, 'r') as fh_annotation:
	anno = []
	for line in fh_annotation:
		line = line.strip('\n').split('\t')
		gene = (line[8].split(";"))[0].split("=")[1]
		chrom, start, end = int(line[0][3:]), int(line[3]), int(line[4])
		gene_info = [gene, chrom, start, end]
		anno.append(gene_info)
	anno = sorted(anno, key = itemgetter(1, 2))

with open(fl_tassel, 'r') as fh_tassel:
	sigsnp = fh_tassel.read()

## give a wing to each snp
def my_SettleWing(sigsnp):
	wing = []
	sigsnp = sigsnp.strip('\n')
	sigsnp = sigsnp.split('\n')
	for line in sigsnp:
		line = line.split('\t')
		trait, snp, p, markerR2, snp_chr, snp_location = line[0], line[1], line[6], line[14], int(line[2]), int(line[3])
		wing_start = int(line[3]) - wing_length
		if wing_start <= 0:
			wing_start = 0
		wing_end = int(line[3]) + wing_length
		## the entry of the output should contain the information below
		## [trait, snp, p, markerR2, snp_chr, snp_location, wing_start, wing_end]
		entry = [trait, snp, p, markerR2, snp_chr, snp_location, wing_start, wing_end]
		wing.append(entry)
	wing = sorted(wing, key = itemgetter(4, 5))
	return wing

## wing merge and QTL Summary
def my_xQTLSummary(wing, anno):
	xQTL = []
	entry = []
	i = 0
	for line in wing:
		wing_chr, wing_start, wing_end = line[4], line[6], line[7]
		if i == 0:
			chr_temp, start_temp, end_temp = wing_chr, wing_start, wing_end
			i += 1
			nameQTL = 'xQTL' + str(i)
			entry = [nameQTL, chr_temp, start_temp, end_temp]
		else:
			if wing_chr == chr_temp and wing_start <= end_temp:
				end_temp = wing_end
				entry = [nameQTL, chr_temp, start_temp, end_temp]
			else:
				entry = [nameQTL, chr_temp, start_temp, end_temp]
				xQTL.append(entry)
				chr_temp, start_temp, end_temp = wing_chr, wing_start, wing_end
				i += 1
				nameQTL = 'xQTL' + str(i)
	entry = [nameQTL, chr_temp, start_temp, end_temp]
	xQTL.append(entry)
	
	xQTL_expand = []
	for line in xQTL:
		gene_list = []
		nameQTL, xQTL_chr, xQTL_start, xQTL_end = line[0], line[1], line[2], line[3]
		for line_anno in anno:
			gene, gene_chr, gene_start, gene_end = line_anno[0], line_anno[1], line_anno[2], line_anno[3]
			if gene_chr == xQTL_chr:
				if gene_start >= xQTL_end or gene_end <= xQTL_start:
					pass
				else:
					gene_list.append(gene)
			else:
				continue
		xQTL_expand.append(line + ['|'.join(gene_list)])
	return xQTL_expand

## put snp into xQTL region
def my_snpSettle(xQTL, wing):
	## the aimed data format of dic_xQTL 
	## {xQTL:{trait:[nameQTL, chrQTL, startQTL, endQTL, lengthQTL, snp_number, lead_snp, lead_p, gene_list, all_snp]}}
	## this step give the information about 
	dic_xQTL, dic_xQTL2, wing_expand= {}, {}, []
	for line in xQTL:
		dic_xQTL[line[0]] = {}
		dic_xQTL2[line[0]] = line[4]

	for line_wing in wing:
		trait, snp, p, chr_wing, location_wing = line_wing[0], line_wing[1], float(line_wing[2]), line_wing[4], line_wing[5]
		for line_xQTL in xQTL:
			name_xQTL, chr_xQTL, start_xQTL, end_xQTL = line_xQTL[0], line_xQTL[1], line_xQTL[2], line_xQTL[3]
			if chr_wing == chr_xQTL and location_wing <= end_xQTL and location_wing >= start_xQTL:
				trait_wing = trait
				entry = [trait_wing, snp, p, name_xQTL, chr_xQTL, start_xQTL, end_xQTL]
				wing_expand.append(entry)

	for key in dic_xQTL.keys():
		dic_entry= {}
		trait_list = []
		gene_list = dic_xQTL2[key]
		for line_wing in wing_expand:
			name_xQTL, chr_xQTL, start_xQTL, end_xQTL, trait, snp, p = line_wing[3], line_wing[4], line_wing[5], line_wing[6], line_wing[0], line_wing[1], line_wing[2]
			if name_xQTL == key and trait not in trait_list:
				dic_entry[trait] = [name_xQTL, chr_xQTL, start_xQTL, end_xQTL, end_xQTL - start_xQTL, 1, snp, p, gene_list, [snp]]
				dic_xQTL[name_xQTL] = dic_entry[trait]
				trait_list.append(trait)
				p_min = p
				snp_temp, p_temp = [snp, '-'], [p_min, 1]
			elif name_xQTL == key and trait in trait_list:
				dic_entry[trait][5] += 1
				dic_entry[trait][9].append(snp)
				snp_temp[1], p_temp[1] = snp, p
				p_min = min(p_temp)
				index_leadsnp = p_temp.index(p_min)
				lead_snp = snp_temp[index_leadsnp]
				dic_entry[trait][6], dic_entry[trait][7] = lead_snp, p_min
				snp_temp, p_temp = [snp, '-'], [p_min, 1]
		for value in dic_entry.values():
			value[9] = '|'.join(value[9])
		dic_xQTL[key] = dic_entry
	return dic_xQTL

### the codes below are main process

wing = my_SettleWing(sigsnp)
xQTL = my_xQTLSummary(wing, anno)
summary_xQTL = my_snpSettle(xQTL, wing)

fl_xQTL = index_export + '_xQTL_index.txt'
fl_summary = index_export + '_xQTL_summary.txt'

with open(fl_xQTL, 'a') as fh_xQTL:
	title = ['trait', 'xQTL', 'chr', 'start', 'end', 'gene']
	fh_xQTL.writelines('\t'.join(title))
	fh_xQTL.writelines('\n')
	for line in xQTL:
		## change the ouput data type
		line = [str(i) for i in line]
		fh_xQTL.writelines('\t'.join(line))
		fh_xQTL.writelines('\n')

with open(fl_summary, 'a') as fh_summary:
	title = ['trait', 'xQTL', 'chr', 'start', 'end', 'length', 'snp_number', 'leadsnp','leadp', 'gene_list', 'allsnp']
	fh_summary.writelines('\t'.join(title))
	fh_summary.writelines('\n')
	for key1 in sorted(summary_xQTL.keys()):
		for key2 in sorted(summary_xQTL[key1].keys()):
			value = [str(i) for i in summary_xQTL[key1][key2]]
			value = [key2] + value
			fh_summary.writelines('\t'.join(value))
			fh_summary.writelines('\n')
