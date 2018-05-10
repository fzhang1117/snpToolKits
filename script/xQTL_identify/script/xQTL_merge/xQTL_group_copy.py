## xQTL_group.py ##
## Zhang Fei <zhangfei-123@foxmail.com> ##
## 2018-04-11 ##
## update: 2018-04-16
## merge xQTL together before set snp by trait

"""
The input file is the significant tassel output stat file,
You can get the file by the line below:
	awk -F '\t' '$7 < (threshold) {print $0}' *sorted.txt(all the tassel output file) > <your output file>

the input file has no title line and the data format should look like below:

Trait	Marker	Chr	Pos	df	F	p	add_effect	add_F	add_p	dom_effect	dom_F	dom_p	errordf	MarkerR2	Genetic_Var	Residual_Var	-2LnLikelihood
GC_001	chr1.S_273786061	1	273786061	1	27.98496	2.2093E-7	NaN	NaN	NaN	NaN	NaN	NaN	341	0.0825	0.27495	0.81302	1055.41494

xQTL_group.py is the first script of the project xQTL_identify

this script grouped significant associated SNPs by physical distance and give a summary of candidate xQTL

we do this job by these steps:

1. Give wings to each SNPs, assume we give wings as 5kb and the mid-step file will look like the case below.

Trait	p	chr	pos	start	end
GC_001	2.2093E-7	1	273786061	273781051	273791061
GC_001	1.8492E-6	1	46626975	46621975	46627975	

2. merge the start and end region as a 'xQTL'

Trait	xQTL	chr	start	end
GC_001	xQTL1	1	10000	10070
GC_001	xQTL2	2	1000	2000
GC_002	xQTL3	1	10000	10070

3. put all snps into the clusters and summary the cluster infomation

Trait SNP chr_snp	location_snp	p_snp	xQTL	chr_QTL	start_QTL	end_QTL

and summary the result
Trait	xQTL	chr	start	end	SNP number	lead SNP	p-value SNPs
GC_001	xQTL1	2	1	10000	10070	chr1.12312	0.0000001	chr1.12312|chr1.15000

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
	wing = sorted(wing, key = itemgetter(0, 4, 5))
	return wing

## wing merge and QTL Summary
def my_xQTLSummary(wing):
	xQTL = []
	entry = []
	i = 0
	for line in wing:
		trait, wing_chr, wing_start, wing_end = line[0], line[4], line[6], line[7]
		if i == 0:
			trait_temp, chr_temp, start_temp, end_temp = trait, wing_chr, wing_start, wing_end
			i += 1
			nameQTL = 'xQTL' + str(i)
			entry = [trait_temp, nameQTL, chr_temp, start_temp, end_temp]
		else:
			if trait == trait_temp and wing_chr == chr_temp and wing_start <= end_temp:
				end_temp = wing_end
				entry = [trait_temp, nameQTL, chr_temp, start_temp, end_temp]
			else:
				entry = [trait_temp, nameQTL, chr_temp, start_temp, end_temp]
				xQTL.append(entry)
				trait_temp, chr_temp, start_temp, end_temp = trait, wing_chr, wing_start, wing_end
				i += 1
				nameQTL = 'xQTL' + str(i)
	entry = [trait_temp, nameQTL, chr_temp, start_temp, end_temp]
	xQTL.append(entry)
	return xQTL

## put snp into xQTL region
def my_snpSettle(xQTL, wing, annotation):
	
	## the aimed data format of dic_xQTL 
	## {xQTL1: [trait, nameQTL, chr, start, end, length, snp_number, lead_snp, lead_p, gene_list, all_snp]}
	dic_xQTL, i = {}, 1
	#xQTL = sorted(xQTL, key = itemgetter(2, 3))
	for line in xQTL:
		dic_xQTL[i] = line + [line[4] - line[3], 0, '', 0, '', '']
		i += 1
	#values = dic_xQTL.values()
	#values = sorted(values, key = itemgetter(2, 3))

	for key in dic_xQTL.keys():
		value = dic_xQTL[key]
		trait_xQTL, chr_xQTL, start_xQTL, end_xQTL = value[0], value[2], value[3], value[4]
		## initialization....
		i, snp_temp, p_temp, gene_list = 1, [], [], []
		for line_wing in wing:
			trait_wing, snp_name, p, markerR2, chr_wing, location_wing = line_wing[0], line_wing[1], float(line_wing[2]), float(line_wing[3]), line_wing[4], line_wing[5]
			if trait_wing == trait_xQTL:
				if chr_wing == chr_xQTL and (location_wing >= start_xQTL and location_wing <= end_xQTL):
					snp_number = i
					snp_temp.append(snp_name)
					p_temp.append(p)
					i += 1
			else:
				continue
		for line_anno in annotation:
			gene, chr_anno, start_anno, end_anno = line_anno[0], line_anno[1], line_anno[2], line_anno[3]
			if chr_anno == chr_xQTL:
				if start_anno >= end_xQTL or end_anno <= start_xQTL:
					pass
				else:
					gene_list.append(gene)
			else:
				continue

		lead_p = min(p_temp)
		## how about we have two equal lead_p
		index_leadsnp = p_temp.index(lead_p)
		lead_snp = snp_temp[index_leadsnp]
		value[6], value[7], value[8], value[9], value[10] = snp_number, lead_snp, lead_p, '|'.join(gene_list), '|'.join(snp_temp)
		dic_xQTL[key] = value
	return dic_xQTL

### the codes below are main process

wing = my_SettleWing(sigsnp)
xQTL = my_xQTLSummary(wing)
summary_xQTL = my_snpSettle(xQTL, wing, anno)

fl_xQTL = index_export + '_xQTL_index.txt'
fl_summary = index_export + '_xQTL_summary.txt'

with open(fl_xQTL, 'a') as fh_xQTL:
	title = ['trait', 'xQTL', 'chr', 'start', 'end']
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
	for key in sorted(summary_xQTL.keys()):
		value = [str(i) for i in summary_xQTL[key]]
		fh_summary.writelines('\t'.join(value))
		fh_summary.writelines('\n')
