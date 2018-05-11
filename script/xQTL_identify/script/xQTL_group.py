## xQTL_group_new.py ##
## zhang fei <zhangfei-123@foxmail.com > ##
## 2018-05-07 ##

"""
first input all sigsnp in all enviroment
then merge all snps in given wings
put sigsnps in in divided enviroment

the sigsnp can get from tassle file by the command below:
	awk -F '\t' '$7 < (threshold) {print $0}' <given_file> > <output>
"""

import sys
from operator import itemgetter, attrgetter

## parameter input ##

# sigsnp in drought condition
# sigsnp in normal condition
# wing length
# index export
# annotation file

fl_drought = sys.argv[1]
fl_normal = sys.argv[2]
wing_length = int(sys.argv[3])
index_export = sys.argv[4]
fl_annotation = sys.argv[5]

## data load and format ##
with open(fl_annotation , 'r') as fh_annotation:
	anno = []
	for line in fh_annotation:
		line = line.strip('\n').split('\t')
		gene = (line[8].split(";"))[0].split("=")[1]
		chrom, start, end = int(line[0][3:]), int(line[3]), int(line[4])
		gene_info = [gene, chrom, start, end]
		anno.append(gene_info)
	anno = sorted(anno, key = itemgetter(1, 2))

with open(fl_drought, 'r') as fh_drought:
	sigsnp_drought = fh_drought.read()

with open(fl_normal, 'r') as fh_normal:
	sigsnp_normal = fh_normal.read()

sigsnp_all = sigsnp_drought + sigsnp_normal

## give a wing to all sig snp
def my_SettleWing(sigsnp, wing_length):
	wing = []
	sigsnp = sigsnp.strip('\n').split('\n')
	for line in sigsnp:
		line = line.split('\t')
		#trait, snp, p, markerR2, snp_chr, snp_location = line[0], line[1], line[6], line[14], int(line[2]), int(line[3])
		snp_chr, snp_location = int(line[2]), int(line[3])
		wing_start = int(line[3]) - wing_length
		if wing_start <= 0:
			wing_start = 0
		wing_end = int(line[3]) + wing_length
		## output : [snp_chr, wing_start, wing_end]
		entry = [snp_chr, wing_start, wing_end]
		wing.append(entry)
	wing = sorted(wing, key = itemgetter(0, 1))
	return wing

## wing merge and gene in this QTL
def my_xQTLSummary(wing, anno):
	xQTL, entry, i = [], [], 0
	for line in wing:
		wing_chr, wing_start, wing_end = line[0], line[1], line[2]
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
	## put last entry into xQTL
	entry = [nameQTL, chr_temp, start_temp, end_temp]
	xQTL.append(entry)

	xQTL_expand = []
	for line in xQTL:
		gene_list, nameQTL, xQTL_chr, xQTL_start, xQTL_end = [], line[0], line[1], line[2], line[3]
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

## put snp in xQTL
def my_snpSettle(snp, xQTL, condition):
	dic_xQTL, dic_xQTL2, list_snp = {}, {}, []
	snp = snp.strip('\n').split('\n')

	for line in xQTL:
		dic_xQTL[line[0]] = {}
		#dic_xQTL2[line[0]] = line[4]

	for line_xQTL in xQTL:
		xQTL_name, xQTL_chr, xQTL_start, xQTL_end, genes = line_xQTL[0], line_xQTL[1], line_xQTL[2], line_xQTL[3], line_xQTL[4]
		for line2 in snp:
			line2 = line2.split('\t')
			snp_trait, snp_name, snp_chr, snp_location, snp_p = line2[0], line2[1], int(line2[2]), int(line2[3]), float(line2[6])
			if snp_chr == xQTL_chr and (snp_location >= xQTL_start and snp_location <= xQTL_end):
				entry = [snp_trait, snp_name, condition, snp_p, xQTL_name, xQTL_chr, xQTL_start, xQTL_end, genes, snp_location]
				list_snp.append(entry)
	for key in dic_xQTL.keys():
		dic_entry, trait_list = {}, []
		for line_snp in list_snp:
			snp_trait, snp_name, condition, snp_p, xQTL_name, xQTL_chr, xQTL_start, xQTL_end, genes, snp_location = line_snp[0], line_snp[1], line_snp[2], line_snp[3], line_snp[4], line_snp[5], line_snp[6], line_snp[7], line_snp[8], line_snp[9]
			if xQTL_name == key and snp_trait not in trait_list:
				dic_entry[snp_trait] = [condition, xQTL_name, xQTL_chr, xQTL_start, xQTL_end, xQTL_end - xQTL_start, 1, snp_name, snp_p, genes, [snp_name], snp_location]
				dic_xQTL[xQTL_name] = dic_entry[snp_trait]
				trait_list.append(snp_trait)
				p_min = snp_p
				snp_temp, p_temp, location_temp  = [snp_name, '-'], [p_min, 1], [snp_location, 0]
			elif xQTL_name == key and snp_trait in trait_list:
				dic_entry[snp_trait][6] += 1
				dic_entry[snp_trait][10].append(snp_name)
				snp_temp[1], p_temp[1], location_temp[1] = snp_name, snp_p, snp_location
				p_min = min(p_temp)
				index_leadsnp = p_temp.index(p_min)
				lead_snp = snp_temp[index_leadsnp]
				index_leadlocation = p_temp.index(p_min)
				lead_location = location_temp[index_leadlocation]
				dic_entry[snp_trait][7], dic_entry[snp_trait][8], dic_entry[snp_trait][11] = lead_snp, p_min, lead_location
				snp_temp, p_temp, location_temp = [snp_name, '-'], [p_min, 1], [snp_location, 0]
		for value in dic_entry.values():
			value[10] = '|'.join(value[10])
		dic_xQTL[key] = dic_entry
	return dic_xQTL


## main process, and output file

wing = my_SettleWing(sigsnp_all, wing_length)
xQTL = my_xQTLSummary(wing, anno)
summary_drought = my_snpSettle(sigsnp_drought, xQTL, 'drought')
summary_normal = my_snpSettle(sigsnp_normal, xQTL, 'normal')

fl_xQTL = index_export + '_xQTL_index.txt'
fl_summary = index_export + '_xQTL_summary.txt'

with open(fl_xQTL, 'a') as fh_xQTL:
	title = ['trait', 'xQTL', 'chr', 'start', 'end', 'gene']
	fh_xQTL.writelines('\t'.join(title))
	fh_xQTL.writelines('\n')
	for line in xQTL:
		line = [str(i) for i in line]
		fh_xQTL.writelines('\t'.join(line))
		fh_xQTL.writelines('\n')

with open(fl_summary, 'a') as fh_summary:
	title = ['trait', 'condition', 'xQTL', 'chr', 'start', 'end', 'length', 'snp_number', 'leadsnp', 'leadp', 'gene_list', 'allsnp', 'lead_location']
	fh_summary.writelines('\t'.join(title))
	fh_summary.writelines('\n')
	for key1 in sorted(summary_drought.keys()):
		for key2 in sorted(summary_drought[key1].keys()):
			value = [str(i) for i in summary_drought[key1][key2]]
			value = [key2] + value
			fh_summary.writelines('\t'.join(value))
			fh_summary.writelines('\n')
	
	for key1 in sorted(summary_normal.keys()):
		for key2 in sorted(summary_normal[key1].keys()):
			value = [str(i) for i in summary_normal[key1][key2]]
			value = [key2] + value
			fh_summary.writelines('\t'.join(value))
			fh_summary.writelines('\n')
