## Significatnt SNP annotation from GFF file ##
## Zhang Fei ##
## 2016-09-02 ##

import sys
from itertools import islice 

fl_sig, fl_gff, fl_output = sys.argv[1], sys.argv[2], sys.argv[3]

def gff_dicbuild(fl_gff):
	dic = {}
	with open(fl_gff, 'r') as fh_gff:
		for line in fh_gff:
			line = line.strip('\n').split('\t')
			key = (line[8].split(';'))[0].split("=")[1]
			value = [key, str(line[0][3: ]), int(line[3]), int(line[4])]
			dic[key] = value
	return dic

gff_dic = gff_dicbuild(fl_gff)
with open(fl_output, 'a') as fh_output:
	title = ["Phenotype", "SNP_name", "Chromosome", "LOCUS", "p-value", "gene"]
	fh_output.writelines('\t'.join(title))
	fh_output.writelines('\n')
	with open(fl_sig, 'r') as fh_sig:
		for line in islice(fh_sig, 1, None):
			line = line.strip('\n').split('\t')
			if line[1] != 'None':
				chr_snp, location_snp = line[2], int(line[3])
				for key in gff_dic.keys():
					gene_name, chr_gene, start_gene, end_gene = gff_dic[key][0], gff_dic[key][1], int(gff_dic[key][2]), int(gff_dic[key][3])
					if (str(chr_snp) == str(chr_gene)) and (location_snp >= start_gene) and (location_snp <= end_gene):
						output = [line[0], line[1], line[2], line[3], line[6], gene_name]
						fh_output.writelines('\t'.join(output))
						fh_output.writelines('\n')
					else:
						continue
