## EntryExtractByLocation.py ##
## Extract entry in given location from Tassel GWAS stat file ##
## Zhang Fei ##
## 2018-04-03 ##

import sys

fl_gwas = sys.argv[1]
chrom = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])
fl_out = sys.argv[5]

with open(fl_out, 'a') as fh_out:
	with open(fl_gwas, 'r') as fh_gwas:
		i = 1
		for line in fh_gwas:
			line = line.strip('\n')
			line = line.split('\t')
			if i == 1:
				fh_out.writelines('\t'.join(line))
				fh_out.writelines('\n')
				i += 1
			elif line[2] != '':
				chrom_snp = int(line[2])
				location_snp = int(line[3])
				if chrom_snp == chrom and (location_snp >= start and location_snp <= end):
					fh_out.writelines('\t'.join(line))
					fh_out.writelines('\n')
