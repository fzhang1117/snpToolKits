## SNP_Extract.py ##
## Extract SNPs by physical location ##
## Zhang Fei ##
## 2018-04-02 ##

"""
This script are used to extract SNPs by physical location
input the hmp file, chrom, start and end range, and export the SNPs in this range
"""

import sys

fl_hmp = sys.argv[1]
chrom = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])
fl_out = sys.argv[5]

with open(fl_out, 'a') as fh_out:
	i = 1
	with open(fl_hmp, 'r') as fh_hmp:
		for line in fh_hmp:
			line = line.strip('\n')
			line = line.split('\t')
			if i == 1:
				fh_out.writelines('\t'.join(line))
				fh_out.writelines('\n')
				i += 1
			else:
				chrom_snp = int(line[2])
				locus_snp = int(line[3])
				if chrom_snp == chrom and (locus_snp <= end and locus_snp >= start):
					fh_out.writelines('\t'.join(line))
					fh_out.writelines('\n')

