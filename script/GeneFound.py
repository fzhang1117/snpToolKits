## Significatnt SNP annotation from GFF file ##
## Zhang Fei ##
## 2016-09-02 ##

import sys
from itertools import islice 

input_sig = open(sys.argv[1], 'r')
input_gff = open(sys.argv[2], 'r')  # must contain only gene in chromosomes in Zmays 5b gff file
output = open(sys.argv[3], 'a')

def gff_dictbuild(FH):
    dic = {}

    for line in FH:
        key = []
        value = []
        line = line.strip('\n')
        line = line.split('\t')
        key = (line[8].split(";"))[0].split("=")[1] # get gene name
        value.append(key)
        value.append((int(line[0][3])))
        #print line[0][3]
        #print type(int(line[0][3]))
        value.append(int(line[3]))
        value.append(int(line[4]))  #The list value contains: Gene name; chromosome; start; end;
        dic[key] = value

    return dic

gff_dic = gff_dictbuild(input_gff)
title = ["Phenotype", "SNP_name", "Chromosome", "LOCUS", "p-value", "gene"]
output.writelines('\t'.join(title))
output.writelines('\n')

print gff_dic['GRMZM2G378836']

for line in islice(input_sig, 1, None):
	line = line.strip('\n')
	line = line.split('\t')
	for key in gff_dic.keys():
		value = gff_dic.get(key)
		if int(line[2]) == value[1] and (value[2] <= int(line[3])) and (value[3] >= int(line[3])):
			out = [line[0], line[1], line[2], line[3], line[6], value[0]]
			output.writelines('\t'.join(out))
			output.writelines('\n')
		else:
			continue

	#print line
"""
for line in islice(input_sig, 1, None):
    line = line.strip('\n')
    line = line.split('\t')
    print line
    for key in gff_dic.keys():
		#print gff_dic.get(key)[1]
        if line[1] == gff_dic.get(key)[1]:
            out = []
            value = gff_dic.get(key)
            if (value[2] <= int(line[2])) and (value[3] >= int(line[2])):
                out = [line[0],line[1],line[2],line[3],line[4],value[0]]
                output.writelines('\t'.join(out))
                output.writelines('\n')
        else:
            continue
"""
input_sig.close()
input_gff.close()
output.close()

