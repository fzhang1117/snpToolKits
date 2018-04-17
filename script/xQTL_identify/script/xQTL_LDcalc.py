## xQTL_LDcalc.py ##
## calculate all the xQTL for a trait one by one ##
## Zhang Fei <zhangfei-123@foxmail.com> ##
## 2018-04-12 ##

from __future__ import division
import sys, copy
import numpy as np
import pandas as pd
from operator import itemgetter, attrgetter
## fl_hmp is the significant snp hmp file
## fl_xQTL is the xQTL_summary file which is the output file of xQTL_group.py
## fl_out_path, outpath , is better same as the path xQTL_group.py
fl_hmp = sys.argv[1]
fl_xQTL = sys.argv[2]
fl_out_path = sys.argv[3]

def my_dicHMPbuild(fl_hmp):
	dicHMP = {}
	i = 1
	with open(fl_hmp, 'r') as fh_hmp:
		for line in fh_hmp:
			line = line.strip('\n')
			line = line.split('\t')
			if i == 1:
				i += 1
			else:
				snp_name = line[0]
				snp_info = line[1: 5]
				hapmap = line[11: ]
				dicHMP[snp_name] = (snp_info, hapmap)
	return dicHMP

dicHMP = my_dicHMPbuild(fl_hmp)
#print dicHMP

def my_LDcalculate(x, y):
	## replace 'N' as nan which pandas and numpy can recognize
	xy = map(lambda a,b: a+b, x,y)
	x = [np.nan if a == 'N' else a for a in x]
	y = [np.nan if a == 'N' else a for a in y]
	## paste x and y one-by-one
	dic = {'x':x, 'y': y, 'xy': xy}
	df1 = pd.DataFrame(dic)
	df1 = df1.dropna()
	x, y, xy = list(df1.x), list(df1.y), list(df1.xy)
	pA = x.count(x[0]) / len(x)
	pB = y.count(y[0]) / len(y)
	pAB = xy.count(xy[0]) / len(xy)
	D = pAB -pA*pB
	if pA*(1-pA)*pB*(1 - pB) != 0:
		R2 = pow(D, 2)/(pA*(1-pA)*pB*(1-pB))
	else:
		R2 = "-Inf"
	## count the factor in each list, solve in a dict
	#summary_x = dict((a, x.count(a)) for a in x)
	#summary_y = dict((a, y.count(a)) for a in y)
	#summary_xy = dict((a, xy.count(a)) for a in xy)
	return R2


with open(fl_xQTL, 'r') as fh_xQTL:
	xQTL = fh_xQTL.read()
	xQTL = xQTL.strip('\n')
	xQTL = xQTL.split('\n')
	xQTL.pop(0)

xQTL_copy = copy.deepcopy(xQTL)

result = []

for line1 in xQTL:
	line1 = line1.split('\t')
	trait1, lead_snp1, xQTL1 = line1[0], line1[7], line1[1]
	for line2 in xQTL_copy:
		line2 = line2.split('\t')
		trait2, lead_snp2, xQTL2 = line2[0], line2[7], line2[1]
		if trait1 == trait2:
			if lead_snp2 != lead_snp1:
				hapmap1 = dicHMP[lead_snp1][1]
				hapmap2 = dicHMP[lead_snp2][1]
				R2 = my_LDcalculate(hapmap1, hapmap2)
				## more test file, imporve the output
				output = [trait1, xQTL1, xQTL2, lead_snp1, lead_snp2, R2]
				result.append(output)
		
	xQTL_copy.pop(0)

fl_out = fl_out_path + '_xQTL_ld.txt'
with open(fl_out, 'a') as fh_out:
	title = ['trait', 'xQTL1', 'xQTL2', 'lead_snp1', 'lead_snp2', 'R2']
	fh_out.writelines('\t'.join(title))
	fh_out.writelines('\n')
	for line in result:
		line = [str(i) for i in line]
		fh_out.writelines('\t'.join(line))
		fh_out.writelines('\n')
