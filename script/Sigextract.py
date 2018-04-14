## Significant SNP extract ##
## Zhang Fei ##
## 2016-09-02 ##

import sys
import math
from itertools import islice

input = open(sys.argv[1],'r')
output = open(sys.argv[2],'a')
SigThreshold = sys.argv[3]

for line in islice(input,2,None):
    line = line.strip('\n')
    line = line.split('\t')
    if (-math.log10(float(line[6]))) >= float(SigThreshold): # The type of variable must be same when compare
        print (-math.log10(float(line[6]))),'\t',SigThreshold
        output.writelines('\t'.join(line))
        output.writelines('\n')

input.close()
output.close()