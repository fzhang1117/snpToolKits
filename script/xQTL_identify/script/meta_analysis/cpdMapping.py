#!/usr/bin/python
# -*- coding: UTF-8 -*-

## meta_analysis.py ##
## zhang fei <zhangfei-123@foxmail.com> ##
## 2018-05-28 ##

import sys, codecs
from itertools import islice

## cpdName cpdID reference ##
fl_cpdID = sys.argv[1]

## literGene list
fl_literGene = sys.argv[2]

## mapping result
fl_myGene = sys.argv[3]

## output
fl_out_candidate = sys.argv[4]
fl_out_myGene = sys.argv[5]

with codecs.open(fl_cpdID, 'r', encoding = 'utf-8') as fh_cpdID:
    dic_cpdID = {}
    for line in islice(fh_cpdID, 1, None):
        line = line.strip('\r\n').split('\t')
        if len(line) >= 3:
            if line[1] == 'NA':
                print line[0], ' is ', 'NA'
            elif dic_cpdID.get(line[2]) is None:
                dic_cpdID[line[2]] = {line[0]: line[1]}
            else:
                dic_cpdID[line[2]][line[0]] = line[1]

def my_cpdMapping(fl, dic):
    result = []
    with codecs.open(fl, 'r', encoding = 'utf-8') as fh:
        for line in islice(fh, 1, None):
            line = line.strip('\r\n').split('\t')
            if dic.get(line[2]) is not None:
                if dic[line[2]].get(line[0]) is not None:
                    entry = [line[0], dic[line[2]][line[0]], line[1], line[2]]
                    result.append(entry)
                else:
                    entry = [line[0], 'NA', line[1], line[2]]
                    result.append(entry)
    return result

cpdMapping_candidate = my_cpdMapping(fl_literGene, dic_cpdID)
cpdMapping_myGene = my_cpdMapping(fl_myGene, dic_cpdID)

with codecs.open(fl_out_candidate, 'w', encoding = 'utf-8') as fh_out_candidate:
    cpdMapping_candidate = ['\t'.join(i) for i in cpdMapping_candidate]
    fh_out_candidate.write('\n'.join(cpdMapping_candidate))

with codecs.open(fl_out_myGene, 'w', encoding = 'utf-8') as fh_out_myGene:
    cpdMapping_myGene = ['\t'.join(i) for i in cpdMapping_myGene]
    fh_out_myGene.write('\n'.join(cpdMapping_myGene))
