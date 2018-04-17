#! /bin/bash

## xQTL_merge.sh ##
## this script is used to merge xQTLs in normal or drought condition

xQTL_summary=$1

awk -v OFS='\t' '{print $3, $4, $5, $2}' $xQTL_summary | sed '1d' | sort -k1,1 -k2,2n > ./xQTL_summary.bed
bedtools merge -i ./xQTL_summary.bed ./xQTL_merged.bed
