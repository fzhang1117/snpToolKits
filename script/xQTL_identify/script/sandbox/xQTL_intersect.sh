#! /bin/bash
xQTL_normal=$1
xQTL_drought=$2

awk -v OFS='\t' '{print $3,$4,$5,$2}' $xQTL_normal | sed '1d' | sort -k1,1n -k2,2n -k3,3n > xQTL_normal.bed
awk -v OFS='\t' '{print $3,$4,$5,$2}' $xQTL_drought | sed '1d'| sort -k1,1n -k2,2n -k3,3n > xQTL_drought.bed

bedtools intersect -a xQTL_normal.bed -b xQTL_drought.bed -wa | wc -l
bedtools intersect -a xQTL_normal.bed -b xQTL_drought.bed -v | wc -l
bedtools intersect -a xQTL_drought.bed -b xQTL_normal.bed -v | wc -l
