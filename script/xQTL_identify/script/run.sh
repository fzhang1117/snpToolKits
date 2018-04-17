#! /bin/bash

## significant snps from gwas
fl_tassel=$1

## hmp file contains all significant snps
fl_hmp=$2

## gff3 annotation file
fl_gff=$3

## window of bin
wing=$4

## the prefix is the share paramater in both xQTL_group.py and xQTL_LDcal.py
prefix=$5

## the threshold of LD
#threshold=$6

## single_filter, should be T/F or TRUE/FALSE
#singlefilter=$7

## for example
## sh run.sh ../data/mlm_GC_drought_select.txt ../data/sorted_mGWAS_sig_snp.all.hmp ../../../data/ZmB73_5b_FGS_GeneOnly.gff3 5000 ../result/mlm_GC_drought_10kb

python ./xQTL_group.py $fl_tassel $wing $prefix $fl_gff
sort -k1,1 -k3,3n -k3,3n ${prefix}_xQTL_summary.txt > ${prefix}_xQTL_summary_sorted.txt
python ./xQTL_LDcalc.py $fl_hmp ${prefix}_xQTL_summary_sorted.txt $prefix
rm -rf ${prefix}_xQTL_summary.txt
#Rscript ./xQTL_merge.R ${prefix}_xQTL_summary.txt ${prefix}_xQTL_ld.txt ${prefix} ${threshold} $singlefilter
