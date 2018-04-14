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

## for example
## sh run.sh ../data/mlm_GC_drought_select.txt ../data/sorted_mGWAS_sig_snp.all.hmp ../../../data/ZmB73_5b_FGS_GeneOnly.gff3 5000 ../result/mlm_GC_drought_10kb

python ./xQTL_group.py $fl_tassel $wing $prefix $fl_gff
python ./xQTL_LDcalc.py $fl_hmp ${prefix}_xQTL_summary.txt $prefix
