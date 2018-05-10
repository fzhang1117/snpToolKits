#! /bin/bash

## significant snps in drought
fl_drought=$1

## significant snps in normal
fl_normal=$2

## hmp file contains all significant snps
fl_hmp=$3

## gff3 annotation file
fl_gff=$4

## window of bin
wing=$5

## the prefix is the share paramater in both xQTL_group.py and xQTL_LDcal.py
prefix=$6

## single_filter, should be T/F or TRUE/FALSE
singlefilter=$7

## the threshold of LD
threshold=$8

## the flanking_region
flanking=$9

## for example
## sh run.sh ../data/mlm_GC_drought_5.00.txt ../data/mlm_GC_normal_5.00.txt ../data/sorted_mGWAS_sig_snp.all.hmp ../../../data/ZmB73_5b_FGS_GeneOnly.gff3 10000 ../result/mlm_GC_5.00_20kb_0.1_100kb T 0.1 100000

python ./xQTL_group.py $fl_drought $fl_normal $wing $prefix $fl_gff
sort -k1,1 -k2,2 -k3,3n ${prefix}_xQTL_summary.txt > ${prefix}_xQTL_summary_sorted.txt
grep "drought" ${prefix}_xQTL_summary_sorted.txt > ${prefix}_xQTL_drought.txt
grep "normal" ${prefix}_xQTL_summary_sorted.txt > ${prefix}_xQTL_normal.txt
python ./xQTL_LDcalc.py $fl_hmp ${prefix}_xQTL_drought.txt ${prefix}_drought
python ./xQTL_LDcalc.py $fl_hmp ${prefix}_xQTL_normal.txt ${prefix}_normal
sed -i '1d' ${prefix}_normal_xQTL_ld.txt
cat ${prefix}_drought_xQTL_ld.txt ${prefix}_normal_xQTL_ld.txt > ${prefix}_xQTL_ld.txt
rm -rf ${prefix}_xQTL_summary_sorted.txt
Rscript ./xQTL_merge.R ${prefix}_xQTL_summary.txt ${prefix}_xQTL_ld.txt ${prefix} ${threshold} $singlefilter
rm -rf ${prefix}_normal_xQTL_ld.txt ${prefix}_drought_xQTL_ld.txt ${prefix}_xQTL_drought.txt ${prefix}_xQTL_normal.txt
python CandidateGene_Search.py ${prefix}_final_summary.txt $fl_gff $flanking ${prefix}_candidate_temp.txt
