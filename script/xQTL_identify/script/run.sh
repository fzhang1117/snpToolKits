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

## the pGWAS result
fl_pGWAS=${10}

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
python CandidateGene_Search.py ${prefix}_xQTL_final_summary.txt $fl_gff $flanking ${prefix}_candidate_temp.txt
awk -F '\t' '{print $9}' ${prefix}_candidate_temp.txt > ${prefix}_gene_temp.txt
python kegg_annotation.py ${prefix}_gene_temp.txt ${prefix}_kegg_annotation_temp.txt
paste ${prefix}_candidate_temp.txt ${prefix}_kegg_annotation_temp.txt > ${prefix}_candidate_anno.txt
python ~/tools/v3tov4_trans/v3tov4_trans.py ~/tools/v3tov4_trans/library/v3_v4_xref.txt ${prefix}_gene_temp.txt ${prefix}_gene_v4.txt
python ~/tools/maize_annotation/auto-annotation.py ${prefix}_gene_temp.txt ${prefix}_gene_annotation.txt
sed '1 itrait\tconditon\txQTL\tchr\tstart\tend\tleadsnp\tleadp\tCandidateGene\tGeneLocation\tDistance\tCandidateGene\tPathway\tPathwayAnnotation' -i ${prefix}_candidate_anno.txt
paste ${prefix}_candidate_anno.txt ${prefix}_gene_v4.txt ${prefix}_gene_annotation.txt > ${prefix}_CandidateGene.txt
## cpd and gene expression file has been fixed here, if change the dir or result address, check and chage here
Rscript ./GeneCpdRegression.R ../data/GC_DN_Normal_GWAS.txt ../data/GC_DN_Drought_GWAS.txt ../data/ephenoMatrixC ../data/ephenoMatrixD ../result/${prefix}_CandidateGene.txt ${prefix}
Rscript ./xQTL_effect.R $fl_hmp ${prefix}_xQTL_final_summary.txt ../data/GC_DN_Drought_GWAS.txt ../data/GC_DN_Normal_GWAS.txt ${prefix}_final
Rscript ./xQTL_effect.R $fl_hmp ${prefix}_xQTL_summary.txt ../data/GC_DN_Drought_GWAS.txt ../data/GC_DN_Normal_GWAS.txt ${prefix}_all

python xQTL_pGWASsigsnpExtract.py $fl_pGWAS ${prefix}_xQTL_final_summary.txt $flanking 0.05 $fl_gff ${prefix}_xQTL_final
## result file has been fixed here, if change the dir or result address, change this line.
#python ./xQTL_multiEnviCampare.py ~/AMP_drought/GC_MS/result/DN_drought_GC/stats/ ~/AMP_drought/GC_MS/result/DN_normal_GC/stats/ ${prefix}_xQTL_summary.txt ${prefix}_xQTL_all_compare.txt
#python ./xQTL_multiEnviCampare.py ~/AMP_drought/GC_MS/result/DN_drought_GC/stats/ ~/AMP_drought/GC_MS/result/DN_normal_GC/stats/ ${prefix}_xQTL_final_summary.txt ${prefix}_xQTL_final_compare.txt
rm -rf ${prefix}_candidate_temp.txt ${prefix}_kegg_annotation_temp.txt ${prefix}_gene_temp.txt ${prefix}_candidate_anno.txt ${prefix}_gene_v4.txt ${prefix}_gene_annotation.txt
