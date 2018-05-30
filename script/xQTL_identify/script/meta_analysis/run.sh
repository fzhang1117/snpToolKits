#! /usr/bin/sh

fl_cpdID=$1
fl_literGene=$2
fl_myGene=$3
prefix=$4

time python cpdMapping.py $fl_cpdID $fl_literGene $fl_myGene ${prefix}_literature.txt ${prefix}_myGene.txt
time Rscript meta_analysis.R ${prefix}_literature.txt ${prefix}_myGene.txt ${prefix}
rm -rf ${prefix}*.log
