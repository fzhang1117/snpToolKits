#! /bin/bash
hmp=$1
chrom=$2
sta=$3
end=$4
out=$5
picName=$6

python SNPExtract.py $hmp $chrom $sta $end $out
sh ld_tassel_pbs.sh $out
Rscript LD_heatmap.R
mv ld_plot.png $picName
rm $out
