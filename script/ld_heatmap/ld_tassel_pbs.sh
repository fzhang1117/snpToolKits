#PBS ld_calculate
#PBS -l nodes=1:ppn=1
#PBS -l mem=20G
#! /bin/bash

hmp=$1
run_pipeline.pl -Xmx20g -fork1 -h $hmp -ld -export ld_output -runfork1
