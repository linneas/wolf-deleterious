#!/bin/bash -l

# Versions used
#vcftools/0.1.16
#htslib/1.10


gzvcf=$1
indlist=$2
chr=$3
outpref=$4

vcftools --gzvcf $gzvcf --keep $indlist --chr $chr --recode --recode-INFO-all --out $outpref && bgzip $outpref.recode.vcf
mv $outpref.recode.vcf.gz $outpref.vcf.gz
tabix -p vcf $outpref.vcf.gz
