#!/bin/bash -l


ind=$1
acclist=$2

echo "Starting download of fastq"
date
err=`awk -v i=$ind '($2==i){print $1}' $acclist`
f1=`echo $err |cut -c1-6`
f2=`echo $err |cut -c10`

path1="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$f1"/00"$f2"/"$err"/"$err"_1.fastq.gz"
path2="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$f1"/00"$f2"/"$err"/"$err"_2.fastq.gz"
wget $path1 -O fastq/$ind"_1.fastq.gz"
echo "Done with first file!"
date
wget $path2 -O fastq/$ind"_2.fastq.gz"
echo "Done with second file!"
date
