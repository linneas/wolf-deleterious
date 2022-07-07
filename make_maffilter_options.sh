#!/bin/bash -l

tree=$1
in=$2
out=$3
log=$4
outmaf=$5

wanted=`sed 's/:/\n/g' $tree |sed 's/,/\n/g' |sed 's/(//g' |sed 's/)//g' |grep ^[a-z] |grep -v "canFam" |tr "\n" "," |sed 's/,$//'`

echo "input.format=Maf" >$out
echo "input.file=$in" >>$out
echo "input.file.compression=gzip" >>$out
echo "output.log=$log" >>$out
echo "maf.filter=Subset(species=($wanted),strict=no,keep=no,remove_duplicates=no),Output(file=$outmaf,compression=none,mask=no)" >>$out
