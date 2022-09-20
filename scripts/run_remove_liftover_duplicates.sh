#!/bin/bash -l


awk '{if(NR==1){prev=$0; ps=$2}else{if($2==ps){prev="0"}else{if(prev!="0"){print prev};prev=$0; ps=$2}}}' $1 >$2
