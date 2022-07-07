#!/bin/bash -l

infile=$1

sort -k2,2n --buffer-size=6G $infile >$infile.sorted
