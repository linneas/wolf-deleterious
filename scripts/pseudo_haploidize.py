#!/usr/bin/python3

################################################################################
# pseudo_haploidize.py
# Author: Linnea Smeds
# Date: 12 May 2021
# This script extract a single allele for each site to be used as a haploid
# genotype. The allele is randomly chosen but using the allele frequency as
# weight (which means that it is more likely that the most common allele is
# chosen). This turned out to be more unbiased than using GATK genotypes, since
# they where biased towards the dog reference (1REF/9ALT was called as hetero-
# zygous much more often than 9REF/1ALT)

################################################################################
##### Import libraries
import argparse
import re
import random
random.seed(a=12345)
#print (random.__file__)


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script takes a vcf file \
from a single individual, and assigns a single most likely allele based on the \
allele frequency from the depth information")
parser.add_argument('-v', '--variationfile', help = 'vcf file with genotypes for one individual', type = str, required = True)
parser.add_argument('-d', '--minDepth', help = 'If depth is lower, the genotype is set as \"N\"', type = str, required = False)
parser.add_argument('-g', '--minGQ', help = 'If GQ is lower, the genotype is set as \"N\"', type = str, required = False)
parser.add_argument('-o', '--oprefix', help = 'Output prefix', type = str, required = True)
args = parser.parse_args()


################################################################################
##### Assign defaults

filtDict={}
if args.minDepth:
	filtDict["DP"]=args.minDepth
else:
	filtDict["DP"]=-1
if args.minGQ:
	filtDict["GQ"]=args.minGQ
	filtDict["RGQ"]=args.minGQ
else:
	filtDict["GQ"]=-1
	filtDict["RGQ"]=-1

################################################################################
##### Main code


# open output and print header
with open(args.oprefix+".bed", 'w') as outbed, open(args.oprefix+".log", 'w') as outlog:
	head1="#CHROM\tSTART\tSTOP\tGT"

	# Go through the variants
	with open(args.variationfile, 'r') as vcf:
		for line in vcf:
			s=re.match("^#", line)
			if not s:
				line = line.rstrip()
				tabs = line.split("\t")
				start=str(int(tabs[1])-1)
				GTdict={}
				format = tabs[8].split(":")
				ind = tabs[9].split(":")
				for i in range(0,len(format)):
					GTdict[format[i]]=ind[i]
				flag="ok"
				err=""
				# Go through all (min) filters and check value exceeds them
				for f in filtDict:
					if f in GTdict and int(GTdict[f])<int(filtDict[f]):
						flag="bad"
						#print("GTdict[f] is "+str(GTdict[f]))
						err=err+f+"="+GTdict[f]+","

				if flag=="ok":
					if GTdict["GT"]=="./.":
						outbed.write(str(tabs[0])+"\t"+start+"\t"+str(tabs[1])+"\tN\n")
					else:
						if "AD" not in GTdict:
							# There is no allele depth -> use reference allele
							gt=tabs[3]
						else:
							# There is allele depth information ->
							# use this as weight when choosing allele
							#print(GTdict["GT"])
							nt=[tabs[3]]+tabs[4].split(",")
							#print(nt)
							ad=GTdict["AD"].split(",")
							ad = [int(x) for x in ad]	#converting read counts to int
							#print(ad)
							# This python module draws k elements from nt,
							# randomly but weighted based on ad.
							gt=random.choices(nt, weights=ad, k=1)[0]
							#print(gt)

						outbed.write(tabs[0]+"\t"+start+"\t"+tabs[1]+"\t"+gt+"\n")

				else:
					outbed.write(tabs[0]+"\t"+start+"\t"+tabs[1]+"\tN\n")
					outlog.write("Removed: "+tabs[0]+"\t"+tabs[1]+"\t"+err+"\n")
