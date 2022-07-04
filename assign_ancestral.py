#!/usr/bin/python3

################################################################################
# assign_ancestral.py
# Author: Linnea Smeds
# Date: 4 May 2021

################################################################################
##### Import libraries
import argparse
import re
import sys

################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script takes a table with \
sites and outgroup genotypes, and assign the most likely ancestral state.")
parser.add_argument('-i', '--infile', help = 'table with sites and genotypes', type = str, required = True)
parser.add_argument('-c', '--colStart', help = 'Genotypes start in this column \
 [default: 5th column]', type = int, required = False)
parser.add_argument('-o', '--oprefix', help = 'Output prefix', type = str, required = True)
args = parser.parse_args()


################################################################################
##### Assign defaults
c=5
if args.colStart:
	c=args.colStart


################################################################################
##### Main code

print("Processing file "+ args.infile + "; Outgroups starting on "+
			str(c)+"th column, (1-based counting)\n")


# open output and print header
with open(args.oprefix+".ancestral.txt", 'w') as out:

	# Go through the sites
	with open(args.infile, 'r') as infile:
		for line in infile:
			s=re.match("^#", line)
			if not s:
				line = line.rstrip()
				tabs = line.split("\t")

				NTdict={"A":0, "C":0, "G":0, "T":0}
				tot=0.0
				colstart=c-1

				for i in range(colstart,len(tabs)):
					if not tabs[i]=="N":
						tot+=1
						al = tabs[i].split("/")
						if(len(al)>1):
							NTdict[al[0]]+=0.5
							NTdict[al[1]]+=0.5
						else:
							NTdict[al[0]]+=1

				ancestral="N"
				support=0.0
				#for k in sorted(NTdict, key=lambda k: NTdict[k], reverse=True):
				best,second=sorted(NTdict, key=NTdict.get, reverse=True)[:2]
				if NTdict[best]==NTdict[second]:
					sys.stderr.write(tabs[0]+":"+tabs[1]+", equal support for "+best+" and "+second+"\n")
				else:
					ancestral=best
					support=NTdict[best]/tot
					#print("support is "+str(NTdict[best])+" divided by "+str(tot)+"="+str(support))

				tot=str(int(tot))
				out.write("\t".join(tabs[0:colstart])+"\t"+ancestral+"\t"+str("%.2f" % support)+"\t"+tot+"\n")
