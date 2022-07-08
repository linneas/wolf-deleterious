#!/usr/bin/python3

################################################################################
# infer_founder_genotypes.py
# Author: Linnéa Smeds
# Date: 05 July 2021

################################################################################
##### Import libraries here
import argparse
import re
import math

################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script takes the asssigned \
haplotypes from the script assign_genotype_to_haplotype.py, and create fake  \
hapotypes for the two unsampled male founders.")
parser.add_argument('-a', '--assignedfile', help = 'vcf file with genotypes for several individual', type = str, required = True)
parser.add_argument('-p', '--haplotypefile', help = 'Bed like file with regions and haplotypes for the individuals', type = str, required = True)
parser.add_argument('-o', '--oprefix', help = 'Output prefix', type = str, required = True)
args = parser.parse_args()

# Fixed parameters
winSize=1000000


################################################################################
##### Main code

# Save regions and haplotypes in a dictionary
hDict={}
header=[]
with open(args.haplotypefile, 'r') as bed:
	for line in bed:
		s=re.match("^#CHROM", line)
		if s:
			line = line.rstrip()
			header = line.split("\t")
		if not s:
			line = line.rstrip()
			tabs = line.split("\t")
			#print("splitting tabs and saving "+tabs[0]+" as key0 and "+tabs[1]+" as key1") ###DEBUG COMMAND
			if not tabs[0] in hDict:
				hDict[tabs[0]]={}
			hDict[tabs[0]][tabs[1]]={}
			#hDict[tabs[0]][tabs[1]]['end']=tabs[2]	#Only needed if window size is unknown
			c=0
			for h in range(3,len(tabs)):
				hDict[tabs[0]][tabs[1]][header[h]]=tabs[h]
				c+=1

				#print("saving element "+tabs[h]+" with key "+header[h]) ###DEBUG COMMAND
				#print("Saving "+str(c)+" individuals in dict") ###DEBUG COMMAND
				# Debug: print the key and values ###DEBUG COMMAND
				#for k, v in hDict[tabs[0]][tabs[1]].items(): ###DEBUG COMMAND
    			#	print(k, v) ###DEBUG COMMAND

# Go through the assigned file
with open(args.oprefix, 'w') as out:
	out.write("#CHROM\tSITE\t"+"\t".join(header[3:])+"\n")

	with open(args.assignedfile, 'r') as infile:
		for line in infile:
			# Save informative header line
			s=re.match("^#CHROM", line)
			if not s:
				line = line.rstrip()
				tabs = line.split("\t")
				#Only look at sites in chromosomes present in the
				# haplotype file!
				if tabs[0] in hDict:
					#print("Chromosome "+tabs[0]+" exists!") ###DEBUG COMMAND
					siteDict={}
					indDict={}
					start=str(math.floor(int(tabs[1])/winSize)*winSize)
					#print(type(start))
					#print("Start is "+str(start)) ###DEBUG COMMAND
					indCount=0
					tempHapDict={}
					for haplo in range(2,len(tabs)):
						hapinfo=tabs[haplo].split(":")
						#print("DEBUG: pos "+tabs[0]+":"+tabs[1]+" adding haplo "+hapinfo[0]+" with value "+hapinfo[1])
						tempHapDict[hapinfo[0]]=hapinfo[1]

					tempout=[]
					for ind in range(3,len(header)):
						(h1,h2)=hDict[tabs[0]][start][header[ind]].split("|")
						newgeno="./."
						if h1 in tempHapDict and h2 in tempHapDict:
							g1=tempHapDict[h1]
							g2=tempHapDict[h2]
							if g1<g2:
								newgeno=g1+"/"+g2
							else:
								newgeno=g2+"/"+g1
						#else:
							#print("DEBUG: pos "+tabs[0]+":"+tabs[1]+" missing haplotype "+h1+" or "+h2)
						tempout.append(newgeno)
					out.write(tabs[0]+"\t"+tabs[1]+"\t"+"\t".join(tempout)+"\n")
