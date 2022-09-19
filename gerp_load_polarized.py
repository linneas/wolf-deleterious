#!/usr/bin/python3

################################################################################
# GerpLoad_polarized.py
# Author: Linnea Smeds
# Date: 28 Jan 2022


################################################################################
##### Import libraries
import argparse
import re


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script takes a vcf file,  \
a list with individuals and populations, and a bed file with ancestral alleles, \
and counts the alleles in different categories")
parser.add_argument('-l', '--populationlist', help = 'Name of the file with individuals and populations', type = str, required = True)
parser.add_argument('-v', '--variationfile', help = 'vcf file with genotypes for all individuals', type = str, required = True)
parser.add_argument('-b', '--ancestralbed', help = 'Name of the bed file with ancestral alleles', type = str, required = True)
parser.add_argument('-g', '--gerpbed', help = 'Name of the bed file with gerp scores in 5th column', type = str, required = True)
parser.add_argument('-t', '--thres', help = 'Gerp threshold for deleterious sites', type = float, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()


################################################################################
##### Main code

# Create a dictionary with the individuals as keys and a list of populations
indDict={}
popList=[]
with open(args.populationlist, 'r') as popListFile:
	for line in popListFile:
		line = line.rstrip()
		(i,p) = line.split("\t")
		indDict[i]={}
		indDict[i]['p']=p
		indDict[i]['sum']=0
		indDict[i]['tot_der']=0
		indDict[i]['hom_der']=0
		indDict[i]['hom_del_sum']=0
		indDict[i]['het_del_sum']=0
		indDict[i]['del']=0
		indDict[i]['hom_del']=0
		indDict[i]['miss']=0
		indDict[i]['calls']=0
		if not p in popList:
			popList.append(p)


print("saved",len(indDict),"individuals as keys")
print("from",len(popList),"populations:")
popList.sort()
print(popList)

# Create a dictionary with ancestral alleles
ancDict={}
with open(args.ancestralbed, 'r') as ancBedFile:
	for line in ancBedFile:
		line = line.rstrip()
		tabs = line.split("\t")
		if tabs[0] not in ancDict:
			ancDict[tabs[0]]={}
		ancDict[tabs[0]][tabs[2]]={}
		ancDict[tabs[0]][tabs[2]]['anc']=tabs[3]
		#print("DEBUG: saving chr"+tabs[0]+" site "+tabs[1]+" with value "+tabs[3])

gthres=args.thres
with open(args.gerpbed, 'r') as gerpBedFile:
	for line in gerpBedFile:
		line = line.rstrip()
		tabs = line.split("\t")
		if tabs[2] in ancDict[tabs[0]]:
			ancDict[tabs[0]][tabs[2]]['gerp']=float(tabs[4])


# open outfile and print header
c_all=0
c_ok=0
c_miss=0

# Saving Header in a list, then go through the variants
header=[]
with open(args.variationfile, 'r') as vcf:
	for line in vcf:
		s=re.match("^#", line)
		if s:
			if "#CHROM" in line:
				line = line.rstrip()
				header = line.split("\t")

		else:

			line = line.rstrip()
			tabs = line.split("\t")
			c_all+=1
			# If we have an ancestral allele assigned and we have a gerp score
			if tabs[0] in ancDict and tabs[1] in ancDict[tabs[0]] and 'gerp' in ancDict[tabs[0]][tabs[1]]:
				anc=ancDict[tabs[0]][tabs[1]]['anc']
				gerp=ancDict[tabs[0]][tabs[1]]['gerp']
				der=""
				if anc==tabs[3]:
					der="alt"
				elif anc==tabs[4]:
					der="ref"
				else:
					print("WARNING: "+tabs[0]+":"+tabs[1]+
					" ancestral allele "+anc+" is neither ref ("
					+tabs[3]+") or alt ("+tabs[4]+")")
					continue

				for i in range(9,len(header)):
					parts=tabs[i].split(":")
					if parts[0]=="0/0":
						indDict[header[i]]['calls']+=1
						if der=="ref":
							indDict[header[i]]['tot_der']+=2
							indDict[header[i]]['hom_der']+=2
							if ancDict[tabs[0]][tabs[1]]['gerp'] > gthres:
								indDict[header[i]]['del']+=2
								indDict[header[i]]['hom_del']+=2
								indDict[header[i]]['sum']+=2*ancDict[tabs[0]][tabs[1]]['gerp']
								indDict[header[i]]['hom_del_sum']+=2*ancDict[tabs[0]][tabs[1]]['gerp']
					elif parts[0]=="0/1":
						indDict[header[i]]['calls']+=1
						indDict[header[i]]['tot_der']+=1
						if ancDict[tabs[0]][tabs[1]]['gerp'] > gthres:
							indDict[header[i]]['del']+=1
							indDict[header[i]]['sum']+=ancDict[tabs[0]][tabs[1]]['gerp']
							indDict[header[i]]['het_del_sum']+=ancDict[tabs[0]][tabs[1]]['gerp']
					elif parts[0]=="1/1":
						indDict[header[i]]['calls']+=1
						if der=="alt":
							indDict[header[i]]['tot_der']+=2
							indDict[header[i]]['hom_der']+=2
							if ancDict[tabs[0]][tabs[1]]['gerp'] > gthres:
								indDict[header[i]]['del']+=2
								indDict[header[i]]['hom_del']+=2
								indDict[header[i]]['sum']+=2*ancDict[tabs[0]][tabs[1]]['gerp']
								indDict[header[i]]['hom_del_sum']+=2*ancDict[tabs[0]][tabs[1]]['gerp']


					else:
						indDict[header[i]]['miss']+=1
				c_ok+=1
			else:
				c_miss+=1
				#print("DEBUG: Site "+tabs[0]+":"+tabs[1]+" is missing ancestral assignment, skipping...")

# Print outfile (one pop at the time)
with open(args.output, 'w') as outfile:

	head="IND\tPOP\tRELATIVE_LOAD\tREALIZED_LOAD\tMASKED_LOAD\tTOT_SUM\tHOMDEL_SUM\tTOT_DERIVED\tHOM_DER\tTOT_DEL\tHOM_DEL\tTOT_CALLS\tMISS"
	print(head, file=outfile)

	for p in popList:
		for k in sorted(indDict.keys()):
			if(indDict[k]['p']==p):
				if(indDict[k]['tot_der']==0):
					load="NA"
					print("There are no derived alleles for ind"+k)
				else:
					reli_load=indDict[k]['sum']/indDict[k]['tot_der']
					real_load=indDict[k]['hom_del_sum']/indDict[k]['calls']
					mask_load=indDict[k]['het_del_sum']/indDict[k]['calls']
				print(k+"\t"+p+"\t"+str(reli_load)+"\t"+str(real_load)+"\t"+str(mask_load)+"\t"+str(indDict[k]['sum'])+"\t"+str(indDict[k]['hom_del_sum'])+"\t"+str(indDict[k]['tot_der'])+"\t"+str(indDict[k]['hom_der'])+"\t"+str(indDict[k]['del'])+"\t"+str(indDict[k]['hom_del'])+"\t"+str(indDict[k]['calls'])+"\t"+str(indDict[k]['miss']), file=outfile)



print("Processed "+str(c_all)+" sites, wrote "+str(c_ok)+" to output. "+str(c_miss)+" sites removed due to missing ancestral state")
