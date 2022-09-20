#!/usr/bin/python3

################################################################################
# assign_genotype_to_haplotype.py
# Author: Linnéa Smeds
# Date: 22 April 2021

################################################################################
##### Import libraries here
import argparse
import re
import math

################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script combines a vcf file \
with genotypes from multiple individuals, and a bed-like file with regions and \
hapotypes for some or all of those individuals, and return a summary of which \
alleles are associated with which haplotypes for each site.")
parser.add_argument('-v', '--variationfile', help = 'vcf file with genotypes for several individual', type = str, required = True)
parser.add_argument('-p', '--haplotypefile', help = 'Bed like file with regions and haplotypes for the individuals', type = str, required = True)
parser.add_argument('-o', '--oprefix', help = 'Output prefix', type = str, required = True)
args = parser.parse_args()

# Fixed parameters
winSize=1000000
supp=0.75

################################################################################
##### Functions

# Function that takes a dictionary with haplotype combinations and genotypes,
# and tries to assign a genotype to each haplotype.
def assignHaplotypeAlleles(d):

	# First, make a new dictionary to save single haplotypes
	hap={}

	#Start with homozygous haplotypes, sorting after the key with most GTs
	for k in sorted(d, key=lambda k: len(d[k]), reverse=True):
		a,b=k.split("|")
		if a==b:
			gt={"0/0":0, "0/1":0, "1/1":0}
			# Go through genotypes
			for g in d[k]:
				if g in gt:
					gt[g]+=1

				else:
					print("Unknown genotype "+g+" detected! Ignored")

			# Take the most common by reverse sorting
			alleles=sorted(gt, key=lambda j: gt[j], reverse=True)[0]
			#print(alleles) ###DEBUG COMMAND
			tot=sumDictElements(gt)
			support=gt[alleles]/tot
			print("Haplocombo "+k+" has genotype "+alleles+" with support "+str(support))

			# only assign if it's a homozygous genotype, with a support higher than x
			if alleles=="0/1" or support<supp:
				print("\tSkipping haplo "+a+" for now..")
			else:

				hap[a]={}
				hap[a]['a']=int(alleles[0])
				hap[a]['s']=support
				hap[a]['t']=tot

	# once the homozygous combinations have been assigned, we can look at
	# heterozygous combinations
	for k in sorted(d, key=lambda k: len(d[k]), reverse=True):
		a,b=k.split("|")
		if a!=b:

			# If both haplotypes are already in hash,
			#we do not need to do an assignment
			if not (a in hap and b in hap):

				gt={"0/0":0, "0/1":0, "1/1":0}
				# Go through genotypes
				for g in d[k]:
					if g in gt:
						gt[g]+=1

					else:
						print("Unknown genotype "+g+" detected! Ignored")

				# Take the most common by reverse sorting
				alleles=sorted(gt, key=lambda j: gt[j], reverse=True)[0]
				#print(alleles)	###DEBUG COMMAND
				tot=sumDictElements(gt)
				support=gt[alleles]/tot
				print("Haplocombo "+k+" has genotype "+alleles+" with support "+str(support))

				## Use already assigned, if there is one,
				# to infer the other
				a1,a2=alleles.split("/")
				a1=int(a1)
				a2=int(a2)
				flag=0

				for i in a,b:
					if i==a:
						o=b
					else:
						o=a

					if i in hap:
						flag+=1
						if hap[i]['a']==a1:
							#print("already assigned haplo "+i+" has allele "+str(a1)+", set "+o+" to "+str(a2)) ###DEBUG COMMAND
							assigned=a2
						elif hap[i]['a']==a2:
							#print("already assigned haplo "+i+" has allele "+str(a2)+", set "+o+" to "+str(a1)) ###DEBUG COMMAND
							assigned=a1
						else:

							#If we end up here, the most common GT is conflicting with previous assignments!
							print("Already assigned haplo "+i+" with allele "+str(hap[i]['a'])+" is conflicting with "+alleles)
							# This should mean that this allele comes from b (we only use biallelic sites,
							# so this only happens if GT is homozygous for not-a, meaning b should have this allele,
							# and this combination should be heterozygous!
							assigned=a1
							support=gt["0/1"]/tot

						hap[o]={}
						hap[o]['a']=assigned
						hap[o]['s']=support
						hap[o]['t']=tot
						break

				# If neither of a and b was saved before, try to assign!
				if flag==0:
					if a1==a2:
						for i in a,b:
							hap[i]={}
							hap[i]['a']=a1
							hap[i]['s']=support
							hap[i]['t']=tot
					else:
						print("Haplocombo "+k+" has no previously assigned haplotypes and gt "+alleles+" is ambigous")






	# After going through all homozygous and heterozygous,
	# go through everything and check overall support
	undef=0
	eq=0
	neq=0
	for k in d:
		a,b=k.split("|")


		#expected genotype
		#If both are assigned
		if a in hap and b in hap:
			exp=str(min(hap[a]['a'],hap[b]['a']))+"/"+str(max(hap[a]['a'],hap[b]['a']))
			#print("Looking at "+k+", expecting gt "+exp) ###DEBUG COMMAND
			for g in d[k]:
				if g==exp:
					eq+=1
				else:
					neq+=1

		# If not both are assigned, we set it as a third category
		else:
			for g in d[k]:
				undef+=1

	totsup=eq/(eq+neq+undef)
	det1=[str("%.2f" % totsup),str(eq),str(neq),str(undef)]

	# And summarize the haplotypes
	det2=[]
	#for k in sorted(hap.keys(), key=lambda t:int(t)):	#This row prints the keys in numberical order - does not work when some haplotypes are letters
	for k in sorted(hap.keys()):
		det2.append(k+":"+str(hap[k]['a'])+":"+str("%.2f" % hap[k]['s'])+":"+str(hap[k]['t']))


	print("Final results: "+":".join(det1)+"\t"+" ".join(det2))
	return det1,det2



# Help function that sum elements in dictionary
def sumDictElements(d):
	s=0
	for i in d:
		s+=d[i]
	return s


################################################################################
##### Main code

# Save regions and haplotypes in a dictionary
hDict={}
header1=[]
with open(args.haplotypefile, 'r') as bed:
	for line in bed:
		s=re.match("^#CHROM", line)
		if s:
			line = line.rstrip()
			header1 = line.split("\t")
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
				hDict[tabs[0]][tabs[1]][header1[h]]=tabs[h]
				c+=1
				#print("saving element "+tabs[h]+" with key "+header1[h]) ###DEBUG COMMAND
				#print("Saving "+str(c)+" individuals in dict") ###DEBUG COMMAND
				# Debug: print the key and values ###DEBUG COMMAND
				#for k, v in hDict[tabs[0]][tabs[1]].items(): ###DEBUG COMMAND
    			#	print(k, v) ###DEBUG COMMAND

# Go through the vcf file
header2=[]
with open(args.oprefix+".site.support.txt", 'w') as out1, open(args.oprefix+".assigned.haplotypes.txt", 'w') as out2:
	out1.write("#CHROM\tSITE\tSUPPORT\tIND_SUPPORT\tIND_MISMATCH\tIND_UNDEF\n")
	out2.write("#CHROM\tSITE\tALL_HAPLOTYPES(TAB_SEPARATED)\n")
	with open(args.variationfile, 'r') as vcf:
		for line in vcf:
			# Save informative header line
			s=re.match("^#", line)
			t=re.match("^#CHROM", line)
			if s and t:
				line = line.rstrip()
				header2 = line.split("\t")
			# For non header lines, convert site to window and save
			# all haplotypes and associated genotypes
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
					for ind in range(9,len(tabs)):
						#print("Looking at ind "+header2[ind]) ###DEBUG COMMAND
						#print("Looking in hash after chr: "+tabs[0]+" and start: "+start) ###DEBUG COMMAND
						#for k, v in hDict.items():  ###DEBUG COMMAND
	    					#	print(k, v, type(v)) ###DEBUG COMMAND
						if header2[ind] in hDict[tabs[0]][start]:
							gt=tabs[ind].split(":")[0]
							if not gt=="./.":
								indCount+=1
								haplo=hDict[tabs[0]][start][header2[ind]]
								if haplo in siteDict:
									siteDict[haplo].append(gt)
									indDict[haplo].append(header2[ind])
								else:
									siteDict[haplo]=[gt]
									indDict[haplo]=[header2[ind]]

					print("#--------")
					print(tabs[0]+"\t"+tabs[1]+", "+str(indCount)+"individuals with genotypes:")
					for h in siteDict:
						print("\t",h,siteDict[h])
						print("\t",h,indDict[h])
					indDict.clear()
					(o1,o2)=assignHaplotypeAlleles(siteDict)
					out1.write(tabs[0]+"\t"+tabs[1]+"\t"+"\t".join(o1)+"\n")
					out2.write(tabs[0]+"\t"+tabs[1]+"\t"+"\t".join(o2)+"\n")
