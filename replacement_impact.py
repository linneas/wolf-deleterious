#!/usr/bin/python3

################################################################################
# replacement_impact.py
# Date: 9 maj 2022
# Modification of the NBIS script simpred.py (written by Marcin Kierzcak)
# This version takes VEP output (instead of snpEff output) and reports the
# absolute scores from the Miyata, Exchangeability and Sneath tables.


################################################################################
##### Import libraries
import argparse
import re
import sys
import re
import numpy as np
import pandas as pd


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script takes the txt output from vep and add replacement scores to each missense variant")
parser.add_argument('-i', '--infile', help = 'Name of input (txt file from vep)', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()

################################################################################
##### Define functions

def readAnnotationMatrices():
	data = {}
	data['miyata'] = pd.read_csv('help_files/table.miyata.csv', index_col = ['AA'])
	data['exchgb'] = pd.read_csv('help_files/table.ex.csv', index_col = ['AA'])
	data['sneath'] = pd.read_csv('help_files/table.sneath.csv', index_col = ['AA'])
	return(data)

def load_aa_data():
	keys = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
	values = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
	aas = dict(zip(keys, values))
	return(aas)

def aa1_to_aa3(aaa):
	global aas
	if aaa != '':
		code = aas[aaa]
	else:
		code = ''
	return(code)

def parseVepFile(f, out):
	for line in f:
		if not line.startswith('#'):						# handle comments
			if 'missense_variant' in line:					# fish out missense variants only
				tmp = line.strip().split('\t')
				(chr,pos) = tmp[1].split(':')
				change = tmp[11]
				(ref,alt)=tmp[10].split('/')
				ref = aa1_to_aa3(ref)
				alt = aa1_to_aa3(alt)
				mi = str(miyata[ref][alt])
				ex1 = str(exchgb[ref][alt])
				ex2 = str(exchgb[alt][ref])
				sn = str(sneath[ref][alt])
				print(chr+"\t"+pos+"\t"+ref+"\t"+alt+"\t"+change+"\t"+mi+"\t"+ex1+"\t"+ex2+"\t"+sn, file=out)


##################### Main ##################
repdata = readAnnotationMatrices()
miyata = repdata['miyata']
exchgb = repdata['exchgb']
sneath = repdata['sneath']
aas = load_aa_data()
miyata_max = np.nanmax(repdata['miyata'].values)
exchgb_min = np.nanmax(repdata['exchgb'].values)
sneath_max = np.nanmax(repdata['sneath'].values)

with open(args.output, 'w') as outfile:
	print("CHROM\tPOS\tRef\tAlt\tChange\tMiyata\tExchgb1\tExchgb2\tSneath", file=outfile)
	with open(args.infile, 'r') as f:
		parseVepFile(f, outfile)
