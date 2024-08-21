import re
import numpy as np
import pandas as pd
import scipy
import sys


def vcf2hapmatrix(vcf):
	hap_matrix_d1 = {} #haplotype 1 of individuals, chromosome number is the key for dict
	hap_matrix_d2 = {} #haplotype 2 of individuals, chromosome number is the key for dict
	variant_names = {}
	variant_positions = {}# Chromosome number is the key for dict
	ref = {} # Chromosome number is the key for dict
	alt = {} # Chromosome number is the key for dict
	chromosome = [] #Chromosome number
	
	with open(vcf,"r") as VCF:
		for line in VCF:
			if re.search("^##",line): ## skip the first annotation lines
				continue
			elif re.search("^#CHROM",line): ## acquire the sample name information
				line = line.strip("\n")
				ind_names = line.split("\t")[9:]
			else:
				line = line.strip("\n")
				items = line.split("\t")
				ch = items[0]
				if items[2] == "\.":
						sys.exit("Found at least one empty variant names, please name variant sites accordingly.")

				if ch not in chromosome:
					chromosome.append(ch)
					variant_names[ch] = [items[2]]
					variant_positions[ch] = [int(items[1])]
					ref[ch] = [items[3]]
					alt[ch] = [items[4]]
					hap_matrix_d1[ch] = []
					hap_matrix_d2[ch] = []
					genotype = items[9:]
					for i in range(len(genotype)):
						m = re.search('([0-9])\|([0-9])',genotype[i])
						hap_matrix_d1[ch].append(int(m.group(1)))
						hap_matrix_d2[ch].append(int(m.group(2)))
				else:
					variant_names[ch].append(items[2])
					variant_positions[ch].append(int(items[1]))
					ref[ch].append(items[3])
					alt[ch].append(items[4])
					genotype = items[9:]
					for i in range(len(genotype)):
						m = re.search('([0-9])\|([0-9])',genotype[i])
						hap_matrix_d1[ch].append(int(m.group(1)))
						hap_matrix_d2[ch].append(int(m.group(2)))

	for ch in chromosome:
		hap_matrix_d1[ch] = np.reshape(np.asarray(hap_matrix_d1[ch],dtype=int),(len(variant_names[ch]),len(ind_names)))
		hap_matrix_d2[ch] = np.reshape(np.asarray(hap_matrix_d2[ch],dtype=int),(len(variant_names[ch]),len(ind_names)))

	return(ind_names,hap_matrix_d1,hap_matrix_d2,variant_names,variant_positions,ref,alt,chromosome)


