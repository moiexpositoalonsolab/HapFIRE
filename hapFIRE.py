
#import modules
from sklearn import preprocessing
import time
import multiprocessing as mp
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import re
import sys
import os
import cvxpy as cp
#import utility scripts

import vcf_processing as VP
import haplotype_generation as HG


parser = argparse.ArgumentParser()
parser.add_argument('-v',type = str, action= 'store',dest='vcf',help='the vcf file')
parser.add_argument('-p',type = str, action= 'store',dest='block')
parser.add_argument('-b',type = str, action= 'store',dest='bam')
parser.add_argument('-f',type = str, action= 'store',dest='reference')
parser.add_argument('-s',type = bool, action= 'store',dest='haplotype_frequency',default = False)
parser.add_argument('-r',type = float, action = 'store', dest = 'corr',default = 0.1)
parser.add_argument('-c',type = float, action = 'store', dest = 'CLQcut',default = 0.5)
parser.add_argument('-w',type = int, action = 'store', dest = 'window',default = 100)
parser.add_argument('-o',type = str, action = 'store', dest = 'output',help = "the prefix of the output files")



args = parser.parse_args()


print("#############################\n program has just started \n#############################")



DIR = os.path.realpath(os.path.dirname(__file__))


"""
Import VCF files, block partition and SNP table generation
"""

########################## create the genotype matrix from vcf file
ind_names,hap_matrix_d1,hap_matrix_d2,variant_names,variant_positions,ref,alt,chromosome = VP.vcf2hapmatrix(vcf=args.vcf) # n*m matrix with n individuals and m snps

geno_matrix = {}
common_geno_matrix = {}
IndepLD_common_breakpoints_index = {}
common_allele_index = {}

common_variant_names = {}
common_variant_positions = {}
alone_common_SNPs_index = {}

gw_independent_breakpoints = {}
gw_fine_breakpoints = {}
snp_frequency = {}
independent_block_haplotype_frequency = {}
ecotype_frequency = {}
ecotype_frequency_selected = {}

block_partitions = {}
length = {}

with open(args.output+"_alone_SNPs.txt", "w") as L:
	for ch in chromosome:
		geno_matrix[ch] = np.transpose(hap_matrix_d1[ch] + hap_matrix_d2[ch])
		r,c = geno_matrix[ch].shape

		# calculating the allele frequency and find common alleles 
		allele_frequency = np.sum(geno_matrix[ch],axis = 0) / (2*r)
		common_allele_index[ch] = [i for i in range(len(allele_frequency)) if allele_frequency[i] > 0.05 and allele_frequency[i] < 0.95 ]
		common_geno_matrix[ch] = geno_matrix[ch][:,common_allele_index[ch]]
		common_variant_names[ch] = [variant_names[ch][i] for i in common_allele_index[ch]]
		common_variant_positions[ch] =[variant_positions[ch][i] for i in common_allele_index[ch]]

		#standardize the genotype matrix
		common_geno_matrix_standard = preprocessing.scale(common_geno_matrix[ch])
		#partition into complete independent LD blocks
		print("start finding complete independent LD blocks using common SNPs with maf > 0.1")

		IndepLD_common_breakpoints_index[ch],alone_common_SNPs_index[ch] = HG.CompleteLDPartition(standardized_genotype_matrix=common_geno_matrix_standard,cutoff=args.corr,window_size=args.window)
		print("hello",IndepLD_common_breakpoints_index[ch])
		print("%d complete independent LD blocks were found" %(len(IndepLD_common_breakpoints_index[ch])))

		if len(alone_common_SNPs_index[ch]) > 0:
			for i in alone_common_SNPs_index[ch]:
				full_index = common_allele_index[ch][i] # conver the common allele index to the index of all the variants 
				L.write(str(ch+"\t"+str(variant_positions[ch][full_index]))+"\n")


## calculate SNP frequencies from large independe blocks

for ch in chromosome:
		common_breakpoints = IndepLD_common_breakpoints_index[ch]
		gw_independent_breakpoints[ch] = HG.convert_independent_genomewide_breakpoints(common_breakpoints,common_allele_index[ch],len(variant_names[ch]))

with open(args.output+"_independent_genomewide_partition.txt", "w") as IND_PARTITION:
	for ch in chromosome:
		for i in range(len(gw_independent_breakpoints[ch])):
			left = gw_independent_breakpoints[ch][i][0]
			right = gw_independent_breakpoints[ch][i][1]
			IND_PARTITION.write(str(ch)+"\t"+ str(variant_positions[ch][left]) + "\t" + str(variant_positions[ch][right]) + "\t"+str(left)+"\t"+str(right)+"\n")


for ch in chromosome:
	snp_frequency[ch],independent_block_haplotype_frequency[ch] = HG.haplotype_frequency_estimation(ch,gw_independent_breakpoints[ch],hap_matrix_d1[ch],hap_matrix_d2[ch],ref[ch],alt[ch],variant_positions[ch],variant_names[ch],args.bam,args.reference,args.output)

OUTPUT_FREQ_NAMES = open(args.output+"_snp_frequency.txt","w")
for ch in chromosome:
	for i in range(len(snp_frequency[ch])):
		print("%s\t%d\t%f" %(ch,variant_positions[ch][i],snp_frequency[ch][i]),file = OUTPUT_FREQ_NAMES)



## calculate ecotype frequency from the haplotype frequency information

# for ch in chromosome:
# 	ecotype_frequency[ch] = HG.ecotype_frequency_estimation(gw_independent_breakpoints[ch],independent_block_haplotype_frequency[ch],hap_matrix_d1[ch],hap_matrix_d2[ch],variant_names[ch],variant_positions[ch])

# ecotype_frequency_avg = np.zeros(len(ind_names))
# for ch in chromosome:
# 	ecotype_frequency_avg += ecotype_frequency[ch]

# ecotype_frequency_avg = ecotype_frequency_avg / len(chromosome)

# with open (args.output+"_ecotype_frequency.txt","w") as ECOTYPE_FREQ:
# 	for i in range(len(ecotype_frequency_avg)):
# 		print("%s\t%f" %(ind_names[i],ecotype_frequency_avg[i]),file = ECOTYPE_FREQ)

## calculate ecotype frequency from the haplotype frequency information using selected blocks

for ch in chromosome:
	ecotype_frequency_selected[ch],length[ch] = HG.ecotype_frequency_estimation_selected(gw_independent_breakpoints[ch],independent_block_haplotype_frequency[ch],hap_matrix_d1[ch],hap_matrix_d2[ch],variant_names[ch],variant_positions[ch])
	# np.savetxt(args.output+"_"+ch+"_ecotype_frequency_selected.txt",ecotype_frequency_selected[ch],delimiter=',')
	# np.savetxt(args.output+"_"+ch+"_ecotype_frequency_selected_length.txt",length[ch],delimiter=',')

length_all = []
ecotype_frequency_all = np.zeros((0,len(ind_names)))
for ch in chromosome:
	length_all.extend(length[ch])
	ecotype_frequency_all = np.concatenate((ecotype_frequency_all,ecotype_frequency_selected[ch]),axis=0)


quantile_9 = np.quantile(length_all,0.8)


index_9 = np.asarray(np.where(length_all >= quantile_9)[0],dtype="i")
weights = np.zeros(len(length_all))
for m in index_9:
	weights[m] = length_all[m]

ecotype_frequency_avg = np.average(ecotype_frequency_all,axis=0,weights = weights)
with open (args.output+"_ecotype_frequency.txt","w") as ECOTYPE_FREQ:
	for i in range(len(ecotype_frequency_avg)):
		print("%s\t%f" %(ind_names[i],ecotype_frequency_avg[i]),file = ECOTYPE_FREQ)


if args.haplotype_frequency == True:
	if args.block == "bigld":
		with open(args.output+"_fine_genomewide_partition.txt", "w") as GW_BREAKPOITS:
			for ch in chromosome:
				common_fine_breakpoints = HG.BigLD_partition(DIR,IndepLD_common_breakpoints_index[ch],common_geno_matrix[ch],common_variant_names[ch],common_variant_positions[ch],args.CLQcut,args.output)
				gw_fine_breakpoints[ch] = HG.convert_fine_genomewide_breakpoints(common_fine_breakpoints,common_allele_index[ch],len(variant_names[ch]),gw_independent_breakpoints[ch])
				for i in range(len(gw_fine_breakpoints[ch])):
					GW_BREAKPOITS.write(str(ch+"\t"+'\t'.join(map(str, gw_fine_breakpoints[ch][i]))+"\n"))	
	else:
		gw_fine_breakpoints = HG.custom_fine_partition(args.block)

	for ch in chromosome:
		block_partitions[ch] = {}
		fine_block_array = np.asarray(gw_fine_breakpoints[ch])
		for i in range(len(gw_independent_breakpoints[ch])):
			left = int(np.where(fine_block_array[:,0] == gw_independent_breakpoints[ch][i][0])[0])
			right = int(np.where(fine_block_array[:,1] == gw_independent_breakpoints[ch][i][1])[0])
		#	left = HG.find_nearest(fine_block_array[:,0],gw_independent_breakpoints[ch][i][0])
		#	right = HG.find_nearest(fine_block_array[:,1],gw_independent_breakpoints[ch][i][1])
			block_partitions[ch][i] = gw_fine_breakpoints[ch][left:right+1]

	with open (args.output+"_unique_haplotype_frequency.txt","w") as UNIQ_HAPLOTYPE_FREQ:
		with open (args.output+"_clustered_haplotype_frequency.txt","w") as CLUSTER_HAPLOTYPE_FREQ:
			for ch in chromosome:
				fine_haplotype_frequency,vcf_haplotype_frequency,fine_haplotype_names,haplotype_cluster_frequency,haplotype_cluster_names,vcf_haplotype_cluster_frequency = HG.fine_haplotype_frequency_calculation(block_partitions[ch],gw_independent_breakpoints[ch],independent_block_haplotype_frequency[ch],gw_fine_breakpoints[ch],hap_matrix_d1[ch],hap_matrix_d2[ch],variant_names[ch],variant_positions[ch],snp_frequency[ch],ch)
				
				for i in range(len(fine_haplotype_frequency)):
					print("%s\t%f\t%f" %(fine_haplotype_names[i],fine_haplotype_frequency[i],vcf_haplotype_frequency[i]),file = UNIQ_HAPLOTYPE_FREQ)

				for i in range(len(haplotype_cluster_frequency)):
					print("%s\t%f\t%f" %(haplotype_cluster_names[i],haplotype_cluster_frequency[i],vcf_haplotype_cluster_frequency[i]),file = CLUSTER_HAPLOTYPE_FREQ)


'''

if args.block == "bigld":

	with open(args.output+"_genomewide_partition.txt", "w") as GW_BREAKPOITS:
		for ch in chromosome:
			common_fine_breakpoints = HG.BigLD_partition(DIR,IndepLD_common_breakpoints_index[ch],common_geno_matrix[ch],common_variant_names[ch],common_variant_positions[ch],args.CLQcut,args.output)
			gw_breakpoints[ch] = HG.convert_genomewide_breakpoints(common_fine_breakpoints,common_allele_index[ch],len(variant_names))
			for i in range(len(gw_breakpoints[ch])):
				GW_BREAKPOITS.write(str(ch+"\t"+'\t'.join(map(str, gw_breakpoints[ch][i]))+"\n"))
elif args.block == None:
	sys.exit("please specify the block partition algoritms or the file contains the partitions with the right format.")

elif args.block == "Naive":
	for ch in chromosome:
		common_fine_breakpoints = IndepLD_common_breakpoints_index[ch]
		gw_breakpoints[ch] = HG.convert_genomewide_breakpoints(common_fine_breakpoints,common_allele_index[ch],len(variant_names))
		print(gw_breakpoints[ch])

else:
	gw_breakpoints = HG.custom_fine_partition(args.block)


######################## generate the haplotype block design matrix from the hap matrix
print("start haplotype design matrix generation.")



for ch in chromosome:
	print(gw_breakpoints[ch])
	snp_frequency[ch] = HG.haplotype_frequency_estimation(ch,gw_breakpoints[ch],hap_matrix_d1[ch],hap_matrix_d2[ch],ref[ch],alt[ch],variant_positions[ch],variant_names[ch],args.bam,args.reference)


OUTPUT_FREQ_NAMES = open(args.output+"_snp_frequency.txt","w")
for ch in chromosome:
	for i in range(len(snp_frequency[ch])):
		print("%s\t%d\t%f" %(ch,variant_positions[ch][i],snp_frequency[ch][i]),file = OUTPUT_FREQ_NAMES)

# HaploBlock_matrix = mp.Manager().dict()
# haplotype_block_name = mp.Manager().dict()
# haplotype_marker_name = mp.Manager().dict()

# processes = []

# # for ch in chromosome:
# # 	p = mp.Process(target = uf.haplotype_frequency_estimation, args=(ch,r,hap_matrix_d1,hap_matrix_d2,geno_matrix,variant_names,variant_positions,fine_breakpoints,HaploBlock_matrix,haplotype_block_name,haplotype_marker_name,args.clustering))
# # 	processes.append(p)
# # 	p.start()

# # for process in processes:
# # 	process.join()

columns = []
H = pd.DataFrame(index=range(r),columns=columns)
	
OUTPUT_BLOCK_NAMES = open(args.output+"_block_names.txt","w")
for ch in chromosome:
	for i in range(len(haplotype_block_name[ch])):
		print("%s" %(haplotype_block_name[ch][i]),file = OUTPUT_BLOCK_NAMES)

OUTPUT_MARKER_NAMES = open(args.output+"_marker_names.txt","w")
for ch in chromosome:
	for i in range(len(haplotype_marker_name[ch])):
		print("%s" %(haplotype_marker_name[ch][i]),file = OUTPUT_MARKER_NAMES)




for ch in chromosome:
	for key in HaploBlock_matrix[ch]:
		H = pd.concat([H,HaploBlock_matrix[ch][key]],axis=1,ignore_index=False)


hap_names = H.columns.values.tolist()
bimbam = uf.format_bimbam(np.transpose(H.values),hap_names)
bimbam.to_csv(args.output+".bimbam",index=False,sep=" ",header=None)

H.to_csv(args.output+"_haplotypeDM.txt",sep="\t",header=True,index=False)


print("finish constructing haplotype design matrix")
'''

