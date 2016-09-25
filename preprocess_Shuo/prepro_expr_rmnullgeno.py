## input:
##	1. GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct
##	2. GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_sample_IDs.txt
## output:
##	1. GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_1_genotype
##
## input:
##	1. phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type
##	2. GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct
## output:
##	1. phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count






##=====================
##==== libraries
##=====================
import math
import numpy as np





##=====================
##==== global variables
##=====================
individual_rep = {}		# hashing all the individuals with genotype information
sample_tissue_map = {}		# mapping all the samples into their tissue types





##==================
##==== sub-routines
##==================
# get the "xxx-yyy" from "xxx-yyy-zzz-aaa-qqq", which is defined as the individual ID of the GTEx samples
def get_individual_id(s):
	## naively find the second '-'
	id = ''
	count = 0
	for i in range(len(s)):
		if s[i] == '-':
			count += 1

		if count == 2:
			break

		id += s[i]

	return id






if __name__ == "__main__":




	""" don't have to do the following if we have that info file
	##========================================================================================
	##==== get the sample-tissue mapping for all the samples
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type
	##========================================================================================
	file = open("../data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt", 'r')
	file1 = open("../data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'w')
	count = 0
	while 1:
		line = file.readline()[:-1]	# can't strip all "\t\t\t\t", as they are place holder
		count += 1
		if count <= 11:  ## 11
			continue

		if not line:
			break

		line = line.split('\t')
		sample = line[1]
		tissue1 = line[12]
		tissue2 = line[14]

		file1.write(sample + '\t' + tissue1 + '\t' + tissue2 + '\n')

	file.close()
	file1.close()
	"""






	##========================================================================================
	##==== remove samples that have no genotype information
	##==== target: ******_gene_rpkm.gct_1_genotype
	##========================================================================================
	file = open("../data_raw/GTEx_Analysis_2015-01-12_OMNI_2.5M_5M_450Indiv_sample_IDs.txt", 'r')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		id = get_individual_id(line)
		individual_rep[id] = 1

	file.close()


	##
	file = open("../data_raw/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct", 'r')
	file1 = open("../data_raw/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_1_genotype", 'w')
	file.readline()
	file.readline()

	# get the effective list for selection
	index_map = {}
	line = (file.readline()).strip()
	line = line.split('\t')
	file1.write(line[0] + '\t')
	for i in range(2, len(line)):		## NOTE: start from 2
		sample = line[i]
		individual = get_individual_id(sample)
		if individual not in individual_rep:
			continue
		else:
			file1.write(line[i] + '\t')
			index_map[i] = 1
	file1.write('\n')

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		file1.write(line[0] + '\t')
		for i in range(2, len(line)):
			if i not in index_map:
				continue
			else:
				file1.write(line[i] + '\t')
		file1.write('\n')

	file.close()
	file1.close()







	##===================================================================================================
	##==== counting samples in each tissue
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count
	##==== target: phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_#size
	##===================================================================================================
	# sample_tissue_map
	file = open("../data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_sample_tissue_type", 'r')
	sample_tissue_map = {}
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')

		if len(line) < 3:
			print line
			continue

		sample = line[0]
		tissue = line[2]

		sample_tissue_map[sample] = tissue
	file.close()

	# counting
	file = open("../data_raw/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_1_genotype", 'r')
	sample_list = (((file.readline()).strip()).split('\t'))[1:]
	file.close()

	print "there are",
	print len(sample_list),
	print "different samples from the rpkm file."

	counting = {}
	for sample in sample_list:
		tissue = sample_tissue_map[sample]
		if tissue not in counting:
			counting[tissue] = 1
		else:
			counting[tissue] += 1

	print "they are distributed among",
	print len(counting),
	print "different tissue types."

	file = open("../data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count", 'w')
	for tissue in counting:
		file.write(tissue + '\t' + str(counting[tissue]) + '\n')
	file.close()



	""" don't have to do the following now
	#==== filtering the counts
	filter_list = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
	for i in range(len(filter_list)):
		filter = filter_list[i]
		file = open("../data_raw/phs000424.v6.pht002743.v6.p1.c1.GTEx_Sample_Attributes.GRU.txt_tissue_type_count_" + str(filter), 'w')
		count = 0
		for tissue in counting:
			if counting[tissue] >= filter:
				count += 1
				file.write(tissue + '\t' + str(counting[tissue]) + '\n')
		print "# of tissue with sample size >= " + str(filter) + ":",
		print count
		file.close()
	"""




