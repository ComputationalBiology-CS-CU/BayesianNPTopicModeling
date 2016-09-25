## I will pick up the SNPs according to a pre-specified SNP list
## that means I preprocess (pruning) the SNP based on dosage info

## the work is under "/ifs/scratch/c2b2/ip_lab/sy2515/GTEx/data.v.6/47024/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/GenotypeFiles/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/genotype_snp"

## dependency:
##	1. genotype_vcf
##	2. list_snp

## where to save data:
##	1. genotype_450_geno_matrix_qc
##	2. genotype_450_geno_matrix_qc/chrX
##	3. genotype_450_geno_matrix_qc/chrX/SNP_info.txt
##	4. genotype_450_geno_matrix_qc/chrX/SNP_geno_"individualID".txt




import numpy as np
import re





###===============
##==== sub-routines
##===============
# get the "xxx-yyy" from "xxx-yyy-zzz-aaa-qqq", which is defined as the individual ID of the GTEx samples
pattern_indiv = re.compile(r'^(\w)+([\-])(\w)+')
def get_individual_id(s):
	match = pattern_indiv.match(s)
	if match:
		return match.group()
	else:
		print "!!! no individual ID is found..."
		return ""






if __name__ == "__main__":





	##==================================================
	##==== extract genotype in a per chromosome fashion
	##==================================================
	for i in range(22):
		chr = i+1
		print "working on chr#",
		print chr


		##== extract all snp geno
		file = open("./genotype_vcf/chr" + str(chr) + ".txt", 'r')
		repo_snp_geno = {}
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split('\t')
			snp = line[2]

			if snp in repo_snp_geno:			## duplicated SNP names, remove the later ones
				continue

			repo_snp_geno[snp] = []

			line = line[9:]
			for pos in range(len(line)):
				# get MAP genotype for this individual
				item = line[pos]
				item = item.split(':')
				item = item[1]
				item = item.split(',')
				prob0 = float(item[0])
				prob1 = float(item[1])
				prob2 = float(item[2])
				if (prob0 >= prob1) and (prob0 >= prob2):
					genotype = 0
					repo_snp_geno[snp].append(genotype)
					continue
				if (prob1 >= prob0) and (prob1 >= prob2):
					genotype = 1
					repo_snp_geno[snp].append(genotype)
					continue
				if (prob2 >= prob0) and (prob2 >= prob1):
					genotype = 2
					repo_snp_geno[snp].append(genotype)
					continue
		file.close()



		##== snp list
		list_snp = []
		file = open("./list_snp/list_snp_chr" + str(chr) + ".txt", 'r')
		file1 = open("./genotype_450_geno_matrix_qc/chr" + str(chr) + "/SNP_info.txt", 'w')
		while 1:
			line = (file.readline()).strip()
			if not line:
				break

			line = line.split(' ')		## NOTE: we have the ' ' as the separator here
			snp = line[0]
			pos = line[1]
			list_snp.append(snp)
			file1.write(snp + ' ' + pos + '\n')
		file.close()
		file1.close()




		##== load individual repo
		repo_individual = {}
		map_index_individual = {}
		file = open("./genotype_vcf/header.txt", 'r')
		n_header = 17				## extra lines
		for i in range(n_header):
			file.readline()

		line = (file.readline()).strip()
		line = line.split('\t')
		line = line[9:]				## NOTE: start from pos#9
		file.close()

		for i in range(len(line)):
			sample = line[i]
			individual = get_individual_id(sample)
			map_index_individual[i] = individual
			repo_individual[individual] = []

		print "there are # of individuals (double times):",
		print len(repo_individual),
		print len(map_index_individual)



		##== fill in repo_individual
		for snp in list_snp:
			for pos in range(len(repo_snp_geno[snp])):
				individual = map_index_individual[pos]
				geno = repo_snp_geno[snp][pos]
				repo_individual[individual].append(geno)
		repo_individual = np.array(repo_individual)
		print "repo shape:",
		print repo_individual.shape()


		##== save
		for individual in repo_individual:
			file = open("./genotype_450_geno_matrix_qc/chr" + str(chr) + "/SNP_geno_" + individual + ".txt", 'w')
			for geno in repo_individual[individual]:
				file.write(str(geno) + '\n')
			file.close()






