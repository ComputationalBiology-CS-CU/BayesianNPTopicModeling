## extract the following from the gene annotation file:
##	1. gene_tss.txt
##	2. gene_xymt.txt




if __name__ == "__main__":


	##=============================
	##========== tss ==============
	##=============================
	## get the chr and tss for all genes (including X, Y and MT genes)
	file = open("../data_raw/gencode.v19.genes.v6p_model.patched_contigs.gtf", 'r')
	file1 = open("../data_raw/gene_tss_gencode.v19.v6p.txt", 'w')
	file.readline()
	file.readline()
	file.readline()
	file.readline()
	file.readline()
	file.readline()		## extra line for this version

	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split()
		chr = line[0]
		tss = line[3]
		gene = line[9][1: -2]

		type = line[2]
		if type == 'gene':
			file1.write(gene + '\t' + chr + '\t' + tss + '\n')

	file.close()
	file1.close()



	##=======================================
	##========== X, Y, MT list ==============
	##=======================================
	## get the list of all X, Y, MT genes
	file = open("../data_raw/gene_tss_gencode.v19.v6p.txt", 'r')
	file1 = open("../data_raw/gene_xymt_gencode.v19.v6p.txt", 'w')
	while 1:
		line = (file.readline()).strip()
		if not line:
			break

		line = line.split('\t')
		gene = line[0]
		chr = line[1]
		tss = line[2]

		if chr == 'X' or chr == 'Y' or chr == 'MT':
			file1.write(gene + '\n')
	file.close()
	file1.close()




