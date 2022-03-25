# Import Function
source("./exclude_overlapped.r")
# total GWAS file, overlapped sample GWAS file 
# output overlap_excluded GWAS file path, sample sizes and phenotypes
# If column names are different, please specify in following order.
# effect size, standard error, MAF, Chromosome, Position, Allele1, Allele2 
# sample code running with files attached
source("./exclude_overlapped.r")
exclude_overlap('input_all.txt', 'input_overlap.txt', 'output.txt',
	col_all = c("beta","SE","maf","chr","pos","A1","A2"),
	col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"),
	n_all = 316253, n_ov = 81538,
	phenotype = "binary")
# For additional information, please refer to annotations on the provided code.