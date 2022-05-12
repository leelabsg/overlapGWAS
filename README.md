# overlapGWAS
In case when validating sample and GWAS building sample have overlaps, overlapGWAS intends to exclude the effect of overlapped samples from original GWAS.
![image](https://user-images.githubusercontent.com/22064612/160033954-88889aaf-b282-487f-bdb4-a2f305d92fc1.png)

## Import Function
```
source("./exclude_overlapped.r")
```
## Input Specifications
### 1. total GWAS file, overlapped sample GWAS file 
### 2. output overlap_excluded GWAS file path, sample sizes and phenotypes
### 3. If column names are different, please specify in following order.
### 4. effect size, standard error, MAF, Chromosome, Position, Allele1, Allele2 
### 5. sample code running with files attached
```
source("./exclude_overlapped.r")
exclude_overlap('input_all.txt', 'input_overlap.txt', 'output.txt',
	col_all = c("beta","SE","maf","chr","pos","A1","A2"),
	col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"),
	n_all = 316253, n_ov = 81538,
	phenotype = "binary")
```
For additional information, please refer to annotations on the provided code.


### 051222 update
1. Changed the joining of all samples and overlapped samples using inner_join()
2. Skipped the CHR_POS generation when given CHR, POS, allele1 and 2 to reduce time
3. Changed write.table to fwrite for faster write
4. Added option whether to drop markers that have standard error to NA's
