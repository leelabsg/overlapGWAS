################################################################## Required Packages
library(data.table)
library(dplyr)
library(stringr)
#library(biomaRt)
################################################################## Function Excluding Overlap Samples 
# Arguments
# 1. Input files consist of total gwas summary and overlapped sample GWAS summary (file_all, file_ov)
# 2. Specify output file path (output_file) 
# 3. Specify column names indicating sample size, effect size, standard error, MAF, Chromosome, Position, Allele1, Allele2 of input, respectively (col1, col2)
#### Or Alternative column names indicating sample size, effect size, standard error, MAF, chromosomal position
# 4. Sample size of each GWAS
# 5. Phenotype can be either binary or continiuous
# 6. Whether to drop markers that has standard error NA's (otherwise use max standard error between overlapped and all)
exclude_overlap = function
(
  file_all,file_ov,output_file, # file paths of Total GWAS summary, overlapped sample GWAS summary, and output 
  col_all = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"), # Column names of Total GWAS 
  col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"), # Column names of overlapped sample GWAS 
  n_all, n_ov, # sample size of Total GWAS and overlapped sample GWAS
  phenotype, # Phenotype if binary : "binary", continuous : "continuous"
  dropna=T # Drop markers that has standard error NA's (otherwise use max standard error between overlapped and all)
) 
{
  # File read
  ss_all = fread(file_all)
  ss_ov = fread(file_ov)

  # Column selection
  ss_all = ss_all%>%select(all_of(col_all))
  ss_ov = ss_ov%>%select(all_of(col_ov))

  # Build chromosomal position and its allele to secure uniqueness
  if (length(col_all)==7){
    colnames(ss_all) = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2")
  } else if (length(col_all)==4){
    colnames(ss_all) = c("BETA","SE","AF_Allele2","CHR_POS")
  }

  if (length(col_ov)==7){
    colnames(ss_ov)= c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2")
  } else if (length(col_ov)==4){
    colnames(ss_ov)= c("BETA","SE","AF_Allele2","CHR_POS")
  }

  # Subset only intersection of chromosal position of each GWAS
  if (length(col_all)==7){ss_tmp = inner_join(ss_all,ss_ov, by = c("CHR","POS","Allele1","Allele2"),suffix = c("_all", "_ov"))}
  if (length(col_all)==4){ss_tmp = inner_join(ss_all,ss_ov, by = c("CHR_POS"),suffix = c("_all", "_ov"))}
  # Sample Size Calculation
  est_all = ss_tmp$BETA_all; se_all = ss_tmp$SE_all ; 
  est_ov = ss_tmp$BETA_ov; se_ov =ss_tmp$SE_ov ; 
  n_exc = n_all-n_ov

  # Get Minor Allele Frequency
  maf = ss_tmp$AF_Allele2_all
  
  # Inverse Variance Weighting
  ss_tmp = ss_tmp%>%mutate(SE_IVW = sqrt(1/((1/SE_all^2)-(1/SE_ov^2))))
  ss_tmp = ss_tmp%>%mutate(BETA_IVW = (BETA_all/SE_all^2-BETA_ov/SE_ov^2)*SE_IVW^2)
  ss_tmp = ss_tmp%>%mutate(P_IVW = exp(pchisq((BETA_IVW^2)/SE_IVW^2,df=1,lower.tail=F,log.p=T)))
  # Partitioning GWAS
  if (phenotype=="binary"){
    xty_all = 0.5*n_all*maf-est_all*n_all*maf/4
    xty_ov = 0.5*n_ov*maf-est_ov*n_ov*maf/4
    xty_ex2  = xty_all - xty_ov
    est_exc2 = (-xty_ex2+0.5*(n_all-n_ov)*maf)*4/((n_all-n_ov)*maf)
    se_exc2 = sqrt((se_all^2*(n_all*maf)^2-se_ov^2*(n_ov*maf)^2)/(n_exc*maf)^2)
  }
  if (phenotype=="continuous"){
    xty_all = est_all*n_all*maf
    xty_ov = est_ov*n_ov*maf
    xty_ex2  = xty_all - xty_ov
    est_exc2 = xty_ex2/((n_exc)*maf)
    se_exc2 = sqrt((se_all^2*(n_all*maf)^2-se_ov^2*(n_ov*maf)^2)/(n_exc*maf)^2)
  }

  # Generate Output GWAS Dataframe
  if (length(col_all)==7){
    new_ss = data.frame(
      CHR = ss_tmp$CHR, POS = ss_tmp$POS,
      Allele1 = ss_tmp$Allele1, Allele2 = ss_tmp$Allele2, AF_Allele = maf, N = n_exc,
      BETA_IVW = ss_tmp$BETA_IVW,SE_IVW = ss_tmp$SE_IVW,P_IVW = ss_tmp$P_IVW,
      BETA_Part = est_exc2,SE_Part = se_exc2,P_Part = exp(pchisq((est_exc2^2)/se_exc2^2,df=1,lower.tail=F,log.p=T))
      )
  } else if (length(col_all)==4){
    new_ss = data.frame(
      SNPID = ss_tmp$CHR_POS, AF_Allele = maf, N = n_exc,
      BETA_IVW = ss_tmp$BETA_exc,SE_IVW = ss_tmp$SE_exc,
      BETA_Part = est_exc2,SE_Part = se_exc2,P_Part = exp(pchisq((est_exc2^2)/se_exc2^2,df=1,lower.tail=F,log.p=T))
      )
  }
  
  if (!dropna){
  new_ss$SE_IVW[is.na(new_ss$SE_IVW)] = max(ss_tmp$SE_all[is.na(new_ss$SE_IVW)],ss_tmp$SE_ov[is.na(new_ss$SE_IVW)])
  new_ss$SE_Part[is.na(new_ss$SE_Part)] = max(ss_tmp$SE_all[is.na(new_ss$SE_Part)],ss_tmp$SE_ov[is.na(new_ss$SE_Part)])
  }
  new_ss = new_ss %>% filter(!is.na(SE_IVW) & !is.na(SE_Part))
  # Check header and size of new overlap excluded GWAS 
  print("New Summary")
  print(head(new_ss))
  print("Dimension of new summary")
  print(dim(new_ss))

  # Write output
  fwrite(new_ss, output_file, col.names = T, row.names = F, quote = F,sep = " ",na = NA)
  #write.table(new_ss, output_file, col.names = T, row.names = F, quote = F)
}
# Sample usage
# source('/media/leelabsg-storage0/seokho/overlap/code/exclude_overlapped.r')
# exclude_overlap('input_all.txt', 'input_overlap.txt', 'output.txt',col_all = c("beta","SE","maf","chr","pos","A1","A2"),col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"),n_all = 316253, n_ov = 81538,phenotype = "binary",dropna = T)
