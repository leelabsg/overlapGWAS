################################################################## Required Packages
library(dplyr)
library(data.table)
library(stringr)

################################################################## Function Excluding Overlap Samples 
# Arguments
# 1. Input files consist of total gwas summary and overlapped sample GWAS summary (file_all, file_ov)
# 2. Specify output file path (output_file) 
# 3. Specify column names indicating sample size, effect size, standard error, MAF, Chromosome, Position, Allele1, Allele2 of input, respectively (col1, col2)
#### Or Alternative column names indicating sample size, effect size, standard error, MAF, chromosomal position
# 4. Sample size of each GWAS
# 5. Phenotype can be either binary or continiuous
exclude_overlap = function
(
  file_all,file_ov,output_file, # file paths of Total GWAS summary, overlapped sample GWAS summary, and output 
  col_all = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"), # Column names of Total GWAS 
  col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"), # Column names of overlapped sample GWAS 
  n_all, n_ov, # sample size of Total GWAS and overlapped sample GWAS
  phenotype # Phenotype if binary : "binary", continuous : "continuous"
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
    ss_all = ss_all %>% mutate(CHR_POS = paste0(CHR,"_",POS,"_",Allele1,"_",Allele2))
  } else if (length(col_all)==4){
    colnames(ss_all) = c("BETA","SE","AF_Allele2","CHR_POS")
  }

  if (length(col_ov)==7){
    colnames(ss_ov)= c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2")
    ss_ov = ss_ov %>% mutate(CHR_POS = paste0(CHR,"_",POS,"_",Allele1,"_",Allele2))
  } else if (length(col_ov)==4){
    colnames(ss_ov)= c("BETA","SE","AF_Allele2","CHR_POS")
  }

  # Subset only intersection of chromosal position of each GWAS
  pos1 = ss_all %>% select(CHR_POS)
  pos2 = ss_ov %>% select(CHR_POS)
  inner_pos = intersect(pos1,pos2)
  ss_all = ss_all %>% filter(CHR_POS%in%inner_pos$CHR_POS)
  ss_all = ss_all[!duplicated(ss_all$CHR_POS),]
  ss_ov = ss_ov %>% filter(CHR_POS%in%inner_pos$CHR_POS)
  ss_ov = ss_ov[!duplicated(ss_ov$CHR_POS),]

  # Sample Size Calculation
  est_all = ss_all$BETA; se_all = ss_all$SE ; 
  est_ov = ss_ov$BETA; se_ov =ss_ov$SE ; 
  n_exc = n_all-n_ov

  # Get Minor Allele Frequency
  maf = ss_all$AF_Allele2
  
  # Inverse Variance Weighting
  se_exc = sqrt(1/((1/se_all^2)-(1/se_ov^2)))
  se_exc[is.na(se_exc)] = max(se_all[is.na(se_exc)],se_ov[is.na(se_exc)])
  est_exc = (est_all/se_all^2-est_ov/se_ov^2)*se_exc^2

  # Partitioning GWAS
  if (phenotype=="binary"){
    xty_all = 0.5*n_all*maf-est_all*n_all*maf/4
    xty_ov = 0.5*n_ov*maf-est_ov*n_ov*maf/4
    xty_ex2  = xty_all - xty_ov
    est_exc2 = (-xty_ex2+0.5*(n_all-n_ov)*maf)*4/((n_all-n_ov)*maf)
    se_exc2 = sqrt((se_all^2*(n_all*maf)^2-se_ov^2*(n_ov*maf)^2)/(n_exc*maf)^2)
    se_exc2[is.na(se_exc2)] = max(se_all[is.na(se_exc2)],se_ov[is.na(se_exc2)])
  }
  if (phenotype=="continuous"){
    xty_all = est_all*n_all*maf
    xty_ov = est_ov*n_ov*maf
    xty_ex2  = xty_all - xty_ov
    est_exc2 = xty_ex2/((n_exc)*maf)
    se_exc2 = sqrt((se_all^2*(n_all*maf)^2-se_ov^2*(n_ov*maf)^2)/(n_exc*maf)^2)
    se_exc2[is.na(se_exc2)] = max(se_all[is.na(se_exc2)],se_ov[is.na(se_exc2)])
  }

  # Generate Output GWAS Dataframe
  if (length(col_all)==7){
    new_ss = data.frame(
      CHR = ss_all$CHR, POS = ss_all$POS,
      Allele1 = ss_all$Allele1, Allele2 = ss_all$Allele2, AF_Allele = maf, N = n_exc,
      BETA_IVW = est_exc,SE_IVW = se_exc,P_IVW = exp(pchisq((est_exc^2)/se_exc^2,df=1,lower.tail=F,log.p=T)),
      BETA_Part = est_exc2,SE_Part = se_exc2,P_Part = exp(pchisq((est_exc2^2)/se_exc2^2,df=1,lower.tail=F,log.p=T))
      )
  } else if (length(col_all)==4){
    new_ss = data.frame(
      SNPID = ss_all$CHR_POS, AF_Allele = maf, N = n_exc,
      BETA_IVW = est_exc,SE_IVW = se_exc,P_IVW = exp(pchisq((est_exc^2)/se_exc^2,df=1,lower.tail=F,log.p=T)),
      BETA_Part = est_exc2,SE_Part = se_exc2,P_Part = exp(pchisq((est_exc2^2)/se_exc2^2,df=1,lower.tail=F,log.p=T))
      )
  }
  # Check header and size of new overlap excluded GWAS 
  print("New Summary")
  print(head(new_ss))
  print("Dimension of new summary")
  print(dim(new_ss))

  # Write output
  write.table(new_ss, output_file, col.names = T, row.names = F, quote = F)
}
# Sample usage
# source('/media/leelabsg-storage0/seokho/overlap/code/exclude_overlapped.r')
# exclude_overlap('input_all.txt', 'input_overlap.txt', 'output.txt',col_all = c("beta","SE","maf","chr","pos","A1","A2"),col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2"),n_all = 316253, n_ov = 81538,phenotype = "binary")
