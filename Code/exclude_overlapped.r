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
# 5. Either Effective sample sizes / case control size needs to be provided
# 6. Phenotype can be either binary or continiuous
# 7. Whether to drop markers that has standard error NA's (otherwise use max standard error between overlapped and all)

# Output
# 0. Summary statistics of All / Overlapped/ IVW/ Partitioned / Z-scoring 
# 1. Effect size, Standard Error, Sample size, P-value
# 2. xty from partitioning Method
# 3. z scores and effective sample sizes

# Transformation of Pvalue into Z-score
PtoZ = function(P,bet){
return(abs(qnorm(P/2))*sign(bet))
}
exclude_overlap = function
(
  file_all,file_ov,output_file, # file paths of Total GWAS summary, overlapped sample GWAS summary, and output 
  col_all = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2","P"), # Column names of Total GWAS 
  col_ov = c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2","P"), # Column names of overlapped sample GWAS 
  n_all, n_ov, # sample size of Total GWAS and overlapped sample GWAS
  eff_all=NULL, eff_ov=NULL, # effective sample size
  n_case_all = NULL, n_ctrl_all =NULL, n_case_ov = NULL, n_ctrl_ov = NULL, # effective sample size with case/control size using Metal 
  phenotype, # Phenotype if binary : "binary", continuous : "continuous"
  dropna=T # Drop markers that has standard error NA's (otherwise use max standard error between overlapped and all)
) 
{
  # File read
  ss_all = fread(file_all,fill=TRUE)
  ss_ov = fread(file_ov,fill=TRUE)
  if (!col_all[3]%in%colnames(ss_ov)){
    ss_ov = ss_ov %>% mutate(af = (control_af*67127+case_af*5083)/(67127+5083))
  }
  # Column selection
  ss_all = ss_all%>%select(all_of(col_all))
  ss_ov = ss_ov%>%select(all_of(col_ov))
  
  # Metal
    if(!is.null(n_case_all) & !is.null(n_ctrl_all)){
        eff_all = 4/(1/n_case_all+1/n_ctrl_all)
    }
    if(!is.null(n_case_ov) & !is.null(n_ctrl_ov)){
        eff_ov = 4/(1/n_case_ov+1/n_ctrl_ov)
    }

  # Build chromosomal position and its allele to secure uniqueness
  colnames(ss_all) = colnames(ss_ov)= c("BETA","SE","AF_Allele2","CHR","POS","Allele1","Allele2","P")

  # Subset only intersection of chromosal position of each GWAS
  ss_tmp = inner_join(ss_all,ss_ov, by = c("CHR","POS","Allele1","Allele2"),suffix = c("_all", "_ov"))
  ss_tmp2 = inner_join(ss_all,ss_ov, by = c("CHR","POS","Allele1"="Allele2","Allele2"="Allele1"),suffix = c("_all", "_ov"))
  ss_tmp2 = ss_tmp2%>%mutate(BETA_ov=(-1)*BETA_ov)
  ss_tmp = rbind(ss_tmp,ss_tmp2)
  print("Number of variant of all, overlapped, subtracted samples")
  print(c(nrow(ss_all),nrow(ss_ov),nrow(ss_tmp)))

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

  # Reversed z-score
  if (!is.null(eff_all)&!is.null(eff_ov)) {
    w_all = sqrt(eff_all); w_ov = sqrt(eff_ov); w_RZ = sqrt(eff_all-eff_ov)
    z_all = PtoZ(ss_tmp$P_all,ss_tmp$BETA_all); z_ov = PtoZ(ss_tmp$P_ov,ss_tmp$BETA_ov)
    z_RZ = z_all*w_all - z_ov*w_ov ; z_RZ = z_RZ/w_RZ
    BETA_RZ =  z_RZ / sqrt(2*maf*(1-maf)*((eff_all-eff_ov) + z_RZ^2)) 
    sign_RZ = sign(z_RZ); P_RZ = exp(pnorm(abs(z_RZ),lower.tail=F,log.p=T))*2 
  }
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
    new_ss = data.frame(
        CHR = ss_tmp$CHR, POS = ss_tmp$POS,
        Allele1 = ss_tmp$Allele1, Allele2 = ss_tmp$Allele2, AF_Allele = maf, N = n_exc, N_all = n_all, N_ov = n_ov, 
        BETA_IVW = ss_tmp$BETA_IVW,SE_IVW = ss_tmp$SE_IVW,P_IVW = ss_tmp$P_IVW,
        BETA_Part = est_exc2,SE_Part = se_exc2,P_Part = exp(pchisq((est_exc2^2)/se_exc2^2,df=1,lower.tail=F,log.p=T)),
        BETA_all = ss_tmp$BETA_all,P_all = ss_tmp$P_all,SE_all = ss_tmp$SE_all, 
        BETA_ov = ss_tmp$BETA_ov,P_ov = ss_tmp$P_ov, SE_ov = ss_tmp$SE_ov,
        xty_all = xty_all, xty_ov = xty_ov
        )
    if (!is.null(eff_all)&!is.null(eff_ov)){
        new_ss = new_ss %>% mutate(sign_RZ = sign_RZ, P_RZ = P_RZ, z_RZ = z_RZ, z_all = z_all, z_ov = z_ov, BETA_RZ = BETA_RZ, ESS_all = eff_all, ESS_ov = eff_ov, ESS_RZ = eff_all-eff_ov)
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

# When case/control sample size is not provided
# exclude_overlap("all.txt","ov.txt","out_n3.txt",n_all = 433540,n_ov = 72210,phenotype = "binary",dropna=F)
#############################################################################################################################################################
# Using case/control sample size and Effective sample size at the same time
# In this case, ESS of all samples are provided and case/control size of overlapped samples are provided
if(FALSE){
exclude_overlap("all_new.txt","/media/leelabsg-storage0/seokho/UKBB/overlapGWAS/Agen_T2D/KoGES_DM.tsv.gz","out_new_0909.txt",n_all = 433540,n_ov = 56918,
col_all = c("Beta","SE","EAF","Chr","Pos","Allele1","Allele2","P"),
eff_all =  211793, n_case_ov = 5083, n_ctrl_ov = 67127,
col_ov = c("beta","sebeta","af","chrom","pos","ref","alt","pval"),phenotype = "binary",dropna=F)
} 
