#Extract the snps based on a given position and window
library(tidyverse)
library(made4)

phenotype = readRDS("../ProcessedData/panicle_count_with_residuals.rds")
phenotype$IID = gsub("PI","PI_",phenotype$genotype)



no_of_PCA = c(60)
n=1
logP_threshold = 10

for(n in 1:length(no_of_PCA)){
  gwas_results = readRDS(paste0("../Results/GWAS_Results/panicle_count/panicle_count_gwas_residuals_PCA",no_of_PCA[n],"PCA.rds"))
  #y = phenotype$residuals_PCA50
  y = phenotype$residuals_PCA60
  chromosome = sapply(rownames(gwas_results),function(x) strsplit(x,"_")[[1]][1])
  bp_position = sapply(rownames(gwas_results),function(x) strsplit(x,"_")[[1]][2])
  gwas_results = data.frame(gwas_results,chromosome,bp_position)#,runMED)
  gwas_results = gwas_results[which(-log10(gwas_results$Pvalue)>=logP_threshold),]
  chrs = as.integer(unique(gwas_results$chromosome))
  for(given_chr in chrs){
    #given_chr = 1
    dirOut=paste0("../Results/Regions_of_interest/panicle_count/PCA",no_of_PCA[n],"/Chr",given_chr)
    if(!dir.exists(dirOut)){
      dir.create(dirOut)
    }
    chr_index = which(gwas_results$chromosome==given_chr)
    tmp_gwas = gwas_results[chr_index,]
    tmp_gwas = tmp_gwas[order(tmp_gwas$Pvalue,decreasing = F),]
    
    for( pos in 1:nrow(tmp_gwas)){
      current_snp = row.names(tmp_gwas)[pos]#$bp_position[pos])
      
      plink_command = paste0("~/Tools/plink --bfile ../ProcessedData/cleaned --snp ",current_snp," --window 10 --recode A --out ../Results/Regions_of_interest/panicle_count/RoI_chr")
      system(plink_command)
      tmp_geno = read.table(paste0("../Results/Regions_of_interest/panicle_count/RoI_chr.raw"),header=T)
      dim(tmp_geno)
      tmp = tmp_geno[,-(1:6)]
      dim(tmp)
      length(which(is.na(tmp)==T))
      tmp <- as.data.frame(apply(tmp, 2, function(col) {
        mean_value <- mean(col, na.rm = TRUE)
        col[is.na(col)] <- round(mean_value)
        return(col)
      }))
      
      
      tmp_geno = data.frame(tmp_geno[,1:6],tmp)
      rm(tmp)
      index = match(phenotype$IID,tmp_geno$IID)
      which(is.na(index)==T)
      current_gwas = data.frame(matrix(1,nrow=ncol(tmp_geno),ncol=2))
      dim(current_gwas)
      rownames(current_gwas)= colnames(tmp_geno)
      
      for(t in 7:ncol(tmp_geno)){
        x=tmp_geno[index,t]
        b=sum((x-mean(x)) * (y-mean(y)))  /  sum((x-mean(x))^2) 
        a=mean(y)-b*mean(x)
        yhat=a+b*x
        # calculate R2, F-value and p-value
        SST = sum((y - mean(y))^2)      #Corrected Sum of Squares Total;  SST = Σi=1n (yi - y)2  This is the sample variance of the y-variable multiplied by n - 1
        SSR = sum((yhat - mean(y))^2)  #Corrected Sum of Squares for Model also called sum of squares for regression
        SSE = sum((y-yhat)^2)   #Sum of Squares for Error also called sum of squares for residuals
        #     
        R2=SSR/SST # r-squared For simple linear regression, R2 is the square of the sample correlation rxy For multiple linear regression with intercept (which includes simple linear regression), it is defined as r2 = SSM / SST
        Fval=SSR /(SSE/(length(y)-2)) # F-value
        Pval=1-pf(Fval,1,(length(y)-2)) # p-value
        current_gwas[t,1] = b
        current_gwas[t,2] = Pval
        
        
      }
      current_gwas = current_gwas[-(1:6),]
      colnames(current_gwas) = c("Effect","P-value")
      gwas_file_name = paste0(dirOut,"/GWAS_",current_snp,".csv")
      write.csv(current_gwas,gwas_file_name)
      
    }
  }
}
  
  





###******************************************************************************************************************************************
# 
# impute_geno=function(x){
#     if(length(which(is.na(x)==T))>0){
#       x[which(is.na(x)==T)] = round(mean(na.omit(x)),0)
#     }
#     return(x)
#   }
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   given_pos = as.integer(tmp_gwas$bp_position[1])
# window_size=15000
# 
# filtered_snps = gwas_results %>%
#   filter(chromosome==given_chr, bp_position<(given_pos+window_size) & bp_position>(given_pos-window_size)  ) %>%
#   arrange(bp_position)
# 
# 
# SNP1 = rownames(filtered_snps)[1]
# SNP2 = rownames(filtered_snps)[nrow(filtered_snps)]
# 
# SNP1 = "9_32597333"
# SNP2 = "9_32626488"
# #Extract the data in PLINK format
# # plink_command = paste0("~/Tools/plink --bfile ProcessedData/cleaned --from ",SNP1," --to ",SNP2," --make-bed --out ProcessedData/filtered_data")
# # system(plink_command)
# 
# 
# #EXTRAC data in VCF format
# plink_command = paste0("~/Tools/plink --bfile ../ProcessedData/cleaned --from ",SNP1," --to ",SNP2," --recode vcf-iid --out ../Results/Regions_of_interest/panicle_count/RoI_chr",given_chr,"_",given_pos)
# system(plink_command)
# #
# # 
# # 
# # #create genotype using VCFtools
#  vcfToolsCommand = paste0("vcftools --vcf ../Results/Regions_of_interest/panicle_count/RoI_chr",given_chr,"_",given_pos,".vcf --012 --out ../Results/Regions_of_interest/panicle_count/filtered_data_chr",given_chr)
#  system(vcfToolsCommand)
# # 
# # 
# genotype = read.table(paste0("../Results/Regions_of_interest/panicle_count/filtered_data_chr",given_chr,".012"),header = F)[,-1]
# id = read.table(paste0("../Results/Regions_of_interest/panicle_count/filtered_data_chr",given_chr,".012.indv"),header = F)[,1]
#  id = gsub("_","",id)
#  row.names(genotype) = id#[,1]
#  pos = read.table(paste0("../Results/Regions_of_interest/panicle_count/filtered_data_chr",given_chr,".012.pos"),header = F)
#  pos=paste0(pos$V1,"_",pos$V2)
# # 
#  colnames(genotype) = pos
# # # saveRDS(genotype,"processed_data/genotype.RDS")
# # # write.csv(genotype,"processed_data/genotype.csv")
# residual_pheno = readRDS("../Results/GWAS_Results/panicle_count/panicle_count_residuals.rds")
# 
# head(residual_pheno)
# dim(residual_pheno)
# 
# # 
# index = match(residual_pheno$genotype,row.names(genotype))
# # #GWAS with residuals for the region of interest
# gwas_results_Region_of_interest = matrix(0,nrow=ncol(genotype),2)
# gwas_results_Region_of_interest[,2]=1
# row.names(gwas_results_Region_of_interest)=colnames(genotype)
# colnames(gwas_results_Region_of_interest) = c("Effect","Pvalue")
# gwas_results_Region_of_interest = data.frame(gwas_results_Region_of_interest)
# # 
# # 
# y = residual_pheno$resuduals
# # system.time({
# #   # for(i in 10000:20000){
#  for(i in 1:ncol(genotype)){
# #     
#      x=genotype[index,i]
#      b=sum((x-mean(x)) * (y-mean(y)))  /  sum((x-mean(x))^2) 
#      a=mean(y)-b*mean(x)
#      yhat=a+b*x
#      # calculate R2, F-value and p-value
#      SST = sum((y - mean(y))^2)      #Corrected Sum of Squares Total;  SST = Σi=1n (yi - y)2  This is the sample variance of the y-variable multiplied by n - 1
#      SSR = sum((yhat - mean(y))^2)  #Corrected Sum of Squares for Model also called sum of squares for regression
#      SSE = sum((y-yhat)^2)   #Sum of Squares for Error also called sum of squares for residuals
# #     
#     R2=SSR/SST # r-squared For simple linear regression, R2 is the square of the sample correlation rxy For multiple linear regression with intercept (which includes simple linear regression), it is defined as r2 = SSM / SST
#     Fval=SSR /(SSE/(length(y)-2)) # F-value
#     Pval=1-pf(Fval,1,(length(y)-2)) # p-value
#     gwas_results_Region_of_interest[i,1] = b
#     gwas_results_Region_of_interest[i,2] = Pval
# 
# 
#   }
#}
#)
# 
# plot(-log10(gwas_results_Region_of_interest$Pvalue))
# #saveRDS(gwas_results_PCA,"Results/GWAS_Results/canopy_height_gwas_180PCA.rds")