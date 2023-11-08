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
