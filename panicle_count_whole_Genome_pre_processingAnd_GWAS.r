genotype = readRDS("ProcessedData/transposed_geno.rds")
genotype=t(genotype)

panicle_count = read.csv("ProcessedData/formatted_Phenotype_files/panicle_count_preprocess.csv")
dim(panicle_count)
plot(density(panicle_count$panicle_count))
#Remove the individuals with panicle_count <=0
panicle_count = panicle_count[which(panicle_count$panicle_count>0),]
dim(panicle_count)
plot(density(panicle_count$panicle_count))

mean_ = mean(panicle_count$panicle_count)
sd_ = sd(panicle_count$panicle_count)

index = which(panicle_count$panicle_count< mean_-3*sd_ | panicle_count$panicle_count> mean_+3*sd_ )
panicle_count[index,]

#remove outlier

panicle_count = panicle_count[-index,]

index = match(panicle_count$genotype,rownames(genotype))

#remove those do not have genotype
panicle_count = panicle_count[which(is.na(index)==F),]
dim(panicle_count) # 5193    7
my_date = as.Date(panicle_count$date)
my_julian <- as.POSIXlt(my_date)$yday    # Convert date to Julian day
head(my_julian)

panicle_count$JulianDay = my_julian

PCA = read.csv("ProcessedData/SVD_scores_PLINK_GRM.csv")
all.equal(PCA$X,rownames(genotype))
index = match(panicle_count$genotype,PCA$X)

panicle_count = data.frame(panicle_count,PCA[index,])

all.equal(panicle_count$genotype,panicle_count$X.1)

panicle_count=panicle_count[,-(which(colnames(panicle_count)=="X.1"))]
#panicle_count =panicle_count[,c(6,2:5,8:ncol(panicle_count))]
PCA = as.matrix(panicle_count[,9:ncol(panicle_count)])

model1 = lm(panicle_count$panicle_count~ panicle_count$range + panicle_count$column + panicle_count$treatment + panicle_count$JulianDay +PCA[,1:60] ) #Adjusted R-squared:  0.889
summary(model1)
plot(density(model1$residuals))


## Get the residuals out from the model and use it as the phenotype in the GWAS
panicle_count$residuals = model1$residuals

head(panicle_count)
# x = panicle_count[which(panicle_count$genotype == "PI570431"),]
# dim(x)
# plot(density(x$residuals))
# mean(x$residuals)
# median(x$residuals)
# plot(density(x$panicle_count))


#GWAS with residuals
gwas_results_PCA = matrix(0,nrow=ncol(genotype),2)
gwas_results_PCA[,2]=1
row.names(gwas_results_PCA)=colnames(genotype)
colnames(gwas_results_PCA) = c("Effect","Pvalue")
gwas_results_PCA = data.frame(gwas_results_PCA)
gwas_results_PCA = data.frame(gwas_results_PCA,Percentage_0=0,Percentage_1=0, Percentage_2=0)

index = match(panicle_count$genotype,row.names(genotype))

y = panicle_count$residuals
system.time({
   for(i in 1:ncol(genotype)){

    x=genotype[index,i]
    b=sum((x-mean(x)) * (y-mean(y)))  /  sum((x-mean(x))^2)
    a=mean(y)-b*mean(x)
    yhat=a+b*x
    # calculate R2, F-value and p-value
    SST = sum((y - mean(y))^2)      #Corrected Sum of Squares Total;  SST = Σi=1n (yi - y)2  This is the sample variance of the y-variable multiplied by n - 1
    SSR = sum((yhat - mean(y))^2)  #Corrected Sum of Squares for Model also called sum of squares for regression
    SSE = sum((y-yhat)^2)   #Sum of Squares for Error also called sum of squares for residuals

    R2=SSR/SST # r-squared For simple linear regression, R2 is the square of the sample correlation rxy For multiple linear regression with intercept (which includes simple linear regression), it is defined as r2 = SSM / SST
    Fval=SSR /(SSE/(length(y)-2)) # F-value
    Pval=1-pf(Fval,1,(length(y)-2)) # p-value
    gwas_results_PCA[i,1] = b
    gwas_results_PCA[i,2] = Pval
    if(i%%50000==0){
      print(i)
      saveRDS(gwas_results_PCA,"Results/GWAS_Results/panicle_count/panicle_count_gwas_60PCA.rds")
    }
  }
})

saveRDS(gwas_results_PCA,"Results/GWAS_Results/panicle_count/panicle_count_gwas_60PCA.rds")
