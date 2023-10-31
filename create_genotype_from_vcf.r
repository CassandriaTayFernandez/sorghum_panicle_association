#########################################################################################

# create PLINK file from vcf

plinkCommand = "~/Tools/plink --vcf BAP_Merged_filteredVCFtools_maf005.vcf --const-fid --make-bed --out processed_data/imp_TERRA_RIL_SNP"
system(plinkCommand)

#Basic missingness check
plinkCommand = "~/Tools/plink --bfile processed_data/imp_TERRA_RIL_SNP  --missing  --out processed_data/imp_TERRA_RIL_SNP"
system(plinkCommand)

plinkCommand = "~/Tools/plink --bfile processed_data/imp_TERRA_RIL_SNP  --maf 0.05 --make-bed  --out processed_data/cleaned"
system(plinkCommand)

#LD prune
plinkCommand = "~/Tools/plink --bfile processed_data/cleaned --indep-pairwise 50 5 0.7 --out processed_data/cleaned_pruned"
system(plinkCommand)

#Extract the filtered SNP PLINK format
plinkCommand = "~/Tools/plink --bfile processed_data/cleaned --extract processed_data/cleaned_pruned.prune.in --make-bed --out processed_data/cleaned_pruned"
system(plinkCommand)


#Extract the filtered SNP in VCF format
plinkCommand = "~/Tools/plink --bfile processed_data/cleaned_pruned --recode vcf-iid --out processed_data/cleaned_pruned"
system(plinkCommand)


#create genotype using VCFtools
vcfToolsCommand = "vcftools --vcf processed_data/cleaned_pruned.vcf --012 --out processed_data/cleaned_pruned"
system(vcfToolsCommand)


genotype = read.table("processed_data/cleaned_pruned.012",header = F)[,-1]
id = read.table("processed_data/cleaned_pruned.012.indv",header = F)
row.names(genotype) = id[,1]
pos = read.table("processed_data/cleaned_pruned.012.pos",header = F)

snp_map = read.table("processed_data/cleaned_pruned.bim",header = F)

#Just check is pos and snp_map have same SNP information
sum(snp_map$V1-pos$V1)
sum(snp_map$V4 - pos$V2)

colnames(genotype) = snp_map$V2
saveRDS(genotype,"processed_data/genotype.RDS")
write.csv(genotype,"processed_data/genotype.csv")

write.csv(geno,"processed_data/genotype_transposed.csv")
