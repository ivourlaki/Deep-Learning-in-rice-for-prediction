#### 1. snps #####################################################

library(data.table)
library(stringr)
library(dplyr)
library(readr)


#intruduce corrected phenotypes for each marke sets
data <- read.csv("all.pheno_snps.csv", header=TRUE)
partitions <- read.csv("partitions_validation.csv", header=TRUE)
names<- names(data)
dim(data)

#introduce markers
load("final_six_nooverlap_datasets.RData")
colnames(final_snps)<-  paste0("SNPs","_",colnames(final_snps))
colnames(final_del)<-  paste0("DEL","_",colnames(final_del))
colnames(final_DUP)<-  paste0("DUP","_",colnames(final_DUP))
colnames(final_INV)<-  paste0("INV","_",colnames(final_INV))
colnames(final_mitedtx)<-  paste0("MITEDTX","_",colnames(final_mitedtx))
colnames(final_rlxrix)<-  paste0("RLXRIX","_",colnames(final_rlxrix))

# 1. first marker set
#rm(rixrlx,mitedtx)
# remove accessions
marker=final_snps[,-1]
dim(marker)

l<- c(1,7,8,10) # for saving with the same names the folds instead of 1,2,3,4

for (j in seq(11)) {
  
  for (i in c(1,2,3,4)) {
    pval.geno1<- NULL
    v=partitions[,j+1] # start from FOLd.1
    
    phenotypes=data[,i]
    phenotypes_train <- phenotypes[which(v=="train")]
    
    # find in marker the train and test based on the partitions 
    markers_train <- marker[which(v=="train"),]
    dim(markers_train)
    
    # remove na from phenotypes
    phenotypes_train_final<- na.omit(phenotypes_train)
    length(phenotypes_train_final)
    
    # remove na based on the initila phenotypes 
    markers_train_final <- markers_train[which(!is.na(phenotypes_train)),]
    dim(markers_train_final)
 
    # do GWAS
    for (k in 1:ncol(markers_train_final)) {
    #for (k in 1:10) {
      m <-  summary(lm(phenotypes_train_final ~ markers_train_final[,k] ) ) 
      r2 <- m$adj.r.squared
      if(!is.null(m$fstatistic)) {
        pvalue <- 1-pf(m$fstatistic[1],m$fstatistic[2], m$fstatistic[3])
      } else {
        pvalue <- NA
      }
      pval.geno1 <- rbind(pval.geno1,c(r2,pvalue))
    }
  #with this way I save all the pvalues in one column and the names of each vairants in the other column
  # when I will have all the results from the six markers I will join them by rowbind
  # Resultsa re saved for each trait and partition separately
 
  pval.trait.partition.with.variant.names<- cbind.data.frame(names=colnames(markers_train_final),pvalues=pval.geno1[,2])
  write_csv(pval.trait.partition.with.variant.names,sprintf("~/GWAS/pvalues/pvalue_trait_%.0f_partitions_%.0f_train_snps.csv",l[i],j))
  
  }
}

#save.image("snps_gwas.RData")








