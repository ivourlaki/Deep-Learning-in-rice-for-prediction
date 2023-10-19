


#### first select best pvalues and markers ######

for (j in seq(11)) {
  percentage.variants.trait<-0
  for (i in c(1,7,8,10)) {
    
    #### 1. select best 10000 pvalues 
    pvalue.snps<- read.csv(sprintf("~/GWAS/pvalues/pvalue_trait_%.0f_partitions_%.0f_train_snps.csv",i,j), header=TRUE)
   
    all.pvalues<- rbind.data.frame(pvalue.snps)
    dim(all.pvalues)
    all.pvalues.order<- all.pvalues[order(all.pvalues$pvalues, decreasing = FALSE),]
    best.pvalues<- all.pvalues.order[1:10000,]
    write_csv(best.pvalues,sprintf("~/GWAS/best_pvalues_partitions_snps/best_pvalues_snps_trait_%.0f_partitions_%.0f_train.csv",i,j))
    
    # how many are snps and how many svs
    detect.SNPs<- grepl('SNPs',best.pvalues[,1])
    how.many.SNPs<- sum(detect.SNPs, na.rm=TRUE)
    
    
    #count the percentages and save them 
    percentage.variants<- rbind.data.frame((how.many.SNPs/10000)*100)
    
    #rownames(percentage.variants)<- c("SNPs","DEL","DUP","INV","MITEDTX","RLXRIX")
    #colnames(percentage.variants)<-c("Percentage")
    percentage.variants.trait<- cbind.data.frame(percentage.variants.trait, percentage.variants)
    
    
    #### 2. select best markers based on pvalues
    load("/Users/ivourlaki/final_SVS_matrix/final_six_nooverlap_datasets.RData")
    colnames(final_snps)<-  paste0("SNPs","_",colnames(final_snps))
    
    
    # remove  accessions 
    Accessions<- final_snps[,1]
    final_snps<-final_snps[,-1]
    
    
    # select the best markers in each set based on the names column
    best.final.snps <- final_snps[,which(colnames(final_snps) %in% best.pvalues[,1])]
    if (dim(best.final.snps)[2] == how.many.SNPs ) { print("TRUE")}
    
    
    # merge all these markers , it has to be 10000
    all.best.markers<- cbind.data.frame(best.final.snps)
    all.best.markers<- cbind.data.frame(Accessions,all.best.markers)
    dim(all.best.markers)
    # save the matrix for all the accessions even it was created based on the pvalues from the training set
    #we will use that in linear models
    write_csv(all.best.markers,sprintf("~/GWAS/best_pvalues_partitions_snps/738_best_markers_snps_trait_%.0f_partitions_%.0f.csv",i,j))
    
  }
  
  # save for each partition. Each csv file includes the values of each trait for one partition
  
  colnames(percentage.variants.trait)<- c("0","trait1","trait7","trait8","trait10")
  variants.types<- c("SNPS")
  percentage.variants.trait<- cbind.data.frame(variants.types,percentage.variants.trait[,-1])
  write_csv(percentage.variants.trait,sprintf("~/GWAS/best_pvalues_partitions_snps/percentage.snps.all.traits_partitions_%.0f.csv",j))
  
}


################################# divide in partiions based on pvaleus from training set ######################################
################################# all ############################################

library(data.table)

all.phenotypes <- read.csv("~/GWAS/all.phenotypes.csv", header=TRUE)
partitions <- read.csv("~/GWAS/partitions_mine.csv", header=TRUE)

data<- cbind(all.phenotypes[,10:20])
names<- names(data)

names[1]<- "culm.diameter.1st.internode"
names[8]<- "grain.weight"
names[9]<- "salt.injury.at.EC12"
names[10]<- "time.to.flowering.from.sowing"
dim(data)

for (j in seq(11)) {
  
  for (i in c(1,7,8,10)) {
    print(names[i])
    marker <- read.csv(sprintf("~/GWAS/best_pvalues_partitions_snps/738_best_markers_snps_trait_%.0f_partitions_%.0f.csv",i,j))
    dim(marker)
    marker<-marker[,-1]
    
    # dont scale here I will do it at the Deep learning and at the Bayessian separately
    #marker<-scale(marker,center=TRUE,scale=TRUE)
    
    v=partitions[,j+1] # start from FOLd.1
    
    phenotypes=data[,i]
    
    
    # find in phenotypes the train and test
    phenotypes=data[,i]
    phenotypes_train <- phenotypes[which(v=="train")]
    phenotypes_test <- phenotypes[which(v=="test")]
    
    # find in marker the train and test based on the partitions 
    markers_train <- marker[which(v=="train"),]
    markers_test <- marker[which(v=="test"),]
    dim(markers_train)
    dim(markers_test)
    
    # remove na from phenotypes
    phenotypes_train_final<- na.omit(phenotypes_train)
    phenotypes_test_final<- na.omit(phenotypes_test)
    
    length(phenotypes_train_final)
    length(phenotypes_test_final)
    
    # remove na based on the initila phenotypes 
    markers_train_final <- markers_train[which(!is.na(phenotypes_train)),]
    markers_test_final <-  markers_test[which(!is.na(phenotypes_test)),]
    
    dim(markers_train_final)
    dim(markers_test_final)
    
    write.table(markers_train_final,sprintf("~/GWAS/best_pvalues_partitions_snps/best_markers_snps_trait_%.0f_partitions_%.0f_train.csv",i,j),sep=",",row.names = FALSE,col.names = FALSE)
    
    write.table(markers_test_final,sprintf("~/GWAS/best_pvalues_partitions_snps/best_markers_snps_trait_%.0f_partitions_%.0f_test.csv",i,j),sep=",",row.names = FALSE,col.names = FALSE)
    
  }
  
}

















