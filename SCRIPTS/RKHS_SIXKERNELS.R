
library(AGHmatrix)
library(data.table)
library(genio)
library(BEDMatrix)
library(BGLR)
library(gridExtra) 




################################# RKHS ######################################################


nIter=100000


all.loss<- NULL

#--- input phenotypes and partitions

all.phenotypes <- read.csv("all.phenotypes.csv", header=TRUE)
partitions <- read.csv("partitions_validation.csv", header=TRUE)



#---- LIST WITH ALL PHENOTYPES FOR LOOP
data<- cbind(all.phenotypes[,10:20])
names<- names(data)

names[1]<- "culm.diameter.1st.internode"
names[8]<- "grain.weight"
names[9]<- "salt.injury.at.EC12"
names[10]<- "time.to.flowering.from.sowing"

list.phen<-list(all.phenotypes[,10:20])


fm= list()

for (j in seq(11)) {

  loss<- NULL
 
  g=0 
  
  
  for (i in c(1,7,8,10)) {
    
    
    #---- SCALE PHENOTYPES
    y<-(unlist(list.phen[[1]][i]))
    print(names[i])
    
   

    
    load("final_six_nooverlap_kernels.RData")
    MITEDTX<- kernel_mitedtx
    RLXRIX<- kernel_rlxrix
    SNP<- kernel_snps
    DEL<- kernel_del
    INV<- kernel_inv
    DUP<- kernel_dup
  
    
    
    if (i %in% c(8,10)) {
      y<-scale(y,center=TRUE,scale=TRUE)
      v=partitions[,j+1] # start from FOLd.1
      yNA=y
      tst <- which(v=="test" & !is.na(y))
      yNA[tst]=NA
      fm1 = BGLR(y=yNA,ETA=list(ETA1=list(K=as.matrix(DEL),model='RKHS'),ETA2=list(K=as.matrix(DUP),model='RKHS'),ETA3=list(K=as.matrix(INV),model='RKHS'),
                                ETA4=list(K=as.matrix(MITEDTX),model='RKHS'), ETA5=list(K=as.matrix(RLXRIX),model='RKHS'),ETA6=list(K=as.matrix(SNP),model='RKHS')), nIter=nIter,verbose=F)
      
      g=g+1
      
      loss1<-mean((fm1$yHat[tst]-y[tst])^2)
      loss<- cbind( loss,as.numeric(loss1))
      colnames(loss)[g]<- paste("average_loss",i)
      
  
      
    } else {
      y<-ifelse(y==2,1,0)
      v=partitions[,j+1] # start from FOLd.1
      yNA=y
      tst <- which(v=="test" & !is.na(y))
      yNA[tst]=NA
      
    

      fm1 = BGLR(y=yNA,response_type="ordinal",ETA=list(ETA1=list(K=as.matrix(DEL),model='RKHS'),ETA2=list(K=as.matrix(DUP),model='RKHS'),ETA3=list(K=as.matrix(INV),model='RKHS'),
                                ETA4=list(K=as.matrix(MITEDTX),model='RKHS'), ETA5=list(K=as.matrix(RLXRIX),model='RKHS'),ETA6=list(K=as.matrix(SNP),model='RKHS')), nIter=nIter,verbose=F)
      
      g=g+1
      
     
      
      loss1= (-1/length(tst))*sum(y[tst]*log(fm1$probs[tst,2])+(1-y[tst])*log(1-fm1$probs[tst,2]))
      loss<- cbind( loss,as.numeric(loss1))
      colnames(loss)[g]<- paste("average_loss",i)
      
    
    }
    
    fm [i] =list(fm1)
  }
  
  all.loss<- rbind.data.frame( all.loss,loss)
  
}

all.metrics <- cbind.data.frame(all.loss, architecture=rep("RKHS",length(all.loss[,1])),
                                marker=rep("SIX_KERNELS",length(all.loss[,1])), Partition=c("Partition_1","Partition_2","Partition_3","Partition_4",
                                                                                    "Partition_5","Partition_6","Partition_7","Partition_8",
                                                                                    "Partition_9","Partition_10","ARO_ADM"))


write.csv(all.metrics,paste0("RKHS_SIXKERNELSp.csv"), row.names = FALSE)

save.image("RKHS_SIXKERNELSp.RData")
rm(list = ls())

