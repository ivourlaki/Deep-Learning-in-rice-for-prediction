
library(AGHmatrix)
library(data.table)
library(genio)
library(BEDMatrix)
library(BGLR)
library(gridExtra)
library(stringr)
library(dplyr)




################################# RKHS ######################################################


nIter=100000



all.loss<- NULL
all.corr<- NULL
all.corr2<- NULL
all.accuracy<- NULL



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
  corr<- NULL
  corr2<- NULL
  loss<- NULL
  accuracy<-NULL
  g=0 
  
  
  for (i in c(1,7,8,10)) {
    
    
    #---- SCALE PHENOTYPES
    y<-(unlist(list.phen[[1]][i]))
    print(names[i])
    
    if (j %in% c(1:10)) {
      marker <- read.csv(paste0("738_best_markers_trait_",i,"_partitions_",j,".csv"),header=FALSE)
      marker<-marker[,-1]
      ## if the following string is detected make the recorde in the TRUE values
      #V1=str_detect(marker, 'DEL|MITE|DTX|RLX|RIX|INV|DUP')
      #V1=str_detect(marker, 'DEL|MITEDTX|RLXRIX|INV|DUP')
      V1=grepl('DEL|MITEDTX|RLXRIX|INV|DUP',marker)
      marker[V1][marker[V1] == 1] <- 2
      marker<-marker[-1,]
      marker<- as.data.frame(sapply(marker,as.numeric))
      G_VanRadenPine <- Gmatrix(SNPmatrix=as.matrix(marker), maf=0.01, method="VanRaden")
      
    }else {
      marker <- read.csv(paste0("738_best_markers_trait_",i,"_partitions_",j,".csv"),header=FALSE)
      marker<-marker[,-1]
      ## if the following string is detected make the recorde in the TRUE values
      V1=grepl('DEL|MITEDTX|RLXRIX|INV|DUP',marker)
      marker[V1][marker[V1] == 1] <- 2
      marker<-marker[-1,]
      marker<- as.data.frame(sapply(marker,as.numeric))
      G_VanRadenPine <- Gmatrix(SNPmatrix=as.matrix(marker), maf=0.01, method="VanRaden")
    }
    
    
    
    
    if (i %in% c(8,10)) {
      y<-scale(y,center=TRUE,scale=TRUE)
      v=partitions[,j+1] # start from FOLd.1
      yNA=y
      tst <- which(v=="test" & !is.na(y))
      yNA[tst]=NA
      fm1 = BGLR(y=yNA,ETA=list(ETA1=list(K=G_VanRadenPine,model='RKHS')), nIter=nIter,verbose=F)
      
      g=g+1
      corel1=cor(fm1$yHat[tst], y[tst])
      corr<- cbind(corr,corel1)
      colnames(corr)[g]<- paste("average_corr",i)
      
      
      loss1<-mean((fm1$yHat[tst]-y[tst])^2)
      loss<- cbind( loss,as.numeric(loss1))
      colnames(loss)[g]<- paste("average_loss",i)
      
      accuracy1=0
      accuracy<- cbind( accuracy,as.numeric(accuracy1))
      colnames(accuracy)[g]<- paste("average_accuracy",i)
      
    } else {
      y<-ifelse(y==2,1,0)
      v=partitions[,j+1] # start from FOLd.1
      yNA=y
      tst <- which(v=="test" & !is.na(y))
      yNA[tst]=NA
      
      fm1 = BGLR(y=yNA,response_type="ordinal",ETA=list(ETA1=list(K=G_VanRadenPine,model='RKHS')), nIter=nIter,verbose=F)
      
      g=g+1
      
      
    
      
      fm1$yHat[tst][which(fm1$probs[tst,2] > 0.5)] <- 1
      fm1$yHat[tst][which(fm1$probs[tst,2]  < 0.5)] <- 0
      
      y2<-scale(y[tst],center=TRUE,scale=TRUE)
      y3<-scale(fm1$yHat[tst],center=TRUE,scale=TRUE)
      corel1=cor(y3, y2)
      corr<- cbind(corr,corel1)
      colnames(corr)[g]<- paste("average_scale_corr",i)
      
      
      loss1= (-1/length(tst))*sum(y[tst]*log(fm1$probs[tst,2])+(1-y[tst])*log(1-fm1$probs[tst,2]))
      loss<- cbind( loss,as.numeric(loss1))
      colnames(loss)[g]<- paste("average_loss",i)
      
      corel2=cor(fm1$yHat[tst], y[tst])
      corr2<- cbind(corr2,corel2)
      colnames(corr2)[g]<- paste("average_corr",i)
      
      if (length(unique(fm1$yHat[tst])) =="2") {
        confusionm<- table(fm1$yHat[tst], y[tst])
        print(confusionm)
        accuracy1<- sum(confusionm[1],confusionm[4])/sum(confusionm[1:4])
        print(accuracy1)
        accuracy<- cbind( accuracy,as.numeric(accuracy1))
        colnames(accuracy)[g]<- paste("average_accuracy",i)
      } else if (fm1$yHat[tst]=="0") {
        
        confusionm<- table(fm1$yHat[tst], y[tst])
        accuracy1<- sum(confusionm[1])/sum(confusionm[1:2])
        print(confusionm)
        print(accuracy1)
        accuracy<- cbind( accuracy,as.numeric(accuracy1))
        colnames(accuracy)[g]<- paste("average_accuracy",i)
      } else { 
        confusionm<- table(fm1$yHat[tst], y[tst])
        print(confusionm)
        accuracy1<- sum(confusionm[2])/sum(confusionm[1:2])
        print(accuracy1)
        accuracy<- cbind( accuracy,as.numeric(accuracy1))
        colnames(accuracy)[g]<- paste("average_accuracy",i)
        
      }
    }
    
    
    
    
    fm [i] =list(fm1)
  }
  
  all.corr<- rbind.data.frame( all.corr,corr)
  all.corr2<- rbind.data.frame( all.corr2,corr2)
  all.loss<- rbind.data.frame( all.loss,loss)
  all.accuracy<- rbind.data.frame( all.accuracy,accuracy)
  

  
}

all.metrics <- cbind.data.frame(all.loss,all.corr,all.corr2,all.accuracy, architecture=rep("RKHS",length(all.loss[,1])),
                                marker=rep("ALL",length(all.loss[,1])), Partition=c("Partition_1","Partition_2","Partition_3","Partition_4",
                                                                                    "Partition_5","Partition_6","Partition_7","Partition_8",
                                                                                    "Partition_9","Partition_10","ARO_ADM"))


#row.names(all.loss.cnn.all)   <-  colnames(all.parameters.as.columns)                                 


write.csv(all.metrics,paste0("RKHS_ALL_10000_cor_corrected.csv"), row.names = FALSE)

save.image("RKHS_ALL_10000_cor_corrected.RData")
rm(list = ls())

