library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)

#########

Pheno<-fread("Phenotype.txt")
## Rename the columns so that they are easier to refer to 
colnames(Pheno)<-c("SGD_ID","Gene_ID","FermDensityFC","RespDensityFC","FermP","RespP","YPG_Phenotype","CLS","CLS%Change","MeanRLS","MaxRLS")
Pheno<-Pheno[-169,c(1,2,3,4,8,10,11)] ## Delete the last line containing wt info. Keep important columns. 
## Reorder the rows according to ordering of the omics dataset
Strain<-read.table("Strain.txt")
Strain<-as.character(Strain$V1)
Pheno<-Pheno %>% slice(match(Strain,Gene_ID))
Pheno<-data.table(Pheno)

write.table(Pheno,"Phenotype_Reordered.txt",sep="\t",quote=F,row.names=F)

##########

load("KNN_Omics.RData")



# PCA<-function(DT) {
#   index<-colSums(is.na(DT))==0
#   pr.out<-prcomp(DT[,..index][,-1],scale=TRUE)
#   percentage<-round(pr.out$sdev^2/sum(pr.out$sdev^2)*100, 2)
#   percentage<-paste(colnames(pr.out$x), "(", as.character(percentage), "%)",sep="" )
#   scores<-as.data.frame(pr.out$x)
#   scores<-cbind(scores, RLS=Pheno$MeanRLS)
#   p<-ggplot(data=scores,aes(x=PC1,y=PC2,col=RLS))+geom_point()+
#     scale_color_gradient2(midpoint=20,low="blue",mid="white",high="red",space="Lab")+
#     xlab(percentage[1])+ylab(percentage[2])+theme_bw()
#   q<-ggplot(data=scores,aes(x=PC1,y=PC3,col=RLS))+geom_point()+
#     scale_color_gradient2(midpoint=20,low="blue",mid="white",high="red",space="Lab")+
#     xlab(percentage[1])+ylab(percentage[3])+theme_bw()
#   r<-p<-ggplot(data=scores,aes(x=PC2,y=PC3,col=RLS))+geom_point()+
#     scale_color_gradient2(midpoint=20,low="blue",mid="white",high="red",space="Lab")+
#     xlab(percentage[2])+ylab(percentage[3])+theme_bw()
#   return(list(p,q,r,pr.out))
# }

### Define function that generates PCA plots with dots colored by phenotype

PCA_Pheno<-function(DataMatrix,Phenotype) {
  PCA<-prcomp(DataMatrix,scale.=TRUE)
  DF<-data.frame(Strain=rownames(PCA$x),PCA$x[,1:3],Phenotype=Phenotype)
  percentage<-round(PCA$sdev^2/sum(PCA$sdev^2)*100, 2)
  percentage<-paste(colnames(PCA$x), "(", as.character(percentage), "%)",sep="" )
  
  r<-ggplot(data=DF, aes(x=PC1,y=PC2,col=Phenotype))+geom_point()+
    scale_color_gradient2(midpoint=median(DF$Phenotype,na.rm=TRUE),low="blue",mid="white",high="red",space="Lab")+
    xlab(percentage[1])+ylab(percentage[2])+theme_bw()
  s<-ggplot(data=DF, aes(x=PC1,y=PC3,col=Phenotype))+geom_point()+
    scale_color_gradient2(midpoint=median(DF$Phenotype,na.rm=TRUE),low="blue",mid="white",high="red",space="Lab")+
    xlab(percentage[1])+ylab(percentage[3])+theme_bw()
  t<-ggplot(data=DF, aes(x=PC2,y=PC3,col=Phenotype))+geom_point()+
    scale_color_gradient2(midpoint=median(DF$Phenotype,na.rm=TRUE),low="blue",mid="white",high="red",space="Lab")+
    xlab(percentage[2])+ylab(percentage[3])+theme_bw()
  grid.arrange(r,s,t,ncol=2)
}

PCA_Pheno(knnRespProt,Pheno$CLS)
PCA_Pheno(knnRespMet,Pheno$CLS)
PCA_Pheno(knnRespLip,Pheno$CLS)

PCA_Pheno(knnFermProt,Pheno$MeanRLS)
PCA_Pheno(knnFermMet,Pheno$MeanRLS)
PCA_Pheno(knnFermLip,Pheno$MeanRLS)


PCA_Pheno(knnRespProt,Pheno$RespDensityFC)
PCA_Pheno(knnRespMet,Pheno$RespDensityFC)
PCA_Pheno(knnRespLip,Pheno$RespDensityFC)

PCA_Pheno(knnFermProt,Pheno$FermDensityFC)
PCA_Pheno(knnFermMet,Pheno$FermDensityFC)
PCA_Pheno(knnFermLip,Pheno$FermDensityFC)


##########################################

## Supervised learning 

###CLS vs RespLip

CLS_RespLip<-data.frame(Pheno[,"CLS"],RespLip[,-1])

library(randomForest)

set.seed(1)
train<-sample(1:nrow(CLS_RespLip),120)
test<-CLS_RespLip[-train,"CLS"]

bag<-randomForest(CLS~.,data=CLS_RespLip,subset=train,mtry=52,importance=TRUE)
bag

yhat.bag<-predict(bag,newdata=CLS_RespLip[-train,])
plot(yhat.bag,test)
abline(0,1)
mean((yhat.bag-test)^2)

###

set.seed(10)
rf<-randomForest(CLS~.,data=CLS_RespLip,subset=train,mtry=17,importance=TRUE)
yhat.rf<-predict(rf,newdata=CLS_RespLip[-train,])
mean((yhat.rf-test)^2)

importance(rf)
varImpPlot(rf)

###

library(gbm)

set.seed(1)
boost<-gbm(CLS~.,data=CLS_RespLip[train,],distribution="gaussian",n.trees=5000,interaction.depth=4)
summary(boost)

yhat.boost<-predict(boost,newdata=CLS_RespLip[-train,],n.trees=5000)
mean((yhat.boost-test)^2)

###





