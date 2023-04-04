# Fuzzy cMeans clustering
# https://2-bitbio.com/post/fuzzy-cmeans-clustering-of-rnaseq-data-using-mfuzz/

library(data.table)
library(dplyr)
library(tidyr)
library(Mfuzz)
library(ggplot2)
library(ggforce)

Pheno<-fread("Phenotype_Reordered.txt")

load("KNN_Omics.RData")

####################

## Reorder Pheno by CLS and RLS

Pheno_CLS<-Pheno[order(CLS,MeanRLS)]
Pheno_RLS<-Pheno[order(MeanRLS,CLS)]

## Reorder data matrix by CLS 

RespProt_CLS<-knnRespProt[Pheno_CLS$Gene_ID,]
RespMet_CLS<-knnRespMet[Pheno_CLS$Gene_ID,]
RespLip_CLS<-knnRespLip[Pheno_CLS$Gene_ID,]

CLS<-Pheno_CLS$CLS
RespProt_Data<-cbind(CLS,RespProt_CLS)
RespProt_Data<-t(RespProt_Data)
RespMet_Data<-cbind(CLS,RespMet_CLS)
RespMet_Data<-t(RespMet_Data)
RespLip_Data<-cbind(CLS,RespLip_CLS)
RespLip_Data<-t(RespLip_Data)

## Reorder data matrix by RLS

FermProt_RLS<-knnFermProt[Pheno_RLS$Gene_ID,]
FermProt_RLS<-FermProt_RLS[1:151,]
FermMet_RLS<-knnFermMet[Pheno_RLS$Gene_ID,]
FermMet_RLS<-FermMet_RLS[1:151,]
FermLip_RLS<-knnFermLip[Pheno_RLS$Gene_ID,]
FermLip_RLS<-FermLip_RLS[1:151,]

RLS<-Pheno_RLS$MeanRLS
RLS<-RLS[1:151]
FermProt_Data<-cbind(RLS,FermProt_RLS)
FermProt_Data<-t(FermProt_Data)
FermMet_Data<-cbind(RLS,FermMet_RLS)
FermMet_Data<-t(FermMet_Data)
FermLip_Data<-cbind(RLS,FermLip_RLS)
FermLip_Data<-t(FermLip_Data)

##################################### 

MfuzzPrep<-function(Data) {
  tmp <- tempfile()
  write.table(Data,file=tmp, sep='\t', quote = F, col.names=NA)
  Dat<-table2eset(file=tmp)
  Dat.s<-standardise(Dat)
  
  m1 <- mestimate(Dat.s)
  
  set.seed(12344)
  Dmin(Dat.s, m=m1, crange=seq(4,40,1), repeats=3, visu=TRUE)
}

MfuzzPrep(RespProt_Data)
MfuzzPrep(RespMet_Data)
MfuzzPrep(RespLip_Data)

MfuzzPrep(FermProt_Data)
MfuzzPrep(FermMet_Data)
MfuzzPrep(FermLip_Data)

## I am looking for the "elbow" after which there is no significant decrease in centroid distance. 
## Seems like clustering is only useful for protein data, as for both metabolomics and lipidomics data, 
## the more clusters we have the smaller the centroid distance. Meaning that it may not make sense to cluster individual molecules. 

###############################################################

### RespProt

tmp <- tempfile()
write.table(RespProt_Data,file=tmp, sep='\t', quote = F, col.names=NA)
RespProt_Dat<-table2eset(file=tmp)
RespProt_Dat.s<-standardise(RespProt_Dat)

m1 <- mestimate(RespProt_Dat.s)
m1

set.seed(12344)
Dmin(RespProt_Dat.s, m=m1, crange=seq(4,40,1), repeats=3, visu=TRUE)

clust<-18

set.seed(123)
c.RespProt <- mfuzz(RespProt_Dat.s,c=clust,m=m1)

# mfuzz.plot(RespProt_Dat.s,cl=c.RespProt,mfrow=c(1,1),new.window=FALSE)

cor(t(c.RespProt[[1]])) 
## Ideally, no two clusters should exhibit a correlation greater than 0.85

center.RespProt<-c.RespProt[[1]]

assignment.RespProt<-c.RespProt[[3]]
table(assignment.RespProt)

RespProt_Long<-data.table(cbind(CLS,t(center.RespProt)),keep.rownames=TRUE)
RespProt_Long<-melt(RespProt_Long,id.vars=c("rn","CLS"),variable.name="Cluster",value.name="Center")
Plot<-ggplot(RespProt_Long,aes(x=CLS,y=Center))+
  geom_point()+theme_bw()+xlab("CLS")+ylab("Cluster Center")+
  geom_smooth(method="lm",colour="red",se=TRUE)
pdf("Cluster_RespProt_CLS.pdf")
for (i in 1:3) {
  print(Plot+facet_wrap_paginate(facets='Cluster',ncol=2,nrow=3,scales='free_y',page=i))
}
dev.off()


P<-rep(NA,nrow(center.RespProt))
Rho<-rep(NA,nrow(center.RespProt))
for (i in 1:nrow(center.RespProt)) {
  Cor<-cor.test(CLS,center.RespProt[i,],method="spearman")
  P[i]<-Cor$p.value
  Rho[i]<-Cor$estimate
}
RespProt_CLS_Cor<-data.frame(Cluster=rownames(center.RespProt),Rho,P,P.adj=p.adjust(P,method="BH"))

score.RespProt<-acore(RespProt_Dat,c.RespProt,min.acore=0)
acore_list.RespProt <- do.call(rbind, lapply(seq_along(score.RespProt), function(i){ data.frame(CLUSTER=i, score.RespProt[[i]])}))
# write.table(acore_list.RespProt,"tmp.txt",sep="\t",quote=F,row.names=F)

############

### FermProt

tmp <- tempfile()
write.table(FermProt_Data,file=tmp, sep='\t', quote = F, col.names=NA)
FermProt_Dat<-table2eset(file=tmp)
FermProt_Dat.s<-standardise(FermProt_Dat)

m1 <- mestimate(FermProt_Dat.s)
m1

set.seed(12344)
Dmin(FermProt_Dat.s, m=m1, crange=seq(4,40,1), repeats=3, visu=TRUE)

clust<-19

set.seed(123)
c.FermProt <- mfuzz(FermProt_Dat.s,c=clust,m=m1)

# mfuzz.plot(FermProt_Dat.s,cl=c.FermProt,mfrow=c(1,1),new.window=FALSE)

cor(t(c.FermProt[[1]])) 
## Ideally, no two clusters should exhibit a correlation greater than 0.85

center.FermProt<-c.FermProt[[1]]

assignment.FermProt<-c.FermProt[[3]]
table(assignment.FermProt)

FermProt_Long<-data.table(cbind(RLS,t(center.FermProt)),keep.rownames=TRUE)
FermProt_Long<-melt(FermProt_Long,id.vars=c("rn","RLS"),variable.name="Cluster",value.name="Center")
Plot<-ggplot(FermProt_Long,aes(x=RLS,y=Center))+
  geom_point()+theme_bw()+xlab("RLS")+ylab("Cluster Center")+
  geom_smooth(method="lm",colour="red",se=TRUE)
pdf("Cluster_FermProt_RLS.pdf")
for (i in 1:4) {
  print(Plot+facet_wrap_paginate(facets='Cluster',ncol=2,nrow=3,scales='free_y',page=i))
}
dev.off()


P<-rep(NA,nrow(center.FermProt))
Rho<-rep(NA,nrow(center.FermProt))
for (i in 1:nrow(center.FermProt)) {
  Cor<-cor.test(RLS,center.FermProt[i,],method="spearman")
  P[i]<-Cor$p.value
  Rho[i]<-Cor$estimate
}
FermProt_RLS_Cor<-data.frame(Cluster=rownames(center.FermProt),Rho,P,P.adj=p.adjust(P,method="BH"))

score.FermProt<-acore(FermProt_Dat,c.FermProt,min.acore=0)
acore_list.FermProt <- do.call(rbind, lapply(seq_along(score.FermProt), function(i){ data.frame(CLUSTER=i, score.FermProt[[i]])}))
# write.table(acore_list.FermProt,"tmp.txt",sep="\t",quote=F,row.names=F)

save(center.RespProt,assignment.RespProt,center.FermProt,assignment.FermProt,RespProt_CLS_Cor,FermProt_RLS_Cor,acore_list.RespProt,acore_list.FermProt,
     file="Cluster.RData")


############################


# ### K-means clustering 
# 
# clusters<-kmeans(t(knnRespProt),18)
# 
# plot(x=Pheno$CLS,y=clusters[[2]][15,],pch=20)
# 
# Cluster<-clusters[[1]]
# 
# Cluster[Cluster==4]

