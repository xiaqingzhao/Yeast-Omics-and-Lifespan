library(data.table)
library(dplyr)
library(caret)
library(ggplot2)
library(gridExtra)
library(DMwR)

Resp<-fread("RespOmics.txt")
Resp<-select(Resp,-c(ATP3,COQ8,SCH9,TIM8,'YAL044W-A','YBR230W-A')) # Delete lines that do not have lifespan data
colnames(Resp)[1]<-'Molecule'
colnames(Resp)[2]<-'Type'

Ferm<-fread("FermOmics.txt")
Ferm<-select(Ferm,-c(ATP3,COQ8,SCH9,TIM8,'YAL044W-A','YBR230W-A')) # Delete lines that do not have lifespan data
colnames(Ferm)[1]<-'Molecule'
colnames(Ferm)[2]<-'Type'

##############################################

## Keep complete data only

CompResp<-Resp[complete.cases(Resp)]
CompFerm<-Ferm[complete.cases(Ferm)]

CompRespMolecule<-CompResp[,Molecule]
CompRespType<-CompResp[,Type]
CompFermMolecule<-CompFerm[,Molecule]
CompFermType<-CompFerm[,Type]

CompResp<-CompResp[,-c(1,2)]
CompFerm<-CompFerm[,-c(1,2)]

CompResp<-t(CompResp)
CompFerm<-t(CompFerm)

colnames(CompResp)<-CompRespMolecule
colnames(CompFerm)<-CompFermMolecule

CompRespProt<-CompResp[,CompRespType=="Protein"]
CompRespMet<-CompResp[,CompRespType=="Metabolite"]
CompRespLip<-CompResp[,CompRespType=="Lipid"]

CompFermProt<-CompFerm[,CompFermType=="Protein"]
CompFermMet<-CompFerm[,CompFermType=="Metabolite"]
CompFermLip<-CompFerm[,CompFermType=="Lipid"]

save(CompRespProt,CompRespMet,CompRespLip,CompFermProt,CompFermMet,CompFermLip,file="Complete_Omics.RData")

PCAPlot<-function(DataMatrix) {
  PCA<-prcomp(DataMatrix,scale.=TRUE)
  DF<-data.frame(Strain=rownames(PCA$x),PCA$x[,1:3])
  percentage<-round(PCA$sdev^2/sum(PCA$sdev^2)*100, 2)
  percentage<-paste(colnames(PCA$x), "(", as.character(percentage), "%)",sep="" )
  ggplot(data=DF, aes(x=PC1, y=PC2))+geom_text(aes(label=Strain),size=2)+
    xlab(percentage[1])+ylab(percentage[2])+theme_bw()
}

PCAPlot(CompRespProt)
PCAPlot(CompRespMet)
PCAPlot(CompRespLip)

PCAPlot(CompFermProt)
PCAPlot(CompFermMet)
PCAPlot(CompFermLip)


# fwrite(CompResp,"CompleteResp.csv")
# write.table(CompRespType,"CompleteRespType.txt",quote=F,row.names=F,col.names=F)
# fwrite(CompFerm,"CompleteFerm.csv")
# write.table(CompFermType,"CompleteFermType.txt",quote=F,row.names=F,col.names=F)

##############################################

## Delete molecules with too many missing values

# RMiss<-rowSums(is.na(Resp))
# cumsum(table(RMiss))
# FMiss<-rowSums(is.na(Ferm))
# cumsum(table(FMiss))
# 
# table(colSums(is.na(Resp)))
# table(colSums(is.na(Ferm)))

## There is a "jump" at rowSums(is.na(Resp/Ferm))==17
# tmp<-Ferm[rowSums(is.na(Ferm))==17,]
# fwrite(tmp,"tmp.csv")
## Would make sense to keep molecules with no more than 16 (10%) missing values

Resp<-Resp[rowSums(is.na(Resp))<17,]
Ferm<-Ferm[rowSums(is.na(Ferm))<17,]

##################################################

## Transpose
RespMolecule<-Resp[,Molecule]
RespType<-Resp[,Type]
FermMolecule<-Ferm[,Molecule]
FermType<-Ferm[,Type]

Strain<-colnames(Resp)[-c(1:2)]

Resp<-Resp[,-c(1,2)]
Ferm<-Ferm[,-c(1,2)]

Resp<-t(Resp)
Ferm<-t(Ferm)

colnames(Resp)<-RespMolecule
colnames(Ferm)<-FermMolecule

# fwrite(Resp,"DeMissingResp.csv")
# write.table(RespType,"DeMissingRespType.txt",quote=F,row.names=F,col.names=F)
# fwrite(Ferm,"DeMissingFerm.csv")
# write.table(FermType,"DeMissingFermType.txt",quote=F,row.names=F,col.names=F)

RespProt<-Resp[,RespType=="Protein"]
RespMet<-Resp[,RespType=="Metabolite"]
RespLip<-Resp[,RespType=="Lipid"]

FermProt<-Ferm[,FermType=="Protein"]
FermMet<-Ferm[,FermType=="Metabolite"]
FermLip<-Ferm[,FermType=="Lipid"]

save(RespProt,RespMet,RespLip,FermProt,FermMet,FermLip,file="DeMissing_Omics.RData")


###################################################

# ## Zero- and near zero-variance predictors
# 
# nzv<-nearZeroVar(Resp,saveMetrics=TRUE)
# table(nzv$zeroVar)
# table(nzv$nzv)
# nzv<-nearZeroVar(Ferm,saveMetrics=TRUE)
# table(nzv$zeroVar)
# table(nzv$nzv)
# ### There are no zero- or near zero-variance predictors. Great! 

#################################################

## KNN Imputation 

knnRespProt<-knnImputation(RespProt,scale=F)
knnRespMet<-knnImputation(RespMet,scale=F)
knnRespLip<-knnImputation(RespLip,scale=F)

knnFermProt<-knnImputation(FermProt,scale=F)
knnFermMet<-knnImputation(FermMet,scale=F)
knnFermLip<-knnImputation(FermLip,scale=F)

PCAPlot(knnRespProt)

save(knnRespProt,knnRespMet,knnRespLip,knnFermProt,knnFermMet,knnFermLip,file="KNN_Omics.RData")

#################################################

## Identify and delete correlated predictors 

FilterCor<-function(DataMatrix) {
  Cor<-cor(DataMatrix)
  HighCor<-findCorrelation(Cor,cutoff=0.7,exact=TRUE)
  return(DataMatrix[,-HighCor])
}

DeCorRespProt<-FilterCor(knnRespProt)
DeCorRespMet<-FilterCor(knnRespMet)
DeCorRespLip<-FilterCor(knnRespLip)

DeCorFermProt<-FilterCor(knnFermProt)
DeCorFermMet<-FilterCor(knnFermMet)
DeCorFermLip<-FilterCor(knnFermLip)

save(DeCorRespProt,DeCorRespMet,DeCorRespLip,DeCorFermProt,DeCorFermMet,DeCorFermLip,file="DeCor_Omics.RData")

# Cor<-cor(Resp,use="pairwise.complete.obs")
# #sum(abs(Cor[upper.tri(Cor)])>0.99)
# summary(Cor[upper.tri(Cor)])
# HighCor<-findCorrelation(Cor,cutoff=0.85, exact=TRUE, verbose=TRUE)
# filteredResp<-Resp[,-HighCor]
# filteredRespType<-RespType[-HighCor]
# Cor2<-cor(filteredResp,use="pairwise.complete.obs")
# summary(Cor2[upper.tri(Cor2)])




########

# ## Filtered separately 
# 
# RespP<-Resp[,RespType=="Protein"]
# RespM<-Resp[,RespType=="Metabolite"]
# RespL<-Resp[,RespType=="Lipid"]
# 
# Cor<-cor(RespL,use="pairwise.complete.obs")
# summary(Cor[upper.tri(Cor)])
# HighCor<-findCorrelation(Cor,cutoff=0.85, exact=TRUE, verbose=TRUE)
# filteredRespL<-RespL[,-HighCor]
# Cor2<-cor(filteredRespL,use="pairwise.complete.obs")
# summary(Cor2[upper.tri(Cor2)])
# 
# FermP<-Ferm[,FermType=="Protein"]
# FermM<-Ferm[,FermType=="Metabolite"]
# FermL<-Ferm[,FermType=="Lipid"]
# 
# Cor<-cor(FermL,use="pairwise.complete.obs")
# summary(Cor[upper.tri(Cor)])
# HighCor<-findCorrelation(Cor,cutoff=0.85, exact=TRUE, verbose=TRUE)
# filteredFermL<-FermL[,-HighCor]
# Cor2<-cor(filteredFermL,use="pairwise.complete.obs")
# summary(Cor2[upper.tri(Cor2)])

##########################################################

# ## Imputation 
# ### Note that KNN imputation with preProcess() will automatically trigger center and scale
# 
# DeCorKNN<-preProcess(filteredResp,method="knnImpute")
# DeCorKNN
# DeCorKnnResp<-predict(DeCorKNN,newdata=filteredResp)
# 
# DeCorKNN<-preProcess(filteredFerm,method="knnImpute")
# DeCorKNN
# DeCorKnnFerm<-predict(DeCorKNN,newdata=filteredFerm)
# 
# fwrite(DeCorKnnResp,"DeCorKNNResp.csv")
# fwrite(DeCorKnnFerm,"DeCorKNNFerm.csv")
# 
# ### See if KNN change the data dramatically by PCA
# 
# library(ggplot2)
# library(ggfortify)
# 
# autoplot(prcomp(CompResp,scale=TRUE),shape=FALSE,label.size=2)+theme_bw()
# autoplot(prcomp(KnnResp),shape=FALSE,label.size=2)+theme_bw()
# autoplot(prcomp(DeCorKnnResp),shape=FALSE,label.size=2)+theme_bw()
# 
# autoplot(prcomp(CompFerm,scale=TRUE),shape=FALSE,label.size=2)+theme_bw()
# autoplot(prcomp(KnnFerm),shape=FALSE,label.size=2)+theme_bw()
# autoplot(prcomp(DeCorKnnFerm),shape=FALSE,label.size=2)+theme_bw()
# 
# CompRespPCA<-prcomp(CompResp,scale=TRUE)
# KnnRespPCA<-prcomp(KnnResp,scale=FALSE) ## scale is set to FALSE because the data is already centered and scaled when performing KNN imputation
# DeCorKnnRespPCA<-prcomp(DeCorKnnResp,scale=FALSE)
# cor.test(CompRespPCA$x[,1],KnnRespPCA$x[,1])
# cor.test(CompRespPCA$x[,1],DeCorKnnRespPCA$x[,1])
# cor.test(KnnRespPCA$x[,1],DeCorKnnRespPCA$x[,1])
# 
# CompFermPCA<-prcomp(CompFerm,scale=TRUE)
# KnnFermPCA<-prcomp(KnnFerm,scale=FALSE) ## scale is set to FALSE because the data is already centered and scaled when performing KNN imputation
# DeCorKnnFermPCA<-prcomp(DeCorKnnFerm,scale=FALSE)
# cor.test(CompFermPCA$x[,1],KnnFermPCA$x[,1])
# cor.test(CompFermPCA$x[,1],DeCorKnnFermPCA$x[,1])
# cor.test(KnnFermPCA$x[,1],DeCorKnnFermPCA$x[,1])
# 
# ### KNN imputation does not change the data dramatically
# ### whether imputation happens before or after removing correlated predictors
# 
