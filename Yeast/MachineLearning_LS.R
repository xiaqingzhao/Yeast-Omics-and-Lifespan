library(data.table)
library(caret)
library(gridExtra)

Pheno<-fread("Phenotype_Reordered.txt")

load("KNN_Omics.RData")
knnRespProt<-scale(knnRespProt)
knnRespMet<-scale(knnRespMet)
knnRespLip<-scale(knnRespLip)

#############

## Note that there are missing data in y (phenotype)
# colSums(is.na(Pheno))
# table(rowSums(is.na(Pheno)))
### All 168 strains have CLS, but 17 of them do not have FermDensityFC, RespDensityFC, MeanRLS and MaxRLS
### In order for ML to run properly, we need to remove the 17 lines with missing phenotype data

Pheno151<-Pheno[complete.cases(Pheno)]

# RespProt151<-knnRespProt[complete.cases(Pheno),]
# RespMet151<-knnRespMet[complete.cases(Pheno),]
# RespLip151<-knnRespLip[complete.cases(Pheno),]
FermProt151<-knnFermProt[complete.cases(Pheno),]
FermMet151<-knnFermMet[complete.cases(Pheno),]
FermLip151<-knnFermLip[complete.cases(Pheno),]


##############

## Partial least squares

PLS<-function(Omics,Pheno) {
  
  ## Split data into training and test sets
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## Tune grid
  tuneGrid <- data.frame(
    .ncomp = c(1,2,3,4,5,10,15,20)
  )
  
  ## Model fitting
  model<-train(
    Omics_train, 
    Pheno_train,
    method = "pls",  ## can be any methods available
    tuneGrid = tuneGrid,
    trControl = trainControl(
      method = "cv", 
      number = 5,  ## 5-fold cross validation
      verboseIter = TRUE
    )
  )
  
  error<-Pheno_test-predict(model,newdata=Omics_test)
  RMSE<-sqrt(mean(error^2)) 
  
  return(list(Model=model,test_RMSE=RMSE))
}


RespProt_CLS.pls<-PLS(knnRespProt,Pheno$CLS)
RespProt_CLS.pls

RespMet_CLS.pls<-PLS(knnRespMet,Pheno$CLS)
RespMet_CLS.pls

RespLip_CLS.pls<-PLS(knnRespLip,Pheno$CLS)
RespLip_CLS.pls


a<-plot(RespProt_CLS.pls[[1]],main="Respiration Protein_CLS")
b<-plot(RespMet_CLS.pls[[1]],main="Respiration Metabolite_CLS")
c<-plot(RespLip_CLS.pls[[1]],main="Respiration Lipid_CLS")
grid.arrange(a,b,c,ncol=2)

## ncomp=1 for protein, metabolite and lipid. 
## Now use all data and ncomp=1 to construct the pls model

PLS.final<-function(Omics,Pheno,ncomp) {
  tuneGrid <- data.frame(.ncomp = c(ncomp))
  model<-train(
    Omics, 
    Pheno,
    method = "pls",
    tuneGrid = tuneGrid,
    trControl = trainControl(
      method = "cv", 
      number = 5,  ## 5-fold cross validation
      verboseIter = TRUE
    )
  )
}


## Rank features by their importance in the model 

FeatureImportance<-function(model) {
  Importance<-varImp(model)
  FeatureImportance<-Importance$importance
  FeatureImportance<-FeatureImportance[order(-FeatureImportance),,drop=FALSE]
  return(FeatureImportance)
}

RespProt_CLS.pls.final<-PLS.final(knnRespProt,Pheno$CLS,1)
RespProt_CLS.pls.FeatureImportance<-FeatureImportance(RespProt_CLS.pls.final)

RespMet_CLS.pls.final<-PLS.final(knnRespMet,Pheno$CLS,1)
RespMet_CLS.pls.FeatureImportance<-FeatureImportance(RespMet_CLS.pls.final)

RespLip_CLS.pls.final<-PLS.final(knnRespLip,Pheno$CLS,1)
RespLip_CLS.pls.FeatureImportance<-FeatureImportance(RespLip_CLS.pls.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(RespProt_CLS.pls.final,newdata=knnRespProt)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespMet_CLS.pls.final,newdata=knnRespMet)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespMet",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespLip_CLS.pls.final,newdata=knnRespLip)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespLip",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

###

# pls.fit<-plsr(Pheno_train~Resp_train, ncomp=2, scale=TRUE,validation="CV")
# summary(pls.fit)

### 

FermProt_RLS.pls<-PLS(FermProt151,Pheno151$MeanRLS)
FermProt_RLS.pls ## ncomp=2

FermMet_RLS.pls<-PLS(FermMet151,Pheno151$MeanRLS)
FermMet_RLS.pls ## ncomp=3

FermLip_RLS.pls<-PLS(FermLip151,Pheno151$MeanRLS)
FermLip_RLS.pls ## ncomp=3


a<-plot(FermProt_RLS.pls[[1]],main="Fermentation Protein_Mean RLS")
b<-plot(FermMet_RLS.pls[[1]],main="Fermentation Metabolite_Mean RLS")
c<-plot(FermLip_RLS.pls[[1]],main="Fermentatino Lipid_Mean RLS")
grid.arrange(a,b,c,ncol=2)

FermProt_RLS.pls.final<-PLS.final(FermProt151,Pheno151$MeanRLS,2)
FermProt_RLS.pls.FeatureImportance<-FeatureImportance(FermProt_RLS.pls.final)

FermMet_RLS.pls.final<-PLS.final(FermMet151,Pheno151$MeanRLS,3)
FermMet_RLS.pls.FeatureImportance<-FeatureImportance(FermMet_RLS.pls.final)

FermLip_RLS.pls.final<-PLS.final(FermLip151,Pheno151$MeanRLS,3)
FermLip_RLS.pls.FeatureImportance<-FeatureImportance(FermLip_RLS.pls.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(FermProt_RLS.pls.final,newdata=FermProt151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermMet_RLS.pls.final,newdata=FermMet151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermLip_RLS.pls.final,newdata=FermLip151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

save(RespProt_CLS.pls.final,RespProt_CLS.pls.FeatureImportance,
     RespMet_CLS.pls.final,RespMet_CLS.pls.FeatureImportance,
     RespLip_CLS.pls.final,RespLip_CLS.pls.FeatureImportance,
     FermProt_RLS.pls.final,FermProt_RLS.pls.FeatureImportance,
     FermMet_RLS.pls.final,FermMet_RLS.pls.FeatureImportance,
     FermLip_RLS.pls.final,FermLip_RLS.pls.FeatureImportance,
     file="PLS.RData")

##########################

### Random forest

RF<-function(Omics,Pheno) {
  
  ## Split data into training and test sets
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## Model fitting 
  model<-train(
    Omics_train, 
    Pheno_train,
    tuneLength = 10,
    method = "ranger", ## A wrapper of rf
    trControl = trainControl(
      method = "cv", 
      number = 5,  ## 5-fold cross validation
      verboseIter = TRUE
    )
  )
  
  error<-Pheno_test-predict(model,newdata=Omics_test)
  RMSE<-sqrt(mean(error^2)) 
  
  return(list(Model=model,test_RMSE=RMSE))
}

RespProt_CLS.rf<-RF(knnRespProt,Pheno$CLS)
RespProt_CLS.rf

RespMet_CLS.rf<-RF(knnRespMet,Pheno$CLS)
RespMet_CLS.rf

RespLip_CLS.rf<-RF(knnRespLip,Pheno$CLS)
RespLip_CLS.rf

a<-plot(RespProt_CLS.rf[[1]],main="Respiration Protein_CLS")
b<-plot(RespMet_CLS.rf[[1]],main="Respiration Metabolite_CLS")
c<-plot(RespLip_CLS.rf[[1]],main="Respiration Lipid_CLS")
grid.arrange(a,b,c,ncol=2)

## Fine tuning

RF.fine<-function(Omics,Pheno,tuneGrid) {
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  model <- train(
    Omics_train, 
    Pheno_train,
    tuneGrid = tuneGrid,
    method = "ranger",
    num.trees = 1000,
    trControl = trainControl(
      method = "cv", 
      number = 5, 
      verboseIter = TRUE
    )
  )
  error<-Pheno_test-predict(model,newdata=Omics_test)
  RMSE<-sqrt(mean(error^2)) 
  return(list(Model=model,test_RMSE=RMSE))
}


tuneGrid <- data.frame(
  .mtry = c(40,50,60,70,80,90,100,200),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2003)
RespProt_CLS.rf.fine<-RF.fine(knnRespProt,Pheno$CLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(10,15,20,30,40,50,60,70,80,90,100,150,200),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2008)
RespMet_CLS.rf.fine<-RF.fine(knnRespMet,Pheno$CLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(3,5,7,10,15,20,30,40),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2015)
RespLip_CLS.rf.fine<-RF.fine(knnRespLip,Pheno$CLS,tuneGrid)


## Final model

RF.final<-function(Omics,Pheno,tuneGrid) {
  model <- train(
    Omics, 
    Pheno,
    tuneGrid = tuneGrid,
    method = "ranger",
    trControl = trainControl(
      method = "cv", 
      number = 5, 
      verboseIter = TRUE
    ),
    importance="impurity"
  )
}

tuneGrid <- data.frame(.mtry = 40,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
RespProt_CLS.rf.final<-RF.final(knnRespProt,Pheno$CLS,tuneGrid)
RespProt_CLS.rf.FeatureImportance<-FeatureImportance(RespProt_CLS.rf.final)

tuneGrid <- data.frame(.mtry = 80,.splitrule = "extratrees",.min.node.size = 5)
set.seed(4321)
RespMet_CLS.rf.final<-RF.final(knnRespMet,Pheno$CLS,tuneGrid)
RespMet_CLS.rf.FeatureImportance<-FeatureImportance(RespMet_CLS.rf.final)

tuneGrid <- data.frame(.mtry = 5,.splitrule = "extratrees",.min.node.size = 5)
set.seed(5678)
RespLip_CLS.rf.final<-RF.final(knnRespLip,Pheno$CLS,tuneGrid)
RespLip_CLS.rf.FeatureImportance<-FeatureImportance(RespLip_CLS.rf.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(RespProt_CLS.rf.final,newdata=knnRespProt)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespMet_CLS.rf.final,newdata=knnRespMet)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespMet",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespLip_CLS.rf.final,newdata=knnRespLip)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespLip",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

### 

FermProt_RLS.rf<-RF(FermProt151,Pheno151$MeanRLS)
FermProt_RLS.rf

FermMet_RLS.rf<-RF(FermMet151,Pheno151$MeanRLS)
FermMet_RLS.rf

FermLip_RLS.rf<-RF(FermLip151,Pheno151$MeanRLS)
FermLip_RLS.rf

a<-plot(FermProt_RLS.rf[[1]],main="Fermentation Protein_RLS")
b<-plot(FermMet_RLS.rf[[1]],main="Fermentation Metabolite_RLS")
c<-plot(FermLip_RLS.rf[[1]],main="Fermentation Lipid_RLS")
grid.arrange(a,b,c,ncol=2)

## Fine tuning 

tuneGrid <- data.frame(
  .mtry = c(20,30,40,50,60,70,80,90,100,200),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2001)
FermProt_RLS.rf.fine<-RF.fine(FermProt151,Pheno151$MeanRLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(15,20,40,50,60,70,80,100,150),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2005)
FermMet_RLS.rf.fine<-RF.fine(FermMet151,Pheno151$MeanRLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(4,5,6,7,8,9,10,15,20),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2000)
FermLip_RLS.rf.fine<-RF.fine(FermLip151,Pheno151$MeanRLS,tuneGrid)

## Final model

tuneGrid <- data.frame(.mtry = 100,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
FermProt_RLS.rf.final<-RF.final(FermProt151,Pheno151$MeanRLS,tuneGrid)
FermProt_RLS.rf.FeatureImportance<-FeatureImportance(FermProt_RLS.rf.final)

tuneGrid <- data.frame(.mtry = 60,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
FermMet_RLS.rf.final<-RF.final(FermMet151,Pheno151$MeanRLS,tuneGrid)
FermMet_RLS.rf.FeatureImportance<-FeatureImportance(FermMet_RLS.rf.final)

tuneGrid <- data.frame(.mtry = 6,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
FermLip_RLS.rf.final<-RF.final(FermLip151,Pheno151$MeanRLS,tuneGrid)
FermLip_RLS.rf.FeatureImportance<-FeatureImportance(FermLip_RLS.rf.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(FermProt_RLS.rf.final,newdata=FermProt151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermMet_RLS.rf.final,newdata=FermMet151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermLip_RLS.rf.final,newdata=FermLip151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

save(RespProt_CLS.rf.final,RespProt_CLS.rf.FeatureImportance,
     RespMet_CLS.rf.final,RespMet_CLS.rf.FeatureImportance,
     RespLip_CLS.rf.final,RespLip_CLS.rf.FeatureImportance,
     FermProt_RLS.rf.final,FermProt_RLS.rf.FeatureImportance,
     FermMet_RLS.rf.final,FermMet_RLS.rf.FeatureImportance,
     FermLip_RLS.rf.final,FermLip_RLS.rf.FeatureImportance,
     file="RF.RData")


##################

### Elastic Net
#### In most cases alpha=0 works better than any alpha>0. So it's ridge regression. 

EN<-function(Omics,Pheno) {
  
  ## Split data into training and test sets
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## Tune grid
  tuneGrid = expand.grid(
    alpha=c(0,0.5,1),
    lambda=c(seq(0.000001,1,length=100),seq(1,1000,length=1001))
  )
  
  ## Model fitting 
  model<-train(
    Omics_train, 
    Pheno_train,
    tuneGrid = tuneGrid,
    method = "glmnet", 
    trControl = trainControl(
      method = "cv", 
      number = 5,  ## 5-fold cross validation
      verboseIter = TRUE
    )
  )
  
  error<-Pheno_test-predict(model,newdata=Omics_test)
  RMSE<-sqrt(mean(error^2)) 
  
  return(list(Model=model,test_RMSE=RMSE))
}

RespProt_CLS.en<-EN(knnRespProt,Pheno$CLS)
RespProt_CLS.en

RespMet_CLS.en<-EN(knnRespMet,Pheno$CLS)
RespMet_CLS.en

RespLip_CLS.en<-EN(knnRespLip,Pheno$CLS)
RespLip_CLS.en

a<-plot(RespProt_CLS.en[[1]],main="Respiration Protein_CLS")
b<-plot(RespMet_CLS.en[[1]],main="Respiration Metabolite_CLS")
c<-plot(RespLip_CLS.en[[1]],main="Respiration Lipid_CLS")
grid.arrange(a,b,c,ncol=2)

## Final model 

EN.final<-function(Omics,Pheno,tuneGrid) {
  model <- train(
    Omics, 
    Pheno,
    tuneGrid = tuneGrid,
    method = "glmnet",
    trControl = trainControl(
      method = "cv", 
      number = 5, 
      verboseIter = TRUE
    )
  )
}

tuneGrid = expand.grid(alpha=0,lambda=584)
RespProt_CLS.en.final<-EN.final(knnRespProt,Pheno$CLS,tuneGrid)
RespProt_CLS.en.FeatureImportance<-FeatureImportance(RespProt_CLS.en.final)

tuneGrid = expand.grid(alpha=1,lambda=2)
RespMet_CLS.en.final<-EN.final(knnRespMet,Pheno$CLS,tuneGrid)
RespMet_CLS.en.FeatureImportance<-FeatureImportance(RespMet_CLS.en.final)

tuneGrid = expand.grid(alpha=0,lambda=38)
RespLip_CLS.en.final<-EN.final(knnRespLip,Pheno$CLS,tuneGrid)
RespLip_CLS.en.FeatureImportance<-FeatureImportance(RespLip_CLS.en.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(RespProt_CLS.en.final,newdata=knnRespProt)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespMet_CLS.en.final,newdata=knnRespMet)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespMet",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespLip_CLS.en.final,newdata=knnRespLip)
real<-Pheno$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespLip",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

### 

FermProt_RLS.en<-EN(FermProt151,Pheno151$MeanRLS)
FermProt_RLS.en

FermMet_RLS.en<-EN(FermMet151,Pheno151$MeanRLS)
FermMet_RLS.en

FermLip_RLS.en<-EN(FermLip151,Pheno151$MeanRLS)
FermLip_RLS.en

a<-plot(FermProt_RLS.en[[1]],main="Fermentation Protein_RLS")
b<-plot(FermMet_RLS.en[[1]],main="Fermentation Metabolite_RLS")
c<-plot(FermLip_RLS.en[[1]],main="Fermentation Lipid_RLS")
grid.arrange(a,b,c,ncol=2)

tuneGrid = expand.grid(alpha=0,lambda=128)
FermProt_RLS.en.final<-EN.final(FermProt151,Pheno151$MeanRLS,tuneGrid)
FermProt_RLS.en.FeatureImportance<-FeatureImportance(FermProt_RLS.en.final)

tuneGrid = expand.grid(alpha=0,lambda=19)
FermMet_RLS.en.final<-EN.final(FermMet151,Pheno151$MeanRLS,tuneGrid)
FermMet_RLS.en.FeatureImportance<-FeatureImportance(FermMet_RLS.en.final)

tuneGrid = expand.grid(alpha=0.5,lambda=1000)
FermLip_RLS.en.final<-EN.final(FermLip151,Pheno151$MeanRLS,tuneGrid)
FermLip_RLS.en.FeatureImportance<-FeatureImportance(FermLip_RLS.en.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(FermProt_RLS.en.final,newdata=FermProt151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermMet_RLS.en.final,newdata=FermMet151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermLip_RLS.en.final,newdata=FermLip151)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

save(RespProt_CLS.en.final,RespProt_CLS.en.FeatureImportance,
     RespMet_CLS.en.final,RespMet_CLS.en.FeatureImportance,
     RespLip_CLS.en.final,RespLip_CLS.en.FeatureImportance,
     FermProt_RLS.en.final,FermProt_RLS.en.FeatureImportance,
     FermMet_RLS.en.final,FermMet_RLS.en.FeatureImportance,
     FermLip_RLS.en.final,FermLip_RLS.en.FeatureImportance,
     file="EN.RData")


#####################

## Compare different models

RespProt_CLS_models<-list(RespProt_CLS.pls.final,RespProt_CLS.rf.final,RespProt_CLS.en.final)
RespProt_CLS_resamples<-resamples(RespProt_CLS_models)
summary(RespProt_CLS_resamples)
bwplot(RespProt_CLS_resamples,metric="RMSE")

FermProt_RLS_models<-list(FermProt_RLS.pls.final,FermProt_RLS.rf.final,FermProt_RLS.en.final)
FermProt_RLS_resamples<-resamples(FermProt_RLS_models)
summary(FermProt_RLS_resamples)
bwplot(FermProt_RLS_resamples,metric="RMSE")








