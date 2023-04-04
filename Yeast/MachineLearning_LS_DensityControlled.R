library(data.table)
library(caret)
library(gridExtra)

## Use density as a covariate in the model, 
## since density is highly correlated with omics

Pheno<-fread("Phenotype_Reordered.txt")
load("KNN_Omics.RData")

## Add density to the omics dataset

RespProt<-cbind(RespDen=Pheno$RespDensityFC,knnRespProt)
RespMet<-cbind(RespDen=Pheno$RespDensityFC,knnRespMet)
RespLip<-cbind(RespDen=Pheno$RespDensityFC,knnRespLip)

FermProt<-cbind(FermDen=Pheno$FermDensityFC,knnFermProt)
FermMet<-cbind(FermDen=Pheno$FermDensityFC,knnFermMet)
FermLip<-cbind(FermDen=Pheno$FermDensityFC,knnFermLip)

#############

## Note that there are missing data in phenotype
### All 168 strains have CLS, but 17 of them do not have FermDensityFC, RespDensityFC, MeanRLS and MaxRLS
### In order for ML to run properly, we need to remove the 17 lines with missing phenotype data

Pheno151<-Pheno[complete.cases(Pheno)]

RespProt<-RespProt[complete.cases(Pheno),]
RespMet<-RespMet[complete.cases(Pheno),]
RespLip<-RespLip[complete.cases(Pheno),]

FermProt<-FermProt[complete.cases(Pheno),]
FermMet<-FermMet[complete.cases(Pheno),]
FermLip<-FermLip[complete.cases(Pheno),]

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


RespProt_CLS.dc.pls<-PLS(RespProt,Pheno151$CLS)
RespProt_CLS.dc.pls

RespMet_CLS.dc.pls<-PLS(RespMet,Pheno151$CLS)
RespMet_CLS.dc.pls

RespLip_CLS.dc.pls<-PLS(RespLip,Pheno151$CLS)
RespLip_CLS.dc.pls


a<-plot(RespProt_CLS.dc.pls[[1]],main="Respiration Protein_CLS")
b<-plot(RespMet_CLS.dc.pls[[1]],main="Respiration Metabolite_CLS")
c<-plot(RespLip_CLS.dc.pls[[1]],main="Respiration Lipid_CLS")
grid.arrange(a,b,c,ncol=2)

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

RespProt_CLS.dc.pls.final<-PLS.final(RespProt,Pheno151$CLS,3)
RespProt_CLS.dc.pls.FeatureImportance<-FeatureImportance(RespProt_CLS.dc.pls.final)

RespMet_CLS.dc.pls.final<-PLS.final(RespMet,Pheno151$CLS,1)
RespMet_CLS.dc.pls.FeatureImportance<-FeatureImportance(RespMet_CLS.dc.pls.final)

RespLip_CLS.dc.pls.final<-PLS.final(RespLip,Pheno151$CLS,1)
RespLip_CLS.dc.pls.FeatureImportance<-FeatureImportance(RespLip_CLS.dc.pls.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(RespProt_CLS.dc.pls.final,newdata=RespProt)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespMet_CLS.dc.pls.final,newdata=RespMet)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespMet",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespLip_CLS.dc.pls.final,newdata=RespLip)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespLip",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

### 

FermProt_RLS.dc.pls<-PLS(FermProt,Pheno151$MeanRLS)
FermProt_RLS.dc.pls ## ncomp=2

FermMet_RLS.dc.pls<-PLS(FermMet,Pheno151$MeanRLS)
FermMet_RLS.dc.pls ## ncomp=1

FermLip_RLS.dc.pls<-PLS(FermLip,Pheno151$MeanRLS)
FermLip_RLS.dc.pls ## ncomp=3


a<-plot(FermProt_RLS.dc.pls[[1]],main="Fermentation Protein_Mean RLS")
b<-plot(FermMet_RLS.dc.pls[[1]],main="Fermentation Metabolite_Mean RLS")
c<-plot(FermLip_RLS.dc.pls[[1]],main="Fermentatino Lipid_Mean RLS")
grid.arrange(a,b,c,ncol=2)

FermProt_RLS.dc.pls.final<-PLS.final(FermProt,Pheno151$MeanRLS,2)
FermProt_RLS.dc.pls.FeatureImportance<-FeatureImportance(FermProt_RLS.dc.pls.final)

FermMet_RLS.dc.pls.final<-PLS.final(FermMet,Pheno151$MeanRLS,3)
FermMet_RLS.dc.pls.FeatureImportance<-FeatureImportance(FermMet_RLS.dc.pls.final)

FermLip_RLS.dc.pls.final<-PLS.final(FermLip,Pheno151$MeanRLS,3)
FermLip_RLS.dc.pls.FeatureImportance<-FeatureImportance(FermLip_RLS.dc.pls.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(FermProt_RLS.dc.pls.final,newdata=FermProt)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermMet_RLS.dc.pls.final,newdata=FermMet)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermLip_RLS.dc.pls.final,newdata=FermLip)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

save(RespProt_CLS.dc.pls.final,RespProt_CLS.dc.pls.FeatureImportance,
     RespMet_CLS.dc.pls.final,RespMet_CLS.dc.pls.FeatureImportance,
     RespLip_CLS.dc.pls.final,RespLip_CLS.dc.pls.FeatureImportance,
     FermProt_RLS.dc.pls.final,FermProt_RLS.dc.pls.FeatureImportance,
     FermMet_RLS.dc.pls.final,FermMet_RLS.dc.pls.FeatureImportance,
     FermLip_RLS.dc.pls.final,FermLip_RLS.dc.pls.FeatureImportance,
     file="PLS_DensityControlled.RData")

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

RespProt_CLS.dc.rf<-RF(RespProt,Pheno151$CLS)
RespProt_CLS.dc.rf

RespMet_CLS.dc.rf<-RF(RespMet,Pheno151$CLS)
RespMet_CLS.dc.rf

RespLip_CLS.dc.rf<-RF(RespLip,Pheno151$CLS)
RespLip_CLS.dc.rf

a<-plot(RespProt_CLS.dc.rf[[1]],main="Respiration Protein_CLS")
b<-plot(RespMet_CLS.dc.rf[[1]],main="Respiration Metabolite_CLS")
c<-plot(RespLip_CLS.dc.rf[[1]],main="Respiration Lipid_CLS")
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
  .mtry = c(500,550,600,650,700),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2000)
RespProt_CLS.dc.rf.fine<-RF.fine(RespProt,Pheno151$CLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(80,85,90,95,100),
  .splitrule = "variance",
  .min.node.size = 5
)
set.seed(2005)
RespMet_CLS.dc.rf.fine<-RF.fine(RespMet,Pheno151$CLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(1,2,5,10,15,20),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2005)
RespLip_CLS.dc.rf.fine<-RF.fine(RespLip,Pheno151$CLS,tuneGrid)



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

tuneGrid <- data.frame(.mtry = 500,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
RespProt_CLS.dc.rf.final<-RF.final(RespProt,Pheno151$CLS,tuneGrid)
RespProt_CLS.dc.rf.FeatureImportance<-FeatureImportance(RespProt_CLS.dc.rf.final)

tuneGrid <- data.frame(.mtry = 90,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
RespMet_CLS.dc.rf.final<-RF.final(RespMet,Pheno151$CLS,tuneGrid)
RespMet_CLS.dc.rf.FeatureImportance<-FeatureImportance(RespMet_CLS.dc.rf.final)

tuneGrid <- data.frame(.mtry = 2,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
RespLip_CLS.dc.rf.final<-RF.final(RespLip,Pheno151$CLS,tuneGrid)
RespLip_CLS.dc.rf.FeatureImportance<-FeatureImportance(RespLip_CLS.dc.rf.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(RespProt_CLS.dc.rf.final,newdata=RespProt)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespMet_CLS.dc.rf.final,newdata=RespMet)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespMet",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespLip_CLS.dc.rf.final,newdata=RespLip)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespLip",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

### 

FermProt_RLS.dc.rf<-RF(FermProt,Pheno151$MeanRLS)
FermProt_RLS.dc.rf

FermMet_RLS.dc.rf<-RF(FermMet,Pheno151$MeanRLS)
FermMet_RLS.dc.rf

FermLip_RLS.dc.rf<-RF(FermLip,Pheno151$MeanRLS)
FermLip_RLS.dc.rf

a<-plot(FermProt_RLS.dc.rf[[1]],main="Fermentation Protein_RLS")
b<-plot(FermMet_RLS.dc.rf[[1]],main="Fermentation Metabolite_RLS")
c<-plot(FermLip_RLS.dc.rf[[1]],main="Fermentation Lipid_RLS")
grid.arrange(a,b,c,ncol=2)

## Fine tune 

tuneGrid <- data.frame(
  .mtry = c(1600,1650,1700,1750,1800),
  .splitrule = "variance",
  .min.node.size = 5
)
set.seed(2000)
FermProt_RLS.dc.rf.fine<-RF.fine(FermProt,Pheno151$MeanRLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(20,40,50,60,70,80,90,100),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2002)
FermMet_RLS.dc.rf.fine<-RF.fine(FermMet,Pheno151$MeanRLS,tuneGrid)

tuneGrid <- data.frame(
  .mtry = c(1,2,3,4,5,7,9,10),
  .splitrule = "extratrees",
  .min.node.size = 5
)
set.seed(2000)
FermLip_RLS.dc.rf.fine<-RF.fine(FermLip,Pheno151$MeanRLS,tuneGrid)

## Final model

tuneGrid <- data.frame(.mtry = 1650,.splitrule = "variance",.min.node.size = 5)
set.seed(1234)
FermProt_RLS.dc.rf.final<-RF.final(FermProt,Pheno151$MeanRLS,tuneGrid)
FermProt_RLS.dc.rf.FeatureImportance<-FeatureImportance(FermProt_RLS.dc.rf.final)

tuneGrid <- data.frame(.mtry = 90,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
FermMet_RLS.dc.rf.final<-RF.final(FermMet,Pheno151$MeanRLS,tuneGrid)
FermMet_RLS.dc.rf.FeatureImportance<-FeatureImportance(FermMet_RLS.dc.rf.final)

tuneGrid <- data.frame(.mtry = 1,.splitrule = "extratrees",.min.node.size = 5)
set.seed(1234)
FermLip_RLS.dc.rf.final<-RF.final(FermLip,Pheno151$MeanRLS,tuneGrid)
FermLip_RLS.dc.rf.FeatureImportance<-FeatureImportance(FermLip_RLS.dc.rf.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(FermProt_RLS.dc.rf.final,newdata=FermProt)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermMet_RLS.dc.rf.final,newdata=FermMet)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermLip_RLS.dc.rf.final,newdata=FermLip)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

save(RespProt_CLS.dc.rf.final,RespProt_CLS.dc.rf.FeatureImportance,
     RespMet_CLS.dc.rf.final,RespMet_CLS.dc.rf.FeatureImportance,
     RespLip_CLS.dc.rf.final,RespLip_CLS.dc.rf.FeatureImportance,
     FermProt_RLS.dc.rf.final,FermProt_RLS.dc.rf.FeatureImportance,
     FermMet_RLS.dc.rf.final,FermMet_RLS.dc.rf.FeatureImportance,
     FermLip_RLS.dc.rf.final,FermLip_RLS.dc.rf.FeatureImportance,
     file="RF_DensityControlled.RData")

##################

### Elastic Net

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

RespProt_CLS.dc.en<-EN(RespProt,Pheno151$CLS)
RespProt_CLS.dc.en

RespMet_CLS.dc.en<-EN(RespMet,Pheno151$CLS)
RespMet_CLS.dc.en

RespLip_CLS.dc.en<-EN(RespLip,Pheno151$CLS)
RespLip_CLS.dc.en

a<-plot(RespProt_CLS.dc.en[[1]],main="Respiration Protein_CLS")
b<-plot(RespMet_CLS.dc.en[[1]],main="Respiration Metabolite_CLS")
c<-plot(RespLip_CLS.dc.en[[1]],main="Respiration Lipid_CLS")
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

tuneGrid = expand.grid(alpha=0,lambda=355)
RespProt_CLS.dc.en.final<-EN.final(RespProt,Pheno151$CLS,tuneGrid)
RespProt_CLS.dc.en.FeatureImportance<-FeatureImportance(RespProt_CLS.dc.en.final)

tuneGrid = expand.grid(alpha=1,lambda=1)
RespMet_CLS.dc.en.final<-EN.final(RespMet,Pheno151$CLS,tuneGrid)
RespMet_CLS.dc.en.FeatureImportance<-FeatureImportance(RespMet_CLS.dc.en.final)

tuneGrid = expand.grid(alpha=0,lambda=28)
RespLip_CLS.dc.en.final<-EN.final(RespLip,Pheno151$CLS,tuneGrid)
RespLip_CLS.dc.en.FeatureImportance<-FeatureImportance(RespLip_CLS.dc.en.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(RespProt_CLS.dc.en.final,newdata=RespProt)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespMet_CLS.dc.en.final,newdata=RespMet)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespMet",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(RespLip_CLS.dc.en.final,newdata=RespLip)
real<-Pheno151$CLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="CLS predicted by RespLip",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

### 

FermProt_RLS.dc.en<-EN(FermProt,Pheno151$MeanRLS)
FermProt_RLS.dc.en

FermMet_RLS.dc.en<-EN(FermMet,Pheno151$MeanRLS)
FermMet_RLS.dc.en

FermLip_RLS.dc.en<-EN(FermLip,Pheno151$MeanRLS)
FermLip_RLS.dc.en

a<-plot(FermProt_RLS.dc.en[[1]],main="Fermentation Protein_RLS")
b<-plot(FermMet_RLS.dc.en[[1]],main="Fermentation Metabolite_RLS")
c<-plot(FermLip_RLS.dc.en[[1]],main="Fermentation Lipid_RLS")
grid.arrange(a,b,c,ncol=2)

tuneGrid = expand.grid(alpha=0,lambda=128)
FermProt_RLS.dc.en.final<-EN.final(FermProt,Pheno151$MeanRLS,tuneGrid)
FermProt_RLS.dc.en.FeatureImportance<-FeatureImportance(FermProt_RLS.dc.en.final)

tuneGrid = expand.grid(alpha=0,lambda=19)
FermMet_RLS.dc.en.final<-EN.final(FermMet,Pheno151$MeanRLS,tuneGrid)
FermMet_RLS.dc.en.FeatureImportance<-FeatureImportance(FermMet_RLS.dc.en.final)

tuneGrid = expand.grid(alpha=0.5,lambda=0.2)
FermLip_RLS.dc.en.final<-EN.final(FermLip,Pheno151$MeanRLS,tuneGrid)
FermLip_RLS.dc.en.FeatureImportance<-FeatureImportance(FermLip_RLS.dc.en.final)

## Visualize the prediction 

par(mfrow=c(2,2))

pred<-predict(FermProt_RLS.dc.en.final,newdata=FermProt)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermMet_RLS.dc.en.final,newdata=FermMet)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

pred<-predict(FermLip_RLS.dc.en.final,newdata=FermLip)
real<-Pheno151$MeanRLS
plot(x=real,y=pred,pch=20,xlab="Real Phenotype Data",ylab="Predicted Value",
     main="RLS predicted by FermProt",ylim=c(min(real),max(real)))
abline(0,1,col="red")

invisible(dev.off())

save(RespProt_CLS.dc.en.final,RespProt_CLS.dc.en.FeatureImportance,
     RespMet_CLS.dc.en.final,RespMet_CLS.dc.en.FeatureImportance,
     RespLip_CLS.dc.en.final,RespLip_CLS.dc.en.FeatureImportance,
     FermProt_RLS.dc.en.final,FermProt_RLS.dc.en.FeatureImportance,
     FermMet_RLS.dc.en.final,FermMet_RLS.dc.en.FeatureImportance,
     FermLip_RLS.dc.en.final,FermLip_RLS.dc.en.FeatureImportance,
     file="EN_DensityControlled.RData")

