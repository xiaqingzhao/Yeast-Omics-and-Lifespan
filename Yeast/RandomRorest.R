library(data.table)
library(caret)
library(gridExtra)
library(DMwR)

Pheno<-fread("Phenotype_Reordered.txt")

load("DeMissing_Omics.RData")

## Note that there are missing data in y (phenotype)
# colSums(is.na(Pheno))
# table(rowSums(is.na(Pheno)))
### All 168 strains have CLS, but 17 of them do not have FermDensityFC, RespDensityFC, MeanRLS and MaxRLS
### In order for ML to run properly, we need to remove the 17 lines with missing phenotype data

Pheno151<-Pheno[complete.cases(Pheno)]
FermProt151<-FermProt[complete.cases(Pheno),]
FermMet151<-FermMet[complete.cases(Pheno),]
FermLip151<-FermLip[complete.cases(Pheno),]

#############

RF<-function(Omics,Pheno) {
  
  ## Split data into training and test sets
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## KNN imputation within training and testing set
  
  knnOmics_train<-knnImputation(Omics_train,scale=F)
  knnOmics_test<-knnImputation(Omics_test,scale=F)
  
  ## Model fitting 
  set.seed(1234)
  model<-train(
    knnOmics_train, 
    Pheno_train,
    tuneLength = 10,
    method = "ranger", ## A wrapper of rf
    num.trees=1000,
    trControl = trainControl(
      method = "cv", 
      number = 5,  ## 5-fold cross validation
      verboseIter = TRUE
    )
  )
  
  Pred<-predict(model,newdata=knnOmics_test)
  error<-Pheno_test-Pred
  RMSE<-sqrt(mean(error^2)) 
  
  return(list(Model=model,test_RMSE=RMSE,TrueVSPred=cbind(Pheno_test,Pred)))
}

RespProt_CLS.rf<-RF(RespProt,Pheno$CLS)
RespProt_CLS.rf

RespMet_CLS.rf<-RF(RespMet,Pheno$CLS)
RespMet_CLS.rf

RespLip_CLS.rf<-RF(RespLip,Pheno$CLS)
RespLip_CLS.rf

plot(RespProt_CLS.rf[[3]][,1],RespProt_CLS.rf[[3]][,2],pch=20)
abline(0,1,col="blue")
abline(lm(RespProt_CLS.rf[[3]][,2]~RespProt_CLS.rf[[3]][,1]),col="red")
cor.test(RespProt_CLS.rf[[3]][,1],RespProt_CLS.rf[[3]][,2])

plot(RespMet_CLS.rf[[3]][,1],RespMet_CLS.rf[[3]][,2],pch=20)
abline(0,1,col="blue")
abline(lm(RespMet_CLS.rf[[3]][,2]~RespMet_CLS.rf[[3]][,1]),col="red")
cor.test(RespMet_CLS.rf[[3]][,1],RespMet_CLS.rf[[3]][,2])

plot(RespLip_CLS.rf[[3]][,1],RespLip_CLS.rf[[3]][,2],pch=20)
abline(0,1,col="blue")
abline(lm(RespLip_CLS.rf[[3]][,2]~RespLip_CLS.rf[[3]][,1]),col="red")
cor.test(RespLip_CLS.rf[[3]][,1],RespLip_CLS.rf[[3]][,2])

a<-plot(RespProt_CLS.rf[[1]],main="Respiration Protein_CLS")
b<-plot(RespMet_CLS.rf[[1]],main="Respiration Metabolite_CLS")
c<-plot(RespLip_CLS.rf[[1]],main="Respiration Lipid_CLS")
grid.arrange(a,b,c,ncol=2)

#######

## Fine tuning 

RF.fine<-function(Omics,Pheno,tuneGrid) {
  
  ## Split data into training and test sets
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## KNN imputation within training and testing set
  
  knnOmics_train<-knnImputation(Omics_train,scale=F)
  knnOmics_test<-knnImputation(Omics_test,scale=F)
  
  ## Model fitting 
  set.seed(4321)
  model<-train(
    knnOmics_train, 
    Pheno_train,
    tuneGrid = tuneGrid,
    method = "ranger", 
    num.trees=2000,
    trControl = trainControl(
      method = "cv", 
      number = 5,  
      verboseIter = TRUE
    ),
    importance="impurity"
  )
  
  Pred_test<-predict(model,newdata=knnOmics_test)
  error<-Pheno_test-Pred_test
  RMSE<-sqrt(mean(error^2)) 
  
  Pred_train<-predict(model,newdata=knnOmics_train)
  
  return(list(Model=model,test_RMSE=RMSE,TrueVSPred_Test=cbind(Pheno_test,Pred_test),
              TrueVSPred_Train=cbind(Pheno_train,Pred_train)))
}

tuneGrid <- data.frame(
  .mtry = c(20,30,40,45,50,60,70,80,90,100,200),
  .splitrule = "extratrees",
  .min.node.size = 5
)
RespProt_CLS.rf.fine<-RF.fine(RespProt,Pheno$CLS,tuneGrid)


tuneGrid <- data.frame(
  .mtry = c(10,15,20,30,40,50,60,70,80,90,100),
  .splitrule = "extratrees",
  .min.node.size = 5
)
RespMet_CLS.rf.fine<-RF.fine(RespMet,Pheno$CLS,tuneGrid)


tuneGrid <- data.frame(
  .mtry = c(3,5,7,10,15,20,30,40),
  .splitrule = "extratrees",
  .min.node.size = 5
)
RespLip_CLS.rf.fine<-RF.fine(RespLip,Pheno$CLS,tuneGrid)

# plot(RespProt_CLS.rf.fine[[1]])
par(mfrow=c(1,3))
plot(RespProt_CLS.rf.fine[[3]][,1],RespProt_CLS.rf.fine[[3]][,2],pch=20,
     xlab="Real Phenotype Data",ylab="Predicted Value",ylim=c(20,48),
     main="CLS predicted by RespProt")
abline(0,1,col="blue")
abline(lm(RespProt_CLS.rf.fine[[3]][,2]~RespProt_CLS.rf.fine[[3]][,1]),col="red")
cor.test(RespProt_CLS.rf.fine[[3]][,1],RespProt_CLS.rf.fine[[3]][,2])

# plot(RespMet_CLS.rf.fine[[1]])
plot(RespMet_CLS.rf.fine[[3]][,1],RespMet_CLS.rf.fine[[3]][,2],pch=20,
     xlab="Real Phenotype Data",ylab="Predicted Value",ylim=c(20,48),
     main="CLS predicted by RespMet")
abline(0,1,col="blue")
abline(lm(RespMet_CLS.rf.fine[[3]][,2]~RespMet_CLS.rf.fine[[3]][,1]),col="red")
cor.test(RespMet_CLS.rf.fine[[3]][,1],RespMet_CLS.rf.fine[[3]][,2])

# plot(RespLip_CLS.rf.fine[[1]])
plot(RespLip_CLS.rf.fine[[3]][,1],RespLip_CLS.rf.fine[[3]][,2],pch=20,
     xlab="Real Phenotype Data",ylab="Predicted Value",ylim=c(20,48),
     main="CLS predicted by RespLip")
abline(0,1,col="blue")
abline(lm(RespLip_CLS.rf.fine[[3]][,2]~RespLip_CLS.rf.fine[[3]][,1]),col="red")
cor.test(RespLip_CLS.rf.fine[[3]][,1],RespLip_CLS.rf.fine[[3]][,2])

invisible(dev.off())

a<-ggplot(data.frame(RespProt_CLS.rf.fine[[3]]),aes(x=Pheno_test,y=Pred_test))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(2, 56))+scale_y_continuous(limits = c(10, 50))+theme_bw()+
  xlab("Real CLS")+ylab("Predicted CLS")+geom_smooth(method="lm",col="blue")+ggtitle("Test set: CLS predicted by RespProt")
b<-ggplot(data.frame(RespMet_CLS.rf.fine[[3]]),aes(x=Pheno_test,y=Pred_test))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(2, 56))+scale_y_continuous(limits = c(10, 50))+theme_bw()+
  xlab("Real CLS")+ylab("Predicted CLS")+geom_smooth(method="lm",col="blue")+ggtitle("Test set: CLS predicted by RespMet")
c<-ggplot(data.frame(RespLip_CLS.rf.fine[[3]]),aes(x=Pheno_test,y=Pred_test))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(2, 56))+scale_y_continuous(limits = c(10, 50))+theme_bw()+
  xlab("Real CLS")+ylab("Predicted CLS")+geom_smooth(method="lm",col="blue")+ggtitle("Test set: CLS predicted by RespLip")
grid.arrange(a,b,c,ncol=3)

d<-ggplot(data.frame(RespProt_CLS.rf.fine[[4]]),aes(x=Pheno_train,y=Pred_train))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(2, 56))+scale_y_continuous(limits = c(10, 50))+theme_bw()+
  xlab("Real CLS")+ylab("Predicted CLS")+geom_smooth(method="lm",col="blue")+ggtitle("Training set: CLS predicted by RespProt")
e<-ggplot(data.frame(RespMet_CLS.rf.fine[[4]]),aes(x=Pheno_train,y=Pred_train))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(2, 56))+scale_y_continuous(limits = c(10, 50))+theme_bw()+
  xlab("Real CLS")+ylab("Predicted CLS")+geom_smooth(method="lm",col="blue")+ggtitle("Training set: CLS predicted by RespMet")
f<-ggplot(data.frame(RespLip_CLS.rf.fine[[4]]),aes(x=Pheno_train,y=Pred_train))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(2, 56))+scale_y_continuous(limits = c(10, 50))+theme_bw()+
  xlab("Real CLS")+ylab("Predicted CLS")+geom_smooth(method="lm",col="blue")+ggtitle("Training set: CLS predicted by RespLip")
grid.arrange(d,e,f,ncol=3)

CLS_Pred<-data.table(Obs=RespProt_CLS.rf.fine[[3]][,1],RespProt_Pred=RespProt_CLS.rf.fine[[3]][,2],
                           RespMet_Pred=RespMet_CLS.rf.fine[[3]][,2],RespLip_Pred=RespLip_CLS.rf.fine[[3]][,2])

pairs(CLS_Pred[,1:4],pch=20)
summary(lm(RespMet_Pred~Obs+RespProt_Pred,data=CLS_Pred))
summary(lm(RespLip_Pred~Obs+RespProt_Pred,data=CLS_Pred))
summary(lm(RespLip_Pred~Obs+RespMet_Pred,data=CLS_Pred))

#######

tuneGrid <- data.frame(
  .mtry = c(20,30,40,45,50,60,70,80,90,100,200),
  .splitrule = "extratrees",
  .min.node.size = 5
)
FermProt_RLS.rf.fine<-RF.fine(FermProt151,Pheno151$MeanRLS,tuneGrid)


## There are not enough complete cases in the test set to perform KNN imputation
## So we do knn imputation with all 151 genotypes first, then do modeling
## This is not ideal, as there may be data leakage between test and training set
## I do what I have to do 
FermMet151_Imp<-knnImputation(FermMet151)

RF.fine.Missing<-function(Omics,Pheno,tuneGrid) {
  
  ## Split data into training and test sets
  set.seed(1000)
  index<-sample(1:nrow(Omics),round(nrow(Omics)*0.8))
  Omics_train<-Omics[index,]
  Omics_test<-Omics[-index,]
  Pheno_train<-Pheno[index]
  Pheno_test<-Pheno[-index]
  
  ## Model fitting 
  set.seed(4321)
  model<-train(
    Omics_train, 
    Pheno_train,
    tuneGrid = tuneGrid,
    method = "ranger", 
    num.trees=2000,
    trControl = trainControl(
      method = "cv", 
      number = 5,  
      verboseIter = TRUE
    ),
    importance="impurity"
  )
  
  Pred_test<-predict(model,newdata=Omics_test)
  error<-Pheno_test-Pred_test
  RMSE<-sqrt(mean(error^2)) 
  
  Pred_train<-predict(model,newdata=Omics_train)
  
  return(list(Model=model,test_RMSE=RMSE,TrueVSPred_Test=cbind(Pheno_test,Pred_test),
              TrueVSPred_Train=cbind(Pheno_train,Pred_train)))
}

tuneGrid <- data.frame(
  .mtry = c(10,15,20,30,40,50,60,70,80,90,100),
  .splitrule = "extratrees",
  .min.node.size = 5
)

FermMet_RLS.rf.fine<-RF.fine.Missing(FermMet151_Imp,Pheno151$MeanRLS,tuneGrid)


tuneGrid <- data.frame(
  .mtry = c(4,5,6,7,8,9,10,15,20),
  .splitrule = "extratrees",
  .min.node.size = 5
)
FermLip_RLS.rf.fine<-RF.fine(FermLip151,Pheno151$MeanRLS,tuneGrid)

par(mfrow=c(1,3))
# plot(FermProt_RLS.rf.fine[[1]])
plot(FermProt_RLS.rf.fine[[3]][,1],FermProt_RLS.rf.fine[[3]][,2],pch=20)
abline(0,1,col="blue")
abline(lm(FermProt_RLS.rf.fine[[3]][,2]~FermProt_RLS.rf.fine[[3]][,1]),col="red")
cor.test(FermProt_RLS.rf.fine[[3]][,1],FermProt_RLS.rf.fine[[3]][,2])

# plot(FermMet_RLS.rf.fine[[1]])
plot(FermMet_RLS.rf.fine[[3]][,1],FermMet_RLS.rf.fine[[3]][,2],pch=20)
abline(0,1,col="blue")
abline(lm(FermMet_RLS.rf.fine[[3]][,2]~FermMet_RLS.rf.fine[[3]][,1]),col="red")
cor.test(FermMet_RLS.rf.fine[[3]][,1],FermMet_RLS.rf.fine[[3]][,2])

# plot(FermLip_RLS.rf.fine[[1]])
plot(FermLip_RLS.rf.fine[[3]][,1],FermLip_RLS.rf.fine[[3]][,2],pch=20)
abline(0,1,col="blue")
abline(lm(FermLip_RLS.rf.fine[[3]][,2]~FermLip_RLS.rf.fine[[3]][,1]),col="red")
cor.test(FermLip_RLS.rf.fine[[3]][,1],FermLip_RLS.rf.fine[[3]][,2])

invisible(dev.off())

h<-ggplot(data.frame(FermProt_RLS.rf.fine[[3]]),aes(x=Pheno_test,y=Pred_test))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(7, 36))+scale_y_continuous(limits = c(12, 32))+theme_bw()+
  xlab("Real RLS")+ylab("Predicted RLS")+geom_smooth(method="lm",col="blue")+ggtitle("Test set: RLS predicted by FermProt")
i<-ggplot(data.frame(FermMet_RLS.rf.fine[[3]]),aes(x=Pheno_test,y=Pred_test))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(7, 36))+scale_y_continuous(limits = c(12, 32))+theme_bw()+
  xlab("Real RLS")+ylab("Predicted RLS")+geom_smooth(method="lm",col="blue")+ggtitle("Test set: RLS predicted by FermMet")
j<-ggplot(data.frame(FermLip_RLS.rf.fine[[3]]),aes(x=Pheno_test,y=Pred_test))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(7, 36))+scale_y_continuous(limits = c(12, 32))+theme_bw()+
  xlab("Real RLS")+ylab("Predicted RLS")+geom_smooth(method="lm",col="blue")+ggtitle("Test set: RLS predicted by FermLip")
grid.arrange(h,i,j,ncol=3)

k<-ggplot(data.frame(FermProt_RLS.rf.fine[[4]]),aes(x=Pheno_train,y=Pred_train))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(7, 36))+scale_y_continuous(limits = c(12, 32))+theme_bw()+
  xlab("Real RLS")+ylab("Predicted RLS")+geom_smooth(method="lm",col="blue")+ggtitle("Training set: RLS predicted by FermProt")
l<-ggplot(data.frame(FermMet_RLS.rf.fine[[4]]),aes(x=Pheno_train,y=Pred_train))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(7, 36))+scale_y_continuous(limits = c(12, 32))+theme_bw()+
  xlab("Real RLS")+ylab("Predicted RLS")+geom_smooth(method="lm",col="blue")+ggtitle("Training set: RLS predicted by FermMet")
m<-ggplot(data.frame(FermLip_RLS.rf.fine[[4]]),aes(x=Pheno_train,y=Pred_train))+geom_point()+geom_abline(slope=1,intercept=0,col="red")+
  scale_x_continuous(limits = c(7, 36))+scale_y_continuous(limits = c(12, 32))+theme_bw()+
  xlab("Real RLS")+ylab("Predicted RLS")+geom_smooth(method="lm",col="blue")+ggtitle("Training set: RLS predicted by FermLip")
grid.arrange(k,l,m,ncol=3)

RLS_Pred<-data.table(Obs=FermProt_RLS.rf.fine[[3]][,1],FermProt_Pred=FermProt_RLS.rf.fine[[3]][,2],
                     FermMet_Pred=FermMet_RLS.rf.fine[[3]][,2],FermLip_Pred=FermLip_RLS.rf.fine[[3]][,2])

pairs(RLS_Pred[,1:4],pch=20)
summary(lm(FermMet_Pred~Obs+FermProt_Pred,data=RLS_Pred))
summary(lm(FermLip_Pred~Obs+FermProt_Pred,data=RLS_Pred))
summary(lm(FermLip_Pred~Obs+FermMet_Pred,data=RLS_Pred))


## Rank features by their importance in the model 

FeatureImportance<-function(model) {
  Importance<-varImp(model)
  FeatureImportance<-Importance$importance
  FeatureImportance<-FeatureImportance[order(-FeatureImportance),,drop=FALSE]
  return(FeatureImportance)
}


RespProt_CLS.rf.fine_FeatureImportance<-FeatureImportance(RespProt_CLS.rf.fine[[1]])
RespMet_CLS.rf.fine_FeatureImportance<-FeatureImportance(RespMet_CLS.rf.fine[[1]])
RespLip_CLS.rf.fine_FeatureImportance<-FeatureImportance(RespLip_CLS.rf.fine[[1]])

FermProt_RLS.rf.fine_FeatureImportance<-FeatureImportance(FermProt_RLS.rf.fine[[1]])
FermMet_RLS.rf.fine_FeatureImportance<-FeatureImportance(FermMet_RLS.rf.fine[[1]])
FermLip_RLS.rf.fine_FeatureImportance<-FeatureImportance(FermLip_RLS.rf.fine[[1]])

