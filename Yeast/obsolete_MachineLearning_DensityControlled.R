library(data.table)
library(caret)
library(splines)

Pheno<-fread("Phenotype_Reordered.txt")
Strain<-scan("Strain.txt",character())

Resp<-fread("KNNResp.csv")
Resp<-as.matrix(Resp)
rownames(Resp)<-Strain
RespType<-scan("DeMissingRespType.txt",character())

Ferm<-fread("KNNFerm.csv")
Ferm<-as.matrix(Ferm)
rownames(Ferm)<-Strain

#############

## Some strains do not have density. Take out their omics measures. 

RespDen<-Pheno$RespDensityFC
RespDen<-RespDen[!is.na(RespDen)]
Resp<-Resp[!is.na(RespDen),]


# RespDenCon<-sapply(1:ncol(Resp), function(x) {lm(Resp[,x]~RespDen)$residual})
### Linear regression does not seem right since molecule level~Respiration density is not always linear

Denslims<-range(RespDen)
Dens.grid<-seq(from=Denslims[1],to=Denslims[2],by=0.01)

# Omic<-Resp[,2160]
# 
# plot(RespDen,Omic,xlim=Denslims,cex=.5,col="darkgrey")
# title("Smoothing Spline")
# fit<-smooth.spline(RespDen,Omic,cv=TRUE)
# fit
# lines(fit,col="red",lwd=2)
# # residuals(fit)
# 
# fit2<-loess(Omic~RespDen)
# fit2
# plot(RespDen,Omic,xlim=Denslims,cex=.5,col="darkgrey")
# title("Local Regression")
# Pred<-predict(fit2,newdata=data.frame(RespDen=Dens.grid),se=T)
# lines(Dens.grid,Pred$fit,col="blue",lwd=2)

### Conclusion: Use local regression and take the residual to control for density

RespDenCon<-sapply(1:ncol(Resp), function(x) {loess(Resp[,x]~RespDen)$residual})
colnames(RespDenCon)<-colnames(Resp)

###################

### Split data into training and test sets
set.seed(1000)
index<-sample(1:151,round(151*0.8))
Resp_train<-RespDenCon[index,]
Resp_test<-RespDenCon[-index,]
Pheno<-Pheno[!is.na(Pheno$RespDensityFC),]
Pheno_train<-Pheno$CLS[index]
Pheno_test<-Pheno$CLS[-index]

#################

### Partial least squares

tuneGrid <- data.frame(
  .ncomp = c(1,2,3,5,10)
)

model<-train(
  Resp_train, 
  Pheno_train,
  method = "pls",  ## can be any methods available
  tuneGrid = tuneGrid,
  trControl = trainControl(
    method = "cv", 
    number = 5,  ## 5-fold cross validation
    verboseIter = TRUE
  )
)

# ncomp  RMSE       Rsquared   MAE     
# 1      9.661625  0.2369583  7.467240
# RMSE was used to select the optimal model using the smallest value.
# The final value used for the model was ncomp = 1.
# Cross validation found that the optimal number of retained dimensions is 1. 

error<-Pheno_test-predict(model,newdata=Resp_test)
sqrt(mean(error^2)) 

###

pls.fit<-plsr(Pheno_train~Resp_train, ncomp=2, scale=TRUE,validation="CV")
summary(pls.fit)

###############

### Random forest

model<-train(
  Resp_train, 
  Pheno_train,
  method = "ranger", 
  trControl = trainControl(
    method = "cv", 
    number = 5,  ## 5-fold cross validation
    verboseIter = TRUE
  )
)

error<-Pheno_test-predict(model,newdata=Resp_test)
sqrt(mean(error^2)) 

#### Custom tuning grid

tuneGrid <- data.frame(
  .mtry = rep(c(50,60,70,80,90,100,200,400),each=2),
  .splitrule = c("variance","extratrees"),
  .min.node.size = 5
)

model <- train(
  Resp_train, 
  Pheno_train,
  tuneGrid = tuneGrid,
  method = "ranger",
  trControl = trainControl(
    method = "cv", 
    number = 5, 
    verboseIter = TRUE
  )
)

plot(model)

error<-Pheno_test-predict(model,newdata=Resp_test)
sqrt(mean(error^2)) 

##################

### Ridge regression and lasso

tuneGrid = expand.grid(
  alpha=0,
  lambda=seq(1,1000,length=101)
)

model <- train(
  Resp_train, 
  Pheno_train,
  tuneGrid = tuneGrid,
  method = "glmnet",
  trControl = trainControl(
    method = "cv", 
    number = 5, 
    verboseIter = TRUE
  )
)

plot(model)

error<-Pheno_test-predict(model,newdata=Resp_test)
sqrt(mean(error^2))










#FermDen<-Pheno$FermDensityFC
#FermDen[is.na(FermDen)]<-mean(Pheno$FermDensityFC,na.rm=TRUE)



