library(data.table)
library(dplyr)
library(ggplot2)
library(ggExtra)

Pheno<-fread("Phenotype.txt")

## Rename the columns so that they are easier to refer to 
colnames(Pheno)<-c("SGD_ID","Gene_ID","FermDensityFC","RespDensityFC","FermP","RespP","YPG_Phenotype","CLS","CLS%Change","MeanRLS","MaxRLS")
Pheno<-Pheno[-169,c(1,2,3,4,8,10,11)]

## Check distributions of density fold change and lifespan
par(mfrow=c(1,2))

plot(density(Pheno$RespDensityFC,na.rm=TRUE),main="Respiration Density Fold Change")  ##Bimodal
qqnorm(Pheno$RespDensityFC,pch=20)
qqline(Pheno$RespDensityFC,col="red")

plot(density(Pheno$FermDensityFC,na.rm=TRUE),main="Fermentation Density Fold Change")
qqnorm(Pheno$FermDensityFC,pch=20)
qqline(Pheno$FermDensityFC,col="red")

plot(density(Pheno$CLS),main="Chronological Lifespan")
qqnorm(Pheno$CLS,pch=20)
qqline(Pheno$CLS,col="red")

plot(density(Pheno$MeanRLS,na.rm=TRUE),main="Mean Respiration Lifespan")
qqnorm(Pheno$MeanRLS,pch=20)
qqline(Pheno$MeanRLS,col="red")

plot(density(Pheno$MaxRLS,na.rm=TRUE),main="Maximal Respiration Lifespan")
qqnorm(Pheno$MaxRLS,pch=20)
qqline(Pheno$MaxRLS,col="red")

## Pairwise correlation between different genotypes
pairs(Pheno[,3:7],pch=20)
cor(Pheno[,3:7],use="pairwise.complete.obs",method="spearman")
# cor(Pheno[,3:7],use="pairwise.complete.obs",method="pearson")

## Check correlation between pairs that might be correlated
### RespDensityFC vs FermDensityFC 
cor.test(Pheno$FermDensityFC,Pheno$RespDensityFC,method="spearman")
p<-ggplot(Pheno,aes(x=FermDensityFC,y=RespDensityFC))+geom_point()+theme_bw()+
  geom_smooth(method='lm')
ggMarginal(p, type ="histogram")

##CLS vs Mean RLS 
cor.test(Pheno$CLS,Pheno$MeanRLS,method="spearman") # not correlated 
q<-ggplot(Pheno,aes(x=CLS,y=MeanRLS))+geom_point()+theme_bw()
ggMarginal(q, type ="histogram")

##Mean RLS vs Max RLS
cor.test(Pheno$MeanRLS,Pheno$MaxRLS,method="spearman")
r<-ggplot(Pheno,aes(x=MeanRLS,y=MaxRLS))+geom_point()+theme_bw()+
  geom_smooth(method='lm')
ggMarginal(r, type ="histogram")

##RespDensity vs. CLS
cor.test(Pheno$RespDensityFC,Pheno$CLS,method="spearman") 
# p-value = 0.0005897 rho=0.2764692
s<-ggplot(Pheno,aes(x=RespDensityFC,y=CLS))+geom_point()+theme_bw()+
  geom_smooth(method='lm')
ggMarginal(s, type ="histogram")

##FermDensity vs. MeanRLS
cor.test(Pheno$FermDensityFC,Pheno$MeanRLS,method="spearman") 
# p-value = 0.001533 rho=0.2556684
t<-ggplot(Pheno,aes(x=FermDensityFC,y=MeanRLS))+geom_point()+theme_bw()+
  geom_smooth(method='lm')
ggMarginal(t, type ="histogram")

##############

## Test of normality

### Shapiro-Wilk test of normality
apply(Pheno[,-c(1,2)],2,shapiro.test)
### Only MeanRLS is normal. 

### QQ plot
par(mfrow=c(2,3))
apply(Pheno[,-c(1,2)],2,qqnorm)

