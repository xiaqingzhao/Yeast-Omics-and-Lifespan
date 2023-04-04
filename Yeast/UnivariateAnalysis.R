library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggplus)

Pheno<-fread("Phenotype_Reordered.txt")

load("DeMissing_Omics.RData")

####################

## Shapiro-Wilk's test of normality on omics data for individual molecules
## Calculate the proportion of molecules that are NOT normally distributed

NormRespProt<-apply(RespProt,2,function(x) shapiro.test(x)$p.value)
sum(NormRespProt<0.05)/length(NormRespProt) # 0.8533679
NormRespMet<-apply(RespMet,2,function(x) shapiro.test(x)$p.value)
sum(NormRespMet<0.05)/length(NormRespMet) # 0.6834862
NormRespLip<-apply(RespLip,2,function(x) shapiro.test(x)$p.value)
sum(NormRespLip<0.05)/length(NormRespLip) # 0.7735849

NormFermProt<-apply(FermProt,2,function(x) shapiro.test(x)$p.value)
sum(NormFermProt<0.05)/length(NormFermProt) # 0.9039488
NormFermMet<-apply(FermMet,2,function(x) shapiro.test(x)$p.value)
sum(NormFermMet<0.05)/length(NormFermMet) # 0.7814208
NormFermLip<-apply(FermLip,2,function(x) shapiro.test(x)$p.value)
sum(NormFermLip<0.05)/length(NormFermLip) # 0.6730769

## Most molecules do not follow normal distribution across the 168 strains

####################

## Function that does Spearman test of correlation, and multiple testing correction

SpearmanCor<-function(Phenotype,Predictor) {
  tmp<-apply(Predictor,2,cor.test,Phenotype,method="spearman")
  Rho<-sapply(1:length(tmp),function(x){tmp[[x]]$estimate})
  P<-sapply(1:length(tmp),function(x){tmp[[x]]$p.value})
  FDR<-p.adjust(P,method="BH")
  Bonferroni<-p.adjust(P,method="bonferroni")
  return(data.table(Molecule=colnames(Predictor),Rho,P,FDR,Bonferroni))
}

#####################

RespProt_CLS<-SpearmanCor(Pheno$CLS,RespProt)
RespProt_CLS[FDR<0.05]

RespMet_CLS<-SpearmanCor(Pheno$CLS,RespMet)
RespMet_CLS[FDR<0.05]

RespLip_CLS<-SpearmanCor(Pheno$CLS,RespLip)
RespLip_CLS[FDR<0.05]

### 

RespProt_RespDensityFC<-SpearmanCor(Pheno$RespDensityFC,RespProt)
RespProt_RespDensityFC[FDR<0.05]

RespMet_RespDensityFC<-SpearmanCor(Pheno$RespDensityFC,RespMet)
RespMet_RespDensityFC[FDR<0.05]

RespLip_RespDensityFC<-SpearmanCor(Pheno$RespDensityFC,RespLip)
RespLip_RespDensityFC[FDR<0.05]

### 

FermProt_MeanRLS<-SpearmanCor(Pheno$MeanRLS,FermProt)
FermProt_MeanRLS[FDR<0.05]

FermMet_MeanRLS<-SpearmanCor(Pheno$MeanRLS,FermMet)
FermMet_MeanRLS[FDR<0.05]

FermLip_MeanRLS<-SpearmanCor(Pheno$MeanRLS,FermLip)
FermLip_MeanRLS[FDR<0.05]

###

FermProt_FermDensityFC<-SpearmanCor(Pheno$FermDensityFC,FermProt)
FermProt_FermDensityFC[FDR<0.05]

FermMet_FermDensityFC<-SpearmanCor(Pheno$FermDensityFC,FermMet)
FermMet_FermDensityFC[FDR<0.05]

FermLip_FermDensityFC<-SpearmanCor(Pheno$FermDensityFC,FermLip)
FermLip_FermDensityFC[FDR<0.05]
# write.table(FermLip_FermDensityFC[FDR<0.05]$Molecule,"temp",sep="\n",quote=F,row.names=F)

save(RespProt_CLS,RespMet_CLS,RespLip_CLS,
     RespProt_RespDensityFC,RespMet_RespDensityFC,RespLip_RespDensityFC,
     FermProt_MeanRLS,FermMet_MeanRLS,FermLip_MeanRLS,
     FermProt_FermDensityFC,FermMet_FermDensityFC,FermLip_FermDensityFC,
     file="SpearmanCor.RData")

#####################

### Is the relationship between omics and LS linear? 

MakeLong<-function(Matrix) {
  tmp<-data.table(Matrix,keep.rownames=TRUE)
  tmp<-data.table(Pheno[,3:7],tmp)
  tmp<-melt(tmp,id.vars=c("FermDensityFC","RespDensityFC","CLS","MeanRLS","MaxRLS","rn"),
            variable.name="Molecule",value.name="FC")
}

RespProtLong<-MakeLong(RespProt)
RespMetLong<-MakeLong(RespMet)
RespLipLong<-MakeLong(RespLip)
FermProtLong<-MakeLong(FermProt)
FermMetLong<-MakeLong(FermMet)
FermLipLong<-MakeLong(FermLip)

Plot<-ggplot(RespProtLong,aes(x=CLS,y=FC))+
  geom_point()+theme_bw()+xlab("CLS")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("RespProt_CLS.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(RespMetLong,aes(x=CLS,y=FC))+
  geom_point()+theme_bw()+xlab("CLS")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("RespMet_CLS.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(RespLipLong,aes(x=CLS,y=FC))+
  geom_point()+theme_bw()+xlab("CLS")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("RespLip_CLS.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(FermProtLong,aes(x=MeanRLS,y=FC))+
  geom_point()+theme_bw()+xlab("RLS")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("FermProt_RLS.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(FermMetLong,aes(x=MeanRLS,y=FC))+
  geom_point()+theme_bw()+xlab("RLS")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("FermMet_RLS.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(FermLipLong,aes(x=MeanRLS,y=FC))+
  geom_point()+theme_bw()+xlab("RLS")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("FermLip_RLS.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()




Plot<-ggplot(RespProtLong,aes(x=RespDensityFC,y=FC))+
  geom_point()+theme_bw()+xlab("DensityFC")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("RespProt_RespDensity.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(RespMetLong,aes(x=RespDensityFC,y=FC))+
  geom_point()+theme_bw()+xlab("DensityFC")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("RespMet_RespDensity.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(RespLipLong,aes(x=RespDensityFC,y=FC))+
  geom_point()+theme_bw()+xlab("DensityFC")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("RespLip_RespDensity.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(FermProtLong,aes(x=FermDensityFC,y=FC))+
  geom_point()+theme_bw()+xlab("DensityFC")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("FermProt_FermDensity.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(FermMetLong,aes(x=FermDensityFC,y=FC))+
  geom_point()+theme_bw()+xlab("DensityFC")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("FermMet_FermDensity.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()

Plot<-ggplot(FermLipLong,aes(x=FermDensityFC,y=FC))+
  geom_point()+theme_bw()+xlab("DensityFC")+ylab("FC")+
  geom_smooth(method="lm",colour="red",se=FALSE)+geom_smooth(method="loess",colour="blue",se=FALSE)
pdf("FermLip_FermDensity.pdf")
facet_multiple(Plot, facets ='Molecule', ncol=2, nrow=3, scales ='free_y')
dev.off()



################### 

## Density as a confounder 
### Note that in all models that has density as covariate, we have only 151 genotypes to work with
### Analysis of variance in regression (The R Book page 396)

### Linear model 

LinearModel<-function(Phenotype,Covariate,DataMatrix) {
  Pheno_Beta<-rep(NA,ncol(DataMatrix))
  Pheno_Beta_P<-rep(NA,ncol(DataMatrix))
  Pheno_ANOVA_P<-rep(NA,ncol(DataMatrix))
  Cov_Beta_P<-rep(NA,ncol(DataMatrix))
  Cov_ANOVA_P<-rep(NA,ncol(DataMatrix))
  
  for (i in 1:ncol(DataMatrix)) {
    mod<-lm(DataMatrix[,i]~Phenotype+Covariate)
    Pheno_Beta[i]<-summary(mod)$coefficients[2,1]
    Pheno_Beta_P[i]<-summary(mod)$coefficients[2,4]
    Pheno_ANOVA_P[i]<-anova(mod)$'Pr(>F)'[1]
    Cov_Beta_P[i]<-summary(mod)$coefficients[3,4]
    Cov_ANOVA_P[i]<-anova(mod)$'Pr(>F)'[2]
  }
  
  return(data.table(Molecule=colnames(DataMatrix),Pheno_Beta,Pheno_Beta_P,Pheno_Beta_P.adj=p.adjust(Pheno_Beta_P),
                    Pheno_ANOVA_P,Pheno_ANOVA_P.adj=p.adjust(Pheno_ANOVA_P),
                    Cov_Beta_P,Cov_ANOVA_P))
}

LM_RespProt_CLS<-LinearModel(Pheno$CLS,Pheno$RespDensityFC,RespProt)
LM_RespProt_CLS[Pheno_Beta_P.adj<0.1,]
LM_RespProt_CLS[Pheno_ANOVA_P.adj<0.1,]

LM_RespMet_CLS<-LinearModel(Pheno$CLS,Pheno$RespDensityFC,RespMet)
LM_RespMet_CLS[Pheno_Beta_P.adj<0.1,]
LM_RespMet_CLS[Pheno_ANOVA_P.adj<0.1,]

LM_RespLip_CLS<-LinearModel(Pheno$CLS,Pheno$RespDensityFC,RespLip)
LM_RespLip_CLS[Pheno_Beta_P.adj<0.1,]
LM_RespLip_CLS[Pheno_ANOVA_P.adj<0.1,]


LM_FermProt_RLS<-LinearModel(Pheno$MeanRLS,Pheno$FermDensityFC,FermProt)
LM_FermProt_RLS[Pheno_Beta_P.adj<0.1,]
LM_FermProt_RLS[Pheno_ANOVA_P.adj<0.1,]

LM_FermMet_RLS<-LinearModel(Pheno$MeanRLS,Pheno$FermDensityFC,FermMet)
LM_FermMet_RLS[Pheno_Beta_P.adj<0.1,]
LM_FermMet_RLS[Pheno_ANOVA_P.adj<0.1,]

LM_FermLip_RLS<-LinearModel(Pheno$MeanRLS,Pheno$FermDensityFC,FermLip)
LM_FermLip_RLS[Pheno_Beta_P.adj<0.1,]
LM_FermLip_RLS[Pheno_ANOVA_P.adj<0.1,]

save(LM_RespProt_CLS,LM_RespMet_CLS,LM_RespLip_CLS,
     LM_FermProt_RLS,LM_FermMet_RLS,LM_FermLip_RLS,file="LinearModel.RData")

# mod<-lm(RespProt[,1]~Pheno$CLS+Pheno$RespDensityFC)
# summary(mod)
# anova(mod)

### Generalized additive model 

library(mgcv)

GAM<-function(Phenotype,Covariate,DataMatrix) {
  Pheno_Beta<-rep(NA,ncol(DataMatrix))
  Pheno_Beta_P<-rep(NA,ncol(DataMatrix))
  Cov_Beta_P<-rep(NA,ncol(DataMatrix))
  
  for (i in 1:ncol(DataMatrix)) {
    model<-gam(DataMatrix[,i]~Phenotype+s(Covariate))
    Pheno_Beta[i]<-summary(model)$p.coeff[2]
    Pheno_Beta_P[i]<-summary(model)$p.pv[2]
    Cov_Beta_P[i]<-summary(model)$s.pv
  }
  
  return(data.table(Molecule=colnames(DataMatrix),Pheno_Beta,Pheno_Beta_P,Pheno_Beta_P.adj=p.adjust(Pheno_Beta_P),
                    Cov_Beta_P))
}

GAM_RespProt_CLS<-GAM(Pheno$CLS,Pheno$RespDensityFC,RespProt)
GAM_RespProt_CLS[Pheno_Beta_P<0.05,]

GAM_RespMet_CLS<-GAM(Pheno$CLS,Pheno$RespDensityFC,RespMet)
GAM_RespMet_CLS[Pheno_Beta_P<0.05,]

GAM_RespLip_CLS<-GAM(Pheno$CLS,Pheno$RespDensityFC,RespLip)
GAM_RespLip_CLS[Pheno_Beta_P<0.05,]


GAM_FermProt_RLS<-GAM(Pheno$MeanRLS,Pheno$FermDensityFC,FermProt)
GAM_FermProt_RLS[Pheno_Beta_P<0.05,]

GAM_FermMet_RLS<-GAM(Pheno$MeanRLS,Pheno$FermDensityFC,FermMet)
GAM_FermMet_RLS[Pheno_Beta_P<0.05,]

GAM_FermLip_RLS<-GAM(Pheno$MeanRLS,Pheno$FermDensityFC,FermLip)
GAM_FermLip_RLS[Pheno_Beta_P<0.05,]

save(GAM_RespProt_CLS,GAM_RespMet_CLS,GAM_RespLip_CLS,
     GAM_FermProt_RLS,GAM_FermMet_RLS,GAM_FermLip_RLS,file="GAM.RData")

# model<-gam(RespProt[,1]~s(Pheno$CLS)+s(Pheno$RespDensityFC))
# summary(model)
# plot(model)
# 
# model<-gam(RespProt[,4]~Pheno$CLS+s(Pheno$RespDensityFC))
# summary(model)
# anova(model)
# 
# library(tree)
# model<-tree(RespProt[,3]~Pheno$CLS+Pheno$RespDensityFC)
# plot(model)
# text(model)



#################################################

## Did not end up using the following code... 

Resp_CLS<-SpearmanCor(Pheno$CLS,Resp)
Resp_CLS
Resp_CLS[FDR<0.05]  ### 966 significant molecules with FDR<0.05
mean(abs(Resp_CLS[FDR<0.05]$Rho)) ### 0.2427931
Resp_CLS[Bonferroni<0.1] ### 39 significant molecules with Bonferroni<0.1
mean(abs(Resp_CLS[Bonferroni<0.05]$Rho)) ### 0.3362457
plot(hist(Resp_CLS[FDR<0.05]$Rho))
Random<-sample(Pheno$CLS,168)
SpearmanCor(Random,Resp)

Resp_RespDensity<-SpearmanCor(Pheno$RespDensityFC,Resp)
Resp_RespDensity    ### The Resp omics seems to be better related to RespDensity, as compared to CLS
Resp_RespDensity[FDR<0.05] ### 1797 significant molecules with FDR<0.05
mean(abs(Resp_RespDensity[FDR<0.05]$Rho)) ### 0.5225872
Resp_RespDensity[Bonferroni<0.1] ### 1375 significant molecules with Bonferroni<0.1. This is striking! 
mean(abs(Resp_RespDensity[Bonferroni<0.1]$Rho)) ### 0.6064454
plot(hist(Resp_RespDensity[FDR<0.05]$Rho))  
Random<-sample(Pheno$RespDensityFC,168)
SpearmanCor(Random,Resp)

Resp_MeanRLS<-SpearmanCor(Pheno$MeanRLS,Resp)
Resp_MeanRLS
Resp_MeanRLS[FDR<0.05]  ### 513 significant molecules with FDR<0.05
mean(abs(Resp_MeanRLS[FDR<0.05]$Rho)) ### 0.2525442
Resp_MeanRLS[Bonferroni<0.1] ### 17 significant molecules with Bonferroni<0.1
mean(abs(Resp_MeanRLS[Bonferroni<0.1]$Rho)) ### 0.3448247
plot(hist(Resp_MeanRLS[FDR<0.05]$Rho))
Random<-sample(Pheno$MeanRLS,168)
SpearmanCor(Random,Resp)

# Resp_MaxRLS<-SpearmanCor(Pheno$MaxRLS,Resp)
# Resp_MaxRLS
# Resp_MaxRLS[FDR<0.05]  ### 6 significant molecules with FDR<0.05
# mean(abs(Resp_MaxRLS[FDR<0.05]$Rho)) ### 0.3225118
# Resp_MaxRLS[Bonferroni<0.1] ### 1 significant molecules with Bonferroni<0.1
# mean(abs(Resp_MaxRLS[Bonferroni<0.1]$Rho)) ### 0.3656386
# plot(hist(Resp_MaxRLS[FDR<0.05]$Rho))
# Random<-sample(Pheno$MaxRLS,168)
# SpearmanCor(Random,Resp)

Resp_FermDensity<-SpearmanCor(Pheno$FermDensityFC,Resp)
Resp_FermDensity    ### The Resp omics seems to be better related to RespDensity, as compared to CLS
Resp_FermDensity[FDR<0.05] ### 1661 significant molecules with FDR<0.05
mean(abs(Resp_FermDensity[FDR<0.05]$Rho)) ### 0.4288944
Resp_FermDensity[Bonferroni<0.1] ### 1119 significant molecules with Bonferroni<0.1. This is striking! 
mean(abs(Resp_FermDensity[Bonferroni<0.1]$Rho)) ### 0.5165099
plot(hist(Resp_FermDensity[FDR<0.05]$Rho))  
Random<-sample(Pheno$FermDensityFC,168)
SpearmanCor(Random,Resp)

#####

Ferm_CLS<-SpearmanCor(Pheno$CLS,Ferm)
Ferm_CLS
Ferm_CLS[FDR<0.05]  ### 76 significant molecules with FDR<0.05
mean(abs(Ferm_CLS[FDR<0.05]$Rho)) ### 0.2705297
Ferm_CLS[Bonferroni<0.1] ### 6 significant molecules with Bonferroni<0.1
mean(abs(Ferm_CLS[Bonferroni<0.1]$Rho)) ### 0.3379317
plot(hist(Ferm_CLS[FDR<0.05]$Rho))
Random<-sample(Pheno$CLS,168)
SpearmanCor(Random,Ferm)

Ferm_RespDensity<-SpearmanCor(Pheno$RespDensityFC,Ferm)
Ferm_RespDensity    ### The Resp omics seems to be better related to RespDensity, as compared to CLS
Ferm_RespDensity[FDR<0.05] ### 1322 significant molecules with FDR<0.05
mean(abs(Ferm_RespDensity[FDR<0.05]$Rho)) ### 0.3645886
Ferm_RespDensity[Bonferroni<0.1] ### 721 significant molecules with Bonferroni<0.1. 
mean(abs(Ferm_RespDensity[Bonferroni<0.1]$Rho)) ### 0.4606165
plot(hist(Ferm_RespDensity[FDR<0.05]$Rho))  
Random<-sample(Pheno$RespDensityFC,168)
SpearmanCor(Random,Ferm)

Ferm_MeanRLS<-SpearmanCor(Pheno$MeanRLS,Ferm)
Ferm_MeanRLS
Ferm_MeanRLS[FDR<0.05]  ### 443 significant molecules with FDR<0.05
mean(abs(Ferm_MeanRLS[FDR<0.05]$Rho)) ### 0.2598153
Ferm_MeanRLS[Bonferroni<0.1] ### 37 significant molecules with Bonferroni<0.1
mean(abs(Ferm_MeanRLS[Bonferroni<0.1]$Rho)) ### 0.3479196
plot(hist(Ferm_MeanRLS[FDR<0.05]$Rho))
Random<-sample(Pheno$MeanRLS,168)
SpearmanCor(Random,Ferm)

# Ferm_MaxRLS<-SpearmanCor(Pheno$MaxRLS,Ferm)
# Ferm_MaxRLS
# Ferm_MaxRLS[FDR<0.05]  ### 10 significant molecules with FDR<0.05
# mean(abs(Ferm_MaxRLS[FDR<0.05]$Rho)) ### 0.3179246
# Ferm_MaxRLS[Bonferroni<0.1] ### 2 significant molecules with Bonferroni<0.1
# mean(abs(Ferm_MaxRLS[Bonferroni<0.1]$Rho)) ### 0.3669584
# plot(hist(Ferm_MaxRLS[FDR<0.05]$Rho))
# Random<-sample(Pheno$MaxRLS,168)
# SpearmanCor(Random,Ferm)

Ferm_FermDensity<-SpearmanCor(Pheno$FermDensityFC,Ferm)
Ferm_FermDensity    
Ferm_FermDensity[FDR<0.05] ### 1368 significant molecules with FDR<0.05
mean(abs(Ferm_FermDensity[FDR<0.05]$Rho)) ### 0.3808724
Ferm_FermDensity[Bonferroni<0.1] ### 796 significant molecules with Bonferroni<0.1. 
mean(abs(Ferm_FermDensity[Bonferroni<0.1]$Rho))
plot(hist(Ferm_FermDensity[FDR<0.05]$Rho))  
Random<-sample(Pheno$FermDensityFC,168)
SpearmanCor(Random,Ferm)

