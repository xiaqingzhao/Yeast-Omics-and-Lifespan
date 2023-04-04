library(data.table)
library(dplyr)
library(ggplot2)

Resp<-fread("RespOmics.txt")
Resp<-select(Resp,-c(ATP3,COQ8,SCH9,TIM8,'YAL044W-A','YBR230W-A')) # Delete lines that do not have lifespan data
colnames(Resp)[1]<-'Molecule'
colnames(Resp)[2]<-'Type'
plot(density(as.matrix(Resp[,c(-1,-2)]),na.rm=T),xlim=c(-2,2))
quantile(Resp[,c(-1,-2)],probs=seq(0,1,0.01),na.rm=T)

Ferm<-fread("FermOmics.txt")
Ferm<-select(Ferm,-c(ATP3,COQ8,SCH9,TIM8,'YAL044W-A','YBR230W-A')) # Delete lines that do not have lifespan data
colnames(Ferm)[1]<-'Molecule'
colnames(Ferm)[2]<-'Type'
plot(density(as.matrix(Ferm[,c(-1,-2)]),na.rm=T),xlim=c(-2,2))
quantile(Ferm[,c(-1,-2)],probs=seq(0,1,0.01),na.rm=T)

RespProt<-Resp[Type=="Protein"][,-2]
RespMet<-Resp[Type=="Metabolite"][,-2]
RespLip<-Resp[Type=="Lipid"][,-2]

FermProt<-Ferm[Type=="Protein"][,-2]
FermMet<-Ferm[Type=="Metabolite"][,-2]
FermLip<-Ferm[Type=="Lipid"][,-2]

####################################

### Distribution of log2(FC) across strains. Drag screen wide to see x axis labels
meltRespProt<-melt(RespProt,id=1)
meltRespMet<-melt(RespMet,id=1)
meltRespLip<-melt(RespLip,id=1)

ggplot(meltRespProt, aes(x=variable, y=value))+geom_boxplot()+theme_bw()+labs(title="Respiration Protein", x="Strain",y="log2(FC)")+ 
  theme(plot.subtitle = element_text(vjust=1), plot.caption = element_text(vjust=1), axis.text = element_text(angle=90))
ggplot(meltRespMet, aes(x=variable, y=value))+geom_boxplot()+theme_bw()+labs(title="Respiration Metabolite", x="Strain",y="log2(FC)")+ 
  theme(plot.subtitle = element_text(vjust=1), plot.caption = element_text(vjust=1), axis.text = element_text(angle=90))
ggplot(meltRespLip, aes(x=variable, y=value))+geom_boxplot()+theme_bw()+labs(title="Respiration Lipid", x="Strain",y="log2(FC)")+ 
  theme(plot.subtitle = element_text(vjust=1), plot.caption = element_text(vjust=1), axis.text = element_text(angle=90))


meltFermProt<-melt(FermProt,id=1)
meltFermMet<-melt(FermMet,id=1)
meltFermLip<-melt(FermLip,id=1)

ggplot(meltFermProt, aes(x=variable, y=value))+geom_boxplot()+theme_bw()+labs(title="Fermentation Protein", x="Strain",y="log2(FC)")+ 
  theme(plot.subtitle = element_text(vjust=1), plot.caption = element_text(vjust=1), axis.text = element_text(angle=90))
ggplot(meltFermMet, aes(x=variable, y=value))+geom_boxplot()+theme_bw()+labs(title="Fermentation", x="Strain",y="log2(FC)")+ 
  theme(plot.subtitle = element_text(vjust=1), plot.caption = element_text(vjust=1), axis.text = element_text(angle=90))
ggplot(meltFermLip, aes(x=variable, y=value))+geom_boxplot()+theme_bw()+labs(title="Fermentation Lipid", x="Strain",y="log2(FC)")+ 
  theme(plot.subtitle = element_text(vjust=1), plot.caption = element_text(vjust=1), axis.text = element_text(angle=90))

##### Medium is around 0 across all types of molecules. Relatively stable across samples.
##### Lots of outliers

#####################################

# ## Strain names vs protein names
# ## Check to see if e.g. strain ACO1 has very low expression level of protein ACO1
# 
# Knockout<-colnames(RespProt[,-1])  ## Strains
# RespProtein<-unlist(lapply(strsplit(RespProt$Molecule,split=" "),'[',1)) ##Proteins with measurements
# FermProtein<-unlist(lapply(strsplit(FermProt$Molecule,split=" "),'[',1)) ##Proteins with measurements
# 
# a<-Knockout[Knockout %in% RespProtein]
# b<-Knockout[Knockout %in% FermProtein]
# 
# RespProtDF<-as.data.frame(RespProt)
# FermProtDF<-as.data.frame(FermProt)
# RespMeasure<-rep(NA,length(a))
# FermMeasure<-rep(NA,length(b))
# for (i in 1:length(a)) {
#   Row<-c(1:length(RespProtein))[(RespProtein==a[i])]
#   Col<-c(1:length(colnames(RespProt)))[(colnames(RespProt)==a[i])]
#   RespMeasure[i]<-RespProtDF[Row,Col]
# }
# 
# for (i in 1:length(b)) {
#   Row<-c(1:length(FermProtein))[(FermProtein==b[i])]
#   Col<-c(1:length(colnames(FermProt)))[(colnames(FermProt)==b[i])]
#   FermMeasure[i]<-FermProtDF[Row,Col]
# }
# 
# 
# data.frame(KnockoutGene=a,RespMeasure)
# data.frame(KnockoutGene=b,FermMeasure)
# 
# DF<-merge(data.frame(KnockoutGene=a,RespMeasure),data.frame(KnockoutGene=b,FermMeasure),all=FALSE)
# 
# sum(is.na(RespProtDF[c(1:length(RespProtein))[(RespProtein=="COQ4")],]))
# sum(is.na(FermProtDF[c(1:length(FermProtein))[(FermProtein=="COQ4")],]))
# 
# # write.table(RespRelevant,"temp.txt",sep="\t",quote=F,row.names=F)
# 
# #############################################

##Transpose
RespProt<-dcast(melt(RespProt, id.vars="Molecule"), variable~Molecule)
RespMet<-dcast(melt(RespMet, id.vars="Molecule"), variable~Molecule)
RespLip<-dcast(melt(RespLip, id.vars="Molecule"), variable~Molecule)

FermProt<-dcast(melt(FermProt, id.vars="Molecule"), variable~Molecule)
FermMet<-dcast(melt(FermMet, id.vars="Molecule"), variable~Molecule)
FermLip<-dcast(melt(FermLip, id.vars="Molecule"), variable~Molecule)

#######################################

# ## Strain names and ordering in omics dataset
# 
# Strain<-as.character(RespProt$variable)
# as.character(FermProt$variable)
# write.table(Strain, "Strain.txt", sep="\t",quote=F,row.names=F,col.names=F)

#######################################

## Correlation within omics
RespProtCor<-cor(RespProt[,-1],use="pairwise.complete.obs",method="pearson")
sum(abs(RespProtCor[upper.tri(RespProtCor)])>0.9,na.rm=T)
sum(abs(RespProtCor[upper.tri(RespProtCor)])>0.95,na.rm=T)
sum(abs(RespProtCor[upper.tri(RespProtCor)])>0.99,na.rm=T)

RespMetCor<-cor(RespMet[,-1],use="pairwise.complete.obs",method="pearson")
sum(abs(RespMetCor[upper.tri(RespMetCor)])>0.9,na.rm=T)

RespLipCor<-cor(RespLip[,-1],use="pairwise.complete.obs",method="pearson")
sum(abs(RespLipCor[upper.tri(RespLipCor)])>0.9,na.rm=T)


FermProtCor<-cor(FermProt[,-1],use="pairwise.complete.obs",method="pearson")
sum(abs(FermProtCor[upper.tri(FermProtCor)])>0.9,na.rm=T)
sum(abs(FermProtCor[upper.tri(FermProtCor)])>0.95,na.rm=T)
sum(abs(FermProtCor[upper.tri(FermProtCor)])>0.99,na.rm=T)

FermMetCor<-cor(FermMet[,-1],use="pairwise.complete.obs",method="pearson")
sum(abs(FermMetCor[upper.tri(FermMetCor)])>0.9,na.rm=T)

FermLipCor<-cor(FermLip[,-1],use="pairwise.complete.obs",method="pearson")
sum(abs(FermLipCor[upper.tri(FermLipCor)])>0.9,na.rm=T)


### Output highly correlated pairs
RespProtCor[upper.tri(RespProtCor, TRUE)] <- NA
i<-which(abs(RespProtCor) >= 0.99, arr.ind = TRUE)
RespProtHighCor<-data.frame(matrix(colnames(RespProtCor)[as.vector(i)], ncol = 2), Cor=RespProtCor[i])
write.table(RespProtHighCor,"HighlyCorrelatedProteinPairs_Resp.txt",sep="\t",quote=F,row.names=F)
# Note that a lot of 1 and -1 are generaeted because there are only two strains that have measurements of both molecules,
# and correlation between two dots is always 1 or -1

RespMetCor[upper.tri(RespMetCor, TRUE)] <- NA
i<-which(abs(RespMetCor) >= 0.9, arr.ind = TRUE)
RespMetHighCor<-data.frame(matrix(colnames(RespMetCor)[as.vector(i)], ncol = 2), Cor=RespMetCor[i])
# write.table(RespMetHighCor,"HighlyCorrelatedMetabolitePairs_Resp.txt",sep="\t",quote=F,row.names=F)

RespLipCor[upper.tri(RespLipCor, TRUE)] <- NA
i<-which(abs(RespLipCor) >= 0.9, arr.ind = TRUE)
RespLipHighCor<-data.frame(matrix(colnames(RespLipCor)[as.vector(i)], ncol = 2), Cor=RespLipCor[i])
# write.table(RespLipHighCor,"HighlyCorrelatedLipidPairs_Resp.txt",sep="\t",quote=F,row.names=F)

FermProtCor[upper.tri(FermProtCor, TRUE)] <- NA
i<-which(abs(FermProtCor) >= 0.99, arr.ind = TRUE)
FermProtHighCor<-data.frame(matrix(colnames(FermProtCor)[as.vector(i)], ncol = 2), Cor=FermProtCor[i])
write.table(FermProtHighCor,"HighlyCorrelatedProteinPairs_Ferm.txt",sep="\t",quote=F,row.names=F)

FermMetCor[upper.tri(FermMetCor, TRUE)] <- NA
i<-which(abs(FermMetCor) >= 0.9, arr.ind = TRUE)
FermMetHighCor<-data.frame(matrix(colnames(FermMetCor)[as.vector(i)], ncol = 2), Cor=FermMetCor[i])
# write.table(FermMetHighCor,"HighlyCorrelatedMetabolitePairs_Ferm.txt",sep="\t",quote=F,row.names=F)

FermLipCor[upper.tri(FermLipCor, TRUE)] <- NA
i<-which(abs(FermLipCor) >= 0.9, arr.ind = TRUE)
FermLipHighCor<-data.frame(matrix(colnames(FermLipCor)[as.vector(i)], ncol = 2), Cor=FermLipCor[i])
# write.table(FermLipHighCor,"HighlyCorrelatedLipidPairs_Ferm.txt",sep="\t",quote=F,row.names=F)

##############

### Check for missing values 

table(sapply(RespProt[,-1],function(x) sum(is.na(x))))
table(sapply(RespMet[,-1],function(x) sum(is.na(x))))
table(sapply(RespLip[,-1],function(x) sum(is.na(x))))

table(sapply(FermProt[,-1],function(x) sum(is.na(x))))
table(sapply(FermMet[,-1],function(x) sum(is.na(x))))
table(sapply(FermLip[,-1],function(x) sum(is.na(x))))

### Are percentage of missingness correlate with expression level? 

par(mfrow=c(1,3))
plot(y=sapply(RespProt[,-1],function(x) sum(is.na(x))),x=sapply(RespProt[,-1],mean,na.rm=TRUE),pch=20,
     xlab="mean fold change",ylab="number of missingness",main="Respiration Proteins")
abline(h=17,col="red")
plot(y=sapply(RespMet[,-1],function(x) sum(is.na(x))),x=sapply(RespMet[,-1],mean,na.rm=TRUE),pch=20,
     xlab="mean fold change",ylab="number of missingness",main="Respiration Metabolites")
abline(h=17,col="red")
plot(y=sapply(RespLip[,-1],function(x) sum(is.na(x))),x=sapply(RespLip[,-1],mean,na.rm=TRUE),pch=20,
     xlab="mean fold change",ylab="number of missingness",main="Respiration Lipids")
abline(h=17,col="red")

plot(y=sapply(FermProt[,-1],function(x) sum(is.na(x))),x=sapply(FermProt[,-1],mean,na.rm=TRUE),pch=20,
     xlab="mean fold change",ylab="number of missingness",main="Fermentation Proteins")
abline(h=17,col="red")
plot(y=sapply(FermMet[,-1],function(x) sum(is.na(x))),x=sapply(FermMet[,-1],mean,na.rm=TRUE),pch=20,
     xlab="mean fold change",ylab="number of missingness",main="Fermentation Metabolites")
abline(h=17,col="red")
plot(y=sapply(FermLip[,-1],function(x) sum(is.na(x))),x=sapply(FermLip[,-1],mean,na.rm=TRUE),pch=20,
     xlab="mean fold change",ylab="number of missingness",main="Fermentation Lipids")
abline(h=17,col="red")

dev.off()

###########

Strain<-as.character(RespProt[[1]])

RespProt<-as.matrix(RespProt[,-1])
rownames(RespProt)<-Strain
RespMet<-as.matrix(RespMet[,-1])
rownames(RespMet)<-Strain
RespLip<-as.matrix(RespLip[,-1])
rownames(RespLip)<-Strain

FermProt<-as.matrix(FermProt[,-1])
rownames(FermProt)<-Strain
FermMet<-as.matrix(FermMet[,-1])
rownames(FermMet)<-Strain
FermLip<-as.matrix(FermLip[,-1])
rownames(FermLip)<-Strain





###########

### Factor analysis 

RespProt.pca<-prcomp(RespProt)
factanal(RespProt,factors=3,rotation="varimax")


