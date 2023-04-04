library(data.table)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

load('Cluster.RData')
load('KNN_Omics.RData')
Pheno<-fread("Phenotype_Reordered.txt")
load('SpearmanCor.RData')

## RespProt_CLS

RespProt_Sorted<-sort(assignment.RespProt[-1])
RespProtName<-names(RespProt_Sorted)
knnRespProt_Sorted<-knnRespProt[,RespProtName]
Rho_RespProt_CLS<-as.vector(RespProt_CLS$Rho)
names(Rho_RespProt_CLS)<-RespProt_CLS$Molecule
Rho_RespProt_CLS_Sorted<-Rho_RespProt_CLS[RespProtName]

h1<-Heatmap(t(knnRespProt_Sorted),name="Respiration Proteins \nFold Change",
            show_row_names=FALSE,show_column_names=FALSE,
            show_row_dend=FALSE,show_column_dend=FALSE,cluster_rows = FALSE,
            row_split=as.numeric(RespProt_Sorted))
h2<-Heatmap(Rho_RespProt_CLS_Sorted,name="Correlation with CLS \nSpearman's Rho",
            show_row_names=FALSE,show_column_names=FALSE,cluster_columns=FALSE,cluster_rows=FALSE)
# h3<-Heatmap(as.character(RespProt_Sorted),name="Cluster",
#            show_row_names=FALSE,show_column_names=FALSE,cluster_columns=FALSE,cluster_rows=FALSE)
h1+h2


FermGenotype<-Pheno[!is.na(Pheno$MeanRLS),Gene_ID]
knnFermProt_Complete<-knnFermProt[FermGenotype,]
FermProt_Sorted<-sort(assignment.FermProt[-1])
FermProtName<-names(FermProt_Sorted)
knnFermProt_Sorted<-knnFermProt_Complete[,FermProtName]
Rho_FermProt_RLS<-as.vector(FermProt_MeanRLS$Rho)
names(Rho_FermProt_RLS)<-FermProt_MeanRLS$Molecule
Rho_FermProt_RLS_Sorted<-Rho_FermProt_RLS[FermProtName]


## FermProt_RLS

h3<-Heatmap(t(knnFermProt_Sorted),name="Fermentation Proteins \nFold Change",
            show_row_names=FALSE,show_column_names=FALSE,
            show_row_dend=FALSE,show_column_dend=FALSE,cluster_rows = FALSE,
            row_split=as.numeric(FermProt_Sorted))
h4<-Heatmap(Rho_FermProt_RLS_Sorted,name="Correlation with RLS \nSpearman's Rho",
            show_row_names=FALSE,show_column_names=FALSE,cluster_columns=FALSE,cluster_rows=FALSE)
h3+h4


## RespMet_CLS

Rho_RespMet_CLS<-as.vector(RespMet_CLS$Rho)
names(Rho_RespMet_CLS)<-RespMet_CLS$Molecule

h5<-Heatmap(t(knnRespMet),name="Respiration Metabolites \nFold Change",
        show_row_names=FALSE,show_column_names=FALSE,
        show_row_dend=FALSE,show_column_dend=FALSE)
h6<-Heatmap(Rho_RespMet_CLS,name="Correlation with CLS \nSpearman's Rho",
            show_row_names=FALSE,show_column_names=FALSE,cluster_columns=FALSE,cluster_rows=TRUE)
h5+h6


## RespLip_CLS

Rho_RespLip_CLS<-as.vector(RespLip_CLS$Rho)
names(Rho_RespLip_CLS)<-RespLip_CLS$Molecule

h7<-Heatmap(t(knnRespLip),name="Respiration Lipids \nFold Change",
        show_row_names=FALSE,show_column_names=FALSE,
        show_row_dend=FALSE,show_column_dend=FALSE)
h8<-Heatmap(Rho_RespLip_CLS,name="Correlation with CLS \nSpearman's Rho",
            show_row_names=FALSE,show_column_names=FALSE,cluster_columns=FALSE,cluster_rows=TRUE)
h7+h8


## FermMet_RLS

Rho_FermMet_RLS<-as.vector(FermMet_MeanRLS$Rho)
names(Rho_FermMet_RLS)<-FermMet_MeanRLS$Molecule

h9<-Heatmap(t(knnFermMet),name="Fermentation Metabolites \nFold Change",
            show_row_names=FALSE,show_column_names=FALSE,
            show_row_dend=FALSE,show_column_dend=FALSE)
h10<-Heatmap(Rho_FermMet_RLS,name="Correlation with RLS \nSpearman's Rho",
            show_row_names=FALSE,show_column_names=FALSE,cluster_columns=FALSE,cluster_rows=TRUE)
h9+h10


## FermLip_RLS

Rho_FermLip_RLS<-as.vector(FermLip_MeanRLS$Rho)
names(Rho_FermLip_RLS)<-FermLip_MeanRLS$Molecule

h11<-Heatmap(t(knnFermLip),name="Fermentation Lipids \nFold Change",
            show_row_names=FALSE,show_column_names=FALSE,
            show_row_dend=FALSE,show_column_dend=FALSE)
h12<-Heatmap(Rho_FermLip_RLS,name="Correlation with RLS \nSpearman's Rho",
             show_row_names=FALSE,show_column_names=FALSE,cluster_columns=FALSE,cluster_rows=TRUE)
h11+h12

