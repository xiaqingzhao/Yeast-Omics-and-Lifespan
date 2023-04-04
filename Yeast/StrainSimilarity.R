# Unsupervised methods to identify strains that are highly simialr in their omics
# If there are strains that highly similar (we would expect ACO1 and ACO2 to have similar omics),
# then dividing strains randomly into training and test sets may be problematic

library(data.table)

Strain<-scan("Strain.txt",character())

load("KNN_Omics.RData")

#######

# ## K-Means clustering
# 
# Resp_kmeans<-kmeans(Resp,5,nstart=20)
# Resp_kmeans
# Resp_kmeans$cluster
# plot(Resp,col=Resp_kmeans$cluster+1,main="Respiration: K-Means Clustering with K=5")

########

## Hierarchical clustering

### Euclidean distance 

RespProtHC.complete<-hclust(dist(knnRespProt),method="complete")
plot(RespProtHC.complete)
RespProtHC.average<-hclust(dist(knnRespProt),method="average")
plot(RespProtHC.average)

RespMetHC.complete<-hclust(dist(knnRespMet),method="complete")
plot(RespMetHC.complete)
RespMetHC.average<-hclust(dist(knnRespMet),method="average")
plot(RespMetHC.average)

RespLipHC.complete<-hclust(dist(knnRespLip),method="complete")
plot(RespLipHC.complete)
RespLipHC.average<-hclust(dist(knnRespLip),method="average")
plot(RespLipHC.average)


FermProtHC.complete<-hclust(dist(knnFermProt),method="complete")
plot(FermProtHC.complete)
FermProtHC.average<-hclust(dist(knnFermProt),method="average")
plot(FermProtHC.average)

FermMetHC.complete<-hclust(dist(knnFermMet),method="complete")
plot(FermMetHC.complete)
FermMetHC.average<-hclust(dist(knnFermMet),method="average")
plot(FermMetHC.average)

FermLipHC.complete<-hclust(dist(knnFermLip),method="complete")
plot(FermLipHC.complete)
FermLipHC.average<-hclust(dist(knnFermLip),method="average")
plot(FermLipHC.average)

### Not surprisingly, knockout of the same type of proteins tend to cluster together 

