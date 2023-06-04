# Multivariate Analysis Project 2
# Laila El Saadawi: 900201723
# Masa Tantawy: 900201312

# ............. Cluster Analysis .............
library(DescTools)
library(Matrix)
library(GGally)

df = read.csv("fake_bills.csv", header= TRUE, sep=";")
dim(df); View(head(df)); str(df)

# --------- Data Preparation --------------
# *Note: Please check R script Data Preparation for details

# 1. Categorical variable is_genuine
df$is_genuine= as.factor(df$is_genuine)
# 2. Missing values
# Replacing missing values with the mean
df[is.na(df[,'margin_low']), 'margin_low']= mean(df[,'margin_low'], na.rm = TRUE)
is_genuine=df[,1]
df=df[,c(2:7)]; 


# Defining Functions
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
  require(clue)
  idsA <- unique(clusteringA)  # distinct cluster ids in a
  idsB <- unique(clusteringB)  # distinct cluster ids in b
  nA <- length(clusteringA)  # number of instances in a
  nB <- length(clusteringB)  # number of instances in b
  if (length(idsA) != length(idsB) || nA != nB) {
    stop("number of clusters or number of instances do not match")
  }
  
  nC <- length(idsA)
  tupel <- c(1:nA)
  
  # computeing the distance matrix
  assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
  for (i in 1:nC) {
    tupelClusterI <- tupel[clusteringA == i]
    solRowI <- sapply(1:nC, function(i, clusterIDsB, tupelA_I) {
      nA_I <- length(tupelA_I)  # number of elements in cluster I
      tupelB_I <- tupel[clusterIDsB == i]
      nB_I <- length(tupelB_I)
      nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
      return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
    }, clusteringB, tupelClusterI)
    assignmentMatrix[i, ] <- solRowI
  }
  
  # optimization
  result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
  attr(result, "assignmentMatrix") <- assignmentMatrix
  return(result)
}
R2=function(x,clusters,k=2){
  n=nrow(x); tss=var(x); tss=(n-1)*sum(diag(tss));
  wss=0
  for(j in 0:k-1){
    cj=x[clusters==j,]; nj=nrow(cj);
    vj=var(cj); wssj=0
    if(is.matrix(cj)) wssj=(nj-1)*sum(diag(vj));
    wss=wss+wssj
  }
  r2=1-wss/tss; 
  cat("R2 = ",r2,"\n")
  return(r2)
}

# -------- 1. Hierarchical Algorithm (Agglomerative) ------
#------- 1.1 ----------
# Method: intracluster dist: euclidean --- intercluster dist: ward distance
d = dist(df, method = "euclidean")
hc = hclust(d, method="ward.D2")
plot(hc, ) # display dendrogram
mtext('Inter: Euclidean, Intra: Ward')
clusters=cutree(hc, k=2) 
# draw dendogram with red borders around them
rect.hclust(hc, k=2, border="red")

# Measuring the goodness of fit
r2=R2(df,clusters,2)
cm= table(clusters,is_genuine); print(cm)
cat("Error Rate:",100 *(1-sum(diag(cm))/sum(cm)),"%\n")
#library(anticlust)
#n=nrow(df); tss=var(df); tss=(n-1)*sum(diag(tss));
#bss= dispersion_objective(df, clusters)
#bss/tss

#------- 1.2 ----------
# Method: intracluster dist: Manhattan --- intercluster dist: ward distance
d = dist(df, method = "manhattan")
hc = hclust(d, method="ward.D2")
plot(hc, ) # display dendrogram
mtext('Inter: Manhattan, Intra: Ward')
clusters=cutree(hc, k=2) 
# draw dendogram with red borders around them
rect.hclust(hc, k=2, border="red")

# Measuring the goodness of fit
r2=R2(df,clusters,2)
cm= table(clusters,is_genuine)[,c(2,1)]; print(cm)
cat("Error Rate:",100 *(1-sum(diag(cm))/sum(cm)),"%\n")

#------- 1.3 ----------
# Method: intracluster dist: Euclidean --- intercluster dist: Average Linkage
d = dist(df, method = "euclidean")
hc = hclust(d, method="average")
plot(hc, ) # display dendrogram
mtext('Inter: Euclidean, Intra: Average Linkage')
clusters=cutree(hc, k=2) 
# draw dendogram with red borders around them
rect.hclust(hc, k=2, border="red")

# Measuring the goodness of fit
r2=R2(df,clusters,2)
cm= table(clusters,is_genuine)[,c(2,1)]; print(cm)
cat("Error Rate:", 100 *(1-sum(diag(cm))/sum(cm)),"%\n")

#------- 1.4 ----------
# Method: intracluster dist: Manhattan --- intercluster dist: Average Linkage
d = dist(df, method = "manhattan")
hc = hclust(d, method="average")
plot(hc, ) # display dendrogram
mtext('Inter: Manhattan, Intra: Average Linkage')
clusters=cutree(hc, k=2) 
# draw dendogram with red borders around them
rect.hclust(hc, k=2, border="red")

# Measuring the goodness of fit
r2=R2(df,clusters,2)
cm= table(clusters,is_genuine)[,c(2,1)]; print(cm)
cat("Error Rate:",100 *(1-sum(diag(cm))/sum(cm)),"%\n")

#------- 1.5 ----------
# Method: intracluster dist: Euclidean --- intercluster dist: Cluster Centroids
d = dist(df, method = "euclidean")
hc = hclust(d, method="centroid")
plot(hc, ) # display dendrogram
mtext('Inter: Euclidean, Intra: Cluster Centroids')
clusters=cutree(hc, k=2) 
# draw dendogram with red borders around them
rect.hclust(hc, k=2, border="red")

# Measuring the goodness of fit
r2=R2(df,clusters,2)
cm= table(clusters,is_genuine)[,c(2,1)]; print(cm)
cat("Error Rate:",100 *(1-sum(diag(cm))/sum(cm)),"%\n")

# Examining the single point in cluster 2
df[clusters==2,] ; is_genuine[clusters==2]
op <- par(mfrow = c(3, 2))
for(i in 1:5) { #Index Plots
  plot(df[,i], ylab= colnames(df[i]))
  points(df[clusters==2,], pch = 4, col = 2, cex = 1.5)
} 
par(op)

#------- 1.6 ----------
# Method: intracluster dist: Manhattan --- intercluster dist: Cluster Centroids
d = dist(df, method = "manhattan")
hc = hclust(d, method="centroid")
plot(hc, ) # display dendrogram
mtext('Inter: Manhattan, Intra: Cluster Centroids')
clusters=cutree(hc, k=2) 
# draw dendogram with red borders around them
rect.hclust(hc, k=2, border="red")

# Measuring the goodness of fit
r2=R2(df,clusters,2)
cm= table(clusters,is_genuine)[,c(2,1)]; print(cm)
cat("Error Rate:",100 *(1-sum(diag(cm))/sum(cm)),"%\n")

# Examining cluster 2 points
df[clusters==2,]
op <- par(mfrow = c(3, 2))
for(i in 1:5) { #Index Plots
  plot(df[,i], ylab= colnames(df[i]))
  points(df[clusters==2,], pch = 4, col = 2, cex = 1.5)
} 
par(op)


# -------------- 2.k-Means ------------
k=2; kmc = kmeans(df, k);
clusters=kmc$cluster

r2=R2(df,clusters,2)
cm= table(clusters,is_genuine); print(cm)
cat("Error Rate:",100 *(1-sum(diag(cm))/sum(cm)),"\n")


# Determining the number of clusters
x=scale(df, center=TRUE, scale=TRUE);
wss = (nrow(x)-1)*sum(apply(x,2,var))
for (i in 2:10) { wss[i] = sum(kmeans(x,centers=i)$withinss)}
plot(wss, type="b", pch=19, xlab="k", ylab="WSS", main="The L-Curve")



