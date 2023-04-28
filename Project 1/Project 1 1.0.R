# Multivariate Analysis Project 1
# Laila El Saadawi: 900201723
# Masa Tantawy: 900201312

library(DescTools)
library(Matrix)
library(GGally)

df = read.csv("babies.csv", header= TRUE)
dim(df); View(head(df)); str(df)

#------------------- Data Preparation------------------------------
# Removing the case/index column
df= df[,-1] 

# 1. missing values
sum(is.na(df)) # the data contains missing values
dim(df[rowSums(is.na(df)) > 0, colSums(is.na(df)) > 0 ])
# 62 rows and 5 columns have missing values
colSums(is.na(df)) # Columns with NAs: gestation, age, height, weight, smoke

# Replacing missing values with the median for continuous columns
for(i in c(2,4,5,6)){
  df[is.na(df[,i]), i] <- median(df[,i], na.rm = TRUE)
}
# For the variable smoke, using 3-NN imputations to replace missing values
df= ImputeKnn(df[,],k = 3, scale = T, meth = "weighAvg", distData = NULL)
df[df$smoke> 0.5, 'smoke']= 1
df[df$smoke< 0.5, 'smoke']= 0
unique(df$smoke)
str(df)
colSums(is.na(df)) 

# 2. Checking column rank
p= dim(df)[2] # p=7
rankMatrix(df)[1]; rankMatrix(df)[1]==p

# 3. Columns
summary(df)
cor(df)
ggpairs(df, title= 'Pairwise Plot of df') 

#-------------------- Step 1: Detection of Outliers -------------------
require(robustX); library(robustbase);
# Using BACON to detect outliers in numeric columns
new_p= p-2
df_numeric= df[,-c(3,7)]
output= mvBACON(as.matrix(df_numeric))
output$limit
y = cbind(1:nrow(df_numeric),output$dis); colnames(y)= c("Index","Distance");
plot(y, pch=19, main = "BACON Distances", ylim=c(0,6))
points(y[ ! output$subset, ], pch = 4, col = 2, cex = 1.5)
abline(h=output$limit, col= "red", lty=2)

#install.packages('heplots');
library(heplots)
op <- par(mfrow = c(2, 2), pty = "s")

# Plot 1: Q-Q Plot of Squared MD vs Quantiles of Chisquare(p)
cqplot(df_numeric,method='classical', pch=19,
       main= expression('Q-Q Plot of Squared MD vs'~ paste(chi[5]^2)~'Quantlile'),
       ref.col='blue')

# Plot 2: Index Plot of Mahalanobis Distances
d = sqrt(mahalanobis(df_numeric, colMeans(df_numeric), var(df_numeric)))
plot(d,pch=19,main='Index Plot of Mahalanobis Distances',
     ylab='Mahalanobis Distance')
abline(h=sqrt(qchisq(0.95,new_p)),col='red')

# Plot 3: Q-Q Plot of Squared BD vs Quantiles of Chisquare(5)
baconsq=(output$dis)^2
qqplot(qchisq(ppoints(nrow(df_numeric)),df=new_p), baconsq, pch=20,
       ylab='Squared Bacon Distances',
       xlab= expression('Q-Q Plot of Squared BD vs'~ paste(chi[5]^2)~'Quantlile'), 
       main= expression('Q-Q Plot of Squared BD vs'~ paste(chi[5]^2)~'Quantlile'));
qqline(output$dis,col='blue')

# Plot 4: Index Plot of Bacon Distances
plot(output$dis,pch=19, main='Index Plot of Bacon Distances',
     ylab='Bacon Distance' , ylim= c(0,6))
abline(h=output$limit,  col= "red",  lty=2)
points(y[ ! output$subset, ], pch = 4, col = 2, cex = 1.5)
par(op)

outliers= y[!output$subset, 1]
#------------------ Step 2: Hotelling's T test --------------------
# Separating the data into groups
first= df[df$parity == 0,-c(3, 7)]
not_first= df[df$parity == 1, -c(3,7)]

smokers= df[df$smoke == 1,-c(3, 7)]
non_smokers= df[df$smoke == 0, -c(3, 7)]

#-------- FIRST: Testing the 2 groups first vs not first pregnancy -------
# 1. Checking the assumption of normality and equal variances
# 1.1: using histogram
# First Pregnancy
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(first)){
  names=c()
  hist(first[,i],main = colnames(first)[i], col= 'deeppink4', 
       xlab = colnames(first)[i])
} 
mtext("First Pregnancy", side = 3, line = -2, outer = TRUE)
par(op)
# Not First Pregnancy
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(not_first)){
  names=c()
  hist(not_first[,i],main = colnames(not_first)[i], col= 'deeppink4',
       xlab = colnames(not_first)[i])
}
mtext("Not First", side = 3, line = -2, outer = TRUE)
par(op)

# 1.2: Transformations
# The variables age and weight do not really follow a normal distribution,
# so we will use log transformation to make them normal
# Age:
hist(not_first[,'age'], xlab= 'Age', main='Before Transformation', col= 'deeppink4')
hist(log(not_first[,'age']), xlab= 'log(Age)', main='After Transformation', col= 'deeppink4')
# Weight:
hist(not_first[,'weight'], xlab= 'Weight', main='Before Transformation', col= 'deeppink4')
hist(log(not_first[,'weight']), xlab= 'log(Weight)', main='After Transformation', col= 'deeppink4')

first$age= log(first$age)
first$weight= log(first$weight)
not_first$age= log(not_first$age)
not_first$weight= log(not_first$weight)

# Histograms after transformation:
# First Pregnancy
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(first)){
  names=c()
  hist(first[,i],main = colnames(first)[i], col= 'deeppink4', 
       xlab = colnames(first)[i])
} 
mtext("First Pregnancy", side = 3, line = -2, outer = TRUE)
par(op)
# Not First Pregnancy
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(not_first)){
  names=c()
  hist(not_first[,i],main = colnames(not_first)[i], col= 'deeppink4',
       xlab = colnames(not_first)[i])
}
mtext("Not First", side = 3, line = -2, outer = TRUE)
par(op)

# 1.3: equality of variances
all.equal(cov(first),cov(not_first),check.attributes=FALSE) 

# 2. Testing the equality of means
# Non-Robust test (we know that the data has outliers so we shouldn't depend on it)
ht2 <- function(X, Y){
  n <- nrow(X); m <- nrow(Y)
  delta <- colMeans(X) - colMeans(Y)
  p <- ncol(X)
  xcov <- cov(X); ycov <- cov(Y)
  S_pooled <- ((n-1)*xcov + (m-1)*ycov)/(n+m-2)
  T2=t(delta)%*%solve(S_pooled)%*%delta
  T2=T2*n*m/(n+m)
  F <- T2 * (n+m-p-1)/(p*(n+m-2))
  p_value <- 1-pf(F,p,n+m-p-1)
  cat("F :");cat(F)
  cat("\n T2:"); cat(T2)
  cat ("\n DOF:");cat({p} ,"and", {n+m-p-1})
  cat("\n P-value:"); cat(p_value)
  if (p_value < 0.05){
    cat('\n Reject H_0')
    cat('\n Population Means are NOT Equal')
    decision=('Reject')
  }else{
    cat('\n Do Not Reject H_0')
    cat('\n Population Means are EQUAL')
    decision=('not reject')}
  return(list(df=c(p,n+m-p-1),F=F,T2=T2 ,p_value=p_value,decisicion=decision))
}
result= ht2(first, not_first)
f_tab = qf(p=0.05, df1=6, df2=1456,lower.tail = F);
result$F > f_tab

# The data has outliers so we should use robust test. 
# Robust  test: 
library(rrcov)
result2=T2.test(first, not_first)
result2

#-------- SECOND: Testing the 2 groups smokers vs non-smokers -------
# 1. Checking the assumption of normality 
# 1.1: using histograms
# Smokers
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(smokers)){
  names=c()
  hist(smokers[,i],main = colnames(smokers)[i], col= 'deeppink4', 
       xlab = colnames(smokers)[i])
} 
mtext("Smokers", side = 3, line = -2, outer = TRUE)
par(op)
# Non-smokers
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(non_smokers)){
  names=c()
  hist(non_smokers[,i],main = colnames(non_smokers)[i], col= 'deeppink4',
       xlab = colnames(non_smokers)[i])
}
mtext("Non-Smokers", side = 3, line = -2, outer = TRUE)
par(op)

# 1.2: Transformations
# The variables age and weight do not really follow a normal distribution,
# so we will use log transformation to make them normal
# Age:
hist(smokers[,'age'], xlab= 'Age', main='Before Transformation', col= 'deeppink4')
hist(log(smokers[,'age']), xlab= 'log(Age)', main='After Transformation', col= 'deeppink4')
# Weight:
hist(smokers[,'weight'], xlab= 'Weight', main='Before Transformation', col= 'deeppink4')
hist(log(smokers[,'weight']), xlab= 'log(Weight)', main='After Transformation', col= 'deeppink4')

smokers$age= log(smokers$age)
smokers$weight= log(smokers$weight)
non_smokers$age= log(non_smokers$age)
non_smokers$weight= log(non_smokers$weight)

# Histograms after transformation:
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(smokers)){
  names=c()
  hist(smokers[,i],main = colnames(smokers)[i], col= 'deeppink4', 
       xlab = colnames(smokers)[i])
} 
mtext("Smokers", side = 3, line = -2, outer = TRUE)
par(op)
# Non-smokers
op <- par(mfrow = c(3, 2), pty = "s")
for (i in 1:ncol(non_smokers)){
  names=c()
  hist(non_smokers[,i],main = colnames(non_smokers)[i], col= 'deeppink4',
       xlab = colnames(non_smokers)[i])
}
mtext("Non-Smokers", side = 3, line = -2, outer = TRUE)
par(op)


# 1.3: equality of variances
all.equal(cov(smokers),cov(non_smokers),check.attributes=FALSE) 

# 3. Testing the equality of means
# Non-Robust test (we know that the data has outliers so we shouldn't depend on it)
ht2 <- function(X, Y){
  n <- nrow(X); m <- nrow(Y)
  delta <- colMeans(X) - colMeans(Y)
  p <- ncol(X)
  xcov <- cov(X); ycov <- cov(Y)
  S_pooled <- ((n-1)*xcov + (m-1)*ycov)/(n+m-2)
  T2=t(delta)%*%solve(S_pooled)%*%delta
  T2=T2*n*m/(n+m)
  F <- T2 * (n+m-p-1)/(p*(n+m-2))
  p_value <- 1-pf(F,p,n+m-p-1)
  cat("F :");cat(F)
  cat("\n T2:"); cat(T2)
  cat ("\n DOF:");cat({p} ,"and", {n+m-p-1})
  cat("\n P-value:"); cat(p_value)
  if (p_value < 0.05){
    cat('\n Reject H_0')
    cat('\n Population Means are NOT Equal')
    decision=('Reject')
  }else{
    cat('\n Do Not Reject H_0')
    cat('\n Population Means are EQUAL')
    decision=('not reject')}
  return(list(df=c(p,n+m-p-1),F=F,T2=T2 ,p_value=p_value,decisicion=decision))
}
result= ht2(smokers, non_smokers)
f_tab = qf(p=0.05, df1=6, df2=1456,lower.tail = F);
result$F > f_tab

# The data has outliers so we should use robust test. 
# Robust  test: 
library(rrcov)
result2=T2.test(smokers, non_smokers)
result2
