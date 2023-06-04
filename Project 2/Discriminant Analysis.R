# Multivariate Analysis Project 2
# Laila El Saadawi: 900201723
# Masa Tantawy: 900201312

# ............. Discriminant Analysis .............
library(DescTools)
library(Matrix)
library(GGally)
library(MASS)
library(nnet)
library(dplyr)

df = read.csv("fake_bills.csv", header= TRUE, sep=";")
dim(df); View(head(df)); str(df)

#------------------- Data Preparation------------------------------
# *Note: Please check R script Data Preparation for details

# 1. Categorical variable is_genuine
df$is_genuine= as.factor(df$is_genuine)
# 2. Missing values
# Replacing missing values with the mean
df[is.na(df[,'margin_low']), 'margin_low']= mean(df[,'margin_low'], na.rm = TRUE)

x=df[,c(2:7)]; class=df[,1]
dfTR <- df[df$is_genuine == "True", ]
dfFL <- df[df$is_genuine != "True", ]

# --- Checking the Assumptions for Fisher Linear Discriminant Analysis-------
# Normality Check:
ggpairs(df[,-1], aes(color=df$is_genuine,alpha=0.5), title= 'Pairwise Plot of df') + 
  scale_colour_manual(values=c('pink','black'))

op= par(mfrow = c(2, 3)) 
for(i in colnames(x)) { #Q-Q Plots for genuine bills
  qqnorm(dfTR[[i]], main= i); qqline(dfTR[[i]], col = 2) 
} 
for(i in colnames(x)) { #Q-Q Plots for fake bills
  qqnorm(dfFL[[i]], main= i); qqline(dfFL[[i]], col = 2) 
} 
par(op)

#Equality of Covariance Matrices Check:
x_TR=dfTR[,2:7]
x_FL=dfFL[,2:7]
cov(x_TR)==cov(x_FL)
library('covequal')
test_covequal(data.matrix(df[df$is_genuine=='True',-1]),
              data.matrix(df[df$is_genuine=='False',-1]))

# ----------------------- 1. FLDA -----------------
# -------------- 1.1: FLDA on all the data --------------
# Defining Functions
flda=function(x,class) {
  cat("Fisher Linear Discriminant:\n")
  a = lda(x, class); d = predict(a) 
  t=table(class, d$class); print(t) 
  er=100*(sum(t)-sum(diag(t)))/nrow(x) 
  cat("Error Rate =",er,"%\n")
  return(d) }
loo=function(x,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
    a = lda(x[-i,], class[-i])
    b = predict(a,x[i,])
    rslt[i]=b$class #[i]==class[i]
  }
  return(rslt)}

int_rslt=flda(x,class) # Internal validation
ext_rslt=loo(x,class) # External validation
ext_table=table(class,ext_rslt)
cat("\nExternal validation:\n"); print(ext_table)
cat("Error Rate=",100 *(1-sum(diag(ext_table))/sum(ext_table)),"%\n")

#Plotting Vector Projections
lda_pr= lda(df$is_genuine~., df)
p_pr <- predict(lda_pr, df)

plot(p_pr$x ,col="purple", ylab="(cm)", lwd=1, font.lab=2, pch=19,
     main="Index Plot of Vector Projections",
     font=2)


#-------- 1.2: FLDA on split data (75%-25% test-train split) -----------
#Splitting the Data into Training & Testing:
set.seed(666)
ind <- sample(2, nrow(df), replace = TRUE, prob = c(0.75, 0.25))
training <- df[ind == 1,]; testing <- df[ind==2,]

# FLDA on Training
flda_split= lda(training$is_genuine~., training); flda_split
p <- predict(flda_split, training)
# hist of two labels from training set:
ldahist(data = p$x[,1], g = training$is_genuine) 
# joint pdf plot
p.df <- data.frame(LD1 = p$x, class = p$class)
ggplot(p.df) + geom_density(aes(LD1, fill = p$class), alpha = 0.2)

#Training Set Output of LDA
table1 <- table(Predicted = p$class, Actual = training$is_genuine)
table1

acc1 = 100* (sum(diag(table1))/sum(table1))
cat("Accuracy =",acc1,"%\n")
er1 = 100*(sum(table1)-sum(diag(table1)))/nrow(training) 
cat("Error Rate =",er1,"%\n")

#Testing Set
p2 <- predict(flda_split, testing)$class
table2 <- table(Predicted = p2, Actual = testing$is_genuine)
table2

acc2 = 100* (sum(diag(table2))/sum(table2))
cat("Accuracy =",acc2,"%\n")
er2 = 100*(sum(table2)-sum(diag(table2)))/nrow(testing) 
cat("Error Rate =",er2,"%\n")

#---------------- 2. FLDA 2 ----------------
# Fisher Linear Discriminant Analysis; p = 2 and k = 2
# Defining Function
flda2=function(x, class) { 
  if(ncol(x)!=2) {cat("Data should be 2-dimensional\n" ); return()}
  t=factor(class);   level=levels(t);   
  if(length(level)!=2) {cat("Data should have only two groups\n" ); return()}
  y=x[class==level[1],];     x=x[class==level[2],]
  n1=nrow(x);   n2=nrow(y); n=n1+n2;  xcenter=colMeans(x);  ycenter=colMeans(y)
  xcov=cov(x); ycov=cov(y);   sp=(n1-1)*xcov+(n2-1)*ycov; sp=sp/(n-2)
  d=xcenter-ycenter; m=(xcenter+ycenter)/2;  a=solve(sp)%*%d;  
  class=c(rep(1,n1),rep(2,n2));      p=1;      
  z=rbind(x,y);       pred=z-matrix(m,ncol=2,nrow=n, byrow=T); pred=as.matrix(pred)
  pred=(pred%*%a)<log(p);    C=(class!=pred+1);     ce=sum(C)
  cat("--------------------------------------------------\n")
  cat("           Correct         Incorrect\n Class  Classification   Classification    Total\n")  
  cat("--------------------------------------------------\n")
  cd1=n1-sum(C[1:n1]);         cat("  1          ",cd1,"             ", n1-cd1,"            ",n1,"\n")
  cd2=n2-sum(C[(n1+1):n]);  cat("  2          ",cd2,"             ", n2-cd2,"            ",n2,"\n")
  cat("--------------------------------------------------\n")
  cat(" Total:     ",cd1+cd2,"             ", n-(cd1+cd2),"           ",n,"\n")  
  cat("--------------------------------------------------\n")
  cat("Error Rate = ",100*(ce/n),"%\n");  const=(sum(a*m)+log(p))/a[2]; slope=-a[1]/a[2]
  z=rbind(x,y);   print(rbind(xcenter[1:2],ycenter[1:2]))
  plot(z[,c(1,2)],col=class,pch=19);  abline(const,slope, col = 'blue'); 
  points(rbind(xcenter[1:2],ycenter[1:2]),pch=19,col=3,cex=1.5); 
  segments(xcenter[1], xcenter[2], ycenter[1], ycenter[2],col=3)
  list(xcenter=xcenter[2:1],ycenter=ycenter[2:1],xcov=xcov,ycov=ycov,sp=sp,a=a,slope=slope,
       const=const,ce=ce,m=m,z=z) }

# 2.1: margin_low and length
pair1 = df[ , c("margin_low", "length", "is_genuine")] 
x1 = pair1[,1:2]; class1=pair1[,3]
f2_first=flda2(x1,class1); f2_first$slope

# 2.2: margin_up and length
pair2 = df[ , c("margin_up", "length", "is_genuine")] 
x2 = pair2[,1:2]; class2=pair2[,3]
f2_second=flda2(x2,class2); f2_second$slope

# 2.3: height_right and length
pair3 = df[ , c("height_right", "length", "is_genuine")] 
x3 = pair3[,1:2]; class3=pair3[,3]
f2_third=flda2(x3,class3); f2_third$slope

# ---------------- 3. Multinomial Discriminant Analysis ----------------
mn = multinom(is_genuine ~ ., data = df) 
mn_internal= predict(mn) 
looMN=function(x_subset,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
    y=x_subset[-i,]
    a = multinom(class[-i]~. , x_subset[-i,])
    b = predict(a,x_subset[i,])
    rslt[i]=b
  }
  return(rslt)
}
mn_external=looMN(x,class)

cat("Multinomial")
cat("\nExternal validation:\n")
table(df$is_genuine, mn_internal)
cat("Error Rate=",100 *(1-sum(diag(table(df$is_genuine, mn_internal)))/
                          sum(table(df$is_genuine, mn_internal))),"%\n")
cat("\nExternal validation:\n")
table(df$is_genuine, mn_external)
cat("Error Rate=",100 *(1-sum(diag(table(df$is_genuine, mn_external)))/
                          sum(table(df$is_genuine, mn_external))),"%\n")
