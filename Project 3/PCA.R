# Multivariate Analysis Project 3
# Laila El Saadawi: 900201723
# Masa Tantawy: 900201312

#------------- Principle Component Analysis---------------
library(DescTools)
library(Matrix)
library(GGally)

df = read.csv("India_rainfall_act_dep_1901_2016_1.csv", header= TRUE)
dim(df); View(head(df)); str(df)

# ------------DATA PREPARATION-----------
#removing the non-numeric column
df= df[,-1] 
# renaming the columns
colnames(df)[1:10]= c('AR_June','AR_July','AR_Aug','AR_Sep', 'AR_June_Sep',
                      'DP_June','DP_July','DP_Aug','DP_Sep', 'DP_June_Sep')

# 1. Checking Missing values
dim(df)
sum(is.na(df)) # the data contains no missing values

# 2. Outlier Detection
require(robustX); library(robustbase);
output= mvBACON(as.matrix(df))
output$limit
y = cbind(1:nrow(df),output$dis); colnames(y)= c("Index","Distance");
plot(y, pch=19, main = "BACON Distances")
abline(h=output$limit, col= "red", lty=2)
outliers= y[!output$subset, 1]
length(outliers)

# 3. Scaling the data
df= as.data.frame(scale(df))

# 4. Correlation 
cor= cor(df) #View(cor)
library(ggcorrplot) #install.packages("ggcorrplot")
ggcorrplot(cor,lab =TRUE,title='Correlation Matrix')

# --------------- PCA --------
# 1. Non-Robust
pc=princomp(df, cor=T)

#summary(pc,loadings=T); library(psych); #tr(cor)/13

v=pc$sd^2; sv=sum(v); cs=cumsum(v)/sv; #m=length(cs[cs<0.95])
rslt=round(cbind(v,v/sum(v),cumsum(v)/sv),3)
colnames(rslt)=c("Variance","Prop. of Var.","Cum. Prop.")
print(rslt)

#plot(pc) # Scree plot;
pcaCharts <- function(x) {
  v <- x$sd ^ 2
  sv=sum(v)
  pv <- v/sv
  vp<- round(pv*100,1)
  
  par(mfrow=c(2,2))
  barplot(vp,xlab="PCs", 
          ylab="Percentage of Variance", 
          ylim=c(0,100),col='steelblue',
          main='Variation Plot')
  barplot(cumsum(vp),xlab="PCs",
          ylab="Cumulative Percentage of Variance", 
          ylim=c(0,100),col='steelblue',
          main='Cumulative Variation Plot')
  screeplot(x,col='steelblue',ylim=c(0,5),main='Scree Plot')
  plot(cumsum(v),pch=19,main='Cumulative Scree Plot',
       ylim=c(0,max(cumsum(v))+1),
       ylab='Cumulative Variance',type='b',xlab='PCs')
  par(mfrow=c(1,1))
}
pcaCharts(pc)

# 2. Robust
pc=princomp(df[-outliers,], cor=T)

v=pc$sd^2; sv=sum(v); cs=cumsum(v)/sv; #m=length(cs[cs<0.95])
rslt=round(cbind(v,v/sum(v),cumsum(v)/sv),3)
colnames(rslt)=c("Variance","Prop. of Var.","Cum. Prop.")
print(rslt)

#plot(pc) # Scree plot;
pcaCharts(pc)

