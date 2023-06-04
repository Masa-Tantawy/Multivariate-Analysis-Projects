# Multivariate Analysis Project 2
# Laila El Saadawi: 900201723
# Masa Tantawy: 900201312

#------------- EDA & Data Preparation---------------
library(DescTools)
library(Matrix)
library(GGally)

df = read.csv("fake_bills.csv", header= TRUE, sep=";")
dim(df); View(head(df)); str(df)

# 1. Categorical variable is_genuine
df$is_genuine= as.factor(df$is_genuine)
levels(df$is_genuine)

# 2. Missing values
sum(is.na(df)) # the data contains missing values
colSums(is.na(df)) > 0 #only margin low has missing values: 37 observations
dim(df[rowSums(is.na(df)) > 0, ])
# Replacing missing values with the mean
df[is.na(df[,'margin_low']), 'margin_low']= mean(df[,'margin_low'], na.rm = TRUE)
sum(is.na(df))

# 3. Checking column rank
str(df)
dim(df)
p= dim(df[,-1])[2] # p=7
rankMatrix(df[,-1])[1]; rankMatrix(df[,-1])[1]==p

# 4. Outlier Detection
require(robustX); library(robustbase);

output= mvBACON(as.matrix(df[,-1]))
output$limit
y = cbind(1:nrow(df),output$dis); colnames(y)= c("Index","Distance");
plot(y, pch=19, main = "BACON Distances", ylim=c(0,6))
abline(h=output$limit, col= "red", lty=2)
# data has no outliers

# Checking for outliers in each group separately
dfTR <- df[df$is_genuine == "True", ]
dfFL <- df[df$is_genuine != "True", ]

output= mvBACON(as.matrix(dfTR[,-1]))
output= mvBACON(as.matrix(dfFL[,-1]))
# data has no outliers

# 5. Variables Relationship
summary(df)
# correlation
cor(df[,-1])
# Pairwise Scatter plots
#plot(df[,-1], col=df$is_genuine, pch=20) 
ggpairs(df[,-1], aes(color=df$is_genuine,alpha=0.5), 
        title= 'Pairwise Plot of df') + 
  scale_colour_manual(values=c('pink','black'))

head(df)





