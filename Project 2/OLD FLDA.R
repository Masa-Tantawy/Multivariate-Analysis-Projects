# Multivariate Analysis Project 2
# Laila El Saadawi: 900201723
# Masa Tantawy: 900201312

library(DescTools)
library(Matrix)
library(GGally)
library(MASS)
library(nnet)

df = read.csv("fakebills.csv", header= TRUE, sep=";")
dim(df); View(head(df)); str(df)

#------------------- Data Preparation------------------------------
# 1. Categorial variable is_genuine
df$is_genuine= as.factor(df$is_genuine)
#df[df$is_genuine== 'True', "is_genuine"]= 1
#df[df$is_genuine== 'False', "is_genuine"]= 0
#df$is_genuine= as.numeric(df$is_genuine)
#View(head(df))

# 2. Missing values
sum(is.na(df)) # the data contains missing values
colSums(is.na(df)) > 0 #only margin low has missing values: 37 observations
dim(df[rowSums(is.na(df)) > 0, ])
# Replacing missing values with the mean
df[is.na(df[,'margin_low']), 'margin_low']= mean(df[,'margin_low'], na.rm = TRUE)
sum(is.na(df))

#----------------- Fisher Linear Discriminant Analysis-------------------------
data=df; x=df[,2:7]; class=df[,1]
pairs(x,pch=19,col=class)

#Normality: Check 1
ggpairs(df, columns = 2:7, aes(color = is_genuine, alpha = 0.5))   

#Normality: Check 2
dfTR <- df[df$is_genuine == "True", ]
dfFL <- df[df$is_genuine != "True", ]

par(mfrow = c(2, 3)) 
for(i in colnames(x)) { 
  qqnorm(dfTR[[i]]); qqline(dfTR[[i]], col = 2) 
} #genuine 

par(mfrow = c(2, 3)) 
for(i in colnames(x)) { 
  qqnorm(dfFL[[i]]); qqline(dfFL[[i]], col = 2) 
} #fake 

#Equality of Cov Matrices Check
x_TR=dfTR[,2:7]
x_FL=dfFL[,2:7]
cov(x_TR)==cov(x_FL) #false CHECK

flda=function(x,class) {
  cat("Fisher Linear Discriminant:\n")
  a = lda(x, class); d = predict(a) 
  t=table(class, d$class); print(t) 
  er=100*(sum(t)-sum(diag(t)))/nrow(x) 
  cat("Error Rate =",er,"%\n")
  return(d) }

int_rslt=flda(df[,2:7],df[,1]) #generates internal validation table & error rate
#Error rate is 1.067% 

#External Validation: Leave One Out CV
loo=function(x,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
    a = lda(x[-i,], class[-i])
    b = predict(a,x[i,])
    rslt[i]=b$class #[i]==class[i]
  }
  return(rslt)}

ext_rslt=loo(x,class)
t2=table(class,ext_rslt)
cat("\nExternal validation:\n")
print(t2)

er=100*(sum(t2)-sum(diag(t2)))/nrow(df[,2:7]) 
cat("Error Rate =",er,"%\n") #same error rate of 1.06%

#Splitting the Data into Training & Testing:

set.seed(666)
ind <- sample(2, nrow(df),
              replace = TRUE,
              prob = c(0.6, 0.4))
training <- df[ind==1,]
testing <- df[ind==2,]
linear <- lda(training$is_genuine~., training)
linear

p <- predict(linear, training)

#hist of two labels from training set:
ldahist(data = p$x[,1], g = training$is_genuine) 

p.df <- data.frame(LD1 = p$x, class = p$class) #convert pred to df
ggplot(p.df) + geom_density(aes(LD1, fill = p$class), alpha = 0.2)


#FLDA on Training Set
p1  <- predict(linear, training)$class
table1 <- table(Predicted = p1, Actual = training$is_genuine)
table1

acc1 = 100* (sum(diag(table1))/sum(table1))
cat("Accuracy =",acc1,"%\n")
er1 = 100*(sum(table1)-sum(diag(table1)))/nrow(training) 
cat("Error Rate =",er1,"%\n")

#EXTERNAL????

#Multinomial Discriminant Analysis
mn = multinom(is_genuine ~ ., data = df) 
mn_results = predict(mn) 
table(df$is_genuine, mn_results)


