library(dplyr)
library(fitdistrplus)
library(ggplot2)
library(pROC)
library(car)
library(PresenceAbsence)
library(leaps)
library(glmnet)
library(bestglm)
library(corrplot)
library(dismo)

Original_Diabetes <- read.csv("~/Desktop/Original_Diabetes.csv")
summary(Original_Diabetes)

Diabetes <- Original_Diabetes %>%
  filter(SkinThickness!=0 & BloodPressure != 0 & BMI!= 0 & Glucose!=0)

summary(Diabetes)

boxplot(Diabetes$Pregnancies~Diabetes$Outcome, col=c(rgb(0.1,0.1,0.7,0.5), rgb(0.8,0.1,0.3,0.6)),
        xlab="Diabetic", ylab="No. of Pregnancy", main="Relationship between Pregnancy and Diabeties")
boxplot(Diabetes$Age~Diabetes$Outcome, col=c(rgb(0.1,0.1,0.7,0.5), rgb(0.8,0.1,0.3,0.6)),
        xlab="Diabetic", ylab="Age", main="Relationship between Age and Diabetes")

names(Diabetes)
pairs(Diabetes[-9])
correlations <- cor(Diabetes[,1:8])
corrplot(correlations, method="circle")

X1 <- Diabetes$Pregnancies
X2 <- Diabetes$Glucose
X3 <- Diabetes$BloodPressure
X4 <- Diabetes$SkinThickness
X5 <- Diabetes$Insulin
X6 <- Diabetes$BMI
X7 <- Diabetes$DiabetesPedigreeFunction
X8 <- Diabetes$Age
Y<- Diabetes$Outcome
#Multiple Linear Regression Model

mlm1 <- lm(Outcome~.,data=Diabetes) 
summary(mlm1)
resid <- residuals(mlm1)

plot(X1, residuals(mlm1), main="Relationship between 
     Pregnancies and residuals", 
     xlab="Pregnancies", ylab="Residuals")
plot(X2, residuals(mlm1), main="Relationship between 
     Glucose and residuals", 
     xlab="Glucose", ylab="Residuals")
plot(X3, residuals(mlm1), main="Relationship between 
     BP and residuals", 
     xlab="BloodPressure", ylab="Residuals")
plot(X4, residuals(mlm1), main="Relationship between 
     SkinThickness and residuals", 
     xlab="SkinThickness", ylab="Residuals")
plot(X5, residuals(mlm1), main="Relationship between 
     Insulin and residuals", 
     xlab="Insulin", ylab="Residuals")
plot(X6, residuals(mlm1), main="Relationship between 
     BMI and residuals", 
     xlab="BMI", ylab="Residuals")
plot(X7, residuals(mlm1), main="Relationship between 
     DiabetesPedigreeFunction and residuals", 
     xlab="DiabetesPedigreeFunction", ylab="Residuals")
plot(X8, residuals(mlm1), main="Relationship between 
     Age and residuals", 
     xlab="Age", ylab="Residuals")

mlm1 <- glm(Outcome~., data=Diabetes)
summary(mlm1)
fnorm <- fitdist(resid,"norm")
summary(fnorm) #reject null hypothesis means regression is linear
plot(fnorm)

#multiple logistic regression
mlog.diabetes <- glm(Outcome~., data=Diabetes, family=binomial)
summary(mlog.diabetes)
        
exp(coef(mlog.diabetes))
exp(cbind(OR=coef(mlog.diabetes),confint(mlog.diabetes)))

pairs(Diabetes, col=Diabetes$Outcome)
vif(mlog.diabetes)

#bestglm - Subset Selection based on AIC, BIC, CV

bmodelAIC <- bestglm(Diabetes,IC="AIC", family = binomial)
summary(bmodelAIC$BestModel)

bmodelBIC <- bestglm(Diabetes,IC="BIC",family = binomial)
summary(bmodelBIC$BestModel)

RNGkind(sample.kind = "Rounding")
set.seed(100)

bmodelCV <- bestglm(Diabetes,IC="CV",CVArgs = list(Method="HTF",K=10,REP=1),family = binomial)
summary(bmodelCV$BestModel)

bmodelAIC$Subsets
bmodelBIC$Subsets

plot1 = plot(x=c(0:8),bmodelAIC$Subsets$AIC, xlab = "Number of Variables",ylab = "AIC",
             main = "Best Selection with AIC", type = "b")
plot2 = plot(x=c(0:8),bmodelBIC$Subsets$BIC,xlab = "Number of Variables",ylab = "BIC",
             main = "Best Selection with BIC", type = "b")

#lasso prediction

RNGkind(sample.kind = "Rounding")
set.seed(100)
grid <- 10^seq(10,-2,length = 100)

train <- sample(1:nrow(Diabetes),nrow(Diabetes)/2)
test <- (-train)
Diabetes.train <- Diabetes[train,]
Diabetes.test <- Diabetes[test,]

X <- model.matrix(Diabetes$Outcome~.,Diabetes)[,-1] 
Y.test <- Y[test]
Y.train <- Y[train]

lasso.mod <- glmnet(X[train,], Y[train], alpha = 1, lambda = grid, family = "binomial")
plot(lasso.mod)

RNGkind(sample.kind = "Rounding")
set.seed(100)

cv.out <- cv.glmnet(X[train,] , Y[train], alpha = 1 ,lambda = grid, family = "binomial")
cv.out
plot(cv.glmnet(X[train,] , Y[train], alpha = 1,family="binomial"))

bestlam <- cv.out$lambda.min
bestlam

#test
bmodellasso <- glmnet(X,Y, alpha = 1, family = "binomial",lambda = bestlam)

X.test <- model.matrix(Outcome ~., Diabetes.test)[,-1]

probabilities <- bmodellasso%>% predict(newx = X.test)
predicted.classes <- ifelse(probabilities > 0.25, 1, 0)

table(Diabetes.test$Outcome,predicted.classes)
mean(predicted.classes==Diabetes.test$Outcome) 

lasso.coef <- predict(bmodellasso, type = "coefficients", s = bestlam)
lasso.coef
lasso.coef[lasso.coef!=0]

prob.lasso <- NULL
for (i in 1:nrow(Diabetes)){
  w <- exp(-8.798973375 + 0.102940318 * X1[i] + 0.032295005 * X2[i] + 0.003706411 * X4[i] + 0.068979711* X6[i] + 1.063010860 * X7[i] + 0.021402651 * X8[i])
  prob.lasso[i]<- w/(1+w)
}


#accuracy

prob1 <- predict(bmodelAIC$BestModel, type="response")
pred1 <- rep("No Diabetes", 532)
pred1[prob1 > 0.25] <- "Diabetes"
table1 <- table(pred1, Diabetes$Outcome)
table1 
accuracy1 <- (table1[1,2]+table1[2,1])/length(Diabetes$Outcome)
accuracy1

prob2 <- predict(bmodelBIC$BestModel, type="response")
pred2 <- rep("No Diabetes", 532)
pred2[prob2 > 0.25] <- "Diabetes"
table2 <- table(pred2, Diabetes$Outcome)
table2
accuracy2 <- (table2[1,2]+table2[2,1])/length(Diabetes$Outcome)
accuracy2

prob3 <- predict(bmodelCV$BestModel, type="response")
pred3 <- rep("No Diabetes", 532)
pred3[prob3 > 0.25] <- "Diabetes"
table3 <- table(pred3, Diabetes$Outcome)
table3
accuracy3 <- (table3[1,2]+table3[2,1])/length(Diabetes$Outcome)
accuracy3

pred4 <- rep("No Diabetes", 532)
pred4[prob.lasso > 0.25] <- "Diabetes"
table4 <- table(pred4, Diabetes$Outcome)
table4
accuracy4 <- (table4[1,2]+table4[2,1])/length(Diabetes$Outcome)
accuracy4

accuracylist <- c(accuracy1, accuracy2, accuracy3, accuracy4)
accuracylist 

bmodelAIC$BestModel
bmodelBIC$BestModel
bmodelCV$BestModel

#ROC
roc(Diabetes$Outcome, bmodelAIC$BestModel$fitted.values, plot = TRUE, percent = TRUE, xlab = "False positive rate", ylab = "False negative rate", legacy.axes = TRUE, col = "Red", print.auc = TRUE)
plot.roc(Diabetes$Outcome, bmodelBIC$BestModel$fitted.values, add = TRUE , col = "Blue",print.auc = TRUE, print.auc.y = 40, percent = TRUE)
plot.roc(Diabetes$Outcome, bmodelCV$BestModel$fitted.values, add = TRUE , col = "Brown",print.auc = TRUE, print.auc.y = 30, percent = TRUE)
plot.roc(Diabetes$Outcome, prob.lasso, add = TRUE , col = "Green",print.auc = TRUE, print.auc.y = 20, percent = TRUE)

