##========  step 0 prelim  ============##
library(data.table);
library(glmnet);
library(randomForest);
library(ggplot2);
library(e1071);
library(ada);
library(xgboost);
library(nnet);
library(klaR);
load("~/Dropbox/prADictor/01_loading_datasets/snp_info.RData")
load("~/Dropbox/prADictor/01_loading_datasets/genotype_with_grp.RData")
load("~/Dropbox/prADictor/01_loading_datasets/data_ready_2.RData")

rm(Raw_Geno)

##----------------
SNP_Info3 <- SNP.Keep3
grep("snp_", x = SNP_Info3, ignore.case = T)
SNP_Info3[c(3,6)] <- c("snp_a1892800", "snp_a2056758")

X_Train <- Data_Train[ ,-c(1:6)]
X_Train <- X_Train[ ,SNP.Keep3]
colnames(X_Train) <- SNP_Info3
y_Train <- as.factor(Data_Train[ ,6])

X_Test <- Data_Test[ ,-c(1:6)]
X_Test <- X_Test[ ,SNP.Keep3]
colnames(X_Test) <- SNP_Info3
y_Test <- as.factor(Data_Test[ ,6])
##----------------


##======== final model ========##




##========  step 1 generalized linear models  ============##
## Final Pick -ElasticNet
Fit1 <- cv.glmnet(x = as.matrix(X_Train), y = y_Train,
                   alpha = 0.2, standardize = F,
                   nfolds = 5, family = "binomial");
yPred1 <- predict(Fit1, newx = as.matrix(X_Test), type="class")[ ,1]
cat("Precision: ", mean(yPred1==y_Test), "\n") # 0.60
cat("Other Metrics are: ", "\n", sep=" ")
mean(yPred1==1&y_Test==1) # 0.27
mean(yPred1==0&y_Test==0) # 0.33
mean(yPred1==1&y_Test==0) # 0.19
mean(yPred1==0&y_Test==1) # 0.21

# Fit1a <- cv.glmnet(x = as.matrix(X_Train), y = y_Train,
#                    alpha = 1, standardize = F,
#                    nfolds = 5, family = "binomial");
# yPred1 <- predict(Fit1a, newx = as.matrix(X_Test), type="class")[ ,1]
# cat("Precision: ", mean(yPred1==y_Test), "\n") # 0.5616
# cat("Other Metrics are: ", "\n", sep=" ")
# mean(yPred1==1&y_Test==1) # 0.26
# mean(yPred1==0&y_Test==0) # 0.30
# mean(yPred1==1&y_Test==0) # 0.22
# mean(yPred1==0&y_Test==1) # 0.22
# 
# Fit1b <- cv.glmnet(x = as.matrix(X_Train), y = y_Train,
#                    alpha = 0, standardize = F,
#                    nfolds = 5, family = "binomial");
# yPred1 <- predict(Fit1b, newx = as.matrix(X_Test), type="class")[ ,1]
# cat("Precision: ", mean(yPred1==y_Test), "\n") # 0.5479
# cat("Other Metrics are: ", "\n", sep=" ")
# mean(yPred1==1&y_Test==1) # 0.32
# mean(yPred1==0&y_Test==0) # 0.37
# mean(yPred1==1&y_Test==0) # 0.15
# mean(yPred1==0&y_Test==1) # 0.16
# 
# Fit1c <- cv.glmnet(x = as.matrix(X_Train), y = y_Train,
#                    alpha = 0.2, standardize = F,
#                    nfolds = 5, family = "binomial");
# yPred1 <- predict(Fit1c, newx = as.matrix(X_Test), type="class")[ ,1]
# cat("Precision: ", mean(yPred1==y_Test), "\n") # 0.60 
# cat("Other Metrics are: ", "\n", sep=" ")
# mean(yPred1==1&y_Test==1) # 0.27
# mean(yPred1==0&y_Test==0) # 0.33
# mean(yPred1==1&y_Test==0) # 0.19
# mean(yPred1==0&y_Test==1) # 0.21


##========  step 2 random forest  ============##
## Final Pick --ntree 1000 --mtry 50 --nodesize=20
Fit2 <- randomForest(y_Train ~ ., data=data.frame(y_Train,X_Train), 
                     ntree=1000, mtry=50, nodesize=20,
                     importance = TRUE, keep.forest = TRUE);
yPred2 <- predict(Fit2, newdata=data.frame(y_Test, X_Test)[ ,-1])
cat("Precision: ", mean(yPred2==y_Test), "\n") # 0.74
cat("Other Metrics are", "\n", sep=" ")
mean(yPred2==1&y_Test==1) # 0.33
mean(yPred2==0&y_Test==0) # 0.41
mean(yPred2==1&y_Test==0) # 0.11
mean(yPred2==0&y_Test==1) # 0.15


##========  step 3 svm  ============##
# Final Pick -Sigmoid Kernel
Fit3 <- svm(y_Train ~ ., data = data.frame(y_Train, X_Train), 
             type = "C-classification", kernel = "sigmoid",
             cross = 5, scale = FALSE, cost = 2, gamma = 0.04, coef0 = 0.02)
yPred3 <- predict(Fit3, newdata = X_Test)
cat("Precision: ", mean(yPred3==y_Test), "\n") # 0.63
cat("Other Metrics are: ", "\n", sep=" ")
mean(yPred3==1&y_Test==1) # 0.30
mean(yPred3==0&y_Test==0) # 0.33
mean(yPred3==1&y_Test==0) # 0.19
mean(yPred3==0&y_Test==1) # 0.18

# Fit3a <- svm(y_Train ~ ., data = data.frame(y_Train, X_Train), 
#              type = "C-classification", kernel = "linear",
#              scale = FALSE, cross = 5, cost = 0.01)
# yPred3 <- predict(Fit3a, newdata = X_Test)
# cat("Precision: ", mean(yPred3==y_Test), "\n") # 0.59
# cat("Other Metrics are: ", "\n", sep=" ")
# mean(yPred3==1&y_Test==1) # 0.23
# mean(yPred3==0&y_Test==0) # 0.36
# mean(yPred3==1&y_Test==0) # 0.16
# mean(yPred3==0&y_Test==1) # 0.25
# 
# Fit3b <- svm(y_Train ~ ., data = data.frame(y_Train, X_Train),
#              type = "C-classification", kernel = "radial", cross = 5,
#              scale = FALSE, cost = 1, gamma = 0.04)
# yPred3 <- predict(Fit3b, newdata = X_Test)
# cat("Precision: ", mean(yPred3==y_Test), "\n") # 
# cat("Other Metrics are: ", "\n", sep=" ") # 0.62
# mean(yPred3==1&y_Test==1) # 0.26
# mean(yPred3==0&y_Test==0) # 0.36
# mean(yPred3==1&y_Test==0) # 0.16
# mean(yPred3==0&y_Test==1) # 0.22
# 
# Fit3c <- svm(y_Train ~ ., data = data.frame(y_Train, X_Train), 
#              type = "C-classification", kernel = "polynomial", degree = 3,
#              cross = 5, scale = FALSE, cost = 0.2, gamma = 0.02, coef0 = 0)
# yPred3 <- predict(Fit3c, newdata = X_Test)
# cat("Precision: ", mean(yPred3==y_Test), "\n") # 0.62
# cat("Other Metrics are: ", "\n", sep=" ")
# mean(yPred3==1&y_Test==1) # 0.19
# mean(yPred3==0&y_Test==0) # 0.42
# mean(yPred3==1&y_Test==0) # 0.10
# mean(yPred3==0&y_Test==1) # 0.29
# 
# Fit3d <- svm(y_Train ~ ., data = data.frame(y_Train, X_Train), 
#              type = "C-classification", kernel = "sigmoid",
#              cross = 5, scale = FALSE, cost = 2, gamma = 0.04, coef0 = 0.02)
# yPred3 <- predict(Fit3d, newdata = X_Test)
# cat("Precision: ", mean(yPred3==y_Test), "\n") # 0.63
# cat("Other Metrics are: ", "\n", sep=" ")
# mean(yPred3==1&y_Test==1) # 0.30
# mean(yPred3==0&y_Test==0) # 0.33
# mean(yPred3==1&y_Test==0) # 0.19
# mean(yPred3==0&y_Test==1) # 0.18


##========  step 4 boost  ============##
# final model -xgb with 2 splits, learning rate 0.02, 70 trees
dtrain <- xgb.DMatrix(data = as.matrix(X_Train), label = (as.integer(y_Train)-1))
dtest <- xgb.DMatrix(data = as.matrix(X_Test), label = (as.integer(y_Test)-1))
wat_ls <- list(train=dtrain, test=dtest)
Fit4 <- xgb.train(data = dtrain, max.depth = 2, eta = 0.02, nrounds = 100,
                  watchlist = wat_ls, objective = "binary:logistic")
yPred4 <- as.numeric(predict(Fit4, newdata = dtest, ntreelimit = 100) >= 0.5)
cat("Precision: ", mean(yPred4==y_Test), "\n") # 0.75
cat("Other Metrics are", "\n", sep=" ")
mean(yPred4==1&y_Test==1) # 0.34
mean(yPred4==0&y_Test==0) # 0.41
mean(yPred4==1&y_Test==0) # 0.11
mean(yPred4==0&y_Test==1) # 0.14

# Fit4a <- ada(x = X_Train, y = y_Train, loss = "logistic", type = "discrete",
#              iter = 200, nu = 0.001, bag.frac = 0.8)
# yPred4 <- predict(Fit4a, newdata = X_Test, type = "vector", n.iter = 200)
# cat("Precision: ", mean(yPred4==y_Test), "\n") # 0.73
# cat("Other Metrics are", "\n", sep=" ")
# mean(yPred4==1&y_Test==1) # 0.34
# mean(yPred4==0&y_Test==0) # 0.38
# mean(yPred4==1&y_Test==0) # 0.14
# mean(yPred4==0&y_Test==1) # 0.14
# 
# Fit4b <- ada(x = X_Train, y = y_Train, loss = "exponential", type = "discrete",
#              iter = 200, nu = 0.001, bag.frac = 0.8)
# yPred4 <- predict(Fit4b, newdata = X_Test, type = "vector", n.iter = 200)
# cat("Precision: ", mean(yPred4==y_Test), "\n") # 0.73
# cat("Other Metrics are", "\n", sep=" ")
# mean(yPred4==1&y_Test==1) # 0.34
# mean(yPred4==0&y_Test==0) # 0.38
# mean(yPred4==1&y_Test==0) # 0.14
# mean(yPred4==0&y_Test==1) # 0.14
# 
# Fit4c <- xgboost(data = as.matrix(X_Train), label = (as.integer(y_Train)-1),
#                  max.depth = 2, eta = 0.01, nrounds = 200, objective = "binary:logistic",
#                  verbose = 0)
# yPred4 <- predict(Fit4c, newdata = as.matrix(X_Test))
# yPred4 <- as.numeric(yPred4 >= 0.1)
# cat("Precision: ", mean(yPred4==y_Test), "\n") # 0.71
# cat("Other Metrics are", "\n", sep=" ")
# mean(yPred4==1&y_Test==1) # 0.33
# mean(yPred4==0&y_Test==0) # 0.38
# mean(yPred4==1&y_Test==0) # 0.14
# mean(yPred4==0&y_Test==1) # 0.15
# 
# dtrain <- xgb.DMatrix(data = as.matrix(X_Train), label = (as.integer(y_Train)-1))
# dtest <- xgb.DMatrix(data = as.matrix(X_Test), label = (as.integer(y_Test)-1))
# wat_ls <- list(train=dtrain, test=dtest)
# Fit4d <- xgb.train(data = dtrain, max.depth = 2, eta = 0.02, nrounds = 100,
#                    watchlist = wat_ls, objective = "binary:logistic")
# yPred4 <- as.numeric(predict(Fit4d, newdata = dtest, ntreelimit = 70) >= 0.5)
# cat("Precision: ", mean(yPred4==y_Test), "\n") # 0.75
# cat("Other Metrics are", "\n", sep=" ")
# mean(yPred4==1&y_Test==1) # 0.34
# mean(yPred4==0&y_Test==0) # 0.38
# mean(yPred4==1&y_Test==0) # 0.14
# mean(yPred4==0&y_Test==1) # 0.14


##========  step 5 neural network  ============##
# --nlayer1 --nodesize4
Best_Iter <- 100
for (i in 1:20) {
  Fit5t <- nnet(y_Train ~ ., data = data.frame(y_Train, X_Train),
                size = 4, maxit = 1000)
  Now_Iter <- min(Best_Iter, Fit5t$value)
  if (Now_Iter < Best_Iter) {
    Best_Iter <- Now_Iter
    Fit5 <- Fit5t
  };
};
yPred5 <- predict(Fit5, newdata = X_Test, type = "class")
cat("Precision: ", mean(yPred5==y_Test), "\n") # 0.52
cat("Other Metrics are", "\n", sep=" ")
mean(yPred5==1&y_Test==1) # 0.27
mean(yPred5==0&y_Test==0) # 0.24
mean(yPred5==1&y_Test==0) # 0.27
mean(yPred5==0&y_Test==1) # 0.21


##========  step 6 naive bayes  ============##
# final model
Fit6 <- naiveBayes(y_Train ~ ., data = data.frame(y_Train, X_Train),
                    laplace = TRUE)
yPred6 <- predict(Fit6, newdata = X_Test)
cat("Precision: ", mean(yPred6==y_Test), "\n") # 0.67
cat("Other Metrics are: ", "\n", sep=" ")
mean(yPred6==1&y_Test==1) # 0.25
mean(yPred6==0&y_Test==0) # 0.42
mean(yPred6==1&y_Test==0) # 0.10
mean(yPred6==0&y_Test==1) # 0.23