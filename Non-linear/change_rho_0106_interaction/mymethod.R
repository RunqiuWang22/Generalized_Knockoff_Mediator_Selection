library(randomForest)
library(caret)

### cross-validation for random forest to select the best tuning parameter
CV.randomforest <- function(data,mtry_values,ntree_values,cv_folds){
  Y=data[,1]
  # Storage for cross-validation results
  cv_results <- expand.grid(mtry = mtry_values, ntree = ntree_values)
  cv_results$MSE <- NA
  
  # Cross-validation loop
  set.seed(123)
  folds <- sample(rep(1:cv_folds, length.out = nrow(data)))
  
  for (i in seq_len(nrow(cv_results))) {
    fold_mse <- numeric(cv_folds)
    
    for (fold in 1:cv_folds) {
      # Split data into training and testing sets
      train_data <- data[folds != fold, ]
      test_data <- data[folds == fold, ]
      
      # Fit random forest model
      model <- randomForest(
        Y ~ ., 
        data = train_data, 
        mtry = cv_results$mtry[i],
        ntree = cv_results$ntree[i],
        importance = TRUE
      )
      
      # Predict on test set
      preds <- predict(model, newdata = test_data)
      
      # Compute MSE for the fold
      fold_mse[fold] <- mean((preds - test_data$Y)^2)
    }
    
    # Store mean MSE for the parameter combination
    cv_results$MSE[i] <- mean(fold_mse)
  }
  
  # Find the best parameter combination
  best_params <- cv_results[which.min(cv_results$MSE), ]
  return(best_params)
}

#X=data$X;Y=data$Y;M=data$M;V=NULL;q=0.2;scale=T;type="non-linear"
mymethod <- function(X,Y,M,V=NULL,q,scale=T,type){
  n <- nrow(M)
  p <- ncol(M)
  
  if (scale==T){
  X <- scale(X)
  M <- scale(M)
  Y <- scale(Y)
  }
  
  if (!is.null(V)) {
    if (scale==T){
    V = scale(V)
    }
    rY = lm(Y~V)$residuals
    rX = lm(X~V)$residuals
    rM = matrix(NA,nrow=nrow(M),ncol=ncol(M))
    
    for (j in 1:p){ 
      rM[,j] = lm(M[,j]~V)$residuals
    }
    
  Y = rY
  X = rX
  M = rM
  }
  ### step 1: generate knockoff  
  MK = matrix(0, nrow=n, ncol=p)
  for (j in 1:p){
    mu <- mean(M[,j])         
    sigma2 <- var(M[,j])      
    MK[,j] <- rnorm(length(M[,j]), mean = mu, sd = sqrt(sigma2))
  }
  
  MM <- as.matrix(cbind(M,MK))
  
  ### step 2: patha 
  Z1=Z1tilde=rep(NA,p)
  
  for (j in 1:p){
    fit1 = lm(X~MM[,c(j,j+p)])
    Z1[j] = abs(coef(fit1))[2]
    Z1tilde[j] = abs(coef(fit1))[3]
  }
  
  if (type=="linear"){
  ### step 3: pathb
  fit2 = lm(Y~X)
  res1 = fit2$residuals
  
  ###get res2
  res2 = res2_tilde = matrix(NA,n,p)
  for (j in 1:p){ 
    res2[,j] = lm(MM[,j]~X)$residuals
  }
  
  ### get res2_tilde
  mu_res2 <- apply(res2,2,mean)
  Sig_res2 <- cov(res2)
  res2_tilde <- knockoff::create.gaussian(res2, mu_res2, Sig_res2)
  rres2 = cbind(res2,res2_tilde)
  
  penalty = rep(0,ncol(rres2))
  penalty[1:(2*p)] = 1
  fit3 = cv.glmnet(rres2, res1, family = "gaussian", standardize = FALSE, intercept = FALSE, penalty.factor = penalty)
  Z2 = abs(coef(fit3,s = "lambda.min")[1 + (1:p)])
  Z2tilde = abs(coef(fit3,s = "lambda.min")[1 + p + (1:p)])
}
  
  if (type=="non-linear"){
    ### step 3: pathb
    ###generate knockoff  
    MX=as.matrix(cbind(M,X))
    mu_MX <- apply(MX,2,mean)
    Sig_MX <- cov(MX)
    MX_tilde <- knockoff::create.gaussian(MX, mu_MX, Sig_MX)
    M_tilde <- MX_tilde[,-ncol(MX_tilde)] #remove the Xtilde
    fdata=as.data.frame(cbind(Y,X,M,M_tilde))
    names(fdata) = c("Y","X",paste("M",1:p,sep=""),paste("M_tilde",1:p,sep=""))
    pp=round(sqrt(ncol(fdata)-1))
    cvfit3 = CV.randomforest(fdata,mtry_values=c((pp-10),pp,(pp+10)),ntree_values=c(500,1000),cv_folds=5)###cross-validations
    fit3 = randomForest(Y ~.,importance=TRUE,data=fdata,family="gaussian",mtry=cvfit3$mtry,ntree=cvfit3$ntree) 
    tmp=fit3$importance
    #varImpPlot(fit3)
    Z2 = tmp[(1:p)+1,1] 
    Z2tilde = tmp[(1:p)+1+p,1]
  }
  ### combine
  W = (Z1-Z1tilde)*(round(Z2,3)-round(Z2tilde,3))
  mythred=knockoff.threshold(W,fdr=q,offset=1)
  myselect=which(W>=mythred)
  return(myselect)
}

#(length(myselect)- sum(myselect %in% pathi))/length(myselect)
#sum(myselect %in% pathi)/length(pathi)
