mymethod <- function(X,Y,M,V=NULL,q,scale=T){
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
  
  MM <- cbind(M,MK)
  
  ### step 2: patha 
  Z1=Z1tilde=rep(NA,p)
  
  for (j in 1:p){
    fit1 = lm(X~MM[,c(j,j+p)])
    Z1[j] = abs(coef(fit1))[2]
    Z1tilde[j] = abs(coef(fit1))[3]
  }
  
  
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
  
  ### combine
  W = (Z1-Z1tilde)*(Z2-Z2tilde)
  mythred=knockoff.threshold(W,fdr=q,offset=1)
  myselect=which(W>=mythred)
  return(myselect)
}
