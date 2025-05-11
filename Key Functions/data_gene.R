### data without confounding
datagen <- function(n, p, a, b,ss,rho,std,patha,pathb,covmat){
  pathi = intersect(patha, pathb)
  alpha = rep(0,p); alpha[patha] = a; alpha[pathi] = a*ss;
  beta = rep(0,p); beta[pathb] = b; beta[pathi] = b*ss;
  gamma=1
  X=rnorm(n,0,1)
  M = X %o% alpha + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  Y = X*gamma + M %*% beta + rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  
  list(X = X, Y = Y, M = M)
}

### data with confounding
datagen.conf <- function(n, p, a, b,ss,rho,std,patha,pathb,covmat){
  pathi = intersect(patha, pathb)
  XV = mvtnorm::rmvnorm(n, mean = rep(0,2), sigma = matrix(c(1,0.5,0.5,1),2,2)) # sigma is covariance matrix
  X = as.matrix(XV[,1])
  V = as.matrix(XV[,2])
  alpha = rep(0,p); alpha[patha] = a; alpha[pathi] = a*ss;
  eta1 = rep(0,p); eta1[p] = 1 ### need to consider
  M = X %*% alpha + V %*% eta1 + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  
  beta = rep(0,p); beta[pathb] = b; beta[pathi] = b*ss;
  gamma = 1
  eta2 = 1
  Y = X*gamma + M %*% beta + eta2*V+rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  
  list(X = X, Y = Y, M = M, V = V)
}

datagen.interaction <- function(n, p, a, b,ss,rho,std,patha,pathb,covmat){
  pathi = intersect(patha, pathb)
  alpha = rep(0,p); alpha[patha] = a/2; alpha[pathi] = a;
  X=rnorm(n,0,1)
  M = X %o% alpha + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  interactions = cbind(M[,10]*M[,11] , M[,12]*M[,13] , M[,14]*M[,15] , M[,16]*M[,17] , M[,18]*M[,19] , M[,20]*M[,21] ,
                       M[,22]*M[,23] , M[,24]*M[,25] , M[,26]*M[,27] , M[,28]*M[,29] , M[,30]*M[,31] , M[,32]*M[,33] ,
                       M[,34]*M[,35] , M[,36]*M[,37] , M[,38]*M[,39])
  beta = rep(b,ncol(interactions))
  gamma=0.1
  Y = X*gamma +interactions %*% beta + rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  list(X = X, Y = Y, M = M)
}

datagen.cos <- function(n, p, a, b,ss,rho,std,patha,pathb,covmat){
  pathi = intersect(patha, pathb)
  alpha = rep(0,p); alpha[patha] = a; alpha[pathi] = a*ss;
  beta = rep(0,p); beta[pathb] = b*ss; beta[pathi] = b;
  gamma=0
  X=rnorm(n,0,1)
  M = X %o% alpha + mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = covmat) ## sigma: covariance matrix
  M_cos <- apply(M, 2, cos)
  Y = X*gamma + M_cos %*% beta + rnorm(n, mean = 0, sd = std)
  colnames(M) = paste("M", as.character(1:p), sep = "")
  
  list(X = X, Y = Y, M = M)
}