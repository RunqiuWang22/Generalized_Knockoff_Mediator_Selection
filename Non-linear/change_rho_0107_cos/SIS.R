#X=data$X;Y=data$Y; M=data$M; COV.XM = NULL; COV.MY = COV.XM; Y.family = "gaussian";
#M.family = "gaussian";topN = NULL; parallel = FALSE; ncore = 1;scale = TRUE;verbose = FALSE

# Internal function: parallel computing check
library(foreach)
checkParallel <- function(program.name, parallel, ncore, verbose) {
  if (parallel & (ncore > 1)) {
    if (ncore > parallel::detectCores()) {
      message("You requested ", ncore, " cores. There are only ", 
              parallel::detectCores(), " in your machine!")
      ncore <- parallel::detectCores()
    }
    if (verbose) 
      message("    Running ", program.name, " with ", ncore, " cores in parallel...   (", 
              format(Sys.time(), "%X"), ")")
    doParallel::registerDoParallel(ncore)
  } else {
    if (verbose) 
      message("    Running ", program.name, " with single core...   (", 
              format(Sys.time(), "%X"), ")")
    registerDoSEQ()
  }
}

## Internal function: doOne code generater
doOneGen <- function(model.text, colind.text) {
  L <- length(eval(parse(text = colind.text)))
  script <- paste0("doOne <- function(i, datarun, Ydat){datarun$Mone <- Ydat[,i]; model <- ", 
                   model.text, ";if('try-error' %in% class(model)) b <- rep(NA, ", 
                   L, ") else { res=summary(model)$coefficients; b <- res[2,", colind.text, 
                   "]};invisible(b)}")
  return(script)
}

## Internal function: create iterator for bulk matrix by column

iblkcol_lag <- function(M, ...) {
  i <- 1
  it <- iterators::idiv(ncol(M), ...)
  
  nextEl <- function() {
    n <- iterators::nextElem(it)
    r <- seq(i, length = n)
    i <<- i + n
    M[, r, drop = FALSE]
  }
  obj <- list(nextElem = nextEl)
  class(obj) <- c("abstractiter", "iter")
  obj
}

# Internal function: Sure Independent Screening for hima
himasis <- function(Y, M, X, COV, glm.family, modelstatement, 
                    parallel, ncore, verbose, tag) {
  L.M <- ncol(M)
  M.names <- colnames(M)
  
  X <- data.frame(X)
  X <- data.frame(model.matrix(~., X))[, -1]
  
  if (is.null(COV)) {
    if (verbose) message("    No covariate is adjusted")
    datarun <- data.frame(Y = Y, Mone = NA, X = X)
    modelstatement <- modelstatement
  } else {
    COV <- data.frame(COV)
    conf.names <- colnames(COV)
    COV <- data.frame(model.matrix(~., COV))[, -1]
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    datarun <- data.frame(Y = Y, Mone = NA, X = X, COV = COV)
    modelstatement <- eval(parse(text = (paste0(modelstatement, "+", 
                                                paste0(conf.names, collapse = "+")))))
  }
  
  if(glm.family == "gaussian")
  {
    doOne <- eval(parse(text = doOneGen(paste0("try(glm(modelstatement, family = ", 
                                               glm.family, ", data = datarun))"), "c(1,4)")))
  } else if(glm.family == "negbin") {
    doOne <- eval(parse(text = doOneGen(paste0("try(MASS::glm.nb(modelstatement, data = datarun))"), "c(1,4)")))
  } else {
    stop(paste0("Screening family ", glm.family, " is not supported."))
  }
  
  
  checkParallel(tag, parallel, ncore, verbose)
  
  results <- foreach(n = iterators::idiv(L.M, chunks = ncore), 
                     M_chunk = iblkcol_lag(M, chunks = ncore), 
                     .combine = "cbind") %dopar% {sapply(seq_len(n), doOne, datarun, M_chunk)}
  
  colnames(results) <- M.names
  return(results)
}

process_var <- function(var, scale) {
  if (!is.null(var)) {
    if (scale) {
      return(scale(var))
    } else {
      return(as.matrix(var))
    }
  } else {
    return(NULL)
  }
}


SIS <- function (X, Y, M, COV.XM = NULL, COV.MY = COV.XM, Y.family = c("gaussian","binomial"), M.family = c("gaussian", "negbin"), topN = NULL, parallel = FALSE, ncore = 1, scale = TRUE, verbose = FALSE) {
  Y.family <- match.arg(Y.family)
  M.family <- match.arg(M.family)
  if (parallel & (ncore == 1)) 
    ncore <- parallel::detectCores()
  if (!parallel & (ncore > 1)) 
    parallel = TRUE
  n <- nrow(M)
  p <- ncol(M)
  if (scale==T) {
    X <- scale(X)
    M <- scale(M)
    if (!is.null(COV.XM)) 
      COV.XM <- scale(COV.XM)
    if (!is.null(COV.MY)) 
      COV.MY <- scale(COV.MY)
    if (verbose) 
      message("Data scaling is completed.")
  }
  if (scale==F) {
    X <- as.matrix(X)
    M <- as.matrix(M)
    if (!is.null(COV.XM)) 
      COV.XM <- as.matrix(COV.XM)
    if (!is.null(COV.MY)) 
      COV.MY <- as.matrix(COV.MY)
  }
  if (is.null(topN)) {
    if (Y.family == "binomial") 
      d <- ceiling(n/(2 * log(n)))
    else
      d <- ceiling(2 * n/log(n))
  } else {
    d <- topN
  }
  d <- min(p, d)
  message("Step 1: Sure Independent Screening ...", "     (", 
          format(Sys.time(), "%X"), ")")
  if (Y.family == "binomial") {
    if (verbose) 
      message("    Screening M using the association between X (independent variable) and M (dependent variable): ", 
              appendLF = FALSE)
    alpha = SIS_Results <- himasis(NA, M, X, COV.XM, glm.family = M.family, 
                                   modelstatement = "Mone ~ X", parallel = parallel, 
                                   ncore = ncore, verbose, tag = paste0("Sure Independent Screening (M ~ X + COV.XM, family: ", 
                                                                        M.family, ")"))
    SIS_Pvalue <- SIS_Results[2, ]
  }
  else if (Y.family == "gaussian") {
    if (verbose) 
      message("    Screening M using the association between M (independent variable) and Y (dependent variable): ", 
              appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY, glm.family = Y.family, 
                           modelstatement = "Y ~ Mone + X", parallel = parallel, 
                           ncore = ncore, verbose, tag = paste0("Sure Independent Screening (Y ~ M + X + COV.MY, family: ", 
                                                                Y.family, ")"))
    SIS_Pvalue <- SIS_Results[2, ]
  }
  else {
    stop(paste0("Family ", Y.family, " is not supported."))
  }
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d])
  if (verbose) 
    message("    Top ", length(ID), " mediators are selected: ", 
            paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ", "))
  M_SIS <- M[, ID]
  M_ID_name <- colnames(M)
  XM <- cbind(M_SIS, X)

  return(as.numeric(ID))
}

###update HIMA 
classicHIMA1 <- function(X, M, Y, COV.XM = NULL, COV.MY = COV.XM,
                        Y.type = c("continuous", "binary"),
                        M.type = c("gaussian", "negbin"),
                        penalty = c("MCP", "SCAD", "lasso"),
                        topN = NULL,
                        parallel = FALSE,
                        ncore = 1,
                        scale = TRUE,
                        Bonfcut = 0.05,
                        verbose = FALSE,
                        ...) {
  Y.type <- match.arg(Y.type)
  Y.family <- switch(Y.type,
                     continuous = "gaussian",
                     binary = "binomial",
                     stop("Invalid Y.type. Expected 'continuous' or 'binary'.")
  )
  
  M.type <- match.arg(M.type)
  penalty <- match.arg(penalty)
  
  if (Y.family == "gaussian") message("Running linear HIMA with ", penalty, " penalty...")
  if (Y.family == "binomial") message("Running logistic HIMA with ", penalty, " penalty...")
  
  if (parallel && (ncore == 1)) ncore <- parallel::detectCores()
  if (!parallel && (ncore > 1)) parallel <- TRUE
  
  n <- nrow(M)
  p <- ncol(M)
  
  # Process required variables
  X <- process_var(X, scale)
  M <- process_var(M, scale)
  
  # Process optional covariates
  COV.XM <- process_var(COV.XM, scale)
  COV.MY <- process_var(COV.MY, scale)
  
  if (scale && verbose) message("Data scaling is completed.")
  
  if (is.null(topN)) {
    if (Y.type == "binary") d <- ceiling(n / (2 * log(n))) else d <- ceiling(2 * n / log(n))
  } else {
    d <- topN # the number of top mediators that associated with exposure (X)
  }
  
  d <- min(p, d) # if d > p select all mediators
  
  #########################################################################
  ################################ STEP 1 #################################
  #########################################################################
  message("Step 1: Sure Independent Screening ...", "     (", format(Sys.time(), "%X"), ")")
  
  if (Y.type == "binary") {
    # Screen M using X given the limited information provided by Y (binary)
    if (verbose) message("    Screening M using the association between X (independent variable) and M (dependent variable): ", appendLF = FALSE)
    alpha <- SIS_Results <- himasis(NA, M, X, COV.XM,
                                    glm.family = M.type, modelstatement = "Mone ~ X",
                                    parallel = parallel, ncore = ncore, verbose, tag = paste0("Sure Independent Screening (M ~ X + COV.XM, family: ", M.type, ")")
    )
    SIS_Pvalue <- SIS_Results[2, ]
  } else if (Y.type == "continuous") {
    # Screen M using Y (continuous)
    if (verbose) message("    Screening M using the association between M (independent variable) and Y (dependent variable): ", appendLF = FALSE)
    SIS_Results <- himasis(Y, M, X, COV.MY,
                           glm.family = "gaussian", modelstatement = "Y ~ Mone + X",
                           parallel = parallel, ncore = ncore, verbose, tag = paste0("Sure Independent Screening (Y ~ M + X + COV.MY, family: ", Y.type, ")")
    )
    SIS_Pvalue <- SIS_Results[2, ]
  } else {
    stop(paste0("A ", Y.type, " data type for outcome is not supported."))
  }
  
  # Note: ranking using p on un-standardized data is equivalent to ranking using beta on standardized data
  SIS_Pvalue_sort <- sort(SIS_Pvalue)
  ID <- which(SIS_Pvalue <= SIS_Pvalue_sort[d]) # the index of top mediators
  if (verbose) message("    Top ", length(ID), " mediators are selected: ", paste0(names(SIS_Pvalue_sort[seq_len(d)]), collapse = ", "))
  
  M_SIS <- M[, ID]
  M_ID_name <- colnames(M)
  XM <- cbind(M_SIS, X)
  
  #########################################################################
  ################################ STEP 2 #################################
  #########################################################################
  message("Step 2: Penalized estimate (", penalty, ") ...", "     (", format(Sys.time(), "%X"), ")")
  
  ## Based on the screening results in step 1. We will find the most influential M on Y.
  if (is.null(COV.MY)) {
    fit <- ncvreg(XM, Y,
                  family = Y.family,
                  penalty = penalty,
                  penalty.factor = c(rep(1, ncol(M_SIS)), 0), ...
    )
  } else {
    COV.MY <- data.frame(COV.MY)
    COV.MY <- data.frame(model.matrix(~., COV.MY))[, -1]
    conf.names <- colnames(COV.MY)
    if (verbose) message("    Adjusting for covariate(s): ", paste0(conf.names, collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
    fit <- ncvreg(XM_COV, Y,
                  family = Y.family,
                  penalty = penalty,
                  penalty.factor = c(rep(1, ncol(M_SIS)), rep(0, 1 + ncol(as.matrix(COV.MY)))), ...
    )
  }
  lam <- fit$lambda[which.min(BIC(fit))]
  if (verbose) message("    Tuning parameter lambda selected: ", lam)
  Coefficients <- coef(fit, lambda = lam)
  est <- Coefficients[2:(d + 1)]
  ID_1_non <- which(est != 0)
  if (length(ID_1_non) == 0) {
    if (verbose) message("    All ", penalty, " beta estimates of the ", length(ID), " mediators are zero.")
    results <- NULL
    return(results)
  } else {
    if (verbose) message("    Non-zero ", penalty, " beta estimate(s) of mediator(s) found: ", paste0(names(ID_1_non), collapse = ","))
    beta_est <- est[ID_1_non] # The non-zero MCP estimators of beta
    ID_test <- ID[ID_1_non] # The index of the ID of non-zero beta in Y ~ M
    ##
    
    if (Y.type == "binary") {
      ## This has been done in step 1 (when Y is binary, alpha is estimated in M ~ X)
      alpha <- alpha[, ID_test, drop = FALSE]
      message("    Using alpha estimated in Step 1 ...   (", format(Sys.time(), "%X"), ")")
    } else if (Y.type == "continuous") {
      if (verbose) message("    Estimating alpha (effect of X on M): ", appendLF = FALSE)
      alpha <- himasis(NA, M[, ID_test, drop = FALSE], X, COV.XM,
                       glm.family = M.type,
                       modelstatement = "Mone ~ X", parallel = FALSE, ncore = ncore,
                       verbose, tag = paste0("site-by-site alpha estimation (M ~ X + COV.XM, family: ", M.type, ")")
      )
    } else {
      stop(paste0("A ", Y.type, " data type for outcome is not supported."))
    }
    
    #########################################################################
    ################################ STEP 3 #################################
    #########################################################################
    message("Step 3: Joint significance test ...", "     (", format(Sys.time(), "%X"), ")")
    
    alpha_est_ID_test <- as.numeric(alpha[1, ]) #  the estimator for alpha
    P_alpha <- alpha[2, ] # the raw p-value for alpha
    alpha_est <- alpha_est_ID_test
    
    ## Post-test based on the oracle property of the MCP penalty
    if (is.null(COV.MY)) {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X)
    } else {
      YMX <- data.frame(Y = Y, M[, ID_test, drop = FALSE], X = X, COV.MY)
    }
    
    if (Y.type == "continuous") res <- summary(glm(Y ~ ., family = "gaussian", data = YMX))$coefficients
    if (Y.type == "binary") res <- summary(glm(Y ~ ., family = "binomial", data = YMX))$coefficients
    
    est <- res[2:(length(ID_test) + 1), 1] # the estimator for beta
    P_beta <- res[2:(length(ID_test) + 1), 4] # the raw p-value for beta
    
    IDE <- alpha_est * beta_est # mediation(indirect) effect
    
    ## Use the maximum value as p value
    Pmax <- apply(cbind(P_alpha, P_beta), 1, max)
    
    ## Bonferroni
    Pmax_Bonf <- Pmax * length(ID_test)
    sig_ind <- which(Pmax_Bonf < Bonfcut)
    
    # FDRA <- rbind(P_fdr_beta, P_fdr_alpha)
    # FDR <- apply(FDRA, 2, max)
    
    # Total effect
    # if(is.null(COV.MY)) {
    #   YX <- data.frame(Y = Y, X = X)
    # } else {
    #   YX <- data.frame(Y = Y, X = X, COV.MY)
    # }
    #
    # gamma_est <- coef(glm(Y ~ ., family = Y.type, data = YX))[2]
    
    results <- data.frame(
      Index = M_ID_name[ID_test][sig_ind],
      alpha_hat = alpha_est[sig_ind],
      beta_hat = beta_est[sig_ind],
      IDE = IDE[sig_ind],
      rimp = (abs(IDE) / sum(abs(IDE)))[sig_ind] * 100,
      pmax = Pmax[sig_ind], check.names = FALSE
    )
    
    message("Done!", "     (", format(Sys.time(), "%X"), ")")
    
    doParallel::stopImplicitCluster()
    
    return(results)
  }
}

.JCCorrect <- function(pval) {
  z <- stats::qnorm(pval, lower.tail = F)
  res <- .nullParaEst(z)
  pval.JC <- stats::pnorm(z, mean = res$mu, sd = res$s, lower.tail = F)
  return(pval.JC)
}


DACT <- function(p_a, p_b) {
  Z_a <- stats::qnorm(p_a, lower.tail = F)
  Z_b <- stats::qnorm(p_b, lower.tail = F)
  pi0a <- 1 - .nonnullPropEst(Z_a, 0, 1)
  pi0b <- 1 - .nonnullPropEst(Z_b, 0, 1)
  if (pi0a > 1) {
    pi0a <- 1
  }
  if (pi0b > 1) {
    pi0b <- 1
  }
  p.mat <- cbind(p_a, p_b)
  p3 <- (apply(p.mat, 1, max))^2
  wg1 <- pi0a * (1 - pi0b)
  wg2 <- (1 - pi0a) * pi0b
  wg3 <- pi0a * pi0b
  wg.sum <- wg1 + wg2 + wg3
  wg.std <- c(wg1, wg2, wg3) / wg.sum
  p_dact <- wg.std[1] * p_a + wg.std[2] * p_b + wg.std[3] * p3
  p_dact <- .JCCorrect(p_dact)
  return(p_dact)
}


.nonnullPropEst <- function(x, u, sigma) {
  z <- (x - u) / sigma
  xi <- c(0:100) / 100
  tmax <- sqrt(log(length(x)))
  tt <- seq(0, tmax, 0.1)
  
  epsest <- NULL
  
  for (j in seq_along(tt)) {
    t <- tt[j]
    f <- t * xi
    f <- exp(f^2 / 2)
    w <- (1 - abs(xi))
    co <- 0 * xi
    
    for (i in 1:101) {
      co[i] <- mean(cos(t * xi[i] * z))
    }
    epshat <- 1 - sum(w * f * co) / sum(w)
    epsest <- c(epsest, epshat)
  }
  return(epsest = max(epsest))
}

.nullParaEst <- function(x, gamma = 0.1) {
  n <- length(x)
  t <- c(1:1000) / 200
  
  gan <- n^(-gamma)
  that <- 0
  shat <- 0
  uhat <- 0
  epshat <- 0
  
  phiplus <- rep(1, 1000)
  phiminus <- rep(1, 1000)
  dphiplus <- rep(1, 1000)
  dphiminus <- rep(1, 1000)
  phi <- rep(1, 1000)
  dphi <- rep(1, 1000)
  
  for (i in 1:1000) {
    s <- t[i]
    phiplus[i] <- mean(cos(s * x))
    phiminus[i] <- mean(sin(s * x))
    dphiplus[i] <- -mean(x * sin(s * x))
    dphiminus[i] <- mean(x * cos(s * x))
    phi[i] <- sqrt(phiplus[i]^2 + phiminus[i]^2)
  }
  
  ind <- min(c(1:1000)[(phi - gan) <= 0])
  tt <- t[ind]
  a <- phiplus[ind]
  b <- phiminus[ind]
  da <- dphiplus[ind]
  db <- dphiminus[ind]
  c <- phi[ind]
  
  that <- tt
  shat <- -(a * da + b * db) / (tt * c * c)
  shat <- sqrt(shat)
  uhat <- -(da * b - db * a) / (c * c)
  epshat <- 1 - c * exp((tt * shat)^2 / 2)
  
  return(musigma = list(mu = uhat, s = shat))
}

###update eHIMA
eHIMA1 <- function(X, M, Y, COV = NULL,
                   topN = NULL,
                   scale = TRUE,
                   FDRcut = 0.05,
                   verbose = FALSE) {
  n <- nrow(M)
  p <- ncol(M)
  
  # Process required variables
  X <- process_var(X, scale)
  M <- process_var(M, scale)
  
  # Process optional covariates
  COV <- process_var(COV, scale)
  
  if (scale && verbose) message("Data scaling is completed.")
  
  if (is.null(COV)) {
    MZX <- cbind(M, X)
    XZ <- X
    q <- 0
  } else {
    MZX <- cbind(M, COV, X)
    XZ <- cbind(X, COV)
    q <- ncol(COV)
  }
  
  if (is.null(topN)) d <- ceiling(2 * n / log(n)) else d <- topN # the number of top mediators that associated with exposure (X)
  d <- min(p, d) # if d > p select all mediators
  
  M_ID_name <- colnames(M)
  if (is.null(M_ID_name)) M_ID_name <- seq_len(p)
  
  #------------- Step 1:  mediator screening ---------------------------
  message("Step 1: Sure Independent Screening + minimax concave penalty (MCP) ...", "     (", format(Sys.time(), "%X"), ")")
  
  beta_SIS <- matrix(0, 1, p)
  for (i in 1:p) {
    ID_S <- c(i, (p + 1):(p + q + 1))
    MZX_SIS <- MZX[, ID_S]
    fit <- lsfit(MZX_SIS, Y, intercept = TRUE)
    beta_SIS[i] <- fit$coefficients[2]
  }
  
  ## est_a for SIS #########
  alpha_SIS <- matrix(0, 1, p)
  for (i in 1:p) {
    fit_a <- lsfit(XZ, M[, i], intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    alpha_SIS[i] <- est_a
  }
  ab_SIS <- alpha_SIS * beta_SIS
  ID_SIS <- which(-abs(ab_SIS) <= sort(-abs(ab_SIS))[d]) # \Omega_1
  d <- length(ID_SIS)
  
  if (verbose) message("        Top ", d, " mediators are selected: ", paste0(M_ID_name[ID_SIS], collapse = ", "))
  
  if (verbose) {
    if (is.null(COV)) {
      message("        No covariate was adjusted.")
    } else {
      message("        Adjusting for covariate(s): ", paste0(colnames(COV), collapse = ", "))
    }
  }
  
  MZX_SIS <- MZX[, c(ID_SIS, (p + 1):(p + q + 1))] # select m_i in \Omega_1 from M
  fit <- ncvreg(MZX_SIS, Y, family = "gaussian", penalty = "MCP")
  lam <- fit$lambda[which.min(BIC(fit))]
  beta_penalty <- coef(fit, lambda = lam)[2:(d + 1)]
  id_non <- ID_SIS[which(beta_penalty != 0)] # the ID of non-zero
  
  #----------- Step 2: Refitted partial regression ----------------------
  message("Step 2: Refitted partial regression ...", "     (", format(Sys.time(), "%X"), ")")
  ## beta_est ########
  MZX_penalty <- MZX[, c(id_non, (p + 1):(p + q + 1))]
  fit <- lsfit(MZX_penalty, Y, intercept = TRUE)
  if (length(c(id_non, (p + 1):(p + q + 1)))==1){
    beta_est_cox <- fit$coefficients[2]
    beta_SE_cox <- ls.diag(fit)$std.err[2]
  } else {
    beta_est_cox <- fit$coefficients[2:(dim(MZX_penalty)[2] + 1)]
    beta_SE_cox <- ls.diag(fit)$std.err[2:(dim(MZX_penalty)[2] + 1)]
  }
  
  # Computes basic statistics, including standard errors, t- and p-values for the regression coefficients.
  beta_est <- fit$coefficients[2:(length(id_non) + 1)] # estimated beta != 0
  beta_SE <- ls.diag(fit)$std.err[2:(length(id_non) + 1)]
  P_beta_penalty <- 2 * (1 - pnorm(abs(beta_est_cox[seq_along(id_non)]) / beta_SE_cox[seq_along(id_non)], 0, 1))
  
  P_oracle_beta <- matrix(0, 1, p) # an empty vector
  beta_est_orc <- matrix(0, 1, p)
  beta_SE_orc <- matrix(0, 1, p)
  
  j <- 1
  for (i in 1:p) {
    if (i %in% id_non) {
      P_oracle_beta[i] <- P_beta_penalty[j]
      beta_est_orc[i] <- beta_est[j]
      beta_SE_orc[i] <- beta_SE[j]
      j <- j + 1
    } else {
      MZX_ora <- MZX[, c(id_non, i, (p + 1):(p + q + 1))]
      fit_ora <- lsfit(MZX_ora, Y, intercept = TRUE)
      beta_est_cox <- fit_ora$coefficients[2:(dim(MZX_ora)[2] + 1)]
      beta_est_orc[i] <- beta_est_cox[length(id_non) + 1]
      beta_SE_cox <- ls.diag(fit_ora)$std.err[2:(dim(MZX_ora)[2] + 1)]
      
      beta_SE_orc[i] <- beta_SE_cox[length(id_non) + 1]
      P_oracle_beta[i] <- 2 * (1 - pnorm(abs(beta_est_cox) / beta_SE_cox, 0, 1))[length(id_non) + 1]
    }
  }
  
  ## ----- P_oracle_alpha ----------------- ##
  alpha_est_penalty <- matrix(0, 1, length(id_non))
  alpha_SE_penalty <- matrix(0, 1, length(id_non))
  P_alpha_penalty <- matrix(0, 1, length(id_non))
  for (i in seq_along(id_non)) {
    fit_a <- lsfit(XZ, M[, id_non[i]], intercept = TRUE)
    est_a <- matrix(coef(fit_a))[2]
    se_a <- ls.diag(fit_a)$std.err[2]
    sd_1 <- abs(est_a) / se_a
    P_alpha_penalty[i] <- 2 * (1 - pnorm(sd_1, 0, 1)) ## the SIS for alpha
    alpha_est_penalty[i] <- est_a
    alpha_SE_penalty[i] <- se_a
  }
  
  ### P_oracle_alpha #########
  P_oracle_alpha <- matrix(0, 1, p) # an empty vector
  alpha_est_orc <- matrix(0, 1, p)
  alpha_SE_orc <- matrix(0, 1, p)
  
  j <- 1
  for (i in 1:p) {
    if (i %in% id_non) {
      P_oracle_alpha[i] <- P_alpha_penalty[j]
      alpha_est_orc[i] <- alpha_est_penalty[j]
      alpha_SE_orc[i] <- alpha_SE_penalty[j]
      j <- j + 1
    } else {
      fit_a_ora <- lsfit(XZ, M[, i], intercept = TRUE)
      est_a_ora <- matrix(coef(fit_a_ora))[2]
      se_a_ora <- ls.diag(fit_a_ora)$std.err[2]
      sd_1_ora <- abs(est_a_ora) / se_a_ora
      P_alpha_penalty_ora <- 2 * (1 - pnorm(sd_1_ora, 0, 1))
      P_oracle_alpha[i] <- P_alpha_penalty_ora
      alpha_est_orc[i] <- est_a_ora
      alpha_SE_orc[i] <- se_a_ora
    }
  }
  
  #---------- Step 3: DACT  -------------------------
  message("Step 3: Divide-aggregate composite-null test (DACT) ...", "     (", format(Sys.time(), "%X"), ")")
  
  # Mediator selection
  P_oracle_alpha[P_oracle_alpha == 0] <- 10^(-17)
  P_oracle_beta[P_oracle_beta == 0] <- 10^(-17)
  P_BH <- (1:p) * (FDRcut / p)
  
  ## DACT
  DACT_ora <- DACT(p_a = t(P_oracle_alpha), p_b = t(P_oracle_beta))
  P_sort_DACT <- sort(DACT_ora)
  SN <- sum(as.numeric(P_sort_DACT <= P_BH))
  ID_BH_DACT <- order(DACT_ora)[1:SN]
  
  # # Total effect
  # if(is.null(COV)) {
  #   YX <- data.frame(Y = Y, X = X)
  # } else {
  #   YX <- data.frame(Y = Y, X = X, COV)
  # }
  #
  # gamma_est <- coef(glm(Y ~ ., family = Y.family, data = YX))[2]
  
  IDE <- alpha_est_orc[ID_BH_DACT] * beta_est_orc[ID_BH_DACT]
  
  if (length(ID_BH_DACT) > 0) {
    out_result <- data.frame(
      Index = M_ID_name[ID_BH_DACT],
      alpha_hat = alpha_est_orc[ID_BH_DACT],
      alpha_se = alpha_SE_orc[ID_BH_DACT],
      beta_hat = beta_est_orc[ID_BH_DACT],
      beta_se = beta_SE_orc[ID_BH_DACT],
      IDE = IDE,
      rimp = abs(IDE) / sum(abs(IDE)) * 100,
      pmax = DACT_ora[ID_BH_DACT], check.names = FALSE
    )
    if (verbose) message(paste0("        ", length(ID_BH_DACT), " significant mediator(s) identified."))
  } else {
    if (verbose) message("        No significant mediator identified.")
    out_result <- NULL
  }
  
  message("Done!", "     (", format(Sys.time(), "%X"), ")")
  return(out_result)
}


