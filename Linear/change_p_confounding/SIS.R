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

###HIMA 
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



