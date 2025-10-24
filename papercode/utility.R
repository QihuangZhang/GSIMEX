# Function to perform spectral multiplication on a matrix. 
# It calculates eigenvalues and eigenvectors and then multiplies 
# the matrix of eigenvectors with the matrix of eigenvalues raised to the power zeta.
Spectral_multi <- function(MisMatrix, zeta, multi = FALSE){
  r <- eigen(MisMatrix)
  new_matrix <- r$vectors %*% (diag(r$values))^zeta %*% solve(r$vectors)
  if (multi) {
    new_matrix[new_matrix < 0 ] <- 0
  }
  return(new_matrix)
}

# Function to perform a Monte Carlo operation on the input vector (x) using the MisMatrix.
# It applies the Spectral_multi function to MisMatrix and generates a multinomial distribution 
# for each element in x based on the rows of MisMatrix.
MCoperator <- function(MisMatrix, zeta, x, multi = F, start0 = TRUE){
  MisMatrix_current <- Spectral_multi(MisMatrix, zeta, multi)
  if (start0) {x <- x + 1}
  xnew <- sapply (x, FUN = function(xi) {
    return(which( rmultinom( 1, 1, prob = MisMatrix_current[xi,]) == 1))
  })
  if (start0) {xnew <- xnew - 1}
  return(xnew)
}

# Function to perform cross-validation and calculate the Mean Squared Error (MSE).
# It splits the data into training and test sets, calculates the regression model, 
# makes predictions, and calculates the MSE.
get_MSE_CV <- function(B, ZZ_sub, est, nfold){
  folds <- cut(seq(1,B), breaks=nfold, labels=FALSE)
  ZZ_sub <- data.frame(ZZ_sub)
  est    <- data.frame(est)
  MSE_CV <- NULL
  
  for (foldtest in 1:(nfold-1)){
    testdata_X <- ZZ_sub[folds == foldtest,]  
    traindata_X <- ZZ_sub[folds %in% (foldtest+1):nfold,] 
    test_est <- est[folds == foldtest,]  
    train_est <- est[folds %in% (foldtest+1):nfold,]
    
    traindata_X <- as.matrix(traindata_X)
    train_est   <- as.matrix(train_est)
    
    beta_model <- as.vector(solve(t(traindata_X) %*% traindata_X, tol=1e-200) %*% t(traindata_X) %*% train_est)
    
    est_hat_i <-  t(t(testdata_X)) %*%  beta_model
    
    ## Calculating SSR
    MSE_CV <- c(MSE_CV, mean((test_est - est_hat_i)^2) )
  }
  
  MSE_CV_mean <- mean(MSE_CV)
  
  return(MSE_CV_mean)
}

# Function to perform cross-validation and calculate the Mean Squared Error (MSE).
# It splits the data into training and test sets, calculates the regression model, 
# makes predictions, and calculates the MSE.
get_MSE_CV_step <- function(B, ZZ_sub, est, nfold, stepsize, mfold = nfold %/% 5 ){
  folds <- cut(seq(1,B),breaks=nfold,labels=FALSE)
  ZZ_sub <- data.frame(ZZ_sub)
  est    <- data.frame(est)
  MSE_CV <- NULL
  
  for (foldtest in 1: mfold){
    testdata_X <- ZZ_sub[folds %in% 1:foldtest,]  
    traindata_X <- ZZ_sub[folds %in% (foldtest+stepsize):nfold,] 
    test_est <- est[folds %in% 1:foldtest,]  
    train_est <- est[folds %in% (foldtest+stepsize):nfold,]
    
    traindata_X <- as.matrix(traindata_X)
    train_est   <- as.matrix(train_est)
    
    beta_model <- as.vector(solve(t(traindata_X) %*% traindata_X, tol=1e-200) %*% t(traindata_X) %*% train_est)

    est_hat_i <-  t(t(testdata_X)) %*%  beta_model
    
    ## Calculating SSR
    MSE_CV <- c(MSE_CV, mean((test_est - est_hat_i)^2) )
  }
  
  MSE_CV_mean <- mean(MSE_CV)
  
  return(MSE_CV_mean)
}


Summarize_Results <- function(results, B, extrapolation = "all"){
  if (!extrapolation %in% c("linear","quadratic", "best subset", "all")){
    print("Extrapolation argument should be one of linear, quadratic, best subset, or all.")}
  
  results <- results %>%
    mutate(AIC = B*log(SSR/B) + 2*(as.numeric(modelm)+2) ) %>%
    mutate(BIC = B*log(SSR/B) + log(B)*(as.numeric(modelm)+2) )
    
  ## linear extrapolation
  results_linear <- results %>%
    filter(model_id == 2) %>%
    mutate(bestmodel= NA, bestorder = NA) %>% ## Keep consistent with the results of best subset model
    mutate(extrapo="linear")
  
  ## quadratic extrapolation
  results_quadra <- results %>%
    filter(model_id == 4) %>%
    mutate(bestmodel= NA, bestorder = NA) %>% ## Keep consistent with the results of best subset model
    mutate(extrapo="quadratic")
  
  
  ## Best subset model
  results_bestsubset <- results %>%
    mutate(modelm = as.character(modelm)) %>%
    mutate(order = as.character(order)) %>%
    group_by(para_id, modelm) %>%
    filter(SSR == min(SSR)) %>%
    mutate(bestmodel = model_id )  %>%
    mutate(model_id = paste("best of m=", modelm)) %>%
    mutate(bestorder = NA) %>%
    mutate(extrapo=paste0("best subset of m=", modelm)) 
  
  ## Optimal CV Methods
  results_bestorder <-  results_bestsubset  %>%
    group_by(para_id) %>%
    filter(MSE_CV_value == min(MSE_CV_value)) %>%
    mutate(bestorder = order ) %>%
    mutate(model_id = "optimal order") %>%
    mutate(order = "CV optimal") %>%
    mutate(extrapo=paste0("CV optimal order"))

  
  ## AIC Methods
  results_AIC <-  results_bestsubset  %>%
    group_by(para_id) %>%
    filter(AIC == min(AIC)) %>%
    mutate(bestorder = order ) %>%
    mutate(model_id = "optimal order") %>%
    mutate(order = "AIC-chosen order") %>%
    mutate(extrapo=paste0("AIC-chosen order"))
  
  ## BIC Methods
  results_BIC <-  results_bestsubset  %>%
    group_by(para_id) %>%
    filter(BIC == min(BIC)) %>%
    mutate(bestorder = order ) %>%
    mutate(model_id = "optimal order") %>%
    mutate(order = "BIC-chosen order") %>%
    mutate(extrapo=paste0("BIC-chosen order"))
  
  ## Model averaging AIC method
  results_MAAIC <- results_bestsubset %>%
    group_by(para_id) %>%
    mutate(w_AIC = exp(-(AIC - min(AIC)) / 2)) %>%                # Step 1: Create w_AIC
    mutate(w_AIC = w_AIC / sum(w_AIC)) %>%          # Step 2: Normalize w_AIC
    summarise(
      betahat = sum(betahat * w_AIC),               # Step 3: Weighted average of betahat
      model_id = NA,
      para_id = first(para_id),
      order = NA,
      SSR = NA,
      MSE_CV_value = NA,
      modelm = NA,
      AIC = NA,
      BIC = NA,
      bestmodel = NA,
      bestorder = NA,
      extrapo = "AIC ModelAvg",
      .groups = 'drop'
    )
  
  
  ## Model averaging AIC method
  results_MABIC <- results_bestsubset %>%
    group_by(para_id) %>%
    mutate(w_BIC = exp(-(BIC - min(BIC)) / 2)) %>%                # Step 1: Create w_BIC
    mutate(w_BIC = w_BIC / sum(w_BIC)) %>%          # Step 2: Normalize w_BIC
    summarise(
      betahat = sum(betahat * w_BIC),               # Step 3: Weighted average of betahat
      model_id = NA,
      para_id = first(para_id),
      order = NA,
      SSR = NA,
      MSE_CV_value = NA,
      modelm = NA,
      AIC = NA,
      BIC = NA,
      bestmodel = NA,
      bestorder = NA,
      extrapo = "BIC ModelAvg",
      .groups = 'drop'
    )
  
  # ## ALL Model averaging AIC method
  # results_MAAIC_ALL <- results %>%
  #   group_by(para_id) %>%
  #   mutate(w_AIC = exp(-(AIC - min(AIC)) / 2)) %>%                # Step 1: Create w_AIC
  #   mutate(w_AIC = w_AIC / sum(w_AIC)) %>%          # Step 2: Normalize w_AIC
  #   summarise(
  #     betahat = sum(betahat * w_AIC),               # Step 3: Weighted average of betahat
  #     model_id = NA,
  #     para_id = first(para_id),
  #     order = NA,
  #     SSR = NA,
  #     MSE_CV_value = NA,
  #     modelm = NA,
  #     AIC = NA,
  #     BIC = NA,
  #     bestmodel = NA,
  #     bestorder = NA,
  #     extrapo = "ALL AIC ModelAvg",
  #     .groups = 'drop'
  #   )
  # 
  # 
  # ## ALL Model averaging AIC method
  # results_MABIC_ALL <- results %>%
  #   group_by(para_id) %>%
  #   mutate(w_BIC = exp(-(BIC - min(BIC)) / 2)) %>%                # Step 1: Create w_BIC
  #   mutate(w_BIC = w_BIC / sum(w_BIC)) %>%          # Step 2: Normalize w_BIC
  #   summarise(
  #     betahat = sum(betahat * w_BIC),               # Step 3: Weighted average of betahat
  #     model_id = NA,
  #     para_id = first(para_id),
  #     order = NA,
  #     SSR = NA,
  #     MSE_CV_value = NA,
  #     modelm = NA,
  #     AIC = NA,
  #     BIC = NA,
  #     bestmodel = NA,
  #     bestorder = NA,
  #     extrapo = "ALL BIC ModelAvg",
  #     .groups = 'drop'
  #   )
  # 
  
  
  if (extrapolation == "linear") {
    return(results_linear)
  }
  
  if (extrapolation == "quadratic") {
    return(results_quadra)
  }
  
  if (extrapolation == "best subset") {
    return(results_bestorder)
  }
  
  if (extrapolation == "all") {
    return(rbind(results_linear, results_quadra, results_bestsubset, results_bestorder, 
                 results_AIC, results_BIC, results_MAAIC,  results_MABIC
                 # , results_MAAIC_ALL, results_MABIC_ALL
                 ))
  }
}

expit <- function(x){
  exp(x)/(1+exp(x))
}
