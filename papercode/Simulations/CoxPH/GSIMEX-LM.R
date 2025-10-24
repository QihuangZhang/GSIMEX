# For online server purpose only
# pass argument
args = commandArgs(trailingOnly=TRUE)

# for(i in 1:length(args)){
  eval(parse(text=args))
# }


range_eval <- range
# range_eval    <- as.integer( args )  # read second argument as integer
# range <- c(range_eval[1], range_eval[2])


# 0. Preparation ----------------------------------------------------------

library(pscl)
library(tidyverse)
library(doParallel)
# library(survival)

# 1. Global Parameters  ----------------------------------------------------------
ncores <- 50

nsample <- 1000

set.seed(2023)
seed_i <- sample(1000000,1000)


source("code/utility.R")


data_generation_LM <- function(SIM_id, beta, Mis_mtr, sd_error){
   # Data Generation
  set.seed(SIM_id)
  
  Covarmain1 <- data.frame(X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))

  mu <- as.matrix(Covarmain1) %*% t(t(beta))

  epsilon <- rnorm(nsample, 0, 1)
  
  Y <- mu + epsilon

  SimData <- data.frame(Covarmain1, Y)

  # Adding Mismesurement to Correct Data
  SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
  SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
  
  return(SimData)
}

SIMEX_LM <- function(SimData, sd_error, Mis_mtr,  maxq = 6, zeta_m = 2, B = 4000, nfold = 10, CVstep = 2, CVmfold = 3){
  coefs_Z <- NULL
  
  for (i in 1:B){
    epsilon <- i/B*zeta_m
    SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = epsilon)
    SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2star)
    model_ME <- lm(Y ~ W1 + W2, data = SimData)
    coefs_Z <- rbind(coefs_Z, c(epsilon, unlist(model_ME$coefficients)))
  }
  
   
  # Extrapolation
  q <- maxq
  l <- rep(list(0:1),q)
  EX <- expand.grid(l, KEEP.OUT.ATTRS =  F)
  Fullmodel_index <- 2^c(0:q)
  
  
  modelorder <- cut(1:2^q,breaks = Fullmodel_index, labels = 0:(q-1), include.lowest = FALSE, right = F)
  modelorder <- factor(modelorder, levels = c(0:q), labels = 0:q)
  modelorder[length(modelorder)] <- q
  
 
  knots <- (1:B)/B*zeta_m
  
  ZZ <- sapply(knots, FUN = function(x){
    return(x^c(0:q))
  })
  
  ZZ <- t(ZZ)
  
  dec_vec <- (-1)^c(0:q)
  
  
  EX <- cbind(1, EX)
  
  
  betahat_collect <- NULL
  EX_model_id_collect <- NULL
  para_index_collect <- NULL
  poly_order_collect <- NULL
  MSE_CV_value <- NULL
  SSR <- NULL
  
  for (j in 1:2){
    for (perm_idx in 1:(dim(EX)[1])){
        ZZ_sub <- ZZ[,as.logical(EX[perm_idx,])]
        dec_vec_sub <- dec_vec[as.logical(EX[perm_idx,])]
        
        est <- coefs_Z[,j+1]
        
        beta_model <- as.vector(solve(t(ZZ_sub)%*%ZZ_sub) %*% t(ZZ_sub)%*% est)
        
        betahat_poly <- sum(beta_model * dec_vec_sub)
        
        ## Get CV metric
        MSE_CV_value <- c(MSE_CV_value, get_MSE_CV_step(B, ZZ_sub, est, nfold, CVstep, CVmfold) )
        
        ## Calculating SSR
        SSR <- c(SSR, sum((est - ZZ_sub %*% t(t(beta_model)))^2) )
        
        
        betahat_collect <- c(betahat_collect, betahat_poly)
        EX_model_id_collect <- c(EX_model_id_collect, perm_idx)
        para_index_collect <- c(para_index_collect, j)
        poly_order_collect <- c(poly_order_collect, modelorder[perm_idx])
    }
  }
  
  results <- data.frame(betahat = betahat_collect,
                        model_id = EX_model_id_collect,
                        para_id = para_index_collect,
                        order = poly_order_collect-1,
                        SSR = SSR,
                        MSE_CV_value = MSE_CV_value,
                        modelm = rep(rowSums(EX)-1,2))
  
  return(results)
}




# 2. Simulation: trend of measurement error  ----------------------------------------------------------
# # 2.1 Continuous Variable x Discrete Variable-------
# i <- 1
# beta <- c(-0.5, 0.5)
# 
# # set.seed(seed_i[i])
# 
# Covarmain1 <- data.frame(X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))
# 
# 
# 
# mu <- as.matrix(Covarmain1) %*% t(t(beta))
# 
# U = runif(nsample,0,1)
# 
# 
# T = sqrt(-1*exp(-mu) * log(1-U))
# 
# C = runif(nsample,0,max(T))
# 
# Y = (T<C) * T + (T>C) * C
# delta = (T<C)*1
# 
# 
# # myrates <- exp(mu) # the risk exp(beta*x), parameters for exp r.v.
# # y <- rexp(nsample, rate = myrates) # generates the r.v.
# # 
# # # Optional: Censoring
# # censoring_time <- rexp(nsample, rate = 0.015) # Just an example rate
# # T <- pmin(y, censoring_time)
# # status <- ifelse(y <= censoring_time, 1, 0)
# 
# 
# SimData <- data.frame(Covarmain1, Y, delta)
# 
# 
# # Oracle method
# 
# coef_orc <- LM(Surv(Y, delta)~X1+X2,data = SimData)$coef
# print(coef_orc)
# 
# # Measurement Error
# sd_error <- 1
# Mis_mtr <- matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T)
# 
# SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
# SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
# 
# 
# coefs <- c(0, coef_orc)
# for (i in 1:2000){
#   epsilon <- i/2000*2
#   SimData$W1 <- SimData$X1 + rnorm(nsample, mean = 0, sd = epsilon)
#   SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2)
#   model_ME <- LM(Surv(Y, delta) ~ W1 + W2, data = SimData)
#   coefs <- rbind(coefs, c(i/2000*2, unlist(model_ME$coefficients)))
# }
# 
# colnames(coefs) <- c("epsilon", "con_beta1", "bin_beta2")
# 
# coefs_long <- data.frame(coefs) %>%
#   pivot_longer(con_beta1:bin_beta2, names_to = "variable")
# 
# ggplot(coefs_long, aes(x = epsilon, y = value)) +
#   geom_line() +
#   facet_wrap(~variable)

# 3. Simulation: GSIMEX with mixed variables in Section 2.5 ----------------------------------------------------------

# initialization parameters
SIM_id <- 1
beta1 <- c(-0.2,0.7)
beta2 <- c(1,0.5)
Mis_mtr <- matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T)
sd_error <- 1

### One iteration function

# Sim_LM_1 <- function(SIM_id, beta1, beta2, Mis_mtr, sd_error, zeta_m = 2, B = 4000){
#   # Data Generation
#   set.seed(SIM_id)
# 
#   Covarmain1 <- data.frame(1, X1 = runif(nsample,-3,4))
#   Covarmain2 <- data.frame(1, X2 = rbinom(nsample,1,0.5))
# 
# 
#   mu1 <- as.matrix(Covarmain1) %*% t(t(beta1))
#   phi <- expit(as.matrix(Covarmain1) %*% t(t(beta1)))
#   mu2 <- as.matrix(Covarmain2) %*% t(t(beta2))
#   lamdba <- exp(mu2)
# 
#   Y <- rbinom(nsample, 1, prob = phi) * rpois(nsample, lamdba)
# 
#   SimData <- data.frame(Y, X1 = Covarmain1$X1, X2 = Covarmain2$X2)
# 
#   # Oracle method
#   model_oracle <- zeroinfl( Y ~  X2 | X1, data = SimData)
#   coefs <- c(0, unlist(model_oracle$coefficients))
# 
# 
#   # Misclassification matrix
#   SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
#   SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
# 
#   coefs_Z <- NULL
# 
#   for (i in 1:B){
#     epsilon <- i/B*zeta_m
#     SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = epsilon)
#     SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2star)
#     model_ME <- zeroinfl( Y ~  W2 | W1, data = SimData)
#     coefs_Z <- rbind(coefs_Z, c(epsilon, unlist(model_ME$coefficients)))
#   }
# 
# 
#   coefs_Z[, 4:5] <- - coefs_Z[, 4:5]
# 
# 
#   # Extrapolation
#   q <- 5
#   l <- rep(list(0:1),q)
#   EX <- expand.grid(l, KEEP.OUT.ATTRS =  F)
#   Fullmodel_index <- 2^c(0:q)
# 
# 
#   modelorder <- cut(1:64,breaks = Fullmodel_index, labels = 0:(q-1), include.lowest = FALSE, right = F)
#   modelorder <- factor(modelorder, levels = c(0:q), labels = 0:q)
#   modelorder[length(modelorder)] <- q
# 
# 
#   knots <- (1:B)/B*zeta_m
# 
#   ZZ <- sapply(knots, FUN = function(x){
#     return(x^c(0:q))
#   })
# 
#   ZZ <- t(ZZ)
# 
#   dec_vec <- (-1)^c(0:q)
# 
# 
#   EX <- cbind(1, EX)
# 
# 
#   betahat_collect <- NULL
#   EX_model_id_collect <- NULL
#   para_index_collect <- NULL
#   poly_order_collect <- NULL
#   SSR <- NULL
# 
#   for (perm_idx in 1:(dim(EX)[1])){
#     for (j in 1:4){
#       ZZ_sub <- ZZ[,as.logical(EX[perm_idx,])]
#       dec_vec_sub <- dec_vec[as.logical(EX[perm_idx,])]
# 
#       est <- coefs_Z[,j+1]
# 
#       beta_model <- as.vector(solve(t(ZZ_sub)%*%ZZ_sub) %*% t(ZZ_sub)%*% est)
# 
#       betahat_poly <- sum(beta_model * dec_vec_sub)
# 
#       ## Calculating SSR
#       SSR <- c(SSR, sum((est - ZZ_sub %*% t(t(beta_model))^2)) )
# 
# 
#       betahat_collect <- c(betahat_collect, betahat_poly)
#       EX_model_id_collect <- c(EX_model_id_collect, perm_idx)
#       para_index_collect <- c(para_index_collect, j)
#       poly_order_collect <- c(poly_order_collect, modelorder[perm_idx])
#     }
#   }
# 
#   results <- data.frame(betahat = betahat_collect,
#                         model_id = EX_model_id_collect,
#                         para_id = para_index_collect,
#                         order = poly_order_collect-1,
#                         SSR = SSR)
# 
#   return(results)
# }


## Version with unit adaptation for high dimension

Sim_LM_CV <- function(SIM_id, beta, Mis_mtr, sd_error, maxq = 6, zeta_m = 2, B = 4000, nfold = 20, CVstep = 10, CVmfold = 3){
  # Data Generation
  set.seed(SIM_id)
  
  Covarmain1 <- data.frame(1,X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))

  mu <- as.matrix(Covarmain1) %*% t(t(beta))

  epsilon <- rnorm(nsample, 0, 1)
  
  Y <- mu + epsilon
  
  SimData <- data.frame(Covarmain1, Y)
  
  # Oracle method
  coef_orc <- lm( Y ~ X1 + X2,data = SimData)$coef
  
  
  # Adding Mismesurement to Correct Data
  SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
  SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
  
  coefs_Z <- NULL
  
  for (i in 1:B){
    epsilon <- i/B*zeta_m
    SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = epsilon)
    SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2star)
    model_ME <- lm( Y ~ W1 + W2, data = SimData)
    coefs_Z <- rbind(coefs_Z, c(epsilon, unlist(model_ME$coefficients)))
  }
  
   
  # Extrapolation
  q <- maxq
  l <- rep(list(0:1),q)
  EX <- expand.grid(l, KEEP.OUT.ATTRS =  F)
  Fullmodel_index <- 2^c(0:q)
  
  
  modelorder <- cut(1:2^q,breaks = Fullmodel_index, labels = 0:(q-1), include.lowest = FALSE, right = F)
  modelorder <- factor(modelorder, levels = c(0:q), labels = 0:q)
  modelorder[length(modelorder)] <- q
  

  # modelm <- rep(rowSums(EX)-1,2)
  
  knots <- (1:B)/B*zeta_m
  
  ZZ <- sapply(knots, FUN = function(x){
    return(x^c(0:q))
  })
  
  ZZ <- t(ZZ)
  
  dec_vec <- (-1)^c(0:q)
  
  
  EX <- cbind(1, EX)
  
  
  betahat_collect <- NULL
  EX_model_id_collect <- NULL
  para_index_collect <- NULL
  poly_order_collect <- NULL
  MSE_CV_value <- NULL
  SSR <- NULL
  
  for (j in 1:3){
    for (perm_idx in 1:(dim(EX)[1])){
        ZZ_sub <- ZZ[,as.logical(EX[perm_idx,])]
        dec_vec_sub <- dec_vec[as.logical(EX[perm_idx,])]
        
        est <- coefs_Z[,j+1]
        
        beta_model <- as.vector(solve(t(ZZ_sub)%*%ZZ_sub) %*% t(ZZ_sub)%*% est)
        
        betahat_poly <- sum(beta_model * dec_vec_sub)
        
        ## Get CV metric
        MSE_CV_value <- c(MSE_CV_value, get_MSE_CV_step(B, ZZ_sub, est, nfold, CVstep, CVmfold) )
        
        ## Calculating SSR
        SSR <- c(SSR, sum((est - ZZ_sub %*% t(t(beta_model)))^2) )
        
        
        betahat_collect <- c(betahat_collect, betahat_poly)
        EX_model_id_collect <- c(EX_model_id_collect, perm_idx)
        para_index_collect <- c(para_index_collect, j)
        poly_order_collect <- c(poly_order_collect, modelorder[perm_idx])
    }
  }
  
  results <- data.frame(betahat = betahat_collect,
                        model_id = EX_model_id_collect,
                        para_id = para_index_collect,
                        order = poly_order_collect-1,
                        SSR = SSR,
                        MSE_CV_value = MSE_CV_value,
                        modelm = rep(rowSums(EX)-1,3))
  
  return(results)
}

Sim_LM_CV_MEonly <- function(SIM_id, beta, Mis_mtr, sd_error, maxq = 6, zeta_m = 2, B = 4000, nfold = 20, CVstep = 10, CVmfold = 3){
  # Data Generation
  set.seed(SIM_id)
  
  Covarmain1 <- data.frame(1, X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))
  
  mu <- as.matrix(Covarmain1) %*% t(t(beta))
  
  epsilon <- rnorm(nsample, 0, 1)
  
  Y <- mu + epsilon
  
  SimData <- data.frame(Covarmain1, Y)
  
  # Oracle method
  coef_orc <- lm( Y ~ X1 + X2,data = SimData)$coef
  
  
  # Adding Mismesurement to Correct Data
  SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
  # SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
  
  coefs_Z <- NULL
  
  for (i in 1:B){
    epsilon <- i/B*zeta_m
    SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = epsilon)
    # SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2star)
    model_ME <- lm( Y ~ W1 + X2, data = SimData)
    coefs_Z <- rbind(coefs_Z, c(epsilon, unlist(model_ME$coefficients)))
  }
  
  
  # Extrapolation
  q <- maxq
  l <- rep(list(0:1),q)
  EX <- expand.grid(l, KEEP.OUT.ATTRS =  F)
  Fullmodel_index <- 2^c(0:q)
  
  
  modelorder <- cut(1:2^q,breaks = Fullmodel_index, labels = 0:(q-1), include.lowest = FALSE, right = F)
  modelorder <- factor(modelorder, levels = c(0:q), labels = 0:q)
  modelorder[length(modelorder)] <- q
  
  
  # modelm <- rep(rowSums(EX)-1,2)
  
  knots <- (1:B)/B*zeta_m
  
  ZZ <- sapply(knots, FUN = function(x){
    return(x^c(0:q))
  })
  
  ZZ <- t(ZZ)
  
  dec_vec <- (-1)^c(0:q)
  
  
  EX <- cbind(1, EX)
  
  
  betahat_collect <- NULL
  EX_model_id_collect <- NULL
  para_index_collect <- NULL
  poly_order_collect <- NULL
  MSE_CV_value <- NULL
  SSR <- NULL
  
  for (j in 1:3){
    for (perm_idx in 1:(dim(EX)[1])){
      ZZ_sub <- ZZ[,as.logical(EX[perm_idx,])]
      dec_vec_sub <- dec_vec[as.logical(EX[perm_idx,])]
      
      est <- coefs_Z[,j+1]
      
      beta_model <- as.vector(solve(t(ZZ_sub)%*%ZZ_sub) %*% t(ZZ_sub)%*% est)
      
      betahat_poly <- sum(beta_model * dec_vec_sub)
      
      ## Get CV metric
      MSE_CV_value <- c(MSE_CV_value, get_MSE_CV_step(B, ZZ_sub, est, nfold, CVstep, CVmfold) )
      
      ## Calculating SSR
      SSR <- c(SSR, sum((est - ZZ_sub %*% t(t(beta_model)))^2) )
      
      
      betahat_collect <- c(betahat_collect, betahat_poly)
      EX_model_id_collect <- c(EX_model_id_collect, perm_idx)
      para_index_collect <- c(para_index_collect, j)
      poly_order_collect <- c(poly_order_collect, modelorder[perm_idx])
    }
  }
  
  results <- data.frame(betahat = betahat_collect,
                        model_id = EX_model_id_collect,
                        para_id = para_index_collect,
                        order = poly_order_collect-1,
                        SSR = SSR,
                        MSE_CV_value = MSE_CV_value,
                        modelm = rep(rowSums(EX)-1,3))
  
  return(results)
}


# Sim_LM_CV_bootstrap <- function(SIM_id, beta1, beta2, Mis_mtr, sd_error, maxq = 6, zeta_m = 2, B = 4000, nfold = 10, CVstep = 2, CVmfold = 3, nboot = 50){
#   
#   SimData <- data_generation_LM(SIM_id, beta1, beta2, Mis_mtr, sd_error)
#   
#   results <- SIMEX_LM (SimData, sd_error, Mis_mtr, maxq, zeta_m, B, nfold, CVstep, CVmfold) 
#   
#   results$bootstrap_id <- 0
#   
#   for (boot_i in 1:nboot){
#     SimData_boot_index <- sample(dim(SimData)[1], replace = T)
#     SimData_boot <- SimData[SimData_boot_index,]
#     results_boot <- SIMEX_LM (SimData_boot, sd_error, Mis_mtr, maxq, zeta_m, B, nfold, CVstep, CVmfold)
#     
#     results_boot$bootstrap_id <- boot_i
#     results <- rbind(results, results_boot)
#     
#   }
#   
#   return(results)
# }





## 3.1 Simulation Start-------
### xx 3.1.1 Simple Model (not selected) --------

# cl <- makeCluster(ncores, outfile=paste0("GSIMEX.txt"))
# registerDoParallel(cl)
# 
# beta1true <- c(-0.2,0.7)
# beta2true <- c(1,0.5)
# 
# Simulation1 <- foreach(var1 = 1:30) %dopar% {
#   start <- proc.time()
# 
#   library(pscl)
#   library(tidyverse)
#   set.seed(2023)
#   seed_i <- sample(1000000,1000)
# 
# 
#   result <- Sim_LM_1(SIM_id = seed_i[var1],
#             beta1 = beta1true,
#             beta2 = beta2true,
#             Mis_mtr = matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T),
#             sd_error = 1,
#             zeta_m = 5)
# 
#   end <- proc.time() - start
#   cat(paste0(var1, ": ",end[3],"\n"))
# 
#   return(result)
# }
# 
# stopCluster(cl)
# 
# save(Simulation1, file=paste0("output/Simulation/GIMEX-LM-Sim1.1.2.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-LM-Sim1.1.2.RData"))
# 
# # str(Simulation1)
# 
# Sim1_results <- bind_rows(Simulation1, .id = "SimID")
# 
# Sim1_summary <- Sim1_results %>%
#   group_by(model_id, para_id) %>%
#   mutate(truth = factor(para_id, levels = 1:4, labels = c(beta2true, beta1true)) ) %>%
#   mutate(truth = as.numeric(as.character(truth))) %>%
#   mutate(bias = abs(betahat  - truth))
# 
# pdf("output/Simulation/GIMEX-LM-Sim1.1.2.pdf", width = 12, height = 6)
# 
# Sim1_summary$order <- factor(Sim1_summary$order)
# Sim1_summary$model_id <- factor(Sim1_summary$model_id)
# 
# boxplot <- ggplot(data = Sim1_summary) +
#   geom_boxplot(aes(x=model_id, y=bias, color = order)) +
#   facet_grid(para_id~.,scales = "free") +
#   coord_cartesian (ylim = c(0, 2))
# 
# boxplot
# 
# dev.off()


### 3.1.1.3 Simple Model (explore a smaller amount of B) (selected)--------

# cl <- makeCluster(ncores, outfile=paste0("GSIMEX.txt"))
# registerDoParallel(cl)
# 
# beta1true <- c(-0.2,0.7)
# beta2true <- c(1,0.5)
# 
# Simulation1.1.3 <- foreach(var1 = 1:500) %dopar% {
#   start <- proc.time()
#   
#   library(pscl)
#   library(tidyverse)
#   set.seed(2023)
#   seed_i <- sample(1000000,1000)
#   
#   
#   result <- Sim_LM_1(SIM_id = seed_i[var1],
#                       beta1 = beta1true,
#                       beta2 = beta2true,
#                       Mis_mtr = matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T),
#                       sd_error = 1,
#                       zeta_m = 5,
#                       B = 500)
#   
#   end <- proc.time() - start
#   cat(paste0(var1, ": ",end[3],"\n"))
#   
#   return(result)
# }
# 
# stopCluster(cl)
# 
# save(Simulation1.1.3, file=paste0("output/Simulation/GIMEX-LM-Sim1.1.3.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-LM-Sim1.1.3.RData"))
# 
# # str(Simulation1)
# 
# Sim1_results <- bind_rows(Simulation1.1.3, .id = "SimID")
# 
# 
# # Summarize the 
# 
# Sim1_summary <- Sim1_results %>%
#   mutate(truth = factor(para_id, levels = 1:4, labels = c(beta2true, beta1true)) ) %>%
#   mutate(truth = as.numeric(as.character(truth))) %>%
#   mutate(bias = abs(betahat  - truth)) %>%    
#   group_by(SimID, model_id) %>%
#   mutate(avg_SSR = mean(SSR)) %>%
#   mutate(bestmodel = NA )
# 
# Sim1_bestsubset <-  Sim1_summary  %>% 
#     group_by(SimID, para_id, order) %>%
#     filter(SSR == min(SSR)) %>%
#     mutate(bestmodel = model_id ) %>%
#     filter(order>0) %>%
#     mutate(model_id = paste("best of", order))
# 
# Sim1_summary$model_id <- as.character(Sim1_summary$model_id)
# Sim1_BS_all <- rbind(Sim1_summary, Sim1_bestsubset)
# 
# 
# 
# pdf("output/Simulation/GIMEX-LM-Sim1.1.3.pdf", width = 12, height = 6)
# 
# Sim1_BS_all$order <- factor(Sim1_BS_all$order)
# Sim1_BS_all$model_id <- factor(Sim1_BS_all$model_id, levels = c(1:64, paste("best of", 1:6)))
# 
# boxplot <- ggplot(data = Sim1_BS_all) +
#   geom_boxplot(aes(x=model_id, y=bias, color = order)) +
#   facet_grid(para_id~.,scales = "free") +
#   coord_cartesian (ylim = c(0, 2))
# 
# boxplot
# 
# dev.off()


### xx 3.1.1.4 The Model 3.1.1.3 with cross validation --------
# Evolved to 3.1.1.5 by introducing more fold and only focus on the 
# prediction for the first few folds


# ### (!exploreing) 3.1.1.5 The Model 3.1.1.4 with cross validation --------
# 
# ## update 2023-10-16:
# ## Increase B from 500 to 5000
# 
# 
# 
# # cl <- makeCluster(ncores, outfile=paste0("GSIMEX_LM_3.1.1.5", range_eval[1], "-", range_eval[2],".txt"))
# # registerDoParallel(cl)
# # 
# betatrue <- c(-0.5, 0.5)
# # 
# # 
# # Simulation1.1.5 <- foreach(var1 = range_eval[1]:range_eval[2]) %dopar% {
# #   start <- proc.time()
# # 
# #   library(pscl)
# #   library(tidyverse)
# #   library(survival)
# #   set.seed(2023)
# #   seed_i <- sample(1000000, 1000)
# # 
# # 
# #   result <- Sim_LM_CV(SIM_id = seed_i[var1],
# #                       beta = betatrue,
# #                       Mis_mtr = matrix(c(0.9,0.1,0.1,0.9), nrow = 2, byrow = T),
# #                       sd_error = 0.5,
# #                       zeta_m = 5,
# #                       B = 10000,
# #                       nfold = 20,
# #                       CVstep = 10,
# #                       CVmfold = 1)
# # 
# #   end <- proc.time() - start
# #   cat(paste0(var1, ": ",end[3],"\n"))
# # 
# #   return(result)
# # }
# # 
# # stopCluster(cl)
# # 
# # save(Simulation1.1.5, file=paste0("output/Simulation/GIMEX-LM-Sim1.1.5", range_eval[1], "-", range_eval[2],".RData"))
# # 
# # 
# # List files matching the specified pattern
# matching_files <- list.files(path = "output/Simulation", pattern = "GIMEX-LM-Sim1.1.5.*\\.RData$", full.names = TRUE)
# 
# Simulation1.1.5_list <- NULL
# 
# for (file in matching_files) {
#   # Load the file
#   loaded_data <- load(file)
# 
#   # The load function returns the name of the loaded object, which we can then use to access the object
#   data_name <- loaded_data[1]
# 
#   # Store the loaded data frame in the data_list
#   Simulation1.1.5_list <- c(Simulation1.1.5_list, get(data_name))
# 
#   # Cleanup: remove the loaded data frame from the environment
#   rm(list = data_name)
# }
# 
# Simulation1.1.5 <- Simulation1.1.5_list
# 
# 
# 
# # str(Simulation1)
# 
# 
# Sim1_results <- bind_rows(Simulation1.1.5, .id = "SimID")
# results_all <- Sim1_results %>%
#               group_by(SimID) %>%
#               group_modify(~Summarize_Results(.x, B = 5000, extrapolation = "all"))
# 
# 
# # Attach the true value and calculate the bias
# 
# results_all_plot <- results_all %>%
#   mutate(truth = factor(para_id, levels = 1:2, labels = c(betatrue)) ) %>%
#   mutate(truth = as.numeric(as.character(truth))) %>%
#   mutate(bias = abs(betahat  - truth)) %>%
#   mutate(variable_greek = factor(para_id, levels = 1:2, labels = c( expression(beta["con"]),expression(beta["dis"]) )))
# 
# 
# 
# 
# # # Summarize the
# pdf("output/Simulation/GIMEX-LM-Sim1.1.5.pdf", width = 14, height = 11)
# 
# # Sim1_BS_all$modelm <- factor(Sim1_BS_all$modelm)
# # Sim1_BS_all$model_id <- factor(Sim1_BS_all$model_id, levels = c(1:64, paste("best of", 1:6), "CV optimal m"))
# 
# results_all_plot$extrapo <- factor(results_all$extrapo,
#                                    levels = c("linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal order"), 
#                                    labels = c("linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal m") )
# 
# boxplot <- results_all_plot %>%
#   mutate(modelm = ifelse(extrapo == "CV optimal m", "optimal", modelm)) %>%
#   filter(!is.na(extrapo)) %>%
#   ggplot() +
#   geom_violin(aes(x=extrapo, y=bias, fill = modelm), trim=FALSE) +
#   geom_boxplot(aes(x=extrapo, y=bias), width=0.1, position= position_nudge(x=.15)) +
#   facet_grid(.~variable_greek,labeller = label_parsed, scales = "free") +
#   # coord_cartesian (ylim = c(0, 1.2))
#   scale_fill_manual(values=c( "#27BDBE", "#0096C9", "#B2D235",  "#FFD400", "#F7941D", "#C768A9", "#9E0918", "#5D6770")) + 
#   theme_bw() +
#   theme(text=element_text(size=22), axis.text = element_text(size = 22),
#         axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
#   theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
#   theme(strip.text = element_text(colour = 'white')) +
#   theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "none")  +
#   labs(x = "Extrapolation", y = "Bias", fill = "m" )
# 
# boxplot
# 
# dev.off()
# # 
# # 
# # pdf("output/Simulation/GIMEX-LM-Sim1.1.5_distribution.pdf", width = 14, height = 11)
# # 
# # # Sim1_BS_all$modelm <- factor(Sim1_BS_all$modelm)
# # # Sim1_BS_all$model_id <- factor(Sim1_BS_all$model_id, levels = c(1:64, paste("best of", 1:6), "CV optimal m"))
# # 
# # results_all_plot$extrapo <- factor(results_all$extrapo,
# #                                    levels = c("linear", "quadratic", paste0("best subset of m=",0:6), "CV optimal m") )
# # 
# # count_dist <- results_all_plot %>%
# #   filter(extrapo == "CV optimal m")  %>%
# #   group_by(para_id, modelm) %>%
# #   summarise(count = n()) %>%
# #   group_by(para_id) %>%
# #   mutate(per =  count/sum(count)) %>%
# #   data.frame()
# # 
# # boxplot <- results_all_plot %>%
# #   filter(extrapo == "CV optimal m")  %>%
# #   ggplot() +
# #   geom_bar( data = count_dist, aes(x=modelm, y=per, group = para_id), color = "#3F4536",stat="identity") +
# #   geom_boxplot(aes(x=modelm, y=bias, color = modelm)) +
# #   facet_grid(para_id~.,scales = "free") +
# #   # coord_cartesian (ylim = c(0, 1.2))
# #   theme_bw() +
# #   theme(text=element_text(size=20), axis.text = element_text(size = 20),
# #         axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
# #   theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
# #   theme(strip.text = element_text(colour = 'white')) +
# #   theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "none")  +
# #   labs(x = "Extrapolation", y = "Bias", fill = "m" )
# # 
# # boxplot
# # 
# # dev.off()
# # 
# 
# 
# 
# 
# 
# ### (!exploreing) 3.1.1.6 The Model 3.1.1.4 with cross validation --------
# 
# ## update 2023-10-16:
# ## Increase B from 500 to 5000
# 
# 
# 
# cl <- makeCluster(ncores, outfile=paste0("GSIMEX_LM_3.1.1.6", range_eval[1], "-", range_eval[2],".txt"))
# registerDoParallel(cl)
# 
# betatrue <- c(-0.5, 0.5)
# 
# 
# Simulation1.1.6 <- foreach(var1 = range_eval[1]:range_eval[2]) %dopar% {
#   start <- proc.time()
#   
#   library(pscl)
#   library(tidyverse)
#   library(survival)
#   set.seed(2023)
#   seed_i <- sample(1000000, 1000)
#   
#   
#   result <- Sim_LM_CV(SIM_id = seed_i[var1],
#                          beta = betatrue,
#                          Mis_mtr = matrix(c(0.7,0.3,0.3,0.7), nrow = 2, byrow = T),
#                          sd_error = 2,
#                          zeta_m = 5,
#                          B = 10000,
#                          nfold = 20,
#                          CVstep = 10,
#                          CVmfold = 1)
#   
#   end <- proc.time() - start
#   cat(paste0(var1, ": ",end[3],"\n"))
#   
#   return(result)
# }
# 
# stopCluster(cl)
# 
# save(Simulation1.1.6, file=paste0("output/Simulation/GIMEX-LM-Sim1.1.6", range_eval[1], "-", range_eval[2],".RData"))
# 
# 
# # List files matching the specified pattern
# matching_files <- list.files(path = "output/Simulation", pattern = "GIMEX-LM-Sim1.1.6.*\\.RData$", full.names = TRUE)
# 
# Simulation1.1.6_list <- NULL
# 
# for (file in matching_files) {
#   # Load the file
#   loaded_data <- load(file)
#   
#   # The load function returns the name of the loaded object, which we can then use to access the object
#   data_name <- loaded_data[1]
#   
#   # Store the loaded data frame in the data_list
#   Simulation1.1.6_list <- c(Simulation1.1.6_list, get(data_name))
#   
#   # Cleanup: remove the loaded data frame from the environment
#   rm(list = data_name)
# }
# 
# Simulation1.1.6 <- Simulation1.1.6_list
# 
# 
# 
# # str(Simulation1)
# 
# 
# Sim1_results <- bind_rows(Simulation1.1.6, .id = "SimID")
# results_all <- Sim1_results %>%
#   group_by(SimID) %>%
#   group_modify(~Summarize_Results(.x, B = 5000, extrapolation = "all"))
# 
# 
# # Attach the true value and calculate the bias
# 
# results_all_plot <- results_all %>%
#   mutate(truth = factor(para_id, levels = 1:2, labels = c(betatrue)) ) %>%
#   mutate(truth = as.numeric(as.character(truth))) %>%
#   mutate(bias = abs(betahat  - truth)) %>%
#   mutate(variable_greek = factor(para_id, levels = 1:2, labels = c( expression(beta["con"]),expression(beta["dis"]) )))
# 
# 
# 
# # # Summarize the
# pdf("output/Simulation/GIMEX-LM-Sim1.1.6.pdf", width = 14, height = 11)
# 
# # Sim1_BS_all$modelm <- factor(Sim1_BS_all$modelm)
# # Sim1_BS_all$model_id <- factor(Sim1_BS_all$model_id, levels = c(1:64, paste("best of", 1:6), "CV optimal m"))
# 
# results_all_plot$extrapo <- factor(results_all$extrapo,
#                                    levels = c("linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal order"), 
#                                    labels = c("linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal m") )
# 
# boxplot <- results_all_plot %>%
#   mutate(modelm = ifelse(extrapo == "CV optimal m", "optimal", modelm)) %>%
#   filter(!is.na(extrapo)) %>%
#   ggplot() +
#   geom_violin(aes(x=extrapo, y=bias, fill = modelm), trim=FALSE) +
#   geom_boxplot(aes(x=extrapo, y=bias), width=0.1, position= position_nudge(x=.15)) +
#   facet_grid(.~variable_greek,labeller = label_parsed, scales = "free") +
#   # coord_cartesian (ylim = c(0, 1.2))
#   scale_fill_manual(values=c( "#27BDBE", "#0096C9", "#B2D235",  "#FFD400", "#F7941D", "#C768A9", "#9E0918", "#5D6770")) + 
#   theme_bw() +
#   theme(text=element_text(size=22), axis.text = element_text(size = 22),
#         axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
#   theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
#   theme(strip.text = element_text(colour = 'white')) +
#   theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "none")  +
#   labs(x = "Extrapolation", y = "Bias", fill = "m" )
# 
# boxplot
# 
# dev.off()
# 
# 
# pdf("output/Simulation/GIMEX-LM-Sim1.1.6_distribution.pdf", width = 14, height = 11)
# 
# # Sim1_BS_all$modelm <- factor(Sim1_BS_all$modelm)
# # Sim1_BS_all$model_id <- factor(Sim1_BS_all$model_id, levels = c(1:64, paste("best of", 1:6), "CV optimal m"))
# 
# results_all_plot$extrapo <- factor(results_all$extrapo,
#                                    levels = c("linear", "quadratic", paste0("best subset of m=",0:6), "CV optimal m") )
# 
# count_dist <- results_all_plot %>%
#   filter(extrapo == "CV optimal m")  %>%
#   group_by(para_id, modelm) %>%
#   summarise(count = n()) %>%
#   group_by(para_id) %>%
#   mutate(per =  count/sum(count)) %>%
#   data.frame()
# 
# boxplot <- results_all_plot %>%
#   filter(extrapo == "CV optimal m")  %>%
#   ggplot() +
#   geom_bar( data = count_dist, aes(x=modelm, y=per, group = para_id), color = "#3F4536",stat="identity") +
#   geom_boxplot(aes(x=modelm, y=bias, color = modelm)) +
#   facet_grid(para_id~.,scales = "free") +
#   # coord_cartesian (ylim = c(0, 1.2))
#   theme_bw() +
#   theme(text=element_text(size=20), axis.text = element_text(size = 20),
#         axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
#   theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
#   theme(strip.text = element_text(colour = 'white')) +
#   theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "none")  +
#   labs(x = "Extrapolation", y = "Bias", fill = "m" )
# 
# boxplot
# 
# dev.off()



### (!exploreing) 3.1.1.7 The Model 3.1.1.6  but with only one variable subject to measurement error --------


cl <- makeCluster(ncores, outfile=paste0("GSIMEX_LM_3.1.1.7", range_eval[1], "-", range_eval[2],".txt"))
registerDoParallel(cl)

betatrue <- c(1, -0.5, 0.5)


Simulation1.1.7 <- foreach(var1 = range_eval[1]:range_eval[2]) %dopar% {
  start <- proc.time()
  
  library(pscl)
  library(tidyverse)
  # library(survival)
  set.seed(2023)
  seed_i <- sample(1000000, 1000)
  
  
  result <- Sim_LM_CV_MEonly(SIM_id = seed_i[var1],
                         beta = betatrue,
                         Mis_mtr = matrix(c(0.7,0.3,0.3,0.7), nrow = 2, byrow = T),
                         sd_error = 2,
                         zeta_m = 5,
                         B = 10000,
                         nfold = 20,
                         CVstep = 10,
                         CVmfold = 1)
  
  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))
  
  return(result)
}

stopCluster(cl)

save(Simulation1.1.7, file=paste0("output/Simulation/GIMEX-LM-Sim1.1.7", range_eval[1], "-", range_eval[2],".RData"))


# List files matching the specified pattern
matching_files <- list.files(path = "output/Simulation", pattern = "GIMEX-LM-Sim1.1.7.*\\.RData$", full.names = TRUE)

Simulation1.1.7_list <- NULL

for (file in matching_files) {
  # Load the file
  loaded_data <- load(file)

  # The load function returns the name of the loaded object, which we can then use to access the object
  data_name <- loaded_data[1]

  # Store the loaded data frame in the data_list
  Simulation1.1.7_list <- c(Simulation1.1.7_list, get(data_name))

  # Cleanup: remove the loaded data frame from the environment
  rm(list = data_name)
}

Simulation1.1.7 <- Simulation1.1.7_list



# str(Simulation1)


Sim1_results <- bind_rows(Simulation1.1.7, .id = "SimID")
results_all <- Sim1_results %>%
  group_by(SimID) %>%
  group_modify(~Summarize_Results(.x, B = 10000, extrapolation = "all"))


# Attach the true value and calculate the bias

results_all_plot <- results_all %>%
  mutate(truth = factor(para_id, levels = 1:3, labels = c(betatrue)) ) %>%
  mutate(truth = as.numeric(as.character(truth))) %>%
  mutate(bias = abs(betahat  - truth)) %>%
  mutate(variable_greek = factor(para_id, levels = 1:3, labels = c( expression(beta[0]), expression(beta[1]), expression(beta[2]) )))


pdf("output/Simulation/GIMEX-LM-Sim1.1.7.pdf", width = 14, height = 11)

# Sim1_BS_all$modelm <- factor(Sim1_BS_all$modelm)
# Sim1_BS_all$model_id <- factor(Sim1_BS_all$model_id, levels = c(1:64, paste("best of", 1:6), "CV optimal m"))

results_all_plot$extrapo <- factor(results_all$extrapo,
                                   levels = c("linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal order"),
                                   labels = c("linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal m") )

boxplot <- results_all_plot %>%
  mutate(modelm = ifelse(extrapo == "CV optimal m", "optimal", modelm)) %>%
  filter(!is.na(extrapo)) %>%
  ggplot() +
  geom_violin(aes(x=extrapo, y=bias, fill = modelm), trim=FALSE) +
  geom_boxplot(aes(x=extrapo, y=bias), width=0.1, position= position_nudge(x=.15)) +
  facet_grid(.~variable_greek,labeller = label_parsed, scales = "free") +
  # coord_cartesian (ylim = c(0, 1.2))
  scale_fill_manual(values=c( "#27BDBE", "#0096C9", "#B2D235",  "#FFD400", "#F7941D", "#C768A9", "#9E0918", "#5D6770")) +
  theme_bw() +
  theme(text=element_text(size=22), axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
  theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "none")  +
  labs(x = "Extrapolation", y = "Bias", fill = "m" )

boxplot

dev.off()


### xxx High-order adjusted Model -------- 

# cl <- makeCluster(ncores, outfile=paste0("GSIMEX_highorder.txt"))
# registerDoParallel(cl)
# 
# beta1true <- c(-0.2,0.7)
# beta2true <- c(1,0.5)
# 
# Simulation1.2 <- foreach(var1 = 1:30) %dopar% { 
#   start <- proc.time()
#   
#   library(pscl)
#   library(tidyverse)
#   set.seed(2023)
#   seed_i <- sample(1000000,1000)
#   
#   
#   result <- Sim_LM_unit(SIM_id = seed_i[var1], 
#                       beta1 = beta1true, 
#                       beta2 = beta2true, 
#                       Mis_mtr = matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T), 
#                       sd_error = 1)
#   
#   end <- proc.time() - start
#   cat(paste0(var1, ": ",end[3],"\n"))
#   
#   return(result)
# }
# 
# stopCluster(cl)
# 
# save(Simulation1.2, file=paste0("output/Simulation/GIMEX-LM-Sim1.2_Unit.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-LM-Sim1.2_Unit.RData"))
# 
# # str(Simulation1)
# 
# Sim1_results <- bind_rows(Simulation1.2, .id = "SimID")
# 
# Sim1_summary <- Sim1_results %>%
#   group_by(model_id, para_id) %>%
#   mutate(truth = factor(para_id, levels = 1:4, labels = c(beta2true, beta1true)) ) %>%
#   mutate(truth = as.numeric(as.character(truth))) %>%
#   mutate(bias = abs(betahat  - truth))
# 
# pdf("output/Simulation/GIMEX-LM-Sim1.2_Unit.pdf", width = 12, height = 6)
# 
# Sim1_summary$order <- factor(Sim1_summary$order)
# Sim1_summary$model_id <- factor(Sim1_summary$model_id)
# 
# boxplot <- ggplot(data = Sim1_summary) +
#   geom_boxplot(aes(x=model_id, y=bias, color = order)) +
#   facet_grid(para_id~.,scales = "free") + 
#   coord_cartesian (ylim = c(0, 2)) 
# 
# boxplot
# 
# dev.off()


### Exploration: 3.1.2.1 High-order adjusted Model (higher zeta_M)

# cl <- makeCluster(ncores, outfile=paste0("GSIMEX_highorder.txt"))
# registerDoParallel(cl)
# 
# beta1true <- c(-0.2,0.7)
# beta2true <- c(1,0.5)
# 
# Simulation1.2.1 <- foreach(var1 = 1:30) %dopar% { 
#   start <- proc.time()
#   
#   library(pscl)
#   library(tidyverse)
#   set.seed(2023)
#   seed_i <- sample(1000000,1000)
#   
#   
#   result <- Sim_LM_unit(SIM_id = seed_i[var1], 
#                          beta1 = beta1true, 
#                          beta2 = beta2true, 
#                          Mis_mtr = matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T), 
#                          sd_error = 1,
#                          zeta_m = 5)
#   
#   end <- proc.time() - start
#   cat(paste0(var1, ": ",end[3],"\n"))
#   
#   return(result)
# }
# 
# stopCluster(cl)
# 
# save(Simulation1.2.1, file=paste0("output/Simulation/GIMEX-LM-Sim1.2_Unit.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-LM-Sim1.2_Unit.RData"))
# 
# # str(Simulation1)
# 
# Sim1_results <- bind_rows(Simulation1.2.1, .id = "SimID")
# 
# Sim1_summary <- Sim1_results %>%
#   group_by(model_id, para_id) %>%
#   mutate(truth = factor(para_id, levels = 1:4, labels = c(beta2true, beta1true)) ) %>%
#   mutate(truth = as.numeric(as.character(truth))) %>%
#   mutate(bias = abs(betahat  - truth))
# 
# pdf("output/Simulation/GIMEX-LM-Sim1.2.1.pdf", width = 12, height = 6)
# 
# Sim1_summary$order <- factor(Sim1_summary$order)
# Sim1_summary$model_id <- factor(Sim1_summary$model_id)
# 
# boxplot <- ggplot(data = Sim1_summary) +
#   geom_boxplot(aes(x=model_id, y=bias, color = order)) +
#   facet_grid(para_id~.,scales = "free") + 
#   coord_cartesian (ylim = c(0, 2)) 
# 
# boxplot
# 
# dev.off()


rm(list = ls())