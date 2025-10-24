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
library(survival)
library(xtable)

# 1. Global Parameters  ----------------------------------------------------------
ncores <- 50

nsample <- 1000

set.seed(2023)
seed_i <- sample(1000000,1000)


source("code/utility.R")


data_generation_coxPH <- function(SIM_id, beta, Mis_mtr, sd_error){
   # Data Generation
  set.seed(SIM_id)
  
  Covarmain1 <- data.frame(X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))

  mu <- as.matrix(Covarmain1) %*% t(t(beta))

  U = runif(nsample,0,1)


  T = sqrt(-1*exp(-mu) * log(1-U))

  C = runif(nsample,0,max(T))

  Y = (T<C) * T + (T>C) * C
  delta = (T<C)*1

  SimData <- data.frame(Covarmain1, Y, delta)

  # Adding Mismesurement to Correct Data
  SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
  SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
  
  return(SimData)
}

SIMEX_coxPH <- function(SimData, sd_error, Mis_mtr,  maxq = 6, zeta_m = 2, B = 4000, nfold = 10, CVstep = 2, CVmfold = 3){
  coefs_Z <- NULL
  
  for (i in 1:B){
    epsilon <- i/B*zeta_m
    SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = epsilon)
    SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2star)
    model_ME <- coxph(Surv(Y, delta) ~ W1 + W2, data = SimData)
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
i <- 1
beta <- c(-0.5, 0.5)

# set.seed(seed_i[i])

Covarmain1 <- data.frame(X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))



mu <- as.matrix(Covarmain1) %*% t(t(beta))

U = runif(nsample,0,1)


T = sqrt(-1*exp(-mu) * log(1-U))

C = runif(nsample,0,max(T))

Y = (T<C) * T + (T>C) * C
delta = (T<C)*1


# myrates <- exp(mu) # the risk exp(beta*x), parameters for exp r.v.
# y <- rexp(nsample, rate = myrates) # generates the r.v.
# 
# # Optional: Censoring
# censoring_time <- rexp(nsample, rate = 0.015) # Just an example rate
# T <- pmin(y, censoring_time)
# status <- ifelse(y <= censoring_time, 1, 0)


SimData <- data.frame(Covarmain1, Y, delta)


# Oracle method

coef_orc <- coxph(Surv(Y, delta)~X1+X2,data = SimData)$coef
print(coef_orc)

# Measurement Error
sd_error <- 1
Mis_mtr <- matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T)

SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)


coefs <- c(0, coef_orc)
for (i in 1:2000){
  epsilon <- i/2000*2
  SimData$W1 <- SimData$X1 + rnorm(nsample, mean = 0, sd = epsilon)
  SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2)
  model_ME <- coxph(Surv(Y, delta) ~ W1 + W2, data = SimData)
  coefs <- rbind(coefs, c(i/2000*2, unlist(model_ME$coefficients)))
}

colnames(coefs) <- c("epsilon", "con_beta1", "bin_beta2")

coefs_long <- data.frame(coefs) %>%
  pivot_longer(con_beta1:bin_beta2, names_to = "variable")

ggplot(coefs_long, aes(x = epsilon, y = value)) +
  geom_line() +
  facet_wrap(~variable)

# 3. Simulation: GSIMEX with mixed variables in Section 2.5 ----------------------------------------------------------

# initialization parameters
SIM_id <- 1
beta1 <- c(-0.2,0.7)
beta2 <- c(1,0.5)
Mis_mtr <- matrix(c(0.7,0.3,0.2,0.8), nrow = 2, byrow = T)
sd_error <- 1

### One iteration function

# Sim_coxPH_1 <- function(SIM_id, beta1, beta2, Mis_mtr, sd_error, zeta_m = 2, B = 4000){
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

Sim_coxPH_CV <- function(SIM_id, beta, Mis_mtr, sd_error, maxq = 6, zeta_m = 2, B = 4000, nfold = 20, CVstep = 10, CVmfold = 3){
  # Data Generation
  set.seed(SIM_id)
  
  Covarmain1 <- data.frame(X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))

  mu <- as.matrix(Covarmain1) %*% t(t(beta))

  U = runif(nsample,0,1)


  T = sqrt(-1*exp(-mu) * log(1-U))

  C = runif(nsample,0,max(T))

  Y = (T<C) * T + (T>C) * C
  delta = (T<C)*1

  
  SimData <- data.frame(Covarmain1, Y, delta)
  
  # Oracle method
  coef_orc <- coxph(Surv(Y, delta)~X1+X2,data = SimData)$coef
  
  
  # Adding Mismesurement to Correct Data
  SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
  SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
  
  coefs_Z <- NULL
  
  for (i in 1:B){
    epsilon <- i/B*zeta_m
    SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = sqrt(epsilon) )
    SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2star)
    model_ME <- coxph(Surv(Y, delta) ~ W1 + W2, data = SimData)
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

Sim_coxPH_naive <- function(SIM_id, beta, Mis_mtr, sd_error){
  # Data Generation
  set.seed(SIM_id)
  
  Covarmain1 <- data.frame(X1 = runif(nsample,-3,4), X2 = rbinom(nsample,1,0.5))
  
  mu <- as.matrix(Covarmain1) %*% t(t(beta))
  
  U = runif(nsample,0,1)
  
  
  T = sqrt(-1*exp(-mu) * log(1-U))
  
  C = runif(nsample,0,max(T))
  
  Y = (T<C) * T + (T>C) * C
  delta = (T<C)*1
  
  
  SimData <- data.frame(Covarmain1, Y, delta)
  
  # Oracle method
  coef_orc <- coxph(Surv(Y, delta)~X1+X2,data = SimData)$coef
  
  
  # Adding Mismesurement to Correct Data
  SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
  SimData$X2star <- MCoperator(Mis_mtr, zeta = 1, SimData$X2)
  
  model_ME <- coxph(Surv(Y, delta) ~ X1star + X2star, data = SimData)
  
  coefs_Z <- model_ME$coefficients
  
  # for (i in 1:B){
  #   epsilon <- i/B*zeta_m
  #   SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = epsilon)
  #   SimData$W2 <- MCoperator(Mis_mtr, zeta = epsilon, SimData$X2star)
  #   model_ME <- coxph(Surv(Y, delta) ~ W1 + W2, data = SimData)
  #   coefs_Z <- rbind(coefs_Z, c(epsilon, unlist(model_ME$coefficients)))
  # }
  
  results <- data.frame(betahat = coefs_Z,
                        model_id = 0,
                        para_id = 1:2,
                        order = NA,
                        SSR = NA,
                        MSE_CV_value = NA,
                        modelm = NA)
  
  return(results)
}



Sim_coxPH_CV_single <- function(SIM_id, beta, Mis_mtr, sd_error, maxq = 6, zeta_m = 2, B = 4000, nfold = 20, CVstep = 10, CVmfold = 3){
  # Data Generation
  set.seed(SIM_id)
  
  Covarmain1 <- data.frame(X1 = runif(nsample,-3,4))
  
  mu <- as.matrix(Covarmain1) * beta
  
  U = runif(nsample,0,1)
  
  
  T = sqrt(-1*exp(-mu) * log(1-U))
  
  C = runif(nsample,0,max(T))
  
  Y = (T<C) * T + (T>C) * C
  delta = (T<C)*1
  
  
  SimData <- data.frame(Covarmain1, Y, delta)
  
  # Oracle method
  coef_orc <- coxph(Surv(Y, delta)~X1,data = SimData)$coef
  
  
  # Adding Mismesurement to Correct Data
  SimData$X1star <- SimData$X1 + rnorm(nsample, mean = 0, sd = sd_error)
  
  
  # Naive model
  
  model_naive <- coxph(Surv(Y, delta) ~ X1star, data = SimData)
  
  #
  
  coefs_Z <- NULL
  
  for (i in 1:B){
    epsilon <- i/B*zeta_m
    SimData$W1 <- SimData$X1star + rnorm(nsample, mean = 0, sd = epsilon)
    model_ME <- coxph(Surv(Y, delta) ~ W1, data = SimData)
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
  
  
  betahat_collect <- model_naive$coefficients
  EX_model_id_collect <- rep(0, length(model_naive$coefficients))
  para_index_collect <- 1:length(model_naive$coefficients)
  poly_order_collect <-  rep(NA, length(model_naive$coefficients))
  MSE_CV_value <- rep(NA, length(model_naive$coefficients))
  SSR <- rep(NA, length(model_naive$coefficients))
  
  for (j in 1){
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
                        modelm = c(NA, rep(rowSums(EX)-1,1)))
  
  return(results)
}







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
#   result <- Sim_coxPH_1(SIM_id = seed_i[var1],
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
# save(Simulation1, file=paste0("output/Simulation/GIMEX-coxPH-Sim1.1.2.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-coxPH-Sim1.1.2.RData"))
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
# pdf("output/Simulation/GIMEX-coxPH-Sim1.1.2.pdf", width = 12, height = 6)
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
#   result <- Sim_coxPH_1(SIM_id = seed_i[var1],
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
# save(Simulation1.1.3, file=paste0("output/Simulation/GIMEX-coxPH-Sim1.1.3.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-coxPH-Sim1.1.3.RData"))
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
# pdf("output/Simulation/GIMEX-coxPH-Sim1.1.3.pdf", width = 12, height = 6)
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


### (!selected) 3.1.1.5 The Model 3.1.1.4 with cross validation --------

## update 2023-10-16:
## Increase B from 500 to 5000



cl <- makeCluster(ncores, outfile=paste0("GSIMEX_COXPH_3.1.1.5", range_eval[1], "-", range_eval[2],".txt"))
registerDoParallel(cl)

betatrue <- c(-0.5, 0.5)


Simulation1.1.5 <- foreach(var1 = range_eval[1]:range_eval[2]) %dopar% {
  start <- proc.time()

  library(pscl)
  library(tidyverse)
  library(survival)
  set.seed(2023)
  seed_i <- sample(1000000, 1000)


  result <- Sim_coxPH_CV(SIM_id = seed_i[var1],
                      beta = betatrue,
                      Mis_mtr = matrix(c(0.9,0.1,0.1,0.9), nrow = 2, byrow = T),
                      sd_error = 1,
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

save(Simulation1.1.5, file=paste0("output/Simulation/GIMEX-coxPH-Sim1.1.5", range_eval[1], "-", range_eval[2],".RData"))




Simulation1.1.5_naive <- lapply(1:500, FUN = Sim_coxPH_naive, beta = betatrue,
                                Mis_mtr = matrix(c(0.9,0.1,0.1,0.9), nrow = 2, byrow = T),
                                sd_error = 1)


save(Simulation1.1.5_naive, file=paste0("output/Simulation/GIMEXnaive-coxPH-Sim1.1.5.RData"))



# List files matching the specified pattern
matching_files <- list.files(path = "output/Simulation", pattern = "GIMEX-coxPH-Sim1.1.5.*\\.RData$", full.names = TRUE)

Simulation1.1.5_list <- NULL

for (file in matching_files) {
  # Load the file
  loaded_data <- load(file)

  # The load function returns the name of the loaded object, which we can then use to access the object
  data_name <- loaded_data[1]

  # Store the loaded data frame in the data_list
  Simulation1.1.5_list <- c(Simulation1.1.5_list, get(data_name))

  # Cleanup: remove the loaded data frame from the environment
  rm(list = data_name)
}

Simulation1.1.5 <- Simulation1.1.5_list



# str(Simulation1)


Sim1_results <- bind_rows(Simulation1.1.5, .id = "SimID")
results_all <- Sim1_results %>%
  group_by(SimID) %>%
  group_modify(~Summarize_Results(.x, B = 5000, extrapolation = "all"))

load("output/Simulation/GIMEXnaive-coxPH-Sim1.1.5.RData")

Sim1_results_naive <- bind_rows(Simulation1.1.5_naive, .id = "SimID") %>%
  mutate(bestmodel = NA, bestorder = NA, AIC  = NA,  BIC  = NA) %>%
  mutate(extrapo = "naive")




# Attach the true value and calculate the bias

results_all_plot <- rbind(Sim1_results_naive,results_all) %>%
  mutate(truth = factor(para_id, levels = 1:2, labels = c(betatrue)) ) %>%
  mutate(truth = as.numeric(as.character(truth))) %>%
  mutate(bias = (betahat  - truth)) %>%
  mutate(variable_greek = factor(para_id, levels = 1:2, labels = c( expression(beta["con"]),expression(beta["dis"]) ))) %>%
  mutate(extrapo = factor(extrapo,
                          levels = c("naive", "linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal order", "AIC-chosen order", "BIC-chosen order", "AIC ModelAvg", "BIC ModelAvg"),
                          labels = c("naive", "linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal m", "AIC-chosen m", "BIC-chosen m", "AIC ModelAvg", "BIC ModelAvg") ))



# # Summarize the
pdf("output/Simulation/GIMEX-coxPH-Sim1.1.5.pdf", width = 19, height = 10)

boxplot <- results_all_plot %>%
  mutate(modelm = ifelse(extrapo %in% c("CV optimal m", "AIC-chosen m", "BIC-chosen m"), "optimal", modelm)) %>%
  mutate(modelm = ifelse(extrapo %in% c("AIC ModelAvg", "BIC ModelAvg"), "model averaging", modelm)) %>%
  mutate(modelm = ifelse(extrapo == "naive", "naive", modelm)) %>%
  filter(!is.na(extrapo)) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype='dotted', col = '#3F4536')+
  geom_violin(aes(x=extrapo, y=bias, fill = modelm), trim=FALSE) +
  geom_boxplot(aes(x=extrapo, y=bias), width=0.1, position= position_nudge(x=.15)) +
  facet_wrap(~variable_greek,labeller = label_parsed, scales = "free") +
  # coord_cartesian (ylim = c(0, 1.2))
  scale_fill_manual(values=c(   "#27BDBE", "#0096C9", "#B2D235",  "#FFD400", "#F7941D",  "#C768A9", "#305534", "#5D6770", "#9E0918")) +
  theme_bw() +
  theme(text=element_text(size=22), axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
  theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "bottom")  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25)) + 
  labs(x = "Extrapolation", y = "Bias", fill = "m" )

boxplot

dev.off()


pdf("output/Simulation/GIMEX-coxPH-Sim1.1.5_distribution.pdf", width = 14, height = 12)

# Sim1_BS_all$modelm <- factor(Sim1_BS_all$modelm)
# Sim1_BS_all$model_id <- factor(Sim1_BS_all$model_id, levels = c(1:64, paste("best of", 1:6), "CV optimal m"))

count_dist <- results_all_plot %>%
  filter(extrapo == "CV optimal m")  %>%
  group_by(para_id, modelm) %>%
  summarise(count = n()) %>%
  group_by(para_id) %>%
  mutate(per =  count/sum(count)) %>%
  data.frame()

boxplot <- results_all_plot %>%
  filter(extrapo == "CV optimal m")  %>%
  ggplot() +
  geom_bar( data = count_dist, aes(x=modelm, y=per, group = para_id), color = "#3F4536",stat="identity") +
  geom_boxplot(aes(x=modelm, y=bias, color = modelm)) +
  facet_grid(para_id~.,scales = "free") +
  # coord_cartesian (ylim = c(0, 1.2))
  theme_bw() +
  theme(text=element_text(size=20), axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
  theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "none")  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25)) + 
  labs(x = "Extrapolation", y = "Bias", fill = "m" )

boxplot

dev.off()



### Calculate the distribution of the choice of m

results_m <- results_all %>%
  filter( extrapo %in% c("AIC-chosen order", "BIC-chosen order", "CV optimal order")) %>%
  mutate(extrapo = factor(extrapo, 
                          levels = c("AIC-chosen order", "BIC-chosen order", "CV optimal order"), 
                          labels = c("AIC-chosen", "BIC-chosen", "CV optimal"))) %>%
  group_by(extrapo, para_id, modelm) %>%
  summarise(frequency = n(), .groups = 'drop') %>%
  group_by(extrapo, para_id) %>%
  mutate(proportion = frequency / sum(frequency)) %>%
  arrange(extrapo, para_id, modelm) %>%
  ungroup() 

# 1) Summarize across parameters for each method × m
sum_m <- results_m %>%
  mutate(
    extrapo = dplyr::recode(as.character(extrapo),
                            "AIC-chosen" = "AIC",
                            "BIC-chosen" = "BIC",
                            "CV"         = "CV optimal"),
    modelm  = as.numeric(as.character(modelm))
  ) %>%
  group_by(extrapo, modelm) %>%
  summarise(
    median_prop = median(proportion, na.rm = TRUE),  # use median (robust); swap to mean if you prefer
    .groups = "drop"
  ) %>%
  mutate(setting = "low error")


sum_m_full1 <- sum_m

save(sum_m_full1, file=paste0("output/Simulation/GIMEXnaive-coxPH-Sim1.1.5_mchoice.RData"))




### Calculation of the MSE for each scenario

data_getMSE <- results_all_plot %>%
  filter(!is.na(extrapo)) %>%
  mutate(MSE = bias^2) %>%
  group_by(extrapo,para_id) %>%
  summarise(RMSE = sqrt(mean(MSE)))  %>%
  pivot_wider(names_from = para_id, values_from = RMSE)


xtable(data_getMSE, include.rownames = F, digits = 3)

### (!selected) 3.1.1.6 The Model 3.1.1.4 with cross validation --------

## update 2023-10-16:
## Increase B from 500 to 5000



cl <- makeCluster(ncores, outfile=paste0("GSIMEX_COXPH_3.1.1.6", range_eval[1], "-", range_eval[2],".txt"))
registerDoParallel(cl)

betatrue <- c(-0.5, 0.5)


Simulation1.1.6 <- foreach(var1 = range_eval[1]:range_eval[2]) %dopar% {
  start <- proc.time()

  library(pscl)
  library(tidyverse)
  library(survival)
  set.seed(2023)
  seed_i <- sample(1000000, 1000)


  result <- Sim_coxPH_CV(SIM_id = seed_i[var1],
                         beta = betatrue,
                         Mis_mtr = matrix(c(0.7,0.3,0.3,0.7), nrow = 2, byrow = T),
                         sd_error = 2,
                         zeta_m = 5,
                         B = 10000,
                         nfold = 15,
                         CVstep = 4,
                         CVmfold = 3)

  end <- proc.time() - start
  cat(paste0(var1, ": ",end[3],"\n"))

  return(result)
}

stopCluster(cl)

save(Simulation1.1.6, file=paste0("output/Simulation/GIMEX-coxPH-Sim1.1.6", range_eval[1], "-", range_eval[2],".RData"))



Simulation1.1.6_naive <- lapply(1:500, FUN = Sim_coxPH_naive, beta = betatrue,
                                Mis_mtr = matrix(c(0.7,0.3,0.3,0.7), nrow = 2, byrow = T),
                                sd_error = 2)


save(Simulation1.1.6_naive, file=paste0("output/Simulation/GIMEXnaive-coxPH-Sim1.1.6.RData"))



# List files matching the specified pattern
matching_files <- list.files(path = "output/Simulation", pattern = "GIMEX-coxPH-Sim1.1.6.*\\.RData$", full.names = TRUE)

Simulation1.1.6_list <- NULL

for (file in matching_files) {
  # Load the file
  loaded_data <- load(file)

  # The load function returns the name of the loaded object, which we can then use to access the object
  data_name <- loaded_data[1]

  # Store the loaded data frame in the data_list
  Simulation1.1.6_list <- c(Simulation1.1.6_list, get(data_name))

  # Cleanup: remove the loaded data frame from the environment
  rm(list = data_name)
}

Simulation1.1.6 <- Simulation1.1.6_list



# str(Simulation1)


Sim1_results <- bind_rows(Simulation1.1.6, .id = "SimID")
results_all <- Sim1_results %>%
  group_by(SimID) %>%
  group_modify(~Summarize_Results(.x, B = 5000, extrapolation = "all"))

load("output/Simulation/GIMEXnaive-coxPH-Sim1.1.6.RData")

Sim1_results_naive <- bind_rows(Simulation1.1.6_naive, .id = "SimID") %>%
  mutate(bestmodel = NA, bestorder = NA, AIC  = NA,  BIC  = NA) %>%
  mutate(extrapo = "naive")




# Attach the true value and calculate the bias

results_all_plot <- rbind(Sim1_results_naive,results_all) %>%
  mutate(truth = factor(para_id, levels = 1:2, labels = c(betatrue)) ) %>%
  mutate(truth = as.numeric(as.character(truth))) %>%
  mutate(bias = (betahat  - truth)) %>%
  mutate(variable_greek = factor(para_id, levels = 1:2, labels = c( expression(beta["con"]),expression(beta["dis"]) ))) %>%
  mutate(extrapo = factor(extrapo,
                          levels = c("naive", "linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal order", "AIC-chosen order", "BIC-chosen order", "AIC ModelAvg", "BIC ModelAvg"),
                          labels = c("naive", "linear", "quadratic", paste0("best subset of m=",1:6), "CV optimal m", "AIC-chosen m", "BIC-chosen m", "AIC ModelAvg", "BIC ModelAvg") ))


# # Summarize the
pdf("output/Simulation/GIMEX-coxPH-Sim1.1.6.pdf", width = 19, height = 10)

boxplot <- results_all_plot %>%
  mutate(modelm = ifelse(extrapo %in% c("CV optimal m", "AIC-chosen m", "BIC-chosen m"), "optimal", modelm)) %>%
  mutate(modelm = ifelse(extrapo %in% c("AIC ModelAvg", "BIC ModelAvg"), "model averaging", modelm)) %>%
  mutate(modelm = ifelse(extrapo == "naive", "naive", modelm)) %>%
  filter(!is.na(extrapo)) %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype='dotted', col = '#3F4536')+
  geom_violin(aes(x=extrapo, y=bias, fill = modelm), trim=FALSE) +
  geom_boxplot(aes(x=extrapo, y=bias), width=0.1, position= position_nudge(x=.15)) +
  facet_wrap(~variable_greek,labeller = label_parsed, scales = "free") +
  # coord_cartesian (ylim = c(0, 1.2))
  scale_fill_manual(values=c(   "#27BDBE", "#0096C9", "#B2D235",  "#FFD400", "#F7941D",  "#C768A9", "#305534","#5D6770", "#9E0918")) +
  theme_bw() +
  theme(text=element_text(size=22), axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
  theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "bottom")  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25)) + 
  labs(x = "Extrapolation", y = "Bias", fill = "m" )

boxplot

dev.off()


### Calculate the distribution of the choice of m

results_m <- results_all %>%
  filter( extrapo %in% c("AIC-chosen order", "BIC-chosen order", "CV optimal order")) %>%
  mutate(extrapo = factor(extrapo, 
                          levels = c("AIC-chosen order", "BIC-chosen order", "CV optimal order"), 
                          labels = c("AIC-chosen", "BIC-chosen", "CV optimal"))) %>%
  group_by(extrapo, para_id, modelm) %>%
  summarise(frequency = n(), .groups = 'drop') %>%
  group_by(extrapo, para_id) %>%
  mutate(proportion = frequency / sum(frequency)) %>%
  arrange(extrapo, para_id, modelm) %>%
  ungroup() 

# 1) Summarize across parameters for each method × m
sum_m <- results_m %>%
  mutate(
    extrapo = dplyr::recode(as.character(extrapo),
                            "AIC-chosen" = "AIC",
                            "BIC-chosen" = "BIC",
                            "CV"         = "CV optimal"),
    modelm  = as.numeric(as.character(modelm))
  ) %>%
  group_by(extrapo, modelm) %>%
  summarise(
    median_prop = median(proportion, na.rm = TRUE),  # use median (robust); swap to mean if you prefer
    .groups = "drop"
  ) %>%
  mutate(setting = "substantial error")


sum_m_full2 <- sum_m

save(sum_m_full2, file=paste0("output/Simulation/GIMEXnaive-coxPH-Sim1.1.6_mchoice.RData"))



load("output/Simulation/GIMEXnaive-coxPH-Sim1.1.5_mchoice.RData")
load("output/Simulation/GIMEXnaive-coxPH-Sim1.1.6_mchoice.RData")

sum_m_full_all <- rbind(sum_m_full1, sum_m_full2) %>%
  # mutate( setting = factor(setting, levels = 1:2, labels = c("low error", "substantial error")) ) %>%
  complete(extrapo, setting, modelm = 1:6, fill = list(median_prop = 0))

pdf("output/Simulation/GIMEX-coxPH-mchoice.pdf", width = 6, height = 3)
ggplot(sum_m_full_all, aes(x = modelm, y = median_prop, color = extrapo, fill = extrapo)) +
  geom_ribbon(aes(ymin = 0, ymax = median_prop),
              alpha = 0.25, colour = NA, position = "identity") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  # scale_x_continuous(breaks = all_m) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = c("AIC" = "#087F8C", "BIC" = "#8BA04E", "CV optimal" = "#9B5678"),
                     name = "Selector") +
  scale_fill_manual(values  = c("AIC" = "#087F8C", "BIC" = "#8BA04E", "CV optimal" = "#9B5678"),
                    name = "Selector") +
  labs(x = "Selected order m", y = "Selection rate across parameters") +
  facet_wrap(~setting) + 
  theme_bw() +
  theme(legend.position  = "bottom",
        panel.border     = element_rect(colour = "#3F4536"),
        strip.background = element_rect(fill = "#3F4536", colour = "#3F4536")) +
  theme(strip.text = element_text(colour = 'white'))
dev.off()







### Calculation of the MSE for each scenario

data_getMSE <- results_all_plot %>%
  filter(!is.na(extrapo)) %>%
  mutate(MSE = bias^2) %>%
  group_by(extrapo,para_id) %>%
  summarise(RMSE = sqrt(mean(MSE)))  %>%
  pivot_wider(names_from = para_id, values_from = RMSE)


xtable(data_getMSE, include.rownames = F, digits = 3)

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
#   result <- Sim_coxPH_unit(SIM_id = seed_i[var1], 
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
# save(Simulation1.2, file=paste0("output/Simulation/GIMEX-coxPH-Sim1.2_Unit.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-coxPH-Sim1.2_Unit.RData"))
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
# pdf("output/Simulation/GIMEX-coxPH-Sim1.2_Unit.pdf", width = 12, height = 6)
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
#   result <- Sim_coxPH_unit(SIM_id = seed_i[var1], 
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
# save(Simulation1.2.1, file=paste0("output/Simulation/GIMEX-coxPH-Sim1.2_Unit.RData"))
# 
# load(file=paste0("output/Simulation/GIMEX-coxPH-Sim1.2_Unit.RData"))
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
# pdf("output/Simulation/GIMEX-coxPH-Sim1.2.1.pdf", width = 12, height = 6)
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