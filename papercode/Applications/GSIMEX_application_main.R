# This file is for implementing the salary algorithm on a data arising from a single cell RNA-seq. 


# For online server purpose only
# pass argument
args = commandArgs(trailingOnly=TRUE)

# for(i in 1:length(args)){
eval(parse(text=args))
# }


range_eval <- replicate
# range_eval    <- as.integer( args )  # read second argument as integer
# range <- c(range_eval[1], range_eval[2])


# 0. Preparation ----------------------------------------------------------

library(tidyverse)
# library(doParallel)
# library(purrr)
library(dplyr)
# library(svMisc)


# 1. Global Parameters  ----------------------------------------------------------
# ncores <- 80

nsample <- 1000

set.seed(2024)
seed_i <- sample(1000000,1000)


source("code/utility.R")


# cl <- makeCluster(ncores, outfile=paste0("GSIMEX_DA_", replicate,".txt"))
# registerDoParallel(cl)

if (replicate>0) {set.seed(seed_i[replicate])} else {set.seed(0)}

# 2. Preparation ----------------------------------------------------------
data_all <- NULL

for (i in 1:20){
 data_i <- read.csv(paste0("data/Query_select_obs_fold_",i,".csv"), header = T)
 data_all <- rbind(data_all, data_i)
}

data_all <- data_all %>%
  mutate(neuron_type  = case_when(
    final_celltype %in% c("Ex", "In") ~ "Neuron",
    TRUE ~ "Non-Neuron"
  )) %>%
  mutate(neuron_type = factor(neuron_type, level = c("Non-Neuron", "Neuron"))) %>%
  mutate( pred_layer_naive = factor(pred_layer, levels = 1:6, labels = c(paste0("Layer ", 1:6))) )  %>%
  mutate( pred_layer_all = factor(pred_layer, levels = 1:7, labels = c(paste0("Layer ", 1:6), "White Matters")) )  %>%
  mutate( atscore = factor(atscore))

if (replicate >0) {
  data_all <- data_all[sample(nrow(data_all), replace = TRUE), ]
}



## 3 Naive Analysis 
  
# # data_control <-   data_all %>% filter(atscore == "A-T-")
# 
# table(data_all$pred_layer_all, data_all$neuron_type, data_all$atscore)
# # table(data_control$pred_layer_str,data_control$neuron_type)
# 
# 
# Model.naive <- glm(neuron_type ~ pred_layer_all*atscore, data = data_all, family = binomial())
# summary(Model.naive)
# 
# table(data_all$atscore)



# 4. GSIMEX method ----------------------------------------------------------
Mis_Cl_Matrix <- read.csv("output/DataAnalysis/misclassification_matrix.csv", row.names="X")



maxq = 6
zeta_m = 2
B = 10000
nfold = 10
CVstep = 2
CVmfold = 3

# coefs_Z <- NULL

data_all_simp <- data_all %>%
  select("atscore", "neuron_type", "pred_layer")

# ** this part modified
coefs_Z <- NULL



for (i in 1:B) {
  # start <- proc.time()
  getseeds <- sqrt(i) * B + replicate + ifelse(replicate > 0, seed_i[replicate],0)
  set.seed(getseeds)
  epsilon <- i/B*zeta_m
  data_all_simp$pred_layer_W <- MCoperator(Mis_Cl_Matrix, zeta = epsilon, data_all_simp$pred_layer, multi = T, start0 = F)
  data_all_simp$pred_layer_W <- as.factor(data_all_simp$pred_layer_W)
  Model_ME <- glm(neuron_type ~ pred_layer_W*atscore, data = data_all_simp, family = binomial())
  # end <- proc.time() - start
  coefs_Z <- rbind(coefs_Z, c(epsilon, unlist(Model_ME$coefficients)))
  # cat(paste0(i, ": seed ", getseeds, "   Time: ", end[3],"s \n"))
  # write.csv(c(),  file = paste0("output/DataAnalysis/time_", i,".csv"))
}


# cat(paste0("Seed: ", getseeds, "   Time: ", end[3],"s \n"))
# coefs_Z <- map_df(SimStep, ~ as_tibble(t(.))) %>% data.frame()


# Extrapolation
q <- maxq
l <- rep(list(0:1),q)
EX <- expand.grid(l, KEEP.OUT.ATTRS =  F)
FullModel_index <- 2^c(0:q)


modelorder <- cut(1:2^q,breaks = FullModel_index, labels = 0:(q-1), include.lowest = FALSE, right = F)
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

for (j in 1:21){
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

EX_results <- data.frame(betahat = betahat_collect,
                      model_id = EX_model_id_collect,
                      para_id = para_index_collect,
                      order = poly_order_collect-1,
                      SSR = SSR,
                      MSE_CV_value = MSE_CV_value,
                      modelm = rep(rowSums(EX)-1,21))

write.csv(EX_results,  file = paste0("output/DataAnalysis/beta_table_", replicate,".csv"))


rm(list = ls())