# 0. Preparation ----------------------------------------------------------

library(tidyverse)
library(dplyr)
library(forcats)


# 1. Global Parameters  ----------------------------------------------------------
nsample <- 1000

set.seed(2024)
seed_i <- sample(1000000,1000)


source("code/utility.R")



# 2. Naive Analysis

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

# if (replicate >0) {
#   data_all <- data_all[sample(nrow(data_all), replace = TRUE), ]
# }

# data_control <-   data_all %>% filter(atscore == "A-T-")

table(data_all$pred_layer_all, data_all$neuron_type, data_all$atscore)
# table(data_control$pred_layer_str,data_control$neuron_type)


Model.naive <- glm(neuron_type ~ pred_layer_all*atscore, data = data_all, family = binomial())
summary(Model.naive)

table(data_all$atscore)

### Try combine the results of the naive methods as in the SIMEX-corrected version to be discussed below

Results_naive <- data.frame(summary(Model.naive)$coefficients) %>%
  select(betahat = Estimate, standard_error = Std..Error) %>%
  mutate( para_id = 1:21) %>%
  mutate( extrapo = "naive")


# 3. Summary of SIMEX  ----------------------------------------------------------


main <- read.csv("output/DataAnalysis/beta_table/beta_table_0.csv") %>%
  select(-X)

main_summary <- Summarize_Results(main, B = 10000) %>%
  arrange(model_id)
dim(main_summary)


boostrap_smr_all <- NULL

for (i in 1:200){
  bootstrap_sample <-  read.csv(paste0("output/DataAnalysis/beta_table/beta_table_", i ,".csv")) %>%
    select(-X)
  bootstrap_summary <- Summarize_Results(bootstrap_sample, B = 10000) %>%
    arrange(model_id) %>%
    mutate(boostrap_id = i)
  boostrap_smr_all <- rbind(boostrap_smr_all, bootstrap_summary)
}


bootstrap_sd <- boostrap_smr_all %>%
  group_by(extrapo, para_id) %>%
  summarize(standard_error = sd(betahat), .groups = 'drop')

# Merge the two dataframes
merged_data <- main_summary %>%
  left_join(bootstrap_sd, by = c("extrapo", "para_id"))


final_data <- bind_rows(merged_data, Results_naive) %>%
  mutate(extrapo = fct_relevel(extrapo, "naive")) %>%
  mutate(para_name = factor(para_id, labels = c("(Intercept)", paste0("L",2:6), "WM", "A+T-", "A+T+", paste0("A+T- x L",2:6), "A+T- x WM", paste0("A+T+ x L",2:6), "A+T+ x WM")))


# Print the first few rows of the merged dataframe
head(final_data)



# Look at optimal choicesen model
final_data %>%
  filter(extrapo  == "AIC-chosen order")

final_data %>%
  filter(extrapo  == "BIC-chosen order")

final_data %>%
  filter(extrapo  == "CV optimal order")

final_data %>%
  filter(extrapo  == "AIC ModelAvg")

### Plot the reesults

# Create the plot
pdf("output/DataAnalysis/GIMEX-DA-results.pdf", width = 18, height = 18)

final_data %>% filter( extrapo %in% c("naive", "linear","quadratic", "CV optimal order", "AIC ModelAvg"))  %>%
  mutate(extrapo = factor(extrapo, levels = c("naive", "linear","quadratic", "CV optimal order", "AIC ModelAvg"), labels = c("Naive", "Linear","Quadratic", "CV optimal order", "AIC ModelAvg")) ) %>%
  mutate( para_group = ifelse(grepl("A\\+T-", para_name), "Interactions with A+T-",
                                      ifelse(grepl("A\\+T\\+", para_name), "Interactions with A+T+", "A-T- (Baseline)"))) %>%
  ggplot(aes(x = para_name, y = betahat)) +
    geom_hline(yintercept = 0, linetype='dotted', col = '#3F4536')+
    geom_errorbar(aes(ymin = betahat - 1.645 * standard_error, ymax = betahat + 1.645 * standard_error), width = 0.2) +
    geom_point(aes(color = para_group), size = 3) +
    facet_grid(extrapo ~ . , scales = "free_y") +
    labs(x = "Parameter",
         y = "Estimate")+
     theme_bw()  +
    scale_color_manual(values=c(   "#27BDBE", "#B2D235","#C768A9")) +
    theme(text=element_text(size=22), axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
    theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
    theme(strip.text = element_text(colour = 'white')) +
    theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "bottom")  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25)) + 
  labs( color = "Parameter Groups" )

dev.off()




pdf("output/DataAnalysis/supplementary_other_order.pdf", width = 28, height = 17)

final_data %>% filter(! extrapo %in% c("naive", "linear","quadratic", "best subset of m=0", "CV optimal order", "AIC ModelAvg"))  %>%
  mutate(extrapo =  factor(extrapo, level =  c(paste0("best subset of m=", rep(1:6)), "AIC-chosen order", "BIC-chosen order", "BIC ModelAvg"))) %>%
  # mutate(extrapo = factor(extrapo, levels = c("naive", "linear","quadratic", "CV optimal order"), labels = c("Naive", "Linear","Quadratic", "CV optimal order")) ) %>%
  mutate( para_group = ifelse(grepl("A\\+T-", para_name), "Interactions with A+T-",
                              ifelse(grepl("A\\+T\\+", para_name), "Interactions with A+T+", "A-T- (Baseline)"))) %>%
  ggplot(aes(x = para_name, y = betahat)) +
  geom_hline(yintercept = 0, linetype='dotted', col = '#3F4536')+
  geom_errorbar(aes(ymin = betahat - 1.645 * standard_error, ymax = betahat + 1.645 * standard_error), width = 0.2) +
  geom_point(aes(color = para_group), size = 3) +
  facet_wrap( ~ extrapo , scales = "free_y") +
  labs(x = "Parameter",
       y = "Estimate")+
  theme_bw()  +
  scale_color_manual(values=c(   "#27BDBE", "#B2D235","#C768A9")) +
  theme(text=element_text(size=22), axis.text = element_text(size = 22),
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
  theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
  theme(strip.text = element_text(colour = 'white')) +
  theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "bottom")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25)) + 
  labs( color = "Parameter Groups" )

dev.off()


### Plot the order chosen
optimal_models <- final_data %>%
  filter(extrapo  %in% c("AIC-chosen order", "BIC-chosen order", "CV optimal order")) %>%
  mutate(extrapo = factor(extrapo, 
                          levels = c("AIC-chosen order", "BIC-chosen order", "CV optimal order"), 
                          labels = c("AIC-chosen", "BIC-chosen", "CV optimal")))
summary(optimal_models)

### Add the bootstrap distribution choice

proportion_table <- boostrap_smr_all %>%
  filter(extrapo  %in% c("AIC-chosen order", "BIC-chosen order", "CV optimal order")) %>%
  mutate(extrapo = factor(extrapo, 
                          levels = c("AIC-chosen order", "BIC-chosen order", "CV optimal order"), 
                          labels = c("AIC-chosen", "BIC-chosen", "CV optimal"))) %>%
  group_by(extrapo, para_id, modelm) %>%
  summarise(frequency = n(), .groups = 'drop') %>%
  group_by(extrapo, para_id) %>%
  mutate(proportion = frequency / sum(frequency)) %>%
  arrange(extrapo, para_id, modelm) %>%
  ungroup() 


merged_data <- proportion_table  %>%
  left_join(optimal_models, by = c("extrapo", "para_id", "modelm")) %>%
  mutate(para_name = factor(para_id, labels = c("(Intercept)", paste0("L",2:6), "WM", "A+T-", "A+T+", paste0("A+T- x L",2:6), "A+T- x WM", paste0("A+T+ x L",2:6), "A+T+ x WM")))




# Create the dot plot
pdf("output/DataAnalysis/GIMEX-dot-plot.pdf", width = 10, height = 8)

# ggplot(optimal_models, aes(x = para_name, y = modelm)) +
#   geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 2) +
#   facet_grid(extrapo ~ . ) +
#   labs(x = "parameter name",
#        y = "choice of m") +
#   theme_bw() +
#   theme(text=element_text(size=22), axis.text = element_text(size = 22),
#         axis.text.x = element_text(angle = 20, vjust = 1, hjust=1), panel.spacing = unit(1, "lines")) +
#   theme(strip.background =element_rect(fill="#3F4536",color="#3F4536"))+ # #535b44
#   theme(strip.text = element_text(colour = 'white')) +
#   theme(panel.border = element_rect(colour = "#3F4536"), legend.position = "bottom")  +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(merged_data, aes(x = para_name)) +
  # Layer for the proportions
  geom_point(aes(y = as.numeric(modelm), size = proportion), alpha = 0.2, color = "#5D6770", fill = "#5D6770") +
  # Layer for the original dot plot
  geom_point(data = optimal_models, aes(y = as.numeric(modelm)), shape = 'x', color = "#9E0918", size = 5) +
  facet_grid(extrapo ~ .) +
  labs(x = "Parameter Name",
       y = "Choice of m") +
  theme_bw() +
  theme(text = element_text(size = 15), 
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = "#3F4536", color = "#3F4536"),
        strip.text = element_text(colour = 'white'),
        panel.border = element_rect(colour = "#3F4536"),
        legend.position = "bottom") +
  scale_shape_identity() +
  scale_size_continuous(range = c(1, 10), guide = FALSE) # Adjust the range for better visibility

dev.off()