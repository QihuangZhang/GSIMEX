# Here it's the misclassification matrix this on the results from the paper of celery


# 0. Preparation ----------------------------------------------------------
library(ggplot2)
library(reshape2)
library(dplyr)

# 1. Data Loading ---------------------------------------------------------

CeLEry_results <- read.csv("data/dataSectionMobs.csv", header = T) %>%
  select( c("Layer", "pred_layer")) %>%
  mutate( Layer = factor(Layer, levels = 1:7, labels = c(paste0("Layer ", 1:6), "White Matter" )) ) %>%
  mutate( pred_layer = factor(pred_layer, levels = 1:7, labels = c(paste0("Layer ", 1:6), "White Matter" )) ) 
  
# names(CeLEry_results) <- c("True Layer", "Predicted Layer", paste0("Layer ", 1:6), "White Matter") 

Twoway_table <- table(CeLEry_results$Layer, CeLEry_results$pred_layer)
row_sums <- rowSums(Twoway_table)
normalized_table <- sweep(Twoway_table, 1, row_sums, "/")
# normalized_table <- Twoway_table / total_sum


### This is a past idea which model the misclassification for each indicator 
# # 2. Obtain the 2x2 miclassification table
# # Function to create 2x2 tables
# generate_2x2_table <- function(A, row, col) {
#   n <- nrow(A)
#   
#   # Element itself
#   element <- A[row, col]
#   
#   # Sum of the remaining elements in the same row
#   row_sum <- sum(A[row, ]) - element
#   
#   # Sum of the remaining elements in the same column
#   col_sum <- sum(A[, col]) - element
#   
#   # Sum of all other elements (excluding the given row and column)
#   total_sum <- sum(A) - sum(A[row, ]) - sum(A[, col]) + element
#   
#   # Constructing the 2x2 matrix
#   table_2x2 <- matrix(c(element, row_sum, col_sum, total_sum), nrow = 2)
#   colnames(table_2x2) <- c("Element", "Row Sum")
#   rownames(table_2x2) <- c("Element", "Column Sum")
#   
#   table_2x2_rotate <- t(table_2x2)
#   
#   return( sweep(table_2x2_rotate, 1, rowSums(table_2x2_rotate), "/") )
# }
# 
# # Example: Generating 2x2 table for the cell (1,1)
# table_1_1 <- generate_2x2_table(normalized_table, 1, 1)
# print("2x2 Table for Cell (1,1):")
# print(table_1_1)
# 
# 
# 
# # Generate and print 2x2 tables for cells (2,2) to (7,7)
# for (i in 2:7) {
#   cat(sprintf("\n2x2 Table for Cell (%d,%d):\n", i, i))
#   print(generate_2x2_table(A, i, i))
# }



### Visualization

# Converting matrix to data frame for ggplot

df <- melt(normalized_table) %>%
  mutate(Var1 = factor(Var1, 7:1, rev(c(paste0("Layer ", 1:6), "White Matter")))) %>%
  mutate(Var2 = factor(Var2, 1:7, c(paste0("Layer ", 1:6), "White Matter"))) 


df <- melt(normalized_table) %>%
  mutate(Var1 = factor(Var1, level = rev(c(paste0("Layer ", 1:6), "White Matter")))) %>%
  mutate(Var2 = factor(Var2)) 

names(df) <- c("Actual", "Predicted", "Count")

pdf("output/DataAnalysis/Miclassification.pdf", width = 7, height = 5)

ggplot(df, aes(x = Predicted, y = Actual, fill = Count)) +
  geom_tile(color = "white") +  # add tiles with white borders
  geom_text(aes(label = sprintf("%.3f", Count)), vjust = 1, color = "white") +  # add text to tiles
  scale_fill_gradient(low = "#087F8C", high = "#D3674A") +  # color gradient
  theme_minimal() +  # minimal theme
  labs(title = "Misclassification Matrix Heatmap", fill = "Frequency") 
dev.off()

misclassification_matrix <- as.matrix(normalized_table)
write.csv(misclassification_matrix, file = "output/DataAnalysis/misclassification_matrix.csv")

