library(terra)
library(usdm)
library(sp)
library(corrplot)   
library(RColorBrewer)

# Loading the data
data_2040 <- rast("./clipped_masked_bioclim/bioclim_2021_2040_clipped_masked.tif") 
data_2060 <- rast("./clipped_masked_bioclim/bioclim_2041_2060_clipped_masked.tif")
elevation <-  rast("./drivers/vnp_elevation.tif")
slope <-  rast("./drivers/vnp_slope.tif")

plot(data_2040[[1]])

names(elevation) <- "elevation"
names(slope)     <- "slope"

stack_2040 <- data_2040
stack_2060 <- data_2060

# Correlation matrix plot function
plot_corr_matrix <- function(raster_stack, title, output_path) {
  
  raster_df  <- as.data.frame(raster_stack, na.rm = TRUE)
  cor_matrix <- cor(raster_df, method = "pearson")
  
  png(output_path, width = 1200, height = 1100, res = 130)
  corrplot(cor_matrix,
           method      = "color",        
           type        = "upper",        
           addCoef.col = "black",        
           number.cex  = 0.65,          
           tl.col      = "black",        
           tl.srt      = 45,            
           tl.cex      = 0.8,           
           col         = colorRampPalette(c("#053061", "#2166AC", "#F7F7F7",
                                            "#D6604D", "#67001F"))(200),
           title       = title,
           mar         = c(0, 0, 2, 0),
           diag        = TRUE)          
  dev.off()
  
  return(cor_matrix)
}

# Plot and save for both periods
cor_2040 <- plot_corr_matrix(
  stack_2040,
  title       = "Pearson r — Bioclim (2021–2040)",
  output_path = "./refined_figures/corr_matrix_2040.png"
)

cor_2060 <- plot_corr_matrix(
  stack_2060,
  title       = "Pearson r — Bioclim (2041–2060)",
  output_path = "./refined_figures/corr_matrix_2060.png"
)

# Flag all pairs with |r| > 0.7
flag_high_correlation <- function(cor_matrix, threshold = 0.7) {
  
  cor_long <- as.data.frame(as.table(cor_matrix))
  names(cor_long) <- c("Var1", "Var2", "r")
  
  # Keep only upper triangle, exclude self-correlations
  high_cor <- subset(cor_long, 
                     as.character(Var1) < as.character(Var2) & abs(r) >= threshold)
  high_cor <- high_cor[order(-abs(high_cor$r)), ]
  rownames(high_cor) <- NULL
  return(high_cor)
}

cat("\n--- Pairs with |r| >= 0.7 in 2040 stack ---\n")
print(flag_high_correlation(cor_2040, threshold = 0.7))

cat("\n--- Pairs with |r| >= 0.7 in 2060 stack ---\n")
print(flag_high_correlation(cor_2060, threshold = 0.7))

# VIF selection
run_vif_thresholds <- function(raster_stack, thresholds) {
  
  raster_df    <- as.data.frame(raster_stack, na.rm = TRUE)
  results_list <- list()
  
  for (th in thresholds) {
    cat("Running VIF with threshold =", th, "\n")
    vif_step_result <- vifstep(raster_df, th = th)
    selected_vars   <- vif_step_result@results$Variables
    results_list[[paste0("VIF_", th)]] <- list(
      threshold          = th,
      selected_variables = selected_vars,
      n_selected         = length(selected_vars)
    )
  }
  return(results_list)
}

thresholds   <- 2:10
results_2040 <- run_vif_thresholds(stack_2040, thresholds)
results_2060 <- run_vif_thresholds(stack_2060, thresholds)

summarize_vif_results <- function(results_list) {
  data.frame(
    Threshold  = sapply(results_list, function(x) x$threshold),
    N_Selected = sapply(results_list, function(x) x$n_selected),
    Variables  = sapply(results_list, function(x)
      paste(x$selected_variables, collapse = ", "))
  )
}

summary_2040 <- summarize_vif_results(results_2040)
summary_2060 <- summarize_vif_results(results_2060)

print(summary_2040)
print(summary_2060)

write.csv(summary_2040, "./clipped_masked_bioclim/vif_results/VIF_threshold_comparison_2021_2040.csv")
write.csv(summary_2060, "./clipped_masked_bioclim/vif_results/VIF_threshold_comparison_2041_2060.csv")

best_threshold     <- 3
selected_vars_2040 <- results_2040[[paste0("VIF_", best_threshold)]]$selected_variables
selected_vars_2060 <- results_2060[[paste0("VIF_", best_threshold)]]$selected_variables

cat("\nSelected variables for 2040:", paste(selected_vars_2040, collapse = ", "), "\n")
cat("Selected variables for 2060:", paste(selected_vars_2060, collapse = ", "), "\n")

cat("\nElevation retained in 2040?", "elevation" %in% selected_vars_2040, "\n")
cat("Slope     retained in 2040?", "slope"     %in% selected_vars_2040, "\n")
cat("Elevation retained in 2060?", "elevation" %in% selected_vars_2060, "\n")
cat("Slope     retained in 2060?", "slope"     %in% selected_vars_2060, "\n")

data_2040_final <- stack_2040[[selected_vars_2040]]
data_2060_final <- stack_2060[[selected_vars_2060]]

# Only remain with bio3, bio4, bio15, bio19 which are the same 
# for both periods and remained constant after VIF selection with threshold 3
data_2040_new <- data_2040_final[[c("bio3", "bio4", "bio15", "bio19")]]
data_2060_new <- data_2060_final[[c("bio3", "bio4", "bio15", "bio19")]]

writeRaster(data_2040_new, "./clipped_masked_bioclim/selected_var/wc2.1_30s_bioc_BCC-CSM2-MR_ssp126_2021-2040_selected_variable_vif3.tif", overwrite = TRUE)
writeRaster(data_2060_new, "./clipped_masked_bioclim/selected_var/wc2.1_30s_bioc_BCC-CSM2-MR_ssp126_2041-2060_selected_variable_vif3.tif", overwrite = TRUE)

