## Test-Retest Correlation
# This script is part of the study "Robustness of EEG functional brain networks 
# associated with fluid intelligence: A multiverse analysis of connectivity and 
# thresholding methods". It calculates test-retest correlations for graph-
# theoretical measures across resting-state measurement time points and 
# conditions.
#
# Written by: Christoph Fruehlinger
# Last edit: July 2026

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, cowplot)

options(scipen = 999)
rm(list = ls())
set.seed(42)

# Create output folder
savepath = "Results/"

# Load data
datapath         <- "Data/Connectivity/Graph_data.csv"
main_path_conn   <- "imcoh"
main_path_thresh <- "dens"
data <- read.csv(datapath)

data <- data %>%
  mutate(SetCond = paste(Run, Condition, sep = "_"))

# Select main path columns
pattern        <- paste('[a-z0-9]+', main_path_conn, '[a-z0-9]+', main_path_thresh, sep = '_')
main_path_cols <- grep(pattern, colnames(data), value = TRUE)
main_path_data <- data[, c('ID', 'SetCond', main_path_cols)]

# Set conditions to compare
main_cond <- "first_eyes_closed"
otherconds <- c("first_eyes_open", "second_eyes_closed", "third_eyes_closed")
condnames <- c("Run 1 EO", "Run 2 EC", "Run 3 EC")

bands <- c("delta", "theta", "alpha1", "alpha2", "beta")
feat_vars <- c("cc", "pathl", "eglob", "eloc", "smallworld")
featurenames <- c("Clustering Coefficient", "Characteristic Path Length",
                  "Global Efficiency", "Local Efficiency", "Small-Worldness")

# Collect correlations
p_vals <- list()
plot_list <- list()
correlations <- tibble()

# ------------------------------------------------------------------------------
# Correlation Analysis
# ------------------------------------------------------------------------------

# loop over conditions and sets
for (cond2 in otherconds) {
  
  for (band in bands) {
  
    # filter
    T1 <- main_path_data %>%
      filter(SetCond == main_cond) %>%
      select(ID, contains(band))
    
    T2 <- main_path_data %>%
      filter(SetCond == cond2) %>%
      select(ID, contains(band))
    
    # sort IDs
    common_ids <- intersect(T1$ID, T2$ID)
    T1 <- T1 %>%
      filter(ID %in% common_ids) %>%
      arrange(ID)
    T2 <- T2 %>%
      filter(ID %in% common_ids) %>%
      arrange(ID)
    
    if (nrow(T1) == 0 || nrow(T2) == 0) next
    
    stopifnot(identical(T1$ID, T2$ID))
    
    # calculate correlations
    rho_vec <- numeric(length(feat_vars))
    p_vec <- numeric(length(feat_vars))
    
    for (i in seq_along(feat_vars)) {
      col_name <- paste(feat_vars[i], main_path_conn, band, main_path_thresh, sep = "_")
      x <- T1[[col_name]]
      y <- T2[[col_name]]
      
      test <- cor.test(x, y, method = "spearman", exact = FALSE)
      rho_vec[i] <- test$estimate
      p_vec[i] <- test$p.value
    }
    
    p_vals[[paste(cond2, band, sep = "_")]] <- p_vec
    
    # save results
    corr_tbl <- tibble(
      Band = band,
      Condition1 = main_cond,
      Condition2 = cond2,
      Feature = feat_vars,
      SpearmanRho = rho_vec,
      pValue = p_vec
    )
    
    correlations <- bind_rows(
      correlations,
      corr_tbl %>%
        mutate(Comparison = paste(cond2, band, sep = "_"))
    )
    
    correlations <- correlations %>%
      mutate(pValue_adj = p.adjust(pValue, method = "holm"))
    
  }
}

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

correlations$Band <- as.factor(correlations$Band)
correlations$Band <- factor(correlations$Band, 
                            levels = c("delta", "theta", "alpha1", "alpha2",
                                       "beta"))
correlations$Condition1 <- as.factor(correlations$Condition1)
correlations$Condition2 <- as.factor(correlations$Condition2)
correlations$Feature <- as.factor(correlations$Feature)
correlations$Feature <- factor(correlations$Feature,
                               levels = c("cc", "pathl", "eglob", "eloc", 
                                          "smallworld"))

heatmap_plot_list <- list()

for (cond2 in otherconds) {
  
  heatmap_data <- correlations %>%
    filter(Condition1 == main_cond,
           Condition2 == cond2) %>%
    mutate(Feature = recode(Feature,
                            cc = 'Clustering Coefficient',
                            pathl = 'Characteristic Path Length',
                            eglob =  'Global Efficiency',
                            eloc = 'Local Efficiency',
                            smallworld = "Small-Worldness"),
           Band = recode(Band,
                         delta = 'Delta',
                         theta = 'Theta',
                         alpha1 = 'Alpha-1',
                         alpha2 = 'Alpha-2',
                         beta = 'Beta'))
  
  limits = c(0, .75)
  
  heatmap_plot_list[[cond2]] <- ggplot(heatmap_data, aes(x = Feature, y = Band, 
                                                              fill = SpearmanRho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(SpearmanRho, 2)), color = "black", size = 4) +
    labs(title = condnames[match(cond2, otherconds)]) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA)) +
    scale_fill_distiller(palette   = "Spectral", direction = -1, 
                         name = "Correlation", limits = limits)
  
}

retest_grid <- plot_grid(plotlist = heatmap_plot_list[1:3], 
                         labels = c("A", "B", "C"), ncol = 3)

retest_filename <- paste0(savepath, 'retest_correlations.png')

ggsave(filename = retest_filename, plot = retest_grid, width = 16, 
       height = 4, dpi = 600, bg = "white")
