## Statistical Analysis
# This is the main statistical analysis script for the study "Robustness of EEG
# functional brain networks associated with fluid intelligence: A multiverse
# analysis of connectivity and thresholding methods"
#
# Written by: Christoph Fruehlinger
# Last edit: August 2026

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Load packages 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, tidyverse, future, future.apply, progressr)

options(scipen = 999)
rm(list = ls())
set.seed(42)

# Load data

# Intelligence data
IST_datapath <- "Data/IST_fulltable.csv"
IST_data <- read.csv(IST_datapath, sep = ',')
IST_data$gf_score <- rowSums(IST_data[2:21])
gf_data <- IST_data %>% 
  dplyr::select(ID, gf_score)
IST_data$gc_score <- rowSums(IST_data[22:105])
gc_data <- IST_data %>% 
  dplyr::select(ID, gc_score)

# BFI data
BFI_coding <- read.csv("Data/Personality_Items_coding.csv", sep = ';')
new_names <- c("ID", BFI_coding$Subskala2)
BFI_scores <- read.csv("Data/Personality_Items.csv", col.names = new_names)

# Agreeableness
BFI_scores$Agreeableness <- BFI_scores %>% 
  dplyr::select(contains(c("Compassion",
                           "Respectfulness",
                           "Trust"))) %>%
  rowMeans()

# Conscientiousness
BFI_scores$Conscientiousness <- BFI_scores %>% 
  dplyr::select(contains(c("Organization",
                           "Productiveness",
                           "Responsibility"))) %>%
  rowMeans()

# Extraversion
BFI_scores$Extraversion <- BFI_scores %>% 
  dplyr::select(contains(c("Sociability",
                           "Assertiveness",
                           "EnergyLevel"))) %>%
  rowMeans()

# Neuroticism
BFI_scores$Neuroticism <- BFI_scores %>% 
  dplyr::select(contains(c("Anxiety",
                           "Depression",
                           "EmotionalVolatility"))) %>%
  rowMeans()

# Openness
BFI_scores$Openness <- BFI_scores %>% 
  dplyr::select(contains(c("AestheticSensitivity",
                           "IntellectualCuriosity",
                           "CreativeImagination"))) %>%
  rowMeans()

# Lab data
lab_datapath <- "Data/Labs.txt"
lab_data <- read_tsv(lab_datapath)
lab_data <- lab_data[1:2]
colnames(lab_data)[2] <- "Lab"

# Graph metrics
datapath         <- "Data/Connectivity/Graph_data.csv"
main_path_conn   <- "imcoh"
main_path_thresh <- "dens"
data <- read.csv(datapath)

eco_smallworld_cols <- grep("smallworld_.*_eco$", colnames(data), value = TRUE)

data <- data %>%
  mutate(across(all_of(eco_smallworld_cols),
                ~ ifelse(.x == 0 | .x > 1, NA, .x)))

# Select main path columns
pattern        <- paste('[a-z0-9]+', main_path_conn, '[a-z0-9]+', main_path_thresh, sep = '_')
main_path_cols <- grep(pattern, colnames(data), value = TRUE)
main_path_data <- data[, c('ID', 'Run', 'Condition', main_path_cols)]

# Hypothesis 1 data
main_path_data <- main_path_data %>%
  filter(Run == 'first' & Condition == 'eyes_closed')
main_path_data <- left_join(main_path_data, gf_data, by = "ID") %>% 
  left_join(., gc_data, by = "ID") %>% 
  left_join(., BFI_scores[,c("ID", "Agreeableness", "Conscientiousness", "Extraversion", 
                             "Neuroticism", "Openness")], by = "ID") %>% 
  left_join(., lab_data, by = "ID")

# Hypothesis 2 data
full_data <- data %>%
  filter(Run == 'first' & Condition == 'eyes_closed') %>%
  left_join(gf_data, by = "ID") %>% 
  left_join(., gc_data, by = "ID") %>% 
  left_join(., BFI_scores[,c("ID", "Agreeableness", "Conscientiousness", "Extraversion", 
                             "Neuroticism", "Openness")], by = "ID") %>% 
  left_join(., lab_data, by = "ID")

# ------------------------------------------------------------------------------
## Hypothesis 1
# ------------------------------------------------------------------------------

# Analysis parameters
bands    <- c("delta", "theta", "alpha1", "alpha2", "beta")
measures <- c("cc", "pathl", "eglob", "eloc", "smallworld")
k        <- 5000
dep_vars <- c("gf_score", "gc_score", "Agreeableness", "Conscientiousness",
              "Extraversion", "Neuroticism", "Openness")

graph_cols_main <- main_path_cols
graph_cols_full <- grep(paste(measures, collapse = "|"), colnames(data), value = TRUE)

# Standardize independent variables
main_path_data_base <- main_path_data %>%
  mutate(across(all_of(graph_cols_main), ~ as.numeric(scale(.x))))

full_data_base <- full_data %>%
  mutate(across(all_of(graph_cols_full), ~ as.numeric(scale(.x))))

analysis_grid <- expand.grid(band = bands, measure = measures,
                             stringsAsFactors = FALSE)

# Parallel processing
n_cores <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_cores)
message("Using ", n_cores, " cores for parallel processing")

# Helper functions
permute_within_lab <- function(data, dep_var) {
  data %>%
    group_by(Lab) %>%
    mutate(across(all_of(dep_var), ~ sample(.x))) %>%
    ungroup()
}

get_t_value <- function(model, predictor) {
  coef(summary(model))[predictor, "t value"]
}

run_permutation_test <- function(dep_var, band, measure, data, conn, thresh, k, p) {
  
  col_name <- paste(measure, conn, band, thresh, sep = "_")
  p(message = paste0(dep_var, " | ", measure, " | ", band))
  
  if (!col_name %in% colnames(data)) {
    warning(paste("Column not found, skipping:", col_name))
    return(data.frame(dep_var = dep_var, band = band, measure = measure,
                      t_obs = NA, p_perm = NA))
  }
  
  formula <- as.formula(paste(dep_var, "~", col_name, "+ (1 | Lab)"))
  
  # Observed model
  fit_obs <- lmer(formula, data = data, REML = TRUE,
                  control = lmerControl(optimizer = "bobyqa"))
  t_obs   <- get_t_value(fit_obs, col_name)
  
  # Permutation distribution
  t_perm <- replicate(k, {
    data_perm <- permute_within_lab(data, dep_var)
    fit_perm  <- lmer(formula, data = data_perm, REML = TRUE,
                      control = lmerControl(optimizer = "bobyqa"))
    get_t_value(fit_perm, col_name)
  })
  
  p_perm <- mean(abs(t_perm) >= abs(t_obs))
  
  data.frame(dep_var = dep_var, band = band, measure = measure,
             t_obs = t_obs, p_perm = p_perm)
}

results_list <- list()

for (dep_var in dep_vars) {
  
  message("\n", strrep("=", 50))
  message("Dependent variable: ", dep_var)
  message(strrep("=", 50), "\n")
  
  # Standardize dependent variable
  main_path_data_scaled <- main_path_data_base %>%
    mutate(across(all_of(dep_var), ~ as.numeric(scale(.x))))
  
  handlers(handler_progress(
    format   = "[:bar] :percent | :current/:total combinations running",
    width    = 70,
    complete = "="
  ))
  
  message("Starting permutation tests (k = ", k, ") across ",
          nrow(analysis_grid), " band x measure combinations...")
  
  results_list[[dep_var]] <- with_progress({
    p <- progressor(steps = nrow(analysis_grid))
    
    future_mapply(
      FUN      = run_permutation_test,
      dep_var  = dep_var,
      band     = analysis_grid$band,
      measure  = analysis_grid$measure,
      MoreArgs = list(
        data   = main_path_data_scaled,
        conn   = main_path_conn,
        thresh = main_path_thresh,
        k      = k,
        p      = p
      ),
      SIMPLIFY    = FALSE,
      future.seed = 42
    )
  })
}

# FDR-Korrektur per dependent variable and Band
results_df <- map_dfr(dep_vars, function(dv) {
  bind_rows(results_list[[dv]]) %>%
    group_by(band) %>%
    mutate(p_fdr = p.adjust(p_perm, method = "fdr")) %>%
    ungroup()
}) %>%
  mutate(
    significant = p_fdr < 0.01,
    dep_var = factor(dep_var, levels = dep_vars),
    band    = factor(band,    levels = bands),
    measure = factor(measure, levels = measures)
  ) %>%
  arrange(dep_var, band, measure)

# Print results
print(results_df, n = nrow(results_df))

results_main <- results_df %>%
  filter(dep_var == "gf_score")

write.csv(results_main, "Results/H1_mainpath_results.csv", row.names = FALSE)

results_exploratory <- results_df %>%
  filter(dep_var != "gf_score")

write.csv(results_exploratory, "Results/H1_exploratory_results.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# Hypothesis 2
# ------------------------------------------------------------------------------

conn_measures  <- c("imcoh", "wpli", "pli", "pcoh", "oaec")
thresh_methods <- c("dens", "omst", "eco", "mcc")

multiverse_grid <- expand.grid(
  conn    = conn_measures,
  thresh  = thresh_methods,
  band    = bands,
  measure = measures,
  stringsAsFactors = FALSE
)

# Fit multiverse models
fit_multiverse_model <- function(dep_var, conn, thresh, band, measure, data) {
  
  col_name <- paste(measure, conn, band, thresh, sep = "_")
  
  if (!col_name %in% colnames(data)) {
    return(data.frame(
      dep_var  = dep_var,
      conn     = conn,
      thresh   = thresh,
      band     = band,
      measure  = measure,
      col_name = col_name,
      beta     = NA,
      ci_lower = NA,
      ci_upper = NA,
      fail_reason = "missing_column"
    ))
  }
  
  fit <- tryCatch(
    lmer(as.formula(paste(dep_var, "~", col_name, "+ (1 | Lab)")),
         data = data, REML = TRUE, 
         control = lmerControl(optimizer = "bobyqa")),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(data.frame(
      dep_var  = dep_var,
      conn     = conn,
      thresh   = thresh,
      band     = band,
      measure  = measure,
      col_name = col_name,
      beta     = NA,
      ci_lower = NA,
      ci_upper = NA,
      fail_reason = paste0("fit_error: ", conditionMessage(fit))
    ))
  }
  
  fit_summary <- coef(summary(fit))
  
  if (!col_name %in% rownames(fit_summary)) {
    return(data.frame(
      dep_var  = dep_var,
      conn     = conn,
      thresh   = thresh,
      band     = band,
      measure  = measure,
      col_name = col_name,
      beta     = NA,
      ci_lower = NA,
      ci_upper = NA,
      fail_reason = "coef_missing"
    ))
  }
  
  beta_std <- fixef(fit)[col_name]
  se       <- fit_summary[col_name, "Std. Error"]
  ci_std   <- c(beta_std - 1.96 * se,
                beta_std + 1.96 * se)
  
  data.frame(
    dep_var  = dep_var,
    conn     = conn,
    thresh   = thresh,
    band     = band,
    measure  = measure,
    col_name = col_name,
    beta     = beta_std,
    ci_lower = ci_std[1],
    ci_upper = ci_std[2],
    fail_reason = NA_character_
  )
}

# ------------------------------------------------------------------------------
# Run multiverse models for all dep_vars
# ------------------------------------------------------------------------------

handlers(handler_progress(
  format   = "[:bar] :percent | :current/:total combinations running",
  width    = 70,
  complete = "="
))

multiverse_list <- list()

for (dep_var in dep_vars) {
  
  message("\n", strrep("=", 50))
  message("Multiverse — Dependent variable: ", dep_var)
  message(strrep("=", 50), "\n")
  
  # Standardize dependent variable
  full_data_scaled <- full_data_base %>%
    mutate(across(all_of(dep_var), ~ as.numeric(scale(.x))))
  
  message("Starting multiverse analysis (", nrow(multiverse_grid), " models)...")
  
  multiverse_list[[dep_var]] <- with_progress({
    p <- progressor(steps = nrow(multiverse_grid))
    
    future_mapply(
      FUN = function(conn, thresh, band, measure) {
        p(message = paste0(measure, " | ", conn, " | ", thresh, " | ", band))
        fit_multiverse_model(
          dep_var = dep_var,
          conn    = conn,
          thresh  = thresh,
          band    = band,
          measure = measure,
          data    = full_data_scaled
        )
      },
      conn    = multiverse_grid$conn,
      thresh  = multiverse_grid$thresh,
      band    = multiverse_grid$band,
      measure = multiverse_grid$measure,
      SIMPLIFY    = FALSE,
      future.seed = 42
    )
  })
}

# Reset to sequential processing
plan(sequential)

# Combine all results
multiverse_df <- map_dfr(dep_vars, function(dv) {
  bind_rows(multiverse_list[[dv]])
}) %>%
  mutate(
    dep_var = factor(dep_var, levels = dep_vars),
    band    = factor(band,    levels = bands),
    measure = factor(measure, levels = measures)
  )

# ------------------------------------------------------------------------------
# Robustness assessment
# ------------------------------------------------------------------------------

# Helper function
ci_overlap_proportion <- function(lower, upper){
  keep <- complete.cases(lower, upper)
  lower <- lower[keep]
  upper <- upper[keep]
  if(length(lower) < 2){
    return(NA_real_)
  }
  overlap_matrix <- outer(
    seq_along(lower),
    seq_along(lower),
    Vectorize(function(i, j){
      lower[i] <= upper[j] &&
        upper[i] >= lower[j]
    })
  )
  diag(overlap_matrix) <- NA
  mean(rowMeans(overlap_matrix, na.rm = TRUE))
}

robustness_df <- multiverse_df %>%
  group_by(dep_var, band, measure) %>%
  summarise(
    n_specs              = sum(!is.na(beta)),
    prop_positive        = mean(beta > 0, na.rm = TRUE),
    direction_consistent = prop_positive >= 0.75 | prop_positive <= 0.25,
    mean_beta            = mean(beta,   na.rm = TRUE),
    median_beta          = median(beta, na.rm = TRUE),
    sd_beta              = sd(beta,     na.rm = TRUE),
    min_beta             = min(beta,    na.rm = TRUE),
    max_beta             = max(beta,    na.rm = TRUE),
    prop_ci_overlap      = ci_overlap_proportion(ci_lower, ci_upper),
    h2_confirmed         = prop_ci_overlap > 0.50,
    .groups = "drop"
  ) %>%
  arrange(dep_var, band, measure)

print(multiverse_df)
print(robustness_df)

write.csv(multiverse_df,  "Results/H2_multiverse_results.csv",  row.names = FALSE)
write.csv(robustness_df,  "Results/H2_robustness_results.csv",  row.names = FALSE)

# ------------------------------------------------------------------------------
# NA Checks
# ------------------------------------------------------------------------------

missing_pathl <- expand.grid(
  conn   = conn_measures,
  thresh = thresh_methods,
  band   = bands,
  stringsAsFactors = FALSE
) %>%
  mutate(col_name = paste("pathl", conn, band, thresh, sep = "_")) %>%
  rowwise() %>%
  mutate(all_na = if (col_name %in% colnames(full_data)) 
    all(is.na(full_data[[col_name]])) else NA) %>%
  ungroup() %>%
  filter(all_na == TRUE)

missing_smallworld <- expand.grid(
  conn   = conn_measures,
  thresh = thresh_methods,
  band   = bands,
  stringsAsFactors = FALSE
) %>%
  mutate(col_name = paste("smallworld", conn, band, thresh, sep = "_")) %>%
  rowwise() %>%
  mutate(all_na = if (col_name %in% colnames(full_data)) 
    all(is.na(full_data[[col_name]])) else NA) %>%
  ungroup() %>%
  filter(all_na == TRUE)

# ------------------------------------------------------------------------------
# Specification Curve Plots
# ------------------------------------------------------------------------------

col_main  <- "#D55E00"
col_other <- "#0072B2"
col_zero  <- "grey40"

band_labels <- c(
  "delta"  = "Delta",
  "theta"  = "Theta",
  "alpha1" = "Alpha-1",
  "alpha2" = "Alpha-2",
  "beta"   = "Beta"
)

measure_labels <- c(
  "cc"         = "Clustering Coefficient",
  "pathl"      = "Path Length",
  "eglob"      = "Global Efficiency",
  "eloc"       = "Local Efficiency",
  "smallworld" = "Small World Index"
)

dep_var_labels <- c(
  "gf_score"        = "Fluid Intelligence",
  "gc_score"        = "Crystallized Intelligence",
  "Agreeableness"   = "Agreeableness",
  "Conscientiousness" = "Conscientiousness",
  "Extraversion"    = "Extraversion",
  "Neuroticism"     = "Neuroticism",
  "Openness"        = "Openness"
)

main_path_label <- "Main path (ImCoh • Density-based Thresholding)"

n_specs_full <- length(conn_measures) * length(thresh_methods)

for (dep_var in dep_vars) {
  
  message("Plotting specification curve for: ", dep_var)
  
  plot_df <- multiverse_df %>%
    filter(dep_var == !!dep_var, !is.na(beta)) %>%
    mutate(
      main_path = conn == main_path_conn & thresh == main_path_thresh,
      band      = factor(band,    levels = bands),
      measure   = factor(measure, levels = measures)
    ) %>%
    group_by(band, measure) %>%
    arrange(beta, .by_group = TRUE) %>%
    mutate(
      spec_index = row_number(),
      n_specs    = n()
    ) %>%
    ungroup()
  
  # Filter incomplete facets
  incomplete_facets <- plot_df %>%
    distinct(band, measure, n_specs) %>%
    filter(n_specs != n_specs_full)
  
  if (nrow(incomplete_facets) > 0) {
    message(nrow(incomplete_facets), " of ", length(bands) * length(measures),
            " facets have less than ", n_specs_full, " specifications.")
  }
  
  # Plot specification curve
  p_spec <- ggplot(plot_df, aes(x = spec_index, y = beta)) +
    geom_linerange(
      aes(ymin = ci_lower, ymax = ci_upper, color = main_path),
      linewidth = 0.75, alpha = 0.7
    ) +
    geom_point(
      aes(color = main_path)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", 
               color = col_zero, linewidth = 0.4) +
    scale_color_manual(
      values = c("TRUE" = col_main, "FALSE" = col_other),
      labels = c("TRUE" = main_path_label, "FALSE" = "Other specifications")
    ) +
    facet_grid(
      rows     = vars(band),
      cols     = vars(measure),
      scales   = "free_x",
      labeller = labeller(band = band_labels, measure = measure_labels)
    ) +
    labs(
      x        = "Specifications (sorted by β)",
      y        = "Standardized β",
      color    = NULL,
      subtitle = dep_var_labels[dep_var]
    ) +
    theme_minimal(base_size = 16) +
    theme(
      strip.text          = element_text(size = 16, face = "bold"),
      legend.position     = "bottom",
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      axis.text.x         = element_blank(),
      axis.ticks.x        = element_blank(),
      plot.subtitle       = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  print(p_spec)
  
  ggsave(filename = paste0("Results/H2_specification_curve_", dep_var, ".png"),
         plot     = p_spec, width = 14, height = 10, dpi = 600, bg = "white")
  
}
