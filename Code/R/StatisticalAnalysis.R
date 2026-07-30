## Statistical Analysis
# This is the main statistical analysis script for the study "Robustness of EEG
# functional brain networks associated with fluid intelligence: A multiverse
# analysis of connectivity and thresholding methods"
#
# Written by: Christoph Fruehlinger
# Last edit: July 2026

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

# gf data
gf_datapath <- "Data/IST_table.csv"
gf_data <- read.csv(gf_datapath, sep = ';')
gf_data$gf_score <- rowSums(gf_data[2:21])
gf_data <- gf_data[c(1, 22)]

# Lab data
lab_datapath <- "Data/Labs.txt"
lab_data <- read_tsv(lab_datapath)
lab_data <- lab_data[1:2]
colnames(lab_data)[2] <- "Lab"

# Combine gf and lab
gf_lab <- left_join(gf_data, lab_data, by = "ID")

# Graph metrics
datapath         <- "Data/Connectivity/Graph_data.csv"
main_path_conn   <- "imcoh"
main_path_thresh <- "dens"
data <- read.csv(datapath)

# Select main path columns
pattern        <- paste('[a-z0-9]+', main_path_conn, '[a-z0-9]+', main_path_thresh, sep = '_')
main_path_cols <- grep(pattern, colnames(data), value = TRUE)
main_path_data <- data[, c('ID', 'Run', 'Condition', main_path_cols)]

# Hypothesis 1 data
main_path_data <- main_path_data %>%
  filter(Run == 'first' & Condition == 'eyes_closed')
main_path_data <- left_join(main_path_data, gf_lab, by = "ID")

# Hypothesis 2 data
full_data <- data %>%
  filter(Run == 'first' & Condition == 'eyes_closed') %>%
  left_join(gf_lab, by = "ID")

# ------------------------------------------------------------------------------
## Hypothesis 1
# ------------------------------------------------------------------------------

# Analysis parameters
bands    <- c("delta", "theta", "alpha1", "alpha2", "beta")
measures <- c("cc", "pathl", "eglob", "eloc", "smallworld")
k        <- 5000

# Standardize
graph_cols_main <- main_path_cols
graph_cols_full <- grep(paste(measures, collapse = "|"), colnames(data), value = TRUE)

main_path_data <- main_path_data %>%
  mutate(across(all_of(c(graph_cols_main, "gf_score")), ~ as.numeric(scale(.x))))

full_data <- full_data %>%
  mutate(across(all_of(c(graph_cols_full, "gf_score")), ~ as.numeric(scale(.x))))

analysis_grid <- expand.grid(band = bands, measure = measures,
                             stringsAsFactors = FALSE)

# Set up parallel processing
n_cores <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_cores)
message("Using ", n_cores, " cores for parallel processing")

# Helper functions
permute_within_lab <- function(data) {
  data %>%
    group_by(Lab) %>%
    mutate(gf_score = sample(gf_score)) %>%
    ungroup()
}

get_t_value <- function(model, predictor) {
  coef(summary(model))[predictor, "t value"]
}

# Run one permutation test for a given band x measure combination
run_permutation_test <- function(band, measure, data, conn, thresh, k, p) {
  
  col_name <- paste(measure, conn, band, thresh, sep = "_")

  p(message = paste0(measure, " | ", band))
  
  if (!col_name %in% colnames(data)) {
    warning(paste("Column not found, skipping:", col_name))
    return(data.frame(band = band, measure = measure,
                      t_obs = NA, p_perm = NA, p_fdr = NA))
  }
  
  formula <- as.formula(paste("gf_score ~", col_name, "+ (1 | Lab)"))
  
  # Observed model
  fit_obs <- lmer(formula, data = data, REML = TRUE,
                  control = lmerControl(optimizer = "bobyqa"))
  t_obs   <- get_t_value(fit_obs, col_name)
  
  # Permutation distribution
  t_perm <- replicate(k, {
    data_perm <- permute_within_lab(data)
    fit_perm  <- lmer(formula, data = data_perm, REML = TRUE,
                      control = lmerControl(optimizer = "bobyqa"))
    get_t_value(fit_perm, col_name)
  })
  
  # Two-sided permutation p-value
  p_perm <- mean(abs(t_perm) >= abs(t_obs))
  
  data.frame(band = band, measure = measure, t_obs = t_obs, p_perm = p_perm)
}

# Mixed-effects models with permutation testing

# progress bar
handlers(handler_progress(
  format   = "[:bar] :percent | :current/:total combinations running",
  width    = 70,
  complete = "="
))

message("Starting permutation tests (k = ", k, ") across ",
        nrow(analysis_grid), " band x measure combinations...")

results_list <- with_progress({
  p <- progressor(steps = nrow(analysis_grid))
  
  future_mapply(
    FUN      = run_permutation_test,
    band     = analysis_grid$band,
    measure  = analysis_grid$measure,
    MoreArgs = list(
      data   = main_path_data,
      conn   = main_path_conn,
      thresh = main_path_thresh,
      k      = k,
      p      = p
    ),
    SIMPLIFY = FALSE,
    future.seed = 42
  )
})

# Combine and apply FDR correction

results_df <- bind_rows(results_list)

# FDR correction across the 5 measures within each frequency band
results_df <- results_df %>%
  group_by(band) %>%
  mutate(p_fdr = p.adjust(p_perm, method = "fdr")) %>%
  ungroup()

# Bonferroni correction for 5 frequency bands
results_df$significant <- results_df$p_fdr < 0.01

# Order rows
results_df <- results_df %>%
  mutate(band = factor(band, levels = bands),
         measure = factor(measure, levels = measures)) %>%
  arrange(band, measure)

# Reset to sequential processing
plan(sequential)

# Print results
print(results_df)

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
fit_multiverse_model <- function(conn, thresh, band, measure, data) {
  
  col_name <- paste(measure, conn, band, thresh, sep = "_")
  
  if (!col_name %in% colnames(data)) {
    return(data.frame(
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
    lmer(as.formula(paste("gf_score ~", col_name, "+ (1 | Lab)")),
         data = data, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(data.frame(conn, thresh, band, measure, col_name,
                      beta = NA, ci_lower = NA, ci_upper = NA,
                      fail_reason = paste0("fit_error: ", conditionMessage(fit))))
  }
  
  if (is.null(fit)) {
    return(data.frame(
      conn     = conn,
      thresh   = thresh,
      band     = band,
      measure  = measure,
      col_name = col_name,
      beta     = NA,
      ci_lower = NA,
      ci_upper = NA
    ))
  }
  
  fit_summary <- coef(summary(fit))
  
  if (!col_name %in% rownames(fit_summary)) {
    return(data.frame(
      conn     = conn,
      thresh   = thresh,
      band     = band,
      measure  = measure,
      col_name = col_name,
      beta     = NA,
      ci_lower = NA,
      ci_upper = NA
    ))
  }
  
  # Standardized regression coefficient and CI
  beta_std <- fixef(fit)[col_name]
  
  se <- fit_summary[col_name, "Std. Error"]
  
  ci_std <- c(
    beta_std - 1.96 * se,
    beta_std + 1.96 * se
  )
  
  data.frame(
    conn     = conn,
    thresh   = thresh,
    band     = band,
    measure  = measure,
    col_name = col_name,
    beta     = beta_std,
    ci_lower = ci_std[1],
    ci_upper = ci_std[2]
  )
}


# Run multiverse models

handlers(handler_progress(
  format   = "[:bar] :percent | :current/:total combinations running",
  width    = 70,
  complete = "="
))

message(
  "Starting multiverse analysis (",
  nrow(multiverse_grid),
  " models)..."
)

multiverse_results <- with_progress({
  
  p <- progressor(steps = nrow(multiverse_grid))
  
  future_mapply(
    FUN = function(conn, thresh, band, measure){
      p(message = paste0(
        measure, " | ",
        conn, " | ",
        thresh, " | ",
        band
      ))
      fit_multiverse_model(
        conn    = conn,
        thresh  = thresh,
        band    = band,
        measure = measure,
        data    = full_data
      )
    },
    conn    = multiverse_grid$conn,
    thresh  = multiverse_grid$thresh,
    band    = multiverse_grid$band,
    measure = multiverse_grid$measure,
    SIMPLIFY   = FALSE,
    future.seed = 42
  )
  
})

multiverse_df <- bind_rows(multiverse_results)


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


# Robustness assessment
robustness_df <- multiverse_df %>%
  group_by(band, measure) %>%
  summarise(
    n_specs = sum(!is.na(beta)),
    # Consistency in direction
    prop_positive = mean(beta > 0, na.rm = TRUE),
    direction_consistent = prop_positive >= 0.75 | prop_positive <= 0.25,
    # Consistency in magnitude
    mean_beta   = mean(beta, na.rm = TRUE),
    median_beta = median(beta, na.rm = TRUE),
    sd_beta     = sd(beta, na.rm = TRUE),
    min_beta    = min(beta, na.rm = TRUE),
    max_beta    = max(beta, na.rm = TRUE),
    # Proportion of overlapping confidence intervals
    prop_ci_overlap = ci_overlap_proportion(
      ci_lower,
      ci_upper
    ),
    # Hypothesis confirmed if majority of confidence intervals overlap
    h2_confirmed = prop_ci_overlap > 0.50,
    .groups = "drop"
  ) %>%
  mutate(
    band = factor(
      band,
      levels = bands
    ),
    measure = factor(
      measure,
      levels = measures
    )
  ) %>%
  arrange(
    band,
    measure
  )

print(multiverse_df)

print(robustness_df)

# NA Checks
missing_pathl <- expand.grid(
  conn = conn_measures, thresh = thresh_methods, band = bands,
  stringsAsFactors = FALSE) %>%
  mutate(col_name = paste("pathl", conn, band, thresh, sep = "_")) %>%
  rowwise() %>%
  mutate(all_na = if (col_name %in% colnames(full_data)) all(is.na(full_data[[col_name]])) else NA) %>%
  ungroup() %>%
    filter(all_na == TRUE)

missing_smallworld <- expand.grid(
  conn = conn_measures, thresh = thresh_methods, band = bands,
  stringsAsFactors = FALSE) %>%
  mutate(col_name = paste("smallworld", conn, band, thresh, sep = "_")) %>%
  rowwise() %>%
  mutate(all_na = if (col_name %in% colnames(full_data)) all(is.na(full_data[[col_name]])) else NA) %>%
  ungroup() %>%
  filter(all_na == TRUE)


# Specification Curve Plot

multiverse_df <- multiverse_df %>%
  mutate(main_path = conn == main_path_conn & thresh == main_path_thresh)

# Prepare df for plotting
plot_df <- multiverse_df %>%
  filter(!is.na(beta)) %>%
  mutate(
    band    = factor(band, levels = c("delta", "theta", "alpha1", "alpha2", "beta")),
    measure = factor(measure, levels = c("cc", "pathl", "eglob", "eloc", "smallworld"))
  ) %>%
  group_by(band, measure) %>%
  arrange(beta, .by_group = TRUE) %>%
  mutate(
    spec_index = row_number(),
    n_specs    = n()
  ) %>%
  ungroup()

n_specs_full <- length(conn_measures) * length(thresh_methods)

incomplete_facets <- plot_df %>%
  distinct(band, measure, n_specs) %>%
  filter(n_specs != n_specs_full)

if (nrow(incomplete_facets) > 0) {
  message(
    nrow(incomplete_facets),
    " von ", length(bands) * length(measures),
    " Facetten haben weniger als ", n_specs_full,
    " Spezifikationen (fehlende Spalten/NA-Modelle)."
  )
}

# Coloring
col_main  <- "#E63946"
col_other <- "#A8DADC"
col_zero  <- "grey40"

# Labeling
band_labels <- c(
  "delta"  = "Delta",
  "theta"  = "Theta",
  "alpha1" = "Alpha 1",
  "alpha2" = "Alpha 2",
  "beta"   = "Beta"
)

measure_labels <- c(
  "cc"         = "Clustering Coefficient",
  "pathl"      = "Path Length",
  "eglob"      = "Global Efficiency",
  "eloc"       = "Local Efficiency",
  "smallworld" = "Small World Index"
)

main_path_label <- paste0("Main path (ImCoh • Density-based Thresholding)")

# Plotting
p_spec <- ggplot(plot_df, aes(x = spec_index, y = beta)) +
  geom_linerange(
    aes(ymin = ci_lower, ymax = ci_upper, color = main_path),
    linewidth = 0.5, alpha = 0.6
  ) +
  geom_point(
    aes(color = main_path, size = main_path, shape = main_path)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = col_zero, linewidth = 0.4) +
  scale_color_manual(
    values = c("TRUE" = col_main, "FALSE" = col_other),
    labels = c("TRUE" = main_path_label, "FALSE" = "Other specifications")
  ) +
  scale_size_manual(
    values = c("TRUE" = 2.5, "FALSE" = 1.2),
    guide  = "none"
  ) +
  scale_shape_manual(
    values = c("TRUE" = 18, "FALSE" = 16),
    guide  = "none"
  ) +
  facet_grid(
    rows = vars(band),
    cols = vars(measure),
    scales = "free_x",
    labeller = labeller(band = band_labels, measure = measure_labels)
  ) +
  labs(
    x     = "Specifications (sorted by β)",
    y     = "Standardized β",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text         = element_text(size = 9, face = "bold"),
    legend.position     = "bottom",
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_blank(),
    axis.text.x         = element_blank(),
    axis.ticks.x        = element_blank(),
  )

print(p_spec)

# Save Plot
ggsave("Plots/H2_specification_curve.png", p_spec,
      width = 14, height = 10, dpi = 600, bg = "white")
