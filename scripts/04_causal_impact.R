# ==============================================================================
# 04_causal_impact.R
# Phase 3: Bayesian Causal Impact Analysis
#
# Uses Google's CausalImpact package (Bayesian structural time series) to
# estimate the causal effect of the January 2023 PBS copayment reduction on
# general patient (19-64) prescription rates, using concessional patients (65+)
# as a synthetic control series.
#
# Models: Overall (all ATC summed), by ATC Level 2 class, by state/territory
#
# Inputs:  data/processed/pbs_national_monthly.csv
#          data/processed/pbs_state_monthly.csv
#          outputs/tables/its_by_atc_results.csv (Phase 2 comparison)
#          outputs/tables/atc_class_summary.csv (ATC rankings)
#
# Outputs: outputs/tables/ci_overall_results.csv
#          outputs/tables/ci_by_atc_results.csv
#          outputs/tables/ci_by_state_results.csv
#          outputs/tables/ci_vs_its_comparison.csv
#          outputs/tables/ci_pre_period_correlations.csv
#          outputs/figures/fig14_ci_overall.png
#          outputs/figures/fig15_ci_top6_pointwise.png
#          outputs/figures/fig16_ci_forest_plot.png
#          outputs/figures/fig17_ci_vs_its_scatter.png
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(lubridate)
library(CausalImpact)
library(zoo)
library(scales)

proj_dir <- here::here()
proc_dir <- file.path(proj_dir, "data", "processed")
fig_dir  <- file.path(proj_dir, "outputs", "figures")
tab_dir  <- file.path(proj_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold"),
                  legend.position = "bottom"))

intervention_date <- as.Date("2023-01-01")
covid_start       <- as.Date("2020-03-01")
covid_end         <- as.Date("2021-12-31")

# Pre- and post-period for CausalImpact (as Date)
pre_period  <- as.Date(c("2002-05-01", "2022-12-01"))
post_period <- as.Date(c("2023-01-01", "2025-12-01"))

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================
message("Loading data...")

national <- read_csv(file.path(proc_dir, "pbs_national_monthly.csv"),
                     show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

state_data <- read_csv(file.path(proc_dir, "pbs_state_monthly.csv"),
                       show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

message("  National: ", format(nrow(national), big.mark = ","), " rows")
message("  State:    ", format(nrow(state_data), big.mark = ","), " rows")

# Load Phase 2 ITS results for comparison
its_atc_path <- file.path(tab_dir, "its_by_atc_results.csv")
if (file.exists(its_atc_path)) {
  its_atc_results <- read_csv(its_atc_path, show_col_types = FALSE)
  message("  Loaded ITS ATC results (", nrow(its_atc_results), " classes)")
} else {
  its_atc_results <- NULL
  message("  WARNING: its_by_atc_results.csv not found — comparison table will be skipped")
}

# Load ATC rankings
atc_ranking_path <- file.path(tab_dir, "atc_class_summary.csv")
if (file.exists(atc_ranking_path)) {
  atc_rankings <- read_csv(atc_ranking_path, show_col_types = FALSE) %>%
    mutate(rank = row_number()) %>%
    select(atc_name, rank, mean_scripts_erp)
  message("  Loaded ATC rankings (", nrow(atc_rankings), " classes)")
} else {
  message("  Computing ATC rankings from data...")
  atc_rankings <- national %>%
    filter(age_group == "Total") %>%
    group_by(atc_name) %>%
    summarise(mean_scripts_erp = mean(scripts_per_erp, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(mean_scripts_erp)) %>%
    mutate(rank = row_number())
}

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

#' Prepare data for CausalImpact: filter to treatment + control age groups,
#' pivot to wide format, convert to zoo with Date index
#' @param data Long-format data with age_group, date, scripts_per_erp columns
#' @param treatment_ag Treatment age group (e.g. "19-64")
#' @param control_ag Control age group (e.g. "65+")
#' @return zoo object with treatment as column 1, control as column 2
prepare_ci_data <- function(data, treatment_ag = "19-64", control_ag = "65+") {
  wide <- data %>%
    filter(age_group %in% c(treatment_ag, control_ag)) %>%
    group_by(date, age_group) %>%
    summarise(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = age_group, values_from = scripts_per_erp,
                values_fill = 0) %>%
    arrange(date)

  # Rename for clarity — CausalImpact uses column 1 as response
  treatment_col <- wide[[treatment_ag]]
  control_col   <- wide[[control_ag]]

  zoo_data <- zoo(cbind(treatment = treatment_col, control = control_col),
                  order.by = wide$date)

  # Drop any rows with NA (shouldn't happen with values_fill, but defensive)
  zoo_data <- zoo_data[complete.cases(coredata(zoo_data)), ]
  return(zoo_data)
}

#' Fit CausalImpact model with standard settings
#' @param zoo_data Zoo object (col 1 = treatment, col 2+ = control)
#' @param pre_period Length-2 Date vector
#' @param post_period Length-2 Date vector
#' @return CausalImpact object
fit_causal_impact <- function(zoo_data, pre_period, post_period) {
  CausalImpact(zoo_data,
               pre.period  = pre_period,
               post.period = post_period,
               model.args  = list(niter = 5000,
                                  nseasons = 12,
                                  season.duration = 1))
}

#' Extract tidy summary from CausalImpact object
#' @param impact CausalImpact object
#' @param label Character label for this model
#' @return tibble with key results
tidy_ci_results <- function(impact, label) {
  s <- impact$summary

  tibble(
    model_label       = label,
    # Average effects (per time point in post-period)
    abs_effect        = s["Average", "AbsEffect"],
    abs_effect_lower  = s["Average", "AbsEffect.lower"],
    abs_effect_upper  = s["Average", "AbsEffect.upper"],
    rel_effect_pct    = s["Average", "RelEffect"] * 100,
    rel_effect_lower  = s["Average", "RelEffect.lower"] * 100,
    rel_effect_upper  = s["Average", "RelEffect.upper"] * 100,
    # Cumulative effects (over entire post-period)
    cum_effect        = s["Cumulative", "AbsEffect"],
    cum_effect_lower  = s["Cumulative", "AbsEffect.lower"],
    cum_effect_upper  = s["Cumulative", "AbsEffect.upper"],
    cum_rel_pct       = s["Cumulative", "RelEffect"] * 100,
    # Bayesian tail-area probability (one-sided)
    p_value           = s["Average", "p"],
    # Actual vs predicted means in post-period
    post_actual       = s["Average", "Actual"],
    post_predicted    = s["Average", "Pred"],
    post_pred_lower   = s["Average", "Pred.lower"],
    post_pred_upper   = s["Average", "Pred.upper"]
  )
}

#' Add intervention line and COVID shading to a ggplot
add_intervention_lines <- function(p) {
  p +
    annotate("rect",
             xmin = covid_start, xmax = covid_end,
             ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.3) +
    geom_vline(xintercept = intervention_date,
               linetype = "dashed", colour = "red", linewidth = 0.6) +
    annotate("text", x = intervention_date + 60, y = Inf,
             label = "Copayment\nreduction", vjust = 1.5,
             colour = "red", size = 3, fontface = "italic") +
    annotate("text", x = covid_start + 150, y = Inf,
             label = "COVID", vjust = 1.5,
             colour = "grey50", size = 3, fontface = "italic")
}

# ==============================================================================
# 4. OVERALL CAUSALIMPACT MODEL
# ==============================================================================
message("\n--- Fitting overall CausalImpact model (19-64 vs 65+) ---")

# Aggregate all ATC classes for each age group
overall_zoo <- prepare_ci_data(national)

message("  Zoo dimensions: ", nrow(overall_zoo), " time points x ",
        ncol(overall_zoo), " series")
message("  Date range: ", min(index(overall_zoo)), " to ", max(index(overall_zoo)))

# Pre-period correlation check
pre_idx <- index(overall_zoo) >= pre_period[1] & index(overall_zoo) <= pre_period[2]
pre_cor <- cor(coredata(overall_zoo[pre_idx, 1]),
               coredata(overall_zoo[pre_idx, 2]))
message("  Pre-period correlation (treatment ~ control): ", sprintf("%.4f", pre_cor))

# Fit model
message("  Fitting CausalImpact (niter=5000, nseasons=12)...")
impact_overall <- fit_causal_impact(overall_zoo, pre_period, post_period)
results_overall <- tidy_ci_results(impact_overall, "Overall")

message("\n  OVERALL CAUSALIMPACT RESULTS:")
message("  Average causal effect:  ", sprintf("%.6f", results_overall$abs_effect),
        " [", sprintf("%.6f", results_overall$abs_effect_lower), ", ",
        sprintf("%.6f", results_overall$abs_effect_upper), "]")
message("  Relative effect:        ", sprintf("%.2f%%", results_overall$rel_effect_pct),
        " [", sprintf("%.2f%%", results_overall$rel_effect_lower), ", ",
        sprintf("%.2f%%", results_overall$rel_effect_upper), "]")
message("  Cumulative effect:      ", sprintf("%.4f", results_overall$cum_effect))
message("  Bayesian p-value:       ", sprintf("%.4f", results_overall$p_value))
message("  Post-period actual:     ", sprintf("%.6f", results_overall$post_actual))
message("  Post-period predicted:  ", sprintf("%.6f", results_overall$post_predicted))

# Print natural-language summary from CausalImpact
message("\n  --- CausalImpact narrative summary ---")
cat(impact_overall$report)
message("")

# ==============================================================================
# 5. ATC-STRATIFIED CAUSALIMPACT MODELS
# ==============================================================================
message("\n--- Fitting CausalImpact models by ATC class ---")

atc_classes <- national %>%
  filter(age_group %in% c("19-64", "65+")) %>%
  distinct(atc_name) %>%
  pull(atc_name)

message("  Total ATC classes: ", length(atc_classes))

# Determine top 6 by volume (retain full impact objects for figure 15)
top6_atc <- atc_rankings %>% head(6) %>% pull(atc_name)

atc_ci_results <- list()
atc_ci_impacts <- list()  # Only for top 6
n_skipped <- 0
n_errors  <- 0
error_log <- list()

pb <- txtProgressBar(min = 0, max = length(atc_classes), style = 3)

for (i in seq_along(atc_classes)) {
  atc <- atc_classes[i]
  setTxtProgressBar(pb, i)

  # Filter to this ATC class
  atc_data <- national %>% filter(atc_name == atc)

  # Check variance — skip if near-zero in treatment or control
  sd_treat <- atc_data %>%
    filter(age_group == "19-64") %>%
    pull(scripts_per_erp) %>%
    sd(na.rm = TRUE)
  sd_ctrl <- atc_data %>%
    filter(age_group == "65+") %>%
    pull(scripts_per_erp) %>%
    sd(na.rm = TRUE)

  if (is.na(sd_treat) || is.na(sd_ctrl) || sd_treat < 1e-6 || sd_ctrl < 1e-6) {
    n_skipped <- n_skipped + 1
    next
  }

  # Prepare zoo data
  zoo_atc <- prepare_ci_data(atc_data)

  # Check we have enough time points
  if (nrow(zoo_atc) < 50) {
    n_skipped <- n_skipped + 1
    next
  }

  # Fit model (tryCatch only around the CausalImpact call)
  result <- tryCatch({
    impact_atc <- fit_causal_impact(zoo_atc, pre_period, post_period)
    res <- tidy_ci_results(impact_atc, atc)

    # Store full impact object only for top 6
    if (atc %in% top6_atc) {
      atc_ci_impacts[[atc]] <- impact_atc
    }

    atc_ci_results[[atc]] <- res
    "OK"
  }, error = function(e) {
    n_errors  <<- n_errors + 1
    error_log[[atc]] <<- e$message
    "FAILED"
  })
}
close(pb)

n_success <- length(atc_ci_results)
message("  Successfully fitted: ", n_success, " / ", length(atc_classes),
        " (", n_skipped, " skipped, ", n_errors, " errors)")
if (n_errors > 0) {
  message("  First 5 errors:")
  for (j in seq_len(min(5, length(error_log)))) {
    message("    ", names(error_log)[j], ": ", error_log[[j]])
  }
}

# Bind results and apply FDR correction
atc_ci_df <- bind_rows(atc_ci_results) %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH")) %>%
  left_join(atc_rankings, by = c("model_label" = "atc_name")) %>%
  arrange(rank)

n_sig <- sum(atc_ci_df$p_fdr < 0.05, na.rm = TRUE)
message("  FDR-significant (p < 0.05): ", n_sig, " / ", n_success)

# ==============================================================================
# 6. PRE-PERIOD VALIDATION
# ==============================================================================
message("\n--- Pre-period validation ---")

# Overall pre-period R²
pre_overall <- as.data.frame(coredata(overall_zoo[pre_idx, ]))
overall_r2 <- summary(lm(treatment ~ control, data = pre_overall))$r.squared
message("  Overall pre-period R² (treatment ~ control): ", sprintf("%.4f", overall_r2))
message("  Overall pre-period correlation: ", sprintf("%.4f", pre_cor))

# Per-ATC correlations for top 20 classes
top20_atc <- atc_rankings %>% head(20) %>% pull(atc_name)

pre_correlations <- list()
for (atc in top20_atc) {
  atc_data <- national %>% filter(atc_name == atc)

  result <- tryCatch({
    zoo_atc <- prepare_ci_data(atc_data)
    pre_atc <- coredata(zoo_atc[index(zoo_atc) >= pre_period[1] &
                                  index(zoo_atc) <= pre_period[2], ])
    r <- cor(pre_atc[, 1], pre_atc[, 2], use = "complete.obs")
    r2 <- r^2
    pre_correlations[[atc]] <- tibble(
      atc_name = atc,
      pre_correlation = r,
      pre_r_squared = r2
    )
  }, error = function(e) NULL)
}

pre_cor_df <- bind_rows(pre_correlations) %>%
  left_join(atc_rankings, by = "atc_name") %>%
  arrange(rank)

message("  Pre-period correlations for top 20 ATC classes:")
for (j in 1:min(10, nrow(pre_cor_df))) {
  message("    ", sprintf("%-45s", pre_cor_df$atc_name[j]),
          " r = ", sprintf("%.3f", pre_cor_df$pre_correlation[j]),
          "  R² = ", sprintf("%.3f", pre_cor_df$pre_r_squared[j]))
}
if (nrow(pre_cor_df) > 10) {
  message("    ... (", nrow(pre_cor_df) - 10, " more in output file)")
}

# ==============================================================================
# 7. STATE-STRATIFIED CAUSALIMPACT MODELS
# ==============================================================================
message("\n--- Fitting CausalImpact models by state/territory ---")

states <- state_data %>% distinct(state) %>% pull(state) %>% sort()
message("  States: ", paste(states, collapse = ", "))

state_ci_results <- list()

for (st in states) {
  message("  Fitting: ", st, " ... ", appendLF = FALSE)

  result <- tryCatch({
    st_data <- state_data %>% filter(state == st)
    zoo_st <- prepare_ci_data(st_data)

    # Check dimensions
    if (nrow(zoo_st) < 50) {
      message("SKIPPED (insufficient data)")
      next
    }

    impact_st <- fit_causal_impact(zoo_st, pre_period, post_period)
    state_ci_results[[st]] <- tidy_ci_results(impact_st, st)
    message("OK")
  }, error = function(e) {
    message("FAILED (", e$message, ")")
  })
}

state_ci_df <- bind_rows(state_ci_results)
message("  Successfully fitted: ", nrow(state_ci_df), " / ", length(states), " states")

# ==============================================================================
# 8. TABLE OUTPUTS
# ==============================================================================
message("\n--- Saving tables ---")

# 8a. Overall results
write_csv(results_overall, file.path(tab_dir, "ci_overall_results.csv"))
message("  Saved ci_overall_results.csv")

# 8b. ATC class results with FDR
write_csv(atc_ci_df, file.path(tab_dir, "ci_by_atc_results.csv"))
message("  Saved ci_by_atc_results.csv (", nrow(atc_ci_df), " classes)")

# 8c. State results
write_csv(state_ci_df, file.path(tab_dir, "ci_by_state_results.csv"))
message("  Saved ci_by_state_results.csv (", nrow(state_ci_df), " states)")

# 8d. CausalImpact vs ITS comparison
if (!is.null(its_atc_results)) {
  ci_vs_its <- atc_ci_df %>%
    select(atc_name = model_label, ci_abs_effect = abs_effect,
           ci_abs_lower = abs_effect_lower, ci_abs_upper = abs_effect_upper,
           ci_rel_pct = rel_effect_pct, ci_p = p_value, ci_p_fdr = p_fdr) %>%
    left_join(
      its_atc_results %>%
        select(atc_name, its_level_estimate = level_estimate,
               its_level_ci_low = level_ci_low, its_level_ci_high = level_ci_high,
               its_level_p = level_p, its_level_p_fdr = level_p_fdr),
      by = "atc_name"
    ) %>%
    left_join(atc_rankings, by = "atc_name") %>%
    arrange(rank)

  write_csv(ci_vs_its, file.path(tab_dir, "ci_vs_its_comparison.csv"))
  message("  Saved ci_vs_its_comparison.csv (", nrow(ci_vs_its), " matched classes)")
} else {
  ci_vs_its <- NULL
  message("  SKIPPED ci_vs_its_comparison.csv (no ITS results available)")
}

# 8e. Pre-period correlations
write_csv(pre_cor_df, file.path(tab_dir, "ci_pre_period_correlations.csv"))
message("  Saved ci_pre_period_correlations.csv (", nrow(pre_cor_df), " classes)")

# ==============================================================================
# 9. FIGURES
# ==============================================================================
message("\n--- Generating figures ---")

# --- 9a. Fig 14: Overall CausalImpact 3-panel plot ---
# CausalImpact has a built-in plot method; we customise with ggplot layers
fig14 <- plot(impact_overall) +
  labs(title = "CausalImpact: PBS Prescriptions per Capita (19-64 vs 65+ Control)",
       subtitle = paste0("Pre-period: May 2002 – Dec 2022 | Post-period: Jan 2023 – Dec 2025 | ",
                         "Bayesian p = ", sprintf("%.4f", results_overall$p_value)))

ggsave(file.path(fig_dir, "fig14_ci_overall.png"),
       fig14, width = 11, height = 8, dpi = 300)
message("  Saved fig14_ci_overall.png")

# --- 9b. Fig 15: Top 6 ATC classes — pointwise effects ---
if (length(atc_ci_impacts) > 0) {

  # Extract pointwise effects from each stored impact object
  pointwise_list <- list()
  for (atc in names(atc_ci_impacts)) {
    imp <- atc_ci_impacts[[atc]]
    pw <- as.data.frame(imp$series)
    pw$date <- as.Date(rownames(pw))
    pw$atc_name <- atc
    pointwise_list[[atc]] <- pw
  }
  pointwise_df <- bind_rows(pointwise_list)

  # Faceted pointwise effect plot
  fig15_data <- pointwise_df %>%
    filter(date >= as.Date("2015-01-01")) %>%
    mutate(atc_name = factor(atc_name,
                             levels = intersect(top6_atc, names(atc_ci_impacts))))

  fig15 <- ggplot(fig15_data, aes(x = date)) +
    geom_ribbon(aes(ymin = point.effect.lower, ymax = point.effect.upper),
                fill = "steelblue", alpha = 0.2) +
    geom_line(aes(y = point.effect), colour = "steelblue", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.3) +
    geom_vline(xintercept = intervention_date, linetype = "dashed",
               colour = "red", linewidth = 0.5) +
    facet_wrap(~ atc_name, scales = "free_y", ncol = 2,
               labeller = labeller(atc_name = label_wrap_gen(35))) +
    labs(title = "CausalImpact Pointwise Effects: Top 6 ATC Classes by Volume",
         subtitle = "Blue = pointwise effect (actual - predicted) with 95% credible intervals",
         x = NULL, y = "Pointwise causal effect (scripts per ERP)") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y")

  ggsave(file.path(fig_dir, "fig15_ci_top6_pointwise.png"),
         fig15, width = 11, height = 9, dpi = 300)
  message("  Saved fig15_ci_top6_pointwise.png")
} else {
  message("  SKIPPED fig15 (no top-6 impact objects stored)")
}

# --- 9c. Fig 16: Forest plot of CausalImpact absolute effects, top 20 ---
top20_ci <- atc_ci_df %>%
  filter(!is.na(rank)) %>%
  head(20) %>%
  mutate(
    sig_label = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01  ~ "**",
      p_fdr < 0.05  ~ "*",
      TRUE           ~ ""
    ),
    atc_label = str_wrap(model_label, width = 30),
    atc_label = fct_reorder(atc_label, rank, .desc = TRUE)
  )

fig16 <- ggplot(top20_ci, aes(x = abs_effect, y = atc_label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = abs_effect_lower, xmax = abs_effect_upper),
                width = 0.3, colour = "steelblue", orientation = "y") +
  geom_point(aes(colour = p_fdr < 0.05), size = 2.5) +
  geom_text(aes(label = sig_label), nudge_y = 0.3, size = 3.5) +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey60"),
                      labels = c("TRUE" = "FDR p < 0.05", "FALSE" = "Not significant"),
                      name = NULL) +
  labs(title = "Forest Plot: CausalImpact Average Causal Effect on Prescriptions",
       subtitle = "Top 20 ATC Level 2 classes by volume, 95% Bayesian credible intervals",
       x = "Average causal effect (scripts per ERP)",
       y = NULL) +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(fig_dir, "fig16_ci_forest_plot.png"),
       fig16, width = 10, height = 8, dpi = 300)
message("  Saved fig16_ci_forest_plot.png")

# --- 9d. Fig 17: CausalImpact vs ITS scatter plot ---
if (!is.null(ci_vs_its)) {

  scatter_data <- ci_vs_its %>%
    filter(!is.na(ci_abs_effect), !is.na(its_level_estimate))

  if (nrow(scatter_data) > 0) {
    # Axis range: symmetric around max absolute value
    max_abs <- max(c(abs(scatter_data$ci_abs_effect),
                     abs(scatter_data$its_level_estimate)), na.rm = TRUE) * 1.15

    fig17 <- ggplot(scatter_data,
                    aes(x = its_level_estimate, y = ci_abs_effect)) +
      # 1:1 reference line
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  colour = "grey60", linewidth = 0.4) +
      # Zero reference lines
      geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
      geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
      # Error bars: ITS horizontal, CI vertical
      geom_errorbar(aes(ymin = ci_abs_lower, ymax = ci_abs_upper),
                    width = 0, colour = "steelblue", alpha = 0.3, linewidth = 0.3) +
      geom_errorbarh(aes(xmin = its_level_ci_low, xmax = its_level_ci_high),
                     height = 0, colour = "firebrick", alpha = 0.3, linewidth = 0.3) +
      # Points
      geom_point(aes(size = mean_scripts_erp), alpha = 0.6, colour = "grey30") +
      scale_size_continuous(name = "Mean scripts/ERP", range = c(1, 5),
                            labels = label_number(accuracy = 0.001)) +
      coord_cartesian(xlim = c(-max_abs, max_abs),
                      ylim = c(-max_abs, max_abs)) +
      labs(title = "CausalImpact vs ITS Estimates by ATC Class",
           subtitle = "Dashed line = 1:1 agreement | Blue bars = CI 95% credible | Red bars = ITS 95% Newey-West CI",
           x = expression(paste("ITS level change (", beta[2], ", scripts per ERP)")),
           y = "CausalImpact average causal effect (scripts per ERP)") +
      theme(legend.position = "right")

    ggsave(file.path(fig_dir, "fig17_ci_vs_its_scatter.png"),
           fig17, width = 9, height = 8, dpi = 300)
    message("  Saved fig17_ci_vs_its_scatter.png")
  } else {
    message("  SKIPPED fig17 (no matched data between CI and ITS)")
  }
} else {
  message("  SKIPPED fig17 (no ITS results for comparison)")
}

# ==============================================================================
# 10. SUMMARY
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PHASE 3: CAUSALIMPACT ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))

message("\nTables saved:")
message("  outputs/tables/ci_overall_results.csv")
message("  outputs/tables/ci_by_atc_results.csv    (", nrow(atc_ci_df), " ATC classes)")
message("  outputs/tables/ci_by_state_results.csv   (", nrow(state_ci_df), " states)")
if (!is.null(ci_vs_its)) {
  message("  outputs/tables/ci_vs_its_comparison.csv  (", nrow(ci_vs_its), " matched)")
}
message("  outputs/tables/ci_pre_period_correlations.csv")

message("\nFigures saved:")
message("  outputs/figures/fig14_ci_overall.png")
message("  outputs/figures/fig15_ci_top6_pointwise.png")
message("  outputs/figures/fig16_ci_forest_plot.png")
message("  outputs/figures/fig17_ci_vs_its_scatter.png")

message("\nKey findings:")
message("  Overall average causal effect: ", sprintf("%.6f", results_overall$abs_effect),
        " [", sprintf("%.6f", results_overall$abs_effect_lower), ", ",
        sprintf("%.6f", results_overall$abs_effect_upper), "]")
message("  Relative effect: ", sprintf("%.2f%%", results_overall$rel_effect_pct),
        " [", sprintf("%.2f%%", results_overall$rel_effect_lower), ", ",
        sprintf("%.2f%%", results_overall$rel_effect_upper), "]")
message("  Bayesian p-value: ", sprintf("%.4f", results_overall$p_value))
message("  Cumulative effect: ", sprintf("%.4f", results_overall$cum_effect))
message("  Pre-period correlation: ", sprintf("%.4f", pre_cor),
        "  |  R² = ", sprintf("%.4f", overall_r2))
message("  FDR-significant ATC classes: ", n_sig, " / ", n_success)

message("\nInterpretation:")
message("  Negative effect = prescriptions *decreased* relative to 65+ counterfactual")
message("  Positive effect = prescriptions *increased* relative to 65+ counterfactual")
message("  Relative effect = % change vs what would have happened without reform")

message("\nState-level variation:")
for (j in 1:nrow(state_ci_df)) {
  message("  ", sprintf("%-4s", state_ci_df$model_label[j]),
          " abs_effect = ", sprintf("%+.6f", state_ci_df$abs_effect[j]),
          " (", sprintf("%.2f%%", state_ci_df$rel_effect_pct[j]),
          ", p = ", sprintf("%.4f", state_ci_df$p_value[j]), ")")
}

message("\nNext: Run 05_equity_analysis.R (Phase 4)")
