# ==============================================================================
# 03_its_analysis.R
# Phase 2: Interrupted Time Series Analysis — Segmented Regression
#
# Fits ITS models with Fourier seasonal adjustment and Newey-West HAC SEs
# to quantify the immediate level change (B2) and trend change (B3) in PBS
# prescription rates following the January 2023 copayment reduction.
#
# Models: Overall (all ATC summed), by ATC Level 2 class, by age group
#
# Inputs:  data/processed/pbs_national_monthly.csv
#          outputs/tables/atc_class_summary.csv (optional, for ATC rankings)
#
# Outputs: outputs/tables/its_overall_results.csv
#          outputs/tables/its_by_atc_results.csv
#          outputs/tables/its_by_age_group.csv
#          outputs/figures/fig10_its_overall.png
#          outputs/figures/fig11_its_top6_atc.png
#          outputs/figures/fig12_its_forest_plot.png
#          outputs/figures/fig13_its_age_comparison.png
#          outputs/figures/diag_acf_overall.png
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(lubridate)
library(lmtest)
library(sandwich)
library(broom)
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

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================
message("Loading data...")

national <- read_csv(file.path(proc_dir, "pbs_national_monthly.csv"),
                     show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

message("  National: ", format(nrow(national), big.mark = ","), " rows, ",
        n_distinct(national$atc_name), " ATC classes")

# --- Add Fourier terms (K=2: annual + semi-annual harmonics) ---
national <- national %>%
  mutate(
    sin1 = sin(2 * pi * 1 * month_num / 12),
    cos1 = cos(2 * pi * 1 * month_num / 12),
    sin2 = sin(2 * pi * 2 * month_num / 12),
    cos2 = cos(2 * pi * 2 * month_num / 12)
  )

# --- Aggregate overall time series (all ATC classes, Total age group) ---
overall_ts <- national %>%
  filter(age_group == "Total") %>%
  group_by(date, time_index, post_reform, time_post_reform, covid_period,
           sin1, cos1, sin2, cos2) %>%
  summarise(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(date)

message("  Overall time series: ", nrow(overall_ts), " months (",
        min(overall_ts$date), " to ", max(overall_ts$date), ")")

# --- Load ATC rankings (from Phase 1 descriptive analysis, or compute) ---
atc_ranking_path <- file.path(tab_dir, "atc_class_summary.csv")
if (file.exists(atc_ranking_path)) {
  atc_rankings <- read_csv(atc_ranking_path, show_col_types = FALSE) %>%
    mutate(rank = row_number()) %>%
    select(atc_name, rank, mean_scripts_erp)
  message("  Loaded ATC rankings from atc_class_summary.csv")
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

#' Fit ITS segmented regression model
#' @param data Data frame with required columns
#' @return lm object
fit_its_model <- function(data) {
  lm(scripts_per_erp ~ time_index + post_reform + time_post_reform +
       sin1 + cos1 + sin2 + cos2 + covid_period,
     data = data)
}

#' Extract Newey-West corrected coefficients into tidy tibble
#' @param model lm object
#' @param label Character label for this model
#' @return tibble with term, estimate, nw_se, ci_low, ci_high, p_value, dw_stat
tidy_its_results <- function(model, label) {
  # Newey-West HAC standard errors (auto bandwidth, no prewhitening)
  nw_vcov <- NeweyWest(model, prewhite = FALSE, adjust = TRUE)
  ct      <- coeftest(model, vcov. = nw_vcov)

  # Durbin-Watson statistic
  dw <- dwtest(model)$statistic

  tibble(
    model_label = label,
    term        = rownames(ct),
    estimate    = ct[, "Estimate"],
    nw_se       = ct[, "Std. Error"],
    t_value     = ct[, "t value"],
    p_value     = ct[, "Pr(>|t|)"],
    ci_low      = estimate - 1.96 * nw_se,
    ci_high     = estimate + 1.96 * nw_se,
    dw_stat     = dw,
    r_squared   = summary(model)$r.squared,
    adj_r_sq    = summary(model)$adj.r.squared,
    n_obs       = nobs(model)
  )
}

#' Generate counterfactual predictions (no reform scenario)
#' @param model lm object
#' @param data Data frame used to fit the model
#' @return Data frame with date, fitted, counterfactual columns
generate_counterfactual <- function(model, data) {
  # Fitted values
  fitted_vals <- predict(model, newdata = data)

  # Counterfactual: set post_reform = 0 and time_post_reform = 0

  cf_data <- data %>%
    mutate(post_reform = 0, time_post_reform = 0)
  cf_vals <- predict(model, newdata = cf_data)

  data %>%
    mutate(fitted = fitted_vals,
           counterfactual = cf_vals)
}

# ==============================================================================
# 4. OVERALL ITS MODEL
# ==============================================================================
message("\n--- Fitting overall ITS model ---")

model_overall <- fit_its_model(overall_ts)
results_overall <- tidy_its_results(model_overall, "Overall")
cf_overall <- generate_counterfactual(model_overall, overall_ts)

# Print key results
message("\nOverall ITS Results (Newey-West HAC SEs):")
message("  Level change (B2):  ", sprintf("%.6f", results_overall$estimate[results_overall$term == "post_reform"]),
        " (SE: ", sprintf("%.6f", results_overall$nw_se[results_overall$term == "post_reform"]),
        ", p = ", sprintf("%.4f", results_overall$p_value[results_overall$term == "post_reform"]), ")")
message("  Trend change (B3):  ", sprintf("%.6f", results_overall$estimate[results_overall$term == "time_post_reform"]),
        " (SE: ", sprintf("%.6f", results_overall$nw_se[results_overall$term == "time_post_reform"]),
        ", p = ", sprintf("%.4f", results_overall$p_value[results_overall$term == "time_post_reform"]), ")")
message("  DW statistic:       ", sprintf("%.4f", results_overall$dw_stat[1]))
message("  R-squared:          ", sprintf("%.4f", results_overall$r_squared[1]))
message("  Observations:       ", results_overall$n_obs[1])

# ==============================================================================
# 5. ITS BY ATC CLASS
# ==============================================================================
message("\n--- Fitting ITS models by ATC class ---")

atc_classes <- national %>%
  filter(age_group == "Total") %>%
  distinct(atc_name) %>%
  pull(atc_name)

message("  Fitting ", length(atc_classes), " ATC-level models...")

atc_models  <- list()
atc_results <- list()
atc_cf      <- list()

for (atc in atc_classes) {
  atc_data <- national %>%
    filter(age_group == "Total", atc_name == atc) %>%
    arrange(date)

  result <- tryCatch({
    mod <- fit_its_model(atc_data)
    res <- tidy_its_results(mod, atc)
    cf  <- generate_counterfactual(mod, atc_data) %>%
      mutate(atc_name = atc)

    atc_models[[atc]]  <- mod
    atc_results[[atc]] <- res
    atc_cf[[atc]]      <- cf
    "OK"
  }, error = function(e) {
    message("  WARNING: Failed for '", atc, "': ", e$message)
    "FAILED"
  })
}

n_success <- length(atc_results)
message("  Successfully fitted: ", n_success, " / ", length(atc_classes))

# Bind all ATC results
atc_results_df <- bind_rows(atc_results)
atc_cf_df      <- bind_rows(atc_cf)

# --- FDR correction (Benjamini-Hochberg) for level change and trend change ---
atc_level_change <- atc_results_df %>%
  filter(term == "post_reform") %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH"))

atc_trend_change <- atc_results_df %>%
  filter(term == "time_post_reform") %>%
  mutate(p_fdr = p.adjust(p_value, method = "BH"))

# Merge FDR-adjusted p-values and rankings
atc_summary_its <- atc_level_change %>%
  select(atc_name = model_label, level_estimate = estimate, level_se = nw_se,
         level_ci_low = ci_low, level_ci_high = ci_high,
         level_p = p_value, level_p_fdr = p_fdr,
         dw_stat, r_squared, n_obs) %>%
  left_join(
    atc_trend_change %>%
      select(atc_name = model_label, trend_estimate = estimate, trend_se = nw_se,
             trend_ci_low = ci_low, trend_ci_high = ci_high,
             trend_p = p_value, trend_p_fdr = p_fdr),
    by = "atc_name"
  ) %>%
  left_join(atc_rankings, by = "atc_name") %>%
  arrange(rank)

n_sig_level <- sum(atc_summary_its$level_p_fdr < 0.05, na.rm = TRUE)
n_sig_trend <- sum(atc_summary_its$trend_p_fdr < 0.05, na.rm = TRUE)
message("  FDR-significant level changes:  ", n_sig_level, " / ", n_success)
message("  FDR-significant trend changes:  ", n_sig_trend, " / ", n_success)

# ==============================================================================
# 6. ITS BY AGE GROUP
# ==============================================================================
message("\n--- Fitting ITS models by age group ---")

age_groups <- c("19-64", "65+")
age_models  <- list()
age_results <- list()
age_cf      <- list()

for (ag in age_groups) {
  ag_data <- national %>%
    filter(age_group == ag) %>%
    group_by(date, time_index, post_reform, time_post_reform, covid_period,
             sin1, cos1, sin2, cos2) %>%
    summarise(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(date)

  mod <- fit_its_model(ag_data)
  age_models[[ag]]  <- mod
  age_results[[ag]] <- tidy_its_results(mod, ag)
  age_cf[[ag]]      <- generate_counterfactual(mod, ag_data) %>%
    mutate(age_group = ag)
}

age_results_df <- bind_rows(age_results)
age_cf_df      <- bind_rows(age_cf)

# Print comparison
for (ag in age_groups) {
  b2 <- age_results_df %>% filter(model_label == ag, term == "post_reform")
  b3 <- age_results_df %>% filter(model_label == ag, term == "time_post_reform")
  message("  ", ag, ":")
  message("    Level change: ", sprintf("%.6f", b2$estimate),
          " (p = ", sprintf("%.4f", b2$p_value), ")")
  message("    Trend change: ", sprintf("%.6f", b3$estimate),
          " (p = ", sprintf("%.4f", b3$p_value), ")")
}

# ==============================================================================
# 7. TABLE OUTPUTS
# ==============================================================================
message("\n--- Saving tables ---")

# 7a. Overall results
write_csv(results_overall, file.path(tab_dir, "its_overall_results.csv"))
message("  Saved its_overall_results.csv")

# 7b. ATC class results (wide format with FDR-adjusted p-values)
write_csv(atc_summary_its, file.path(tab_dir, "its_by_atc_results.csv"))
message("  Saved its_by_atc_results.csv (", nrow(atc_summary_its), " classes)")

# 7c. Age group results
write_csv(age_results_df, file.path(tab_dir, "its_by_age_group.csv"))
message("  Saved its_by_age_group.csv")

# ==============================================================================
# 8. FIGURES
# ==============================================================================
message("\n--- Generating figures ---")

# Helper: add intervention + COVID shading
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

# --- 8a. Fig 10: Overall ITS (observed + fitted + counterfactual) ---
fig10 <- ggplot(cf_overall, aes(x = date)) +
  geom_line(aes(y = scripts_per_erp, colour = "Observed"), linewidth = 0.4) +
  geom_line(aes(y = fitted, colour = "Fitted"), linewidth = 0.5) +
  geom_line(aes(y = counterfactual, colour = "Counterfactual"),
            linewidth = 0.5, linetype = "dashed") +
  scale_colour_manual(
    values = c("Observed" = "grey50", "Fitted" = "steelblue",
               "Counterfactual" = "red"),
    breaks = c("Observed", "Fitted", "Counterfactual")
  ) +
  labs(title = "Interrupted Time Series: Overall PBS Prescriptions per Capita",
       subtitle = "Segmented regression with Fourier seasonal terms (K=2) and COVID covariate",
       x = NULL, y = "Scripts per ERP (all ATC classes summed)", colour = NULL) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")

fig10 <- add_intervention_lines(fig10)
ggsave(file.path(fig_dir, "fig10_its_overall.png"),
       fig10, width = 11, height = 5.5, dpi = 300)
message("  Saved fig10_its_overall.png")

# --- 8b. Fig 11: Top 6 ATC classes faceted ITS ---
top6_atc <- atc_rankings %>% head(6) %>% pull(atc_name)

fig11_data <- atc_cf_df %>%
  filter(atc_name %in% top6_atc) %>%
  mutate(atc_name = factor(atc_name, levels = top6_atc))

fig11 <- ggplot(fig11_data, aes(x = date)) +
  geom_line(aes(y = scripts_per_erp), colour = "grey50", linewidth = 0.3) +
  geom_line(aes(y = fitted), colour = "steelblue", linewidth = 0.4) +
  geom_line(aes(y = counterfactual), colour = "red",
            linewidth = 0.4, linetype = "dashed") +
  facet_wrap(~ atc_name, scales = "free_y", ncol = 2,
             labeller = labeller(atc_name = label_wrap_gen(35))) +
  geom_vline(xintercept = intervention_date,
             linetype = "dashed", colour = "red", linewidth = 0.4, alpha = 0.5) +
  annotate("rect", xmin = covid_start, xmax = covid_end,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  labs(title = "ITS Models for Top 6 ATC Level 2 Classes by Volume",
       subtitle = "Grey = observed, Blue = fitted, Red dashed = counterfactual (no reform)",
       x = NULL, y = "Scripts per ERP") +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y")

ggsave(file.path(fig_dir, "fig11_its_top6_atc.png"),
       fig11, width = 11, height = 9, dpi = 300)
message("  Saved fig11_its_top6_atc.png")

# --- 8c. Fig 12: Forest plot of level change (B2), top 20 ATC classes ---
top20_for_forest <- atc_summary_its %>%
  filter(!is.na(rank)) %>%
  head(20) %>%
  mutate(
    sig_label = case_when(
      level_p_fdr < 0.001 ~ "***",
      level_p_fdr < 0.01  ~ "**",
      level_p_fdr < 0.05  ~ "*",
      TRUE                ~ ""
    ),
    atc_label = str_wrap(atc_name, width = 30),
    atc_label = fct_reorder(atc_label, rank, .desc = TRUE)
  )

fig12 <- ggplot(top20_for_forest,
                aes(x = level_estimate, y = atc_label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = level_ci_low, xmax = level_ci_high),
                width = 0.3, colour = "steelblue", orientation = "y") +
  geom_point(aes(colour = level_p_fdr < 0.05), size = 2.5) +
  geom_text(aes(label = sig_label), nudge_y = 0.3, size = 3.5) +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey60"),
                      labels = c("TRUE" = "FDR p < 0.05", "FALSE" = "Not significant"),
                      name = NULL) +
  labs(title = expression(paste("Forest Plot: Immediate Level Change (", beta[2], ") at Copayment Reform")),
       subtitle = "Top 20 ATC Level 2 classes by prescription volume, Newey-West 95% CIs",
       x = expression(paste("Level change in scripts per ERP (", beta[2], ")")),
       y = NULL) +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(fig_dir, "fig12_its_forest_plot.png"),
       fig12, width = 10, height = 8, dpi = 300)
message("  Saved fig12_its_forest_plot.png")

# --- 8d. Fig 13: Age group comparison ITS ---
fig13 <- ggplot(age_cf_df, aes(x = date)) +
  geom_line(aes(y = scripts_per_erp), colour = "grey50", linewidth = 0.3) +
  geom_line(aes(y = fitted), colour = "steelblue", linewidth = 0.5) +
  geom_line(aes(y = counterfactual), colour = "red",
            linewidth = 0.5, linetype = "dashed") +
  facet_wrap(~ age_group, scales = "free_y", ncol = 1,
             labeller = labeller(age_group = c("19-64" = "19-64 (General Proxy)",
                                               "65+" = "65+ (Concessional Proxy)"))) +
  geom_vline(xintercept = intervention_date,
             linetype = "dashed", colour = "red", linewidth = 0.4, alpha = 0.5) +
  annotate("rect", xmin = covid_start, xmax = covid_end,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  labs(title = "ITS Comparison: General (19-64) vs Concessional (65+) Proxy",
       subtitle = "Grey = observed, Blue = fitted, Red dashed = counterfactual. Reform should affect 19-64 more.",
       x = NULL, y = "Scripts per ERP (all ATC classes summed)") +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y")

ggsave(file.path(fig_dir, "fig13_its_age_comparison.png"),
       fig13, width = 10, height = 7, dpi = 300)
message("  Saved fig13_its_age_comparison.png")

# ==============================================================================
# 9. DIAGNOSTICS
# ==============================================================================
message("\n--- Diagnostics ---")

# Durbin-Watson test
dw_test <- dwtest(model_overall)
message("  Durbin-Watson test (overall model):")
message("    DW = ", sprintf("%.4f", dw_test$statistic),
        ", p = ", sprintf("%.4f", dw_test$p.value))

# ACF/PACF diagnostic plot
png(file.path(fig_dir, "diag_acf_overall.png"),
    width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
acf(residuals(model_overall), main = "ACF of Overall ITS Residuals",
    lag.max = 36)
pacf(residuals(model_overall), main = "PACF of Overall ITS Residuals",
     lag.max = 36)
dev.off()
message("  Saved diag_acf_overall.png")

# Model summary to console
message("\n--- Overall model summary ---")
print(summary(model_overall))

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PHASE 2: ITS ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))

message("\nTables saved:")
message("  outputs/tables/its_overall_results.csv")
message("  outputs/tables/its_by_atc_results.csv  (", nrow(atc_summary_its), " ATC classes)")
message("  outputs/tables/its_by_age_group.csv")

message("\nFigures saved:")
message("  outputs/figures/fig10_its_overall.png")
message("  outputs/figures/fig11_its_top6_atc.png")
message("  outputs/figures/fig12_its_forest_plot.png")
message("  outputs/figures/fig13_its_age_comparison.png")
message("  outputs/figures/diag_acf_overall.png")

message("\nKey findings:")
b2_overall <- results_overall %>% filter(term == "post_reform")
b3_overall <- results_overall %>% filter(term == "time_post_reform")
message("  Overall level change (B2): ", sprintf("%.6f", b2_overall$estimate),
        " [", sprintf("%.6f", b2_overall$ci_low), ", ",
        sprintf("%.6f", b2_overall$ci_high), "]",
        " (p = ", sprintf("%.4f", b2_overall$p_value), ")")
message("  Overall trend change (B3): ", sprintf("%.6f", b3_overall$estimate),
        " [", sprintf("%.6f", b3_overall$ci_low), ", ",
        sprintf("%.6f", b3_overall$ci_high), "]",
        " (p = ", sprintf("%.4f", b3_overall$p_value), ")")
message("  FDR-significant level changes across ATC classes: ", n_sig_level, " / ", n_success)
message("  FDR-significant trend changes across ATC classes: ", n_sig_trend, " / ", n_success)

b2_1964 <- age_results_df %>% filter(model_label == "19-64", term == "post_reform")
b2_65   <- age_results_df %>% filter(model_label == "65+", term == "post_reform")
message("  19-64 level change: ", sprintf("%.6f", b2_1964$estimate),
        " | 65+ level change: ", sprintf("%.6f", b2_65$estimate))

message("\nNext: Run 04_causal_impact.R (Phase 3)")
