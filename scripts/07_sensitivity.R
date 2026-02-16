# ==============================================================================
# 07_sensitivity.R
# Phase 5: Sensitivity and Robustness Analyses
#
# Tests robustness of primary findings from Phases 2-4:
#   A. 60-day dispensing adjustment (ITS + DiD)
#   B. COVID period exclusion (ITS + CI + DiD)
#   C. Placebo intervention dates (ITS only)
#   D. Shorter pre-period (ITS + CI)
#   E. CausalImpact correlation filtering
#   F. Equity window robustness (DiD only)
#
# Inputs:  data/processed/pbs_national_monthly.csv
#          data/processed/pbs_lga_monthly.csv
#          data/processed/lga_characteristics.csv
#          outputs/tables/ (Phase 2-4 base results)
#
# Outputs: outputs/tables/sensitivity_*.csv (10 tables)
#          outputs/figures/fig26_placebo_tests.png
#          outputs/figures/fig27_ci_correlation_filter.png
#          outputs/figures/fig28_sensitivity_summary.png
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(lubridate)
library(lmtest)
library(sandwich)
library(CausalImpact)
library(zoo)
library(data.table)
library(fixest)
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
sixty_day_date    <- as.Date("2023-09-01")

# CausalImpact periods (base)
pre_period  <- as.Date(c("2002-05-01", "2022-12-01"))
post_period <- as.Date(c("2023-01-01", "2025-12-01"))

# ==============================================================================
# 2. LOAD DATA AND DEFINE HELPERS
# ==============================================================================
message("Loading data...")

# --- National data (88K rows) ---
national <- read_csv(file.path(proc_dir, "pbs_national_monthly.csv"),
                     show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

# Add Fourier terms and sixty_day covariate
national <- national %>%
  mutate(
    sin1 = sin(2 * pi * 1 * month_num / 12),
    cos1 = cos(2 * pi * 1 * month_num / 12),
    sin2 = sin(2 * pi * 2 * month_num / 12),
    cos2 = cos(2 * pi * 2 * month_num / 12),
    sixty_day = as.integer(date >= sixty_day_date)
  )

message("  National: ", format(nrow(national), big.mark = ","), " rows")

# Aggregate overall time series (all ATC, Total age group)
overall_ts <- national %>%
  filter(age_group == "Total") %>%
  group_by(date, time_index, post_reform, time_post_reform, covid_period,
           sin1, cos1, sin2, cos2, sixty_day) %>%
  summarise(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(date)

message("  Overall TS: ", nrow(overall_ts), " months")

# Age-group-specific time series
age_ts <- list()
for (ag in c("19-64", "65+")) {
  age_ts[[ag]] <- national %>%
    filter(age_group == ag) %>%
    group_by(date, time_index, post_reform, time_post_reform, covid_period,
             sin1, cos1, sin2, cos2, sixty_day) %>%
    summarise(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(date)
}

# Prepare CausalImpact zoo data (overall: 19-64 vs 65+)
overall_zoo <- {
  wide <- national %>%
    filter(age_group %in% c("19-64", "65+")) %>%
    group_by(date, age_group) %>%
    summarise(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = age_group, values_from = scripts_per_erp,
                values_fill = 0) %>%
    arrange(date)
  zoo(cbind(treatment = wide[["19-64"]], control = wide[["65+"]]),
      order.by = wide$date)
}

message("  Zoo: ", nrow(overall_zoo), " time points")

# --- LGA data (all years, 19-64 only, aggregated across ATC) ---
message("Loading LGA panel (this may take a minute)...")

lga_raw <- fread(
  file.path(proc_dir, "pbs_lga_monthly.csv"),
  select = c("lga_code", "date", "age_group", "scripts_per_erp",
             "time_index", "post_reform", "covid_period", "month_num")
)
lga_raw <- lga_raw[age_group == "19-64"]

lga_agg <- lga_raw[,
  .(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE)),
  by = .(lga_code, date, time_index, post_reform, covid_period, month_num)
]
rm(lga_raw); gc()

lga_full_panel <- as_tibble(lga_agg) %>%
  mutate(
    date = as.Date(date),
    lga_code = as.integer(lga_code),
    sin1 = sin(2 * pi * 1 * month_num / 12),
    cos1 = cos(2 * pi * 1 * month_num / 12),
    sin2 = sin(2 * pi * 2 * month_num / 12),
    cos2 = cos(2 * pi * 2 * month_num / 12),
    sixty_day = as.integer(date >= sixty_day_date)
  )
rm(lga_agg)

# Load LGA characteristics and merge
lga_chars <- read_csv(file.path(proc_dir, "lga_characteristics.csv"),
                      show_col_types = FALSE) %>%
  mutate(lga_code = as.integer(lga_code))

lga_full_panel <- lga_full_panel %>%
  left_join(lga_chars %>% select(lga_code, seifa_quintile, remoteness_area),
            by = "lga_code") %>%
  filter(!is.na(seifa_quintile), !is.na(remoteness_area)) %>%
  mutate(
    seifa_quintile = factor(seifa_quintile, levels = c(5, 1, 2, 3, 4)),
    remoteness_area = factor(remoteness_area,
                             levels = c("Major Cities of Australia",
                                        "Inner Regional Australia",
                                        "Outer Regional Australia",
                                        "Remote Australia",
                                        "Very Remote Australia"))
  )

# Base DiD window: 2018-2025
lga_panel <- lga_full_panel %>%
  filter(date >= as.Date("2018-01-01"), date <= as.Date("2025-12-01"))

message("  LGA full panel: ", format(nrow(lga_full_panel), big.mark = ","),
        " rows (", n_distinct(lga_full_panel$lga_code), " LGAs)")
message("  LGA base panel (2018-2025): ", format(nrow(lga_panel), big.mark = ","),
        " rows")

# --- Helper functions (from scripts 03/04) ---

fit_its_model <- function(data) {
  lm(scripts_per_erp ~ time_index + post_reform + time_post_reform +
       sin1 + cos1 + sin2 + cos2 + covid_period,
     data = data)
}

fit_its_model_60d <- function(data) {
  lm(scripts_per_erp ~ time_index + post_reform + time_post_reform +
       sin1 + cos1 + sin2 + cos2 + covid_period + sixty_day,
     data = data)
}

fit_its_model_nocovid <- function(data) {
  lm(scripts_per_erp ~ time_index + post_reform + time_post_reform +
       sin1 + cos1 + sin2 + cos2,
     data = data)
}

tidy_its_results <- function(model, label) {
  nw_vcov <- NeweyWest(model, prewhite = FALSE, adjust = TRUE)
  ct      <- coeftest(model, vcov. = nw_vcov)
  dw      <- dwtest(model)$statistic
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

fit_causal_impact <- function(zoo_data, pre_period, post_period) {
  CausalImpact(zoo_data,
               pre.period  = pre_period,
               post.period = post_period,
               model.args  = list(niter = 5000,
                                  nseasons = 12,
                                  season.duration = 1))
}

tidy_ci_results <- function(impact, label) {
  s <- impact$summary
  tibble(
    model_label       = label,
    abs_effect        = s["Average", "AbsEffect"],
    abs_effect_lower  = s["Average", "AbsEffect.lower"],
    abs_effect_upper  = s["Average", "AbsEffect.upper"],
    rel_effect_pct    = s["Average", "RelEffect"] * 100,
    rel_effect_lower  = s["Average", "RelEffect.lower"] * 100,
    rel_effect_upper  = s["Average", "RelEffect.upper"] * 100,
    cum_effect        = s["Cumulative", "AbsEffect"],
    cum_effect_lower  = s["Cumulative", "AbsEffect.lower"],
    cum_effect_upper  = s["Cumulative", "AbsEffect.upper"],
    cum_rel_pct       = s["Cumulative", "RelEffect"] * 100,
    p_value           = s["Average", "p"],
    post_actual       = s["Average", "Actual"],
    post_predicted    = s["Average", "Pred"],
    post_pred_lower   = s["Average", "Pred.lower"],
    post_pred_upper   = s["Average", "Pred.upper"]
  )
}

# Summary row helper
add_summary <- function(scenario, param, est, lo, hi, pval) {
  tibble(scenario = scenario, parameter = param,
         estimate = est, ci_low = lo, ci_high = hi,
         p_value = pval, significant = pval < 0.05)
}

# Collector for summary table
summary_rows <- list()

# --- Load base results from Phases 2-4 ---
message("Loading base results...")

base_its <- read_csv(file.path(tab_dir, "its_overall_results.csv"),
                     show_col_types = FALSE)
base_its_b2 <- base_its %>% filter(term == "post_reform")
summary_rows[["base_its_b2"]] <- add_summary(
  "Base", "ITS Overall B2",
  base_its_b2$estimate, base_its_b2$ci_low, base_its_b2$ci_high,
  base_its_b2$p_value)

base_its_age <- read_csv(file.path(tab_dir, "its_by_age_group.csv"),
                         show_col_types = FALSE)
base_its_1964 <- base_its_age %>%
  filter(model_label == "19-64", term == "post_reform")
summary_rows[["base_its_1964"]] <- add_summary(
  "Base", "ITS 19-64 B2",
  base_its_1964$estimate, base_its_1964$ci_low, base_its_1964$ci_high,
  base_its_1964$p_value)

base_ci <- read_csv(file.path(tab_dir, "ci_overall_results.csv"),
                    show_col_types = FALSE)
summary_rows[["base_ci"]] <- add_summary(
  "Base", "CI Overall Rel%",
  base_ci$rel_effect_pct, base_ci$rel_effect_lower,
  base_ci$rel_effect_upper, base_ci$p_value)

base_did_seifa <- read_csv(file.path(tab_dir, "did_seifa_results.csv"),
                           show_col_types = FALSE)
base_did_q1 <- base_did_seifa %>%
  filter(term == "post_reform:seifa_quintile1")
summary_rows[["base_did_q1"]] <- add_summary(
  "Base", "DiD SEIFA Q1",
  base_did_q1$estimate, base_did_q1$conf.low, base_did_q1$conf.high,
  base_did_q1$p.value)

base_did_remote <- read_csv(file.path(tab_dir, "did_remoteness_results.csv"),
                            show_col_types = FALSE)
base_did_remote_row <- base_did_remote %>%
  filter(term == "post_reform:remoteness_areaRemote Australia")
summary_rows[["base_did_remote"]] <- add_summary(
  "Base", "DiD Remote",
  base_did_remote_row$estimate, base_did_remote_row$conf.low,
  base_did_remote_row$conf.high, base_did_remote_row$p.value)

message("  Base results loaded")

# ==============================================================================
# 3. SCENARIO A — 60-DAY DISPENSING ADJUSTMENT
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SCENARIO A: 60-day dispensing adjustment")
message(paste(rep("=", 60), collapse = ""))

# --- A.1: ITS overall with sixty_day ---
message("\n  Fitting ITS overall + sixty_day...")
model_a_overall <- fit_its_model_60d(overall_ts)
results_a_overall <- tidy_its_results(model_a_overall, "A: Overall + 60-day")

b2_a <- results_a_overall %>% filter(term == "post_reform")
sixty_a <- results_a_overall %>% filter(term == "sixty_day")
message("    B2 (post_reform): ", sprintf("%+.6f", b2_a$estimate),
        " (p = ", sprintf("%.4f", b2_a$p_value), ")")
message("    sixty_day coeff:  ", sprintf("%+.6f", sixty_a$estimate),
        " (p = ", sprintf("%.4f", sixty_a$p_value), ")")

summary_rows[["a_its_b2"]] <- add_summary(
  "A: 60-day", "ITS Overall B2",
  b2_a$estimate, b2_a$ci_low, b2_a$ci_high, b2_a$p_value)

# --- A.2: ITS by age group with sixty_day ---
message("  Fitting ITS by age group + sixty_day...")
results_a_age <- list()
for (ag in c("19-64", "65+")) {
  mod <- fit_its_model_60d(age_ts[[ag]])
  results_a_age[[ag]] <- tidy_its_results(mod, paste0("A: ", ag, " + 60-day"))
  b2 <- results_a_age[[ag]] %>% filter(term == "post_reform")
  message("    ", ag, " B2: ", sprintf("%+.6f", b2$estimate),
          " (p = ", sprintf("%.4f", b2$p_value), ")")
}

b2_a_1964 <- results_a_age[["19-64"]] %>% filter(term == "post_reform")
summary_rows[["a_its_1964"]] <- add_summary(
  "A: 60-day", "ITS 19-64 B2",
  b2_a_1964$estimate, b2_a_1964$ci_low, b2_a_1964$ci_high, b2_a_1964$p_value)

# Save ITS results
its_a <- bind_rows(results_a_overall, bind_rows(results_a_age))
write_csv(its_a, file.path(tab_dir, "sensitivity_sixty_day_its.csv"))
message("  Saved sensitivity_sixty_day_its.csv")

# --- A.3: DiD SEIFA with sixty_day ---
message("  Fitting DiD SEIFA + sixty_day...")
did_a_seifa <- feols(
  scripts_per_erp ~ post_reform * seifa_quintile +
    time_index + covid_period + sixty_day +
    sin1 + cos1 + sin2 + cos2 | lga_code,
  data = lga_panel,
  vcov = ~lga_code
)

did_a_seifa_tidy <- broom::tidy(did_a_seifa, conf.int = TRUE)
q1_a <- did_a_seifa_tidy %>% filter(term == "post_reform:seifa_quintile1")
message("    Q1 interaction: ", sprintf("%+.6f", q1_a$estimate),
        " (p = ", sprintf("%.4f", q1_a$p.value), ")")

summary_rows[["a_did_q1"]] <- add_summary(
  "A: 60-day", "DiD SEIFA Q1",
  q1_a$estimate, q1_a$conf.low, q1_a$conf.high, q1_a$p.value)

# --- A.4: DiD Remoteness with sixty_day ---
message("  Fitting DiD Remoteness + sixty_day...")
did_a_remote <- feols(
  scripts_per_erp ~ post_reform * remoteness_area +
    time_index + covid_period + sixty_day +
    sin1 + cos1 + sin2 + cos2 | lga_code,
  data = lga_panel,
  vcov = ~lga_code
)

did_a_remote_tidy <- broom::tidy(did_a_remote, conf.int = TRUE)
remote_a <- did_a_remote_tidy %>%
  filter(term == "post_reform:remoteness_areaRemote Australia")
message("    Remote interaction: ", sprintf("%+.6f", remote_a$estimate),
        " (p = ", sprintf("%.4f", remote_a$p.value), ")")

summary_rows[["a_did_remote"]] <- add_summary(
  "A: 60-day", "DiD Remote",
  remote_a$estimate, remote_a$conf.low, remote_a$conf.high, remote_a$p.value)

# Save DiD results
did_a <- bind_rows(
  did_a_seifa_tidy %>% mutate(model = "SEIFA + sixty_day"),
  did_a_remote_tidy %>% mutate(model = "Remoteness + sixty_day")
)
write_csv(did_a, file.path(tab_dir, "sensitivity_sixty_day_did.csv"))
message("  Saved sensitivity_sixty_day_did.csv")

# ==============================================================================
# 4. SCENARIO B — COVID PERIOD EXCLUSION
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SCENARIO B: COVID period exclusion")
message(paste(rep("=", 60), collapse = ""))

# --- B.1: ITS overall excluding COVID months ---
message("\n  Fitting ITS overall (COVID months removed)...")
overall_nocovid <- overall_ts %>% filter(covid_period == 0)
model_b_overall <- fit_its_model_nocovid(overall_nocovid)
results_b_overall <- tidy_its_results(model_b_overall, "B: Overall no COVID")

b2_b <- results_b_overall %>% filter(term == "post_reform")
message("    B2: ", sprintf("%+.6f", b2_b$estimate),
        " (p = ", sprintf("%.4f", b2_b$p_value), ")")

summary_rows[["b_its_b2"]] <- add_summary(
  "B: COVID excl", "ITS Overall B2",
  b2_b$estimate, b2_b$ci_low, b2_b$ci_high, b2_b$p_value)

# --- B.2: ITS by age group excluding COVID ---
message("  Fitting ITS by age group (COVID months removed)...")
results_b_age <- list()
for (ag in c("19-64", "65+")) {
  ag_nocovid <- age_ts[[ag]] %>% filter(covid_period == 0)
  mod <- fit_its_model_nocovid(ag_nocovid)
  results_b_age[[ag]] <- tidy_its_results(mod, paste0("B: ", ag, " no COVID"))
  b2 <- results_b_age[[ag]] %>% filter(term == "post_reform")
  message("    ", ag, " B2: ", sprintf("%+.6f", b2$estimate),
          " (p = ", sprintf("%.4f", b2$p_value), ")")
}

b2_b_1964 <- results_b_age[["19-64"]] %>% filter(term == "post_reform")
summary_rows[["b_its_1964"]] <- add_summary(
  "B: COVID excl", "ITS 19-64 B2",
  b2_b_1964$estimate, b2_b_1964$ci_low, b2_b_1964$ci_high, b2_b_1964$p_value)

# Save ITS results
its_b <- bind_rows(results_b_overall, bind_rows(results_b_age))
write_csv(its_b, file.path(tab_dir, "sensitivity_covid_excluded_its.csv"))
message("  Saved sensitivity_covid_excluded_its.csv")

# --- B.3: CausalImpact with pre-period ending before COVID ---
message("  Fitting CausalImpact (pre-period: May 2002 - Feb 2020)...")
pre_period_b <- as.Date(c("2002-05-01", "2020-02-01"))
impact_b <- fit_causal_impact(overall_zoo, pre_period_b, post_period)
results_b_ci <- tidy_ci_results(impact_b, "B: CI no COVID")

message("    Relative effect: ", sprintf("%.2f%%", results_b_ci$rel_effect_pct),
        " [", sprintf("%.2f%%", results_b_ci$rel_effect_lower), ", ",
        sprintf("%.2f%%", results_b_ci$rel_effect_upper), "]",
        " (p = ", sprintf("%.4f", results_b_ci$p_value), ")")

summary_rows[["b_ci"]] <- add_summary(
  "B: COVID excl", "CI Overall Rel%",
  results_b_ci$rel_effect_pct, results_b_ci$rel_effect_lower,
  results_b_ci$rel_effect_upper, results_b_ci$p_value)

write_csv(results_b_ci, file.path(tab_dir, "sensitivity_covid_excluded_ci.csv"))
message("  Saved sensitivity_covid_excluded_ci.csv")

# --- B.4: DiD SEIFA excluding COVID months ---
message("  Fitting DiD SEIFA (COVID months removed)...")
lga_nocovid <- lga_panel %>% filter(covid_period == 0)

did_b_seifa <- feols(
  scripts_per_erp ~ post_reform * seifa_quintile +
    time_index + sin1 + cos1 + sin2 + cos2 | lga_code,
  data = lga_nocovid,
  vcov = ~lga_code
)

did_b_seifa_tidy <- broom::tidy(did_b_seifa, conf.int = TRUE)
q1_b <- did_b_seifa_tidy %>% filter(term == "post_reform:seifa_quintile1")
message("    Q1 interaction: ", sprintf("%+.6f", q1_b$estimate),
        " (p = ", sprintf("%.4f", q1_b$p.value), ")")

summary_rows[["b_did_q1"]] <- add_summary(
  "B: COVID excl", "DiD SEIFA Q1",
  q1_b$estimate, q1_b$conf.low, q1_b$conf.high, q1_b$p.value)

write_csv(did_b_seifa_tidy %>% mutate(model = "SEIFA no COVID"),
          file.path(tab_dir, "sensitivity_covid_excluded_did.csv"))
message("  Saved sensitivity_covid_excluded_did.csv")

# ==============================================================================
# 5. SCENARIO C — PLACEBO INTERVENTION DATES
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SCENARIO C: Placebo intervention dates")
message(paste(rep("=", 60), collapse = ""))

placebo_dates <- as.Date(c("2019-01-01", "2020-01-01", "2021-01-01",
                            "2022-01-01", "2023-01-01"))
placebo_labels <- c("Jan 2019", "Jan 2020", "Jan 2021",
                     "Jan 2022", "Jan 2023 (real)")

placebo_results <- list()

for (i in seq_along(placebo_dates)) {
  pd <- placebo_dates[i]
  label <- placebo_labels[i]
  message("  Testing: ", label, "...")

  # Recompute ITS variables for this placebo date
  # Matches Python 02_data_cleaning.py: round(days_since_intervention / 30.44)
  placebo_ts <- overall_ts %>%
    mutate(
      post_reform = as.integer(date >= pd),
      reform_days = as.numeric(difftime(date, pd, units = "days")),
      time_post_reform = if_else(post_reform == 1L,
                                  as.integer(round(reform_days / 30.44)), 0L)
    ) %>%
    select(-reform_days)

  mod <- fit_its_model(placebo_ts)
  res <- tidy_its_results(mod, label)

  b2 <- res %>% filter(term == "post_reform")
  b3 <- res %>% filter(term == "time_post_reform")

  message("    B2: ", sprintf("%+.6f", b2$estimate),
          " (p = ", sprintf("%.4f", b2$p_value), ")",
          if_else(b2$p_value < 0.05, " *", ""))

  placebo_results[[label]] <- tibble(
    placebo_date  = pd,
    placebo_label = label,
    is_real       = (pd == intervention_date),
    b2_estimate   = b2$estimate,
    b2_ci_low     = b2$ci_low,
    b2_ci_high    = b2$ci_high,
    b2_p_value    = b2$p_value,
    b3_estimate   = b3$estimate,
    b3_ci_low     = b3$ci_low,
    b3_ci_high    = b3$ci_high,
    b3_p_value    = b3$p_value,
    r_squared     = b2$r_squared,
    n_obs         = b2$n_obs
  )
}

placebo_df <- bind_rows(placebo_results)
write_csv(placebo_df, file.path(tab_dir, "sensitivity_placebo_tests.csv"))
message("  Saved sensitivity_placebo_tests.csv")

# --- Fig 26: Placebo test B2 estimates ---
fig26_data <- placebo_df %>%
  mutate(placebo_label = fct_inorder(placebo_label))

fig26 <- ggplot(fig26_data, aes(x = b2_estimate, y = placebo_label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = b2_ci_low, xmax = b2_ci_high, colour = is_real),
                width = 0.25, linewidth = 0.6) +
  geom_point(aes(colour = is_real, shape = b2_p_value < 0.05), size = 3) +
  scale_colour_manual(
    values = c("TRUE" = "red", "FALSE" = "steelblue"),
    labels = c("TRUE" = "Real intervention", "FALSE" = "Placebo"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "Not significant"),
    name = NULL
  ) +
  labs(
    title = expression(paste(
      "Placebo Tests: ITS Level Change (", beta[2],
      ") at Alternative Intervention Dates")),
    subtitle = paste0(
      "Significant B2 at placebo dates suggests pre-existing ",
      "structural breaks, not reform-specific"),
    x = expression(paste("Level change (", beta[2],
                          ", scripts per ERP)")),
    y = "Intervention date tested"
  ) +
  theme(axis.text.y = element_text(size = 11))

ggsave(file.path(fig_dir, "fig26_placebo_tests.png"),
       fig26, width = 10, height = 5, dpi = 300)
message("  Saved fig26_placebo_tests.png")

# ==============================================================================
# 6. SCENARIO D — SHORTER PRE-PERIOD (2018+)
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SCENARIO D: Shorter pre-period (2018+)")
message(paste(rep("=", 60), collapse = ""))

# --- D.1: ITS overall on 2018+ data ---
message("\n  Fitting ITS overall (2018+ only)...")
overall_short <- overall_ts %>% filter(date >= as.Date("2018-01-01"))

model_d_overall <- fit_its_model(overall_short)
results_d_overall <- tidy_its_results(model_d_overall, "D: Overall 2018+")

b2_d <- results_d_overall %>% filter(term == "post_reform")
message("    B2: ", sprintf("%+.6f", b2_d$estimate),
        " (p = ", sprintf("%.4f", b2_d$p_value), ")")
message("    Observations: ", b2_d$n_obs, " (vs ", base_its_b2$n_obs, " base)")

summary_rows[["d_its_b2"]] <- add_summary(
  "D: Short pre", "ITS Overall B2",
  b2_d$estimate, b2_d$ci_low, b2_d$ci_high, b2_d$p_value)

# --- D.2: CausalImpact with shorter pre-period ---
message("  Fitting CausalImpact (pre-period: Jan 2018 - Dec 2022)...")
pre_period_d <- as.Date(c("2018-01-01", "2022-12-01"))

impact_d <- fit_causal_impact(overall_zoo, pre_period_d, post_period)
results_d_ci <- tidy_ci_results(impact_d, "D: CI 2018+ pre")

message("    Relative effect: ", sprintf("%.2f%%", results_d_ci$rel_effect_pct),
        " [", sprintf("%.2f%%", results_d_ci$rel_effect_lower), ", ",
        sprintf("%.2f%%", results_d_ci$rel_effect_upper), "]",
        " (p = ", sprintf("%.4f", results_d_ci$p_value), ")")

summary_rows[["d_ci"]] <- add_summary(
  "D: Short pre", "CI Overall Rel%",
  results_d_ci$rel_effect_pct, results_d_ci$rel_effect_lower,
  results_d_ci$rel_effect_upper, results_d_ci$p_value)

# Save combined ITS + CI results
sensitivity_d <- bind_rows(
  results_d_overall %>%
    filter(term == "post_reform") %>%
    transmute(model_label, metric = "ITS_B2",
              estimate, ci_low, ci_high, p_value),
  tibble(model_label = results_d_ci$model_label, metric = "CI_rel_pct",
         estimate = results_d_ci$rel_effect_pct,
         ci_low = results_d_ci$rel_effect_lower,
         ci_high = results_d_ci$rel_effect_upper,
         p_value = results_d_ci$p_value)
)
write_csv(sensitivity_d, file.path(tab_dir, "sensitivity_short_preperiod.csv"))
message("  Saved sensitivity_short_preperiod.csv")

# ==============================================================================
# 7. SCENARIO E — CAUSALIMPACT CORRELATION FILTERING
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SCENARIO E: CausalImpact correlation filtering")
message(paste(rep("=", 60), collapse = ""))

# Compute pre-period correlations for ALL ATC classes
message("\n  Computing pre-period correlations for all ATC classes...")

safe_cor <- function(x, y) {
  complete <- complete.cases(x, y)
  if (sum(complete) < 3) return(NA_real_)
  cor(x[complete], y[complete])
}

all_pre_cors <- national %>%
  filter(age_group %in% c("19-64", "65+"), date < intervention_date) %>%
  select(atc_name, date, age_group, scripts_per_erp) %>%
  pivot_wider(names_from = age_group, values_from = scripts_per_erp) %>%
  group_by(atc_name) %>%
  summarise(
    pre_correlation = safe_cor(`19-64`, `65+`),
    .groups = "drop"
  ) %>%
  filter(!is.na(pre_correlation))

message("  Computed correlations for ", nrow(all_pre_cors), " ATC classes")

# Load ATC-level CausalImpact results from Phase 3
ci_atc <- read_csv(file.path(tab_dir, "ci_by_atc_results.csv"),
                   show_col_types = FALSE)
message("  Loaded CI results for ", nrow(ci_atc), " ATC classes")

# Merge correlations with CI results
ci_with_cors <- ci_atc %>%
  left_join(all_pre_cors, by = c("model_label" = "atc_name")) %>%
  filter(!is.na(pre_correlation))
message("  Matched: ", nrow(ci_with_cors), " classes with correlations")

# Filter at thresholds and summarise
thresholds <- c(0, 0.3, 0.5, 0.7)
threshold_labels <- c("All (no filter)", "r >= 0.3", "r >= 0.5", "r >= 0.7")

corr_filter_results <- map2_dfr(thresholds, threshold_labels,
                                 function(th, label) {
  filtered <- ci_with_cors %>% filter(pre_correlation >= th)

  if (nrow(filtered) == 0) {
    return(tibble(
      threshold = label, min_correlation = th,
      n_classes = 0L, n_fdr_significant = 0L,
      mean_abs_effect = NA_real_, median_abs_effect = NA_real_,
      mean_rel_effect_pct = NA_real_, median_rel_effect_pct = NA_real_,
      weighted_avg_rel_pct = NA_real_, mean_p_value = NA_real_
    ))
  }

  has_weights <- sum(!is.na(filtered$mean_scripts_erp)) > 0
  tibble(
    threshold             = label,
    min_correlation       = th,
    n_classes             = nrow(filtered),
    n_fdr_significant     = sum(filtered$p_fdr < 0.05, na.rm = TRUE),
    mean_abs_effect       = mean(filtered$abs_effect, na.rm = TRUE),
    median_abs_effect     = median(filtered$abs_effect, na.rm = TRUE),
    mean_rel_effect_pct   = mean(filtered$rel_effect_pct, na.rm = TRUE),
    median_rel_effect_pct = median(filtered$rel_effect_pct, na.rm = TRUE),
    weighted_avg_rel_pct  = if (has_weights)
      weighted.mean(filtered$rel_effect_pct, filtered$mean_scripts_erp,
                    na.rm = TRUE) else NA_real_,
    mean_p_value          = mean(filtered$p_value, na.rm = TRUE)
  )
})

write_csv(corr_filter_results,
          file.path(tab_dir, "sensitivity_ci_correlation_filter.csv"))
message("  Saved sensitivity_ci_correlation_filter.csv")

# Print summary
message("\n  Correlation filtering summary:")
for (i in seq_len(nrow(corr_filter_results))) {
  row <- corr_filter_results[i, ]
  message("    ", sprintf("%-18s", row$threshold),
          " n=", sprintf("%-3d", row$n_classes),
          " sig=", sprintf("%-2d", row$n_fdr_significant),
          " mean_rel=", sprintf("%+.2f%%", row$mean_rel_effect_pct),
          " wtd_avg=", sprintf("%+.2f%%", row$weighted_avg_rel_pct))
}

# --- Fig 27: Effect distribution at each correlation threshold ---
fig27_data <- map2_dfr(thresholds, threshold_labels, function(th, label) {
  ci_with_cors %>%
    filter(pre_correlation >= th) %>%
    mutate(threshold = label)
}) %>%
  mutate(threshold = factor(threshold, levels = threshold_labels))

fig27 <- ggplot(fig27_data, aes(x = threshold, y = rel_effect_pct)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(fill = "steelblue", alpha = 0.15, outlier.shape = NA,
               width = 0.5) +
  geom_jitter(aes(colour = p_fdr < 0.05), width = 0.15, size = 1.5,
              alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
               colour = "red") +
  scale_colour_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "grey60"),
    labels = c("TRUE" = "FDR p < 0.05", "FALSE" = "Not significant"),
    name = NULL
  ) +
  labs(
    title = "CausalImpact Relative Effects by Pre-Period Correlation Threshold",
    subtitle = paste0("Each dot = one ATC class | Red diamond = mean | ",
                      "Filtering removes poorly-matched control series"),
    x = "Pre-period correlation threshold",
    y = "Relative causal effect (%)"
  ) +
  theme(axis.text.x = element_text(size = 10))

ggsave(file.path(fig_dir, "fig27_ci_correlation_filter.png"),
       fig27, width = 10, height = 6, dpi = 300)
message("  Saved fig27_ci_correlation_filter.png")

# ==============================================================================
# 8. SCENARIO F — EQUITY WINDOW ROBUSTNESS
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SCENARIO F: Equity window robustness")
message(paste(rep("=", 60), collapse = ""))

# --- F.a: Full series 2002-2025 ---
message("\n  Fitting DiD SEIFA on full series (2002-2025)...")

did_f_full <- feols(
  scripts_per_erp ~ post_reform * seifa_quintile +
    time_index + covid_period + sin1 + cos1 + sin2 + cos2 | lga_code,
  data = lga_full_panel,
  vcov = ~lga_code
)

did_f_full_tidy <- broom::tidy(did_f_full, conf.int = TRUE)
q1_f_full <- did_f_full_tidy %>% filter(term == "post_reform:seifa_quintile1")
message("    Full series Q1: ", sprintf("%+.6f", q1_f_full$estimate),
        " (p = ", sprintf("%.4f", q1_f_full$p.value), ")")
message("    Observations: ", format(nobs(did_f_full), big.mark = ","))

summary_rows[["f_did_q1_full"]] <- add_summary(
  "F: Full 2002+", "DiD SEIFA Q1",
  q1_f_full$estimate, q1_f_full$conf.low, q1_f_full$conf.high,
  q1_f_full$p.value)

# --- F.b: Short window 2021-2025 ---
message("  Fitting DiD SEIFA on short window (2021-2025)...")

lga_short <- lga_full_panel %>%
  filter(date >= as.Date("2021-01-01"), date <= as.Date("2025-12-01"))

did_f_short <- feols(
  scripts_per_erp ~ post_reform * seifa_quintile +
    time_index + covid_period + sin1 + cos1 + sin2 + cos2 | lga_code,
  data = lga_short,
  vcov = ~lga_code
)

did_f_short_tidy <- broom::tidy(did_f_short, conf.int = TRUE)
q1_f_short <- did_f_short_tidy %>%
  filter(term == "post_reform:seifa_quintile1")
message("    Short window Q1: ", sprintf("%+.6f", q1_f_short$estimate),
        " (p = ", sprintf("%.4f", q1_f_short$p.value), ")")
message("    Observations: ", format(nobs(did_f_short), big.mark = ","))

summary_rows[["f_did_q1_short"]] <- add_summary(
  "F: Short 2021+", "DiD SEIFA Q1",
  q1_f_short$estimate, q1_f_short$conf.low, q1_f_short$conf.high,
  q1_f_short$p.value)

# Save combined equity window results
equity_windows <- bind_rows(
  did_f_full_tidy %>% mutate(window = "Full 2002-2025"),
  did_f_short_tidy %>% mutate(window = "Short 2021-2025")
)
write_csv(equity_windows, file.path(tab_dir, "sensitivity_equity_windows.csv"))
message("  Saved sensitivity_equity_windows.csv")

# ==============================================================================
# 9. SUMMARY COMPARISON TABLE + FIGURES
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY COMPARISON TABLE")
message(paste(rep("=", 60), collapse = ""))

# Build master table
summary_table <- bind_rows(summary_rows) %>%
  mutate(
    scenario = factor(scenario,
      levels = c("Base", "A: 60-day", "B: COVID excl",
                 "D: Short pre", "F: Full 2002+", "F: Short 2021+")),
    parameter = factor(parameter,
      levels = c("ITS Overall B2", "ITS 19-64 B2",
                 "CI Overall Rel%", "DiD SEIFA Q1", "DiD Remote"))
  )

write_csv(summary_table, file.path(tab_dir, "sensitivity_summary_table.csv"))
message("  Saved sensitivity_summary_table.csv")

# --- Fig 28: Parameter stability across scenarios ---
fig28 <- ggplot(summary_table, aes(x = estimate, y = fct_rev(scenario))) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.2, colour = "grey50", linewidth = 0.4) +
  geom_point(aes(colour = significant, shape = significant), size = 2.5) +
  scale_colour_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "grey60"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "Not significant"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 1),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "Not significant"),
    name = NULL
  ) +
  facet_wrap(~ parameter, scales = "free_x", ncol = 1) +
  labs(
    title = "Sensitivity Analysis: Parameter Stability Across Robustness Scenarios",
    subtitle = paste0("Estimates with 95% CIs | Consistent direction and ",
                      "significance = robust finding"),
    x = "Estimate (parameter-specific units)",
    y = NULL
  ) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 9))

ggsave(file.path(fig_dir, "fig28_sensitivity_summary.png"),
       fig28, width = 10, height = 12, dpi = 300)
message("  Saved fig28_sensitivity_summary.png")

# ==============================================================================
# 10. CONSOLE SUMMARY
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("PHASE 5: SENSITIVITY & ROBUSTNESS ANALYSES COMPLETE")
message(paste(rep("=", 70), collapse = ""))

message("\nTables saved (10):")
message("  outputs/tables/sensitivity_sixty_day_its.csv")
message("  outputs/tables/sensitivity_sixty_day_did.csv")
message("  outputs/tables/sensitivity_covid_excluded_its.csv")
message("  outputs/tables/sensitivity_covid_excluded_ci.csv")
message("  outputs/tables/sensitivity_covid_excluded_did.csv")
message("  outputs/tables/sensitivity_placebo_tests.csv")
message("  outputs/tables/sensitivity_short_preperiod.csv")
message("  outputs/tables/sensitivity_ci_correlation_filter.csv")
message("  outputs/tables/sensitivity_equity_windows.csv")
message("  outputs/tables/sensitivity_summary_table.csv")

message("\nFigures saved (3):")
message("  outputs/figures/fig26_placebo_tests.png")
message("  outputs/figures/fig27_ci_correlation_filter.png")
message("  outputs/figures/fig28_sensitivity_summary.png")

# --- Print master comparison table ---
message("\n--- Master comparison table ---")
print(
  summary_table %>%
    mutate(across(where(is.numeric), ~ round(., 6))) %>%
    select(scenario, parameter, estimate, ci_low, ci_high, p_value, significant) %>%
    as.data.frame(),
  right = FALSE
)

# --- Key robustness findings ---
message("\n--- Key robustness findings ---")

# 60-day adjustment
message("\n  1. ITS Overall B2 (60-day adjustment):")
message("     Base:      ", sprintf("%+.6f", base_its_b2$estimate),
        " (p = ", sprintf("%.4f", base_its_b2$p_value), ")")
message("     + 60-day:  ", sprintf("%+.6f", b2_a$estimate),
        " (p = ", sprintf("%.4f", b2_a$p_value), ")")
message("     sixty_day: ", sprintf("%+.6f", sixty_a$estimate),
        " (p = ", sprintf("%.4f", sixty_a$p_value), ")")
if (sign(base_its_b2$estimate) != sign(b2_a$estimate)) {
  message("     -> Direction CHANGED after 60-day adjustment")
} else if (base_its_b2$p_value < 0.05 & b2_a$p_value >= 0.05) {
  message("     -> Significance LOST after 60-day adjustment")
} else {
  message("     -> Finding STABLE after 60-day adjustment")
}

# Placebo tests
message("\n  2. Placebo tests (ITS B2 significance):")
for (i in seq_len(nrow(placebo_df))) {
  row <- placebo_df[i, ]
  sig_flag <- if_else(row$b2_p_value < 0.05,
                      " * SIGNIFICANT", "   not significant")
  message("     ", sprintf("%-20s", row$placebo_label),
          " B2 = ", sprintf("%+.6f", row$b2_estimate), sig_flag)
}
n_placebo_sig <- sum(placebo_df$b2_p_value[!placebo_df$is_real] < 0.05)
if (n_placebo_sig == 0) {
  message("     -> PASSED: No significant B2 at any placebo date ",
          "(temporal specificity confirmed)")
} else {
  message("     -> WARNING: B2 significant at ", n_placebo_sig,
          " placebo date(s) - pre-existing structural breaks detected")
}

# CausalImpact COVID exclusion
message("\n  3. CausalImpact COVID exclusion:")
message("     Base:      ", sprintf("%.2f%%", base_ci$rel_effect_pct),
        " (p = ", sprintf("%.4f", base_ci$p_value), ")")
message("     No COVID:  ", sprintf("%.2f%%", results_b_ci$rel_effect_pct),
        " (p = ", sprintf("%.4f", results_b_ci$p_value), ")")
message("     Short pre: ", sprintf("%.2f%%", results_d_ci$rel_effect_pct),
        " (p = ", sprintf("%.4f", results_d_ci$p_value), ")")

# SEIFA Q1 across equity windows
message("\n  4. DiD SEIFA Q1 across windows:")
message("     Base (2018-2025):  ", sprintf("%+.6f", base_did_q1$estimate),
        " (p = ", sprintf("%.4f", base_did_q1$p.value), ")")
message("     Full (2002-2025):  ", sprintf("%+.6f", q1_f_full$estimate),
        " (p = ", sprintf("%.4f", q1_f_full$p.value), ")")
message("     Short (2021-2025): ", sprintf("%+.6f", q1_f_short$estimate),
        " (p = ", sprintf("%.4f", q1_f_short$p.value), ")")

# DiD Remote with 60-day
message("\n  5. DiD Remote (60-day adjustment):")
message("     Base:     ", sprintf("%+.6f", base_did_remote_row$estimate),
        " (p = ", sprintf("%.4f", base_did_remote_row$p.value), ")")
message("     + 60-day: ", sprintf("%+.6f", remote_a$estimate),
        " (p = ", sprintf("%.4f", remote_a$p.value), ")")

# Overall robustness assessment
message("\n  Overall robustness assessment:")
n_stable <- sum(summary_table$significant[summary_table$scenario == "Base"] ==
                  summary_table$significant[summary_table$scenario != "Base"])
message("     Parameters tested across ", n_distinct(summary_table$scenario) - 1,
        " sensitivity scenarios")
message("     See sensitivity_summary_table.csv for full comparison")

message("\nNext: Update PROGRESS.md, then proceed to Phase 6 (manuscript)")
