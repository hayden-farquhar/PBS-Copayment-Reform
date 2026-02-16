# ==============================================================================
# 05_equity_analysis.R
# Phase 4: Equity Analysis — SEIFA Disadvantage and Remoteness
#
# Tests whether the January 2023 PBS copayment reduction ($42.50 → $30.00) had
# distributional effects — benefiting disadvantaged or remote areas even if the
# overall average causal effect was null (+0.3%, p = 0.49 from Phase 3).
#
# Approach: Difference-in-differences panel regressions using fixest::feols()
# with LGA fixed effects and clustered standard errors.
#
# Models:
#   1. DiD by SEIFA quintile (19-64 only, Q5 = reference)
#   2. DiD by remoteness (19-64 only, Major Cities = reference)
#   3. Triple-difference (19-64 vs 65+ × SEIFA, within-LGA control)
#   4. SEIFA × remoteness interaction (exploratory)
#   5. Equity impact ratios with bootstrapped CIs
#   6. Choropleth maps and spatial autocorrelation
#
# Inputs:  data/processed/pbs_lga_monthly.csv (2.9 GB)
#          data/processed/lga_characteristics.csv
#          data/spatial/LGA_2021/LGA_2021_AUST_GDA2020.shp
#
# Outputs: outputs/tables/descriptive_seifa_pre_post.csv
#          outputs/tables/descriptive_remoteness_pre_post.csv
#          outputs/tables/did_seifa_results.csv
#          outputs/tables/did_remoteness_results.csv
#          outputs/tables/did_triple_diff_results.csv
#          outputs/tables/did_seifa_remoteness_interaction.csv
#          outputs/tables/equity_impact_ratios.csv
#          outputs/tables/lga_reform_effects.csv
#          outputs/tables/morans_i_results.csv
#          outputs/figures/fig18_seifa_time_series.png
#          outputs/figures/fig19_remoteness_time_series.png
#          outputs/figures/fig20_did_seifa_coefplot.png
#          outputs/figures/fig21_did_remoteness_coefplot.png
#          outputs/figures/fig22_triple_diff_heatmap.png
#          outputs/figures/fig23_equity_impact_ratios.png
#          outputs/figures/fig24_choropleth_reform_effects.png
#          outputs/figures/fig25_choropleth_seifa.png
# ==============================================================================

# ==============================================================================
# 1. SETUP
# ==============================================================================

library(tidyverse)
library(lubridate)
library(data.table)
library(fixest)
library(broom)
library(scales)
library(sf)
library(spdep)

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

# Analysis window: 4 years pre + 3 years post (avoids pre-2018 structural breaks)
window_start <- as.Date("2018-01-01")
window_end   <- as.Date("2025-12-01")

# ==============================================================================
# 2. LOAD AND AGGREGATE DATA
# ==============================================================================
message("Loading LGA panel data (2.9 GB — this may take a minute)...")

lga_raw <- fread(
  file.path(proc_dir, "pbs_lga_monthly.csv"),
  select = c("lga_code", "date", "age_group", "scripts_per_erp",
             "time_index", "post_reform", "covid_period", "month_num")
)

message("  Raw rows loaded: ", format(nrow(lga_raw), big.mark = ","))

# Filter to analysis window and relevant age groups
lga_raw <- lga_raw[
  date >= as.character(window_start) &
  date <= as.character(window_end) &
  age_group %in% c("19-64", "65+")
]
message("  After filtering to 2018-2025, ages 19-64 & 65+: ",
        format(nrow(lga_raw), big.mark = ","), " rows")

# Aggregate across all ATC classes: sum scripts_per_erp within LGA × month × age group
lga_agg <- lga_raw[,
  .(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE)),
  by = .(lga_code, date, age_group, time_index, post_reform, covid_period, month_num)
]
message("  After aggregation across ATC classes: ",
        format(nrow(lga_agg), big.mark = ","), " rows")

# Free raw data
rm(lga_raw)
gc()

# Convert to tibble and parse date
lga_panel <- as_tibble(lga_agg) %>%
  mutate(date = as.Date(date),
         lga_code = as.integer(lga_code))
rm(lga_agg)

# --- Add Fourier terms (K=2) ---
lga_panel <- lga_panel %>%
  mutate(
    sin1 = sin(2 * pi * 1 * month_num / 12),
    cos1 = cos(2 * pi * 1 * month_num / 12),
    sin2 = sin(2 * pi * 2 * month_num / 12),
    cos2 = cos(2 * pi * 2 * month_num / 12)
  )

# --- Add age group interaction variables ---
lga_panel <- lga_panel %>%
  mutate(
    is_general       = as.integer(age_group == "19-64"),
    post_x_general   = post_reform * is_general
  )

# --- Load LGA characteristics ---
message("Loading LGA characteristics...")
lga_chars <- read_csv(file.path(proc_dir, "lga_characteristics.csv"),
                      show_col_types = FALSE) %>%
  mutate(lga_code = as.integer(lga_code))

# Merge SEIFA quintile and remoteness onto panel
lga_panel <- lga_panel %>%
  left_join(lga_chars %>% select(lga_code, seifa_quintile, remoteness_area,
                                  remoteness_code),
            by = "lga_code")

# Drop rows without SEIFA/remoteness classification
n_before <- nrow(lga_panel)
lga_panel <- lga_panel %>% filter(!is.na(seifa_quintile), !is.na(remoteness_area))
message("  Dropped ", n_before - nrow(lga_panel), " rows without SEIFA/remoteness match")

# Convert SEIFA quintile to factor (Q5 = least disadvantaged = reference)
lga_panel <- lga_panel %>%
  mutate(
    seifa_quintile = factor(seifa_quintile, levels = c(5, 1, 2, 3, 4)),
    remoteness_area = factor(remoteness_area,
                             levels = c("Major Cities of Australia",
                                        "Inner Regional Australia",
                                        "Outer Regional Australia",
                                        "Remote Australia",
                                        "Very Remote Australia"))
  )

# Collapsed remoteness (3 categories) — avoids sparse cells
lga_panel <- lga_panel %>%
  mutate(
    remoteness_3cat = fct_collapse(
      remoteness_area,
      "Major Cities"  = "Major Cities of Australia",
      "Regional"      = c("Inner Regional Australia", "Outer Regional Australia"),
      "Remote"        = c("Remote Australia", "Very Remote Australia")
    )
  )

message("  Final panel: ", format(nrow(lga_panel), big.mark = ","), " rows, ",
        n_distinct(lga_panel$lga_code), " LGAs, ",
        n_distinct(lga_panel$date), " months")

# --- Load shapefile ---
message("Loading LGA shapefile...")
lga_sf <- st_read(file.path(proj_dir, "data", "spatial", "LGA_2021",
                             "LGA_2021_AUST_GDA2020.shp"),
                  quiet = TRUE) %>%
  mutate(lga_code = as.integer(LGA_CODE21))

message("  Shapefile: ", nrow(lga_sf), " features")

# ==============================================================================
# 3. DESCRIPTIVE STATISTICS
# ==============================================================================
message("\n--- Descriptive statistics ---")

# --- 3a. Pre/post means by SEIFA quintile × age group ---
desc_seifa <- lga_panel %>%
  mutate(period = if_else(post_reform == 1, "Post", "Pre")) %>%
  group_by(seifa_quintile, age_group, period) %>%
  summarise(
    mean_scripts   = mean(scripts_per_erp, na.rm = TRUE),
    sd_scripts     = sd(scripts_per_erp, na.rm = TRUE),
    median_scripts = median(scripts_per_erp, na.rm = TRUE),
    n_lga_months   = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = period,
    values_from = c(mean_scripts, sd_scripts, median_scripts, n_lga_months),
    names_sep = "_"
  ) %>%
  mutate(
    pct_change = (mean_scripts_Post - mean_scripts_Pre) / mean_scripts_Pre * 100
  ) %>%
  arrange(age_group, seifa_quintile)

write_csv(desc_seifa, file.path(tab_dir, "descriptive_seifa_pre_post.csv"))
message("  Saved descriptive_seifa_pre_post.csv")

# --- 3b. Pre/post means by remoteness × age group ---
desc_remote <- lga_panel %>%
  mutate(period = if_else(post_reform == 1, "Post", "Pre")) %>%
  group_by(remoteness_area, age_group, period) %>%
  summarise(
    mean_scripts   = mean(scripts_per_erp, na.rm = TRUE),
    sd_scripts     = sd(scripts_per_erp, na.rm = TRUE),
    median_scripts = median(scripts_per_erp, na.rm = TRUE),
    n_lga_months   = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = period,
    values_from = c(mean_scripts, sd_scripts, median_scripts, n_lga_months),
    names_sep = "_"
  ) %>%
  mutate(
    pct_change = (mean_scripts_Post - mean_scripts_Pre) / mean_scripts_Pre * 100
  ) %>%
  arrange(age_group, remoteness_area)

write_csv(desc_remote, file.path(tab_dir, "descriptive_remoteness_pre_post.csv"))
message("  Saved descriptive_remoteness_pre_post.csv")

# --- 3c. Fig 18: Time series by SEIFA quintile (19-64 only) ---
seifa_ts <- lga_panel %>%
  filter(age_group == "19-64") %>%
  group_by(date, seifa_quintile) %>%
  summarise(scripts_per_erp = mean(scripts_per_erp, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(quintile_label = paste0("Q", seifa_quintile,
                                  if_else(seifa_quintile == "1", " (Most disadvantaged)",
                                          if_else(seifa_quintile == "5", " (Least disadvantaged)", ""))))

fig18 <- ggplot(seifa_ts, aes(x = date, y = scripts_per_erp,
                               colour = quintile_label)) +
  geom_line(linewidth = 0.5) +
  annotate("rect", xmin = covid_start, xmax = covid_end,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  geom_vline(xintercept = intervention_date,
             linetype = "dashed", colour = "red", linewidth = 0.6) +
  annotate("text", x = intervention_date + 60, y = Inf,
           label = "Copayment\nreduction", vjust = 1.5,
           colour = "red", size = 3, fontface = "italic") +
  annotate("text", x = covid_start + 150, y = Inf,
           label = "COVID", vjust = 1.5,
           colour = "grey50", size = 3, fontface = "italic") +
  scale_colour_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "PBS Prescription Rates by SEIFA Disadvantage Quintile (Ages 19-64)",
       subtitle = "Mean scripts per ERP across LGAs, Jan 2018 \u2013 Dec 2025",
       x = NULL, y = "Scripts per ERP (all ATC classes)", colour = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

ggsave(file.path(fig_dir, "fig18_seifa_time_series.png"),
       fig18, width = 11, height = 6, dpi = 300)
message("  Saved fig18_seifa_time_series.png")

# --- 3d. Fig 19: Time series by remoteness (19-64 only) ---
remote_ts <- lga_panel %>%
  filter(age_group == "19-64") %>%
  group_by(date, remoteness_area) %>%
  summarise(scripts_per_erp = mean(scripts_per_erp, na.rm = TRUE),
            .groups = "drop")

fig19 <- ggplot(remote_ts, aes(x = date, y = scripts_per_erp,
                                colour = remoteness_area)) +
  geom_line(linewidth = 0.5) +
  annotate("rect", xmin = covid_start, xmax = covid_end,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  geom_vline(xintercept = intervention_date,
             linetype = "dashed", colour = "red", linewidth = 0.6) +
  annotate("text", x = intervention_date + 60, y = Inf,
           label = "Copayment\nreduction", vjust = 1.5,
           colour = "red", size = 3, fontface = "italic") +
  annotate("text", x = covid_start + 150, y = Inf,
           label = "COVID", vjust = 1.5,
           colour = "grey50", size = 3, fontface = "italic") +
  scale_colour_brewer(palette = "YlOrRd") +
  labs(title = "PBS Prescription Rates by Remoteness Area (Ages 19-64)",
       subtitle = "Mean scripts per ERP across LGAs, Jan 2018 \u2013 Dec 2025",
       x = NULL, y = "Scripts per ERP (all ATC classes)", colour = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

ggsave(file.path(fig_dir, "fig19_remoteness_time_series.png"),
       fig19, width = 11, height = 6, dpi = 300)
message("  Saved fig19_remoteness_time_series.png")

# ==============================================================================
# 4. DiD MODEL 1 — SEIFA QUINTILE (PRIMARY)
# ==============================================================================
message("\n--- DiD Model 1: SEIFA quintile (19-64 only) ---")

panel_1964 <- lga_panel %>% filter(age_group == "19-64")

did_seifa <- feols(
  scripts_per_erp ~ post_reform * seifa_quintile +
    time_index + covid_period + sin1 + cos1 + sin2 + cos2 | lga_code,
  data = panel_1964,
  vcov = ~lga_code
)

message("  Model fitted: ", format(nobs(did_seifa), big.mark = ","), " obs, ",
        did_seifa$nparams, " params")

# Tidy results
did_seifa_tidy <- broom::tidy(did_seifa, conf.int = TRUE) %>%
  mutate(sig_label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.10  ~ ".",
    TRUE             ~ ""
  ))

write_csv(did_seifa_tidy, file.path(tab_dir, "did_seifa_results.csv"))
message("  Saved did_seifa_results.csv")

# Print key interaction terms
message("\n  SEIFA DiD interaction terms (post_reform × quintile, vs Q5):")
seifa_interactions <- did_seifa_tidy %>%
  filter(str_detect(term, "^post_reform:seifa_quintile"))
for (i in seq_len(nrow(seifa_interactions))) {
  row <- seifa_interactions[i, ]
  message("    ", sprintf("%-30s", row$term),
          " est = ", sprintf("%+.6f", row$estimate),
          " [", sprintf("%.6f", row$conf.low), ", ",
          sprintf("%.6f", row$conf.high), "]",
          " p = ", sprintf("%.4f", row$p.value), " ", row$sig_label)
}

# --- Fig 20: Coefficient plot of SEIFA interaction terms ---
seifa_coef_data <- seifa_interactions %>%
  mutate(
    quintile = str_extract(term, "\\d$"),
    quintile_label = paste0("Q", quintile, " vs Q5"),
    quintile_label = fct_reorder(quintile_label, as.numeric(quintile))
  )

fig20 <- ggplot(seifa_coef_data, aes(x = estimate, y = quintile_label)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.2, colour = "steelblue", linewidth = 0.7) +
  geom_point(aes(colour = p.value < 0.05), size = 3) +
  geom_text(aes(label = sig_label), nudge_y = 0.25, size = 4) +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey60"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "Not significant"),
                      name = NULL) +
  labs(title = "DiD: Differential Reform Effect by SEIFA Quintile (Ages 19-64)",
       subtitle = "Interaction: post_reform x seifa_quintile | Reference = Q5 (least disadvantaged) | LGA FE, clustered SEs",
       x = "Differential effect vs Q5 (scripts per ERP)",
       y = NULL) +
  theme(axis.text.y = element_text(size = 11))

ggsave(file.path(fig_dir, "fig20_did_seifa_coefplot.png"),
       fig20, width = 9, height = 5, dpi = 300)
message("  Saved fig20_did_seifa_coefplot.png")

# ==============================================================================
# 5. DiD MODEL 2 — REMOTENESS
# ==============================================================================
message("\n--- DiD Model 2: Remoteness (19-64 only) ---")

did_remote <- feols(
  scripts_per_erp ~ post_reform * remoteness_area +
    time_index + covid_period + sin1 + cos1 + sin2 + cos2 | lga_code,
  data = panel_1964,
  vcov = ~lga_code
)

message("  Model fitted: ", format(nobs(did_remote), big.mark = ","), " obs")

did_remote_tidy <- broom::tidy(did_remote, conf.int = TRUE) %>%
  mutate(sig_label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.10  ~ ".",
    TRUE             ~ ""
  ))

write_csv(did_remote_tidy, file.path(tab_dir, "did_remoteness_results.csv"))
message("  Saved did_remoteness_results.csv")

# Print key interaction terms
message("\n  Remoteness DiD interaction terms (vs Major Cities):")
remote_interactions <- did_remote_tidy %>%
  filter(str_detect(term, "^post_reform:remoteness_area"))
for (i in seq_len(nrow(remote_interactions))) {
  row <- remote_interactions[i, ]
  message("    ", sprintf("%-55s", row$term),
          " est = ", sprintf("%+.6f", row$estimate),
          " p = ", sprintf("%.4f", row$p.value), " ", row$sig_label)
}

# --- Fig 21: Coefficient plot of remoteness interaction terms ---
remote_coef_data <- remote_interactions %>%
  mutate(
    area = str_remove(term, "^post_reform:remoteness_area"),
    area = fct_inorder(area)
  )

fig21 <- ggplot(remote_coef_data, aes(x = estimate, y = area)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high),
                width = 0.2, colour = "darkorange", linewidth = 0.7) +
  geom_point(aes(colour = p.value < 0.05), size = 3) +
  geom_text(aes(label = sig_label), nudge_y = 0.25, size = 4) +
  scale_colour_manual(values = c("TRUE" = "darkorange", "FALSE" = "grey60"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "Not significant"),
                      name = NULL) +
  labs(title = "DiD: Differential Reform Effect by Remoteness (Ages 19-64)",
       subtitle = "Interaction: post_reform x remoteness | Reference = Major Cities | LGA FE, clustered SEs",
       x = "Differential effect vs Major Cities (scripts per ERP)",
       y = NULL) +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(fig_dir, "fig21_did_remoteness_coefplot.png"),
       fig21, width = 10, height = 5, dpi = 300)
message("  Saved fig21_did_remoteness_coefplot.png")

# ==============================================================================
# 6. DiD MODEL 3 — TRIPLE DIFFERENCE (GOLD STANDARD)
# ==============================================================================
message("\n--- DiD Model 3: Triple difference (19-64 vs 65+ x SEIFA) ---")

# Uses both age groups; 65+ within-LGA as control
# LGA x age_group FE absorbs all time-invariant within-group differences

did_triple <- feols(
  scripts_per_erp ~ post_reform * is_general * seifa_quintile +
    time_index * is_general + covid_period * is_general +
    sin1 + cos1 + sin2 + cos2 | lga_code^age_group,
  data = lga_panel,
  vcov = ~lga_code
)

message("  Model fitted: ", format(nobs(did_triple), big.mark = ","), " obs")

did_triple_tidy <- broom::tidy(did_triple, conf.int = TRUE) %>%
  mutate(sig_label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.10  ~ ".",
    TRUE             ~ ""
  ))

write_csv(did_triple_tidy, file.path(tab_dir, "did_triple_diff_results.csv"))
message("  Saved did_triple_diff_results.csv")

# Print triple-diff key coefficients
message("\n  Triple-diff key coefficients (post_reform x is_general x seifa):")
triple_key <- did_triple_tidy %>%
  filter(str_detect(term, "post_reform.*is_general.*seifa_quintile|post_reform:is_general$"))
for (i in seq_len(nrow(triple_key))) {
  row <- triple_key[i, ]
  message("    ", sprintf("%-50s", row$term),
          " est = ", sprintf("%+.6f", row$estimate),
          " p = ", sprintf("%.4f", row$p.value), " ", row$sig_label)
}

# --- Fig 22: Heatmap of reform effects by SEIFA quintile x age group ---
# Extract post_reform coefficients for each quintile x age_group combination
# We need to compute marginal effects from the triple-diff model
triple_coefs <- coef(did_triple)

# Build effects matrix
quintiles <- c("5", "1", "2", "3", "4")
age_labels <- c("19-64 (General)", "65+ (Concessional)")

heatmap_data <- expand_grid(
  quintile = quintiles,
  age_label = age_labels
) %>%
  mutate(
    is_gen = as.integer(age_label == "19-64 (General)"),
    # Base post_reform effect (for Q5, 65+)
    effect = triple_coefs["post_reform"],
    # Add is_general main effect if 19-64
    effect = effect + if_else(is_gen == 1,
                              triple_coefs["post_reform:is_general"], 0),
    # Add quintile interaction if not Q5
    effect = effect + case_when(
      quintile == "1" ~ triple_coefs["post_reform:seifa_quintile1"],
      quintile == "2" ~ triple_coefs["post_reform:seifa_quintile2"],
      quintile == "3" ~ triple_coefs["post_reform:seifa_quintile3"],
      quintile == "4" ~ triple_coefs["post_reform:seifa_quintile4"],
      TRUE ~ 0
    ),
    # Add triple interaction if 19-64 and not Q5
    effect = effect + case_when(
      is_gen == 1 & quintile == "1" ~ triple_coefs["post_reform:is_general:seifa_quintile1"],
      is_gen == 1 & quintile == "2" ~ triple_coefs["post_reform:is_general:seifa_quintile2"],
      is_gen == 1 & quintile == "3" ~ triple_coefs["post_reform:is_general:seifa_quintile3"],
      is_gen == 1 & quintile == "4" ~ triple_coefs["post_reform:is_general:seifa_quintile4"],
      TRUE ~ 0
    ),
    quintile_label = paste0("Q", quintile),
    quintile_label = fct_relevel(quintile_label, "Q1", "Q2", "Q3", "Q4", "Q5")
  )

fig22 <- ggplot(heatmap_data, aes(x = quintile_label, y = age_label, fill = effect)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.4f", effect)), size = 4) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                       midpoint = 0, name = "Effect\n(scripts/ERP)",
                       labels = label_number(accuracy = 0.01)) +
  labs(title = "Triple-Difference: Reform Effect by SEIFA Quintile and Age Group",
       subtitle = "Combined marginal effects from triple-diff model | LGA x age_group FE, clustered SEs",
       x = "SEIFA Quintile (Q1 = most disadvantaged)",
       y = NULL) +
  theme(axis.text = element_text(size = 11),
        panel.grid = element_blank())

ggsave(file.path(fig_dir, "fig22_triple_diff_heatmap.png"),
       fig22, width = 8, height = 4.5, dpi = 300)
message("  Saved fig22_triple_diff_heatmap.png")

# ==============================================================================
# 7. DiD MODEL 4 — SEIFA x REMOTENESS INTERACTION (EXPLORATORY)
# ==============================================================================
message("\n--- DiD Model 4: SEIFA x Remoteness interaction (19-64, exploratory) ---")

did_seifa_remote <- feols(
  scripts_per_erp ~ post_reform * seifa_quintile * remoteness_3cat +
    time_index + covid_period + sin1 + cos1 + sin2 + cos2 | lga_code,
  data = panel_1964,
  vcov = ~lga_code
)

message("  Model fitted: ", format(nobs(did_seifa_remote), big.mark = ","), " obs")

did_seifa_remote_tidy <- broom::tidy(did_seifa_remote, conf.int = TRUE) %>%
  mutate(sig_label = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    p.value < 0.10  ~ ".",
    TRUE             ~ ""
  ))

# Flag sparse cells
cell_counts <- panel_1964 %>%
  count(seifa_quintile, remoteness_3cat) %>%
  mutate(sparse = n < 100)
n_sparse <- sum(cell_counts$sparse)
if (n_sparse > 0) {
  message("  WARNING: ", n_sparse, " SEIFA x remoteness cells have < 100 obs")
}

write_csv(did_seifa_remote_tidy,
          file.path(tab_dir, "did_seifa_remoteness_interaction.csv"))
message("  Saved did_seifa_remoteness_interaction.csv")

# ==============================================================================
# 8. EQUITY IMPACT RATIOS
# ==============================================================================
message("\n--- Equity impact ratios ---")

# EIR = (post_rate_group / pre_rate_group) / (post_rate_ref / pre_rate_ref)
# Bootstrap 95% CIs by resampling LGAs

set.seed(42)
n_boot <- 1000

# --- 8a. SEIFA EIRs (19-64 only) ---
# Compute LGA-level pre/post means
lga_pre_post <- panel_1964 %>%
  mutate(period = if_else(post_reform == 1, "post", "pre")) %>%
  group_by(lga_code, seifa_quintile, period) %>%
  summarise(mean_scripts = mean(scripts_per_erp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_scripts,
              names_prefix = "scripts_")

# Reference group: Q5 aggregate mean
ref_pre  <- mean(lga_pre_post$scripts_pre[lga_pre_post$seifa_quintile == "5"], na.rm = TRUE)
ref_post <- mean(lga_pre_post$scripts_post[lga_pre_post$seifa_quintile == "5"], na.rm = TRUE)
ref_ratio <- ref_post / ref_pre

# Point estimates
seifa_eir <- lga_pre_post %>%
  group_by(seifa_quintile) %>%
  summarise(
    pre_mean  = mean(scripts_pre, na.rm = TRUE),
    post_mean = mean(scripts_post, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    group_ratio = post_mean / pre_mean,
    eir = group_ratio / ref_ratio,
    comparison = paste0("Q", seifa_quintile, " vs Q5"),
    dimension = "SEIFA"
  )

# Bootstrap CIs for SEIFA
boot_seifa_eir <- function() {
  # Resample LGAs with replacement
  lga_sample <- sample(unique(lga_pre_post$lga_code), replace = TRUE)
  boot_data <- lga_pre_post[lga_pre_post$lga_code %in% lga_sample, ]

  ref_pre_b  <- mean(boot_data$scripts_pre[boot_data$seifa_quintile == "5"], na.rm = TRUE)
  ref_post_b <- mean(boot_data$scripts_post[boot_data$seifa_quintile == "5"], na.rm = TRUE)
  ref_ratio_b <- ref_post_b / ref_pre_b

  boot_data %>%
    group_by(seifa_quintile) %>%
    summarise(pre = mean(scripts_pre, na.rm = TRUE),
              post = mean(scripts_post, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(eir = (post / pre) / ref_ratio_b) %>%
    select(seifa_quintile, eir)
}

message("  Bootstrapping SEIFA EIRs (", n_boot, " iterations)...")
boot_results_seifa <- map_dfr(1:n_boot, ~ boot_seifa_eir(), .id = "iter")

seifa_eir_ci <- boot_results_seifa %>%
  group_by(seifa_quintile) %>%
  summarise(
    eir_lower = quantile(eir, 0.025, na.rm = TRUE),
    eir_upper = quantile(eir, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

seifa_eir <- seifa_eir %>%
  left_join(seifa_eir_ci, by = "seifa_quintile")

# --- 8b. Remoteness EIRs ---
lga_pre_post_remote <- panel_1964 %>%
  mutate(period = if_else(post_reform == 1, "post", "pre")) %>%
  group_by(lga_code, remoteness_area, period) %>%
  summarise(mean_scripts = mean(scripts_per_erp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_scripts,
              names_prefix = "scripts_")

ref_pre_r  <- mean(lga_pre_post_remote$scripts_pre[
  lga_pre_post_remote$remoteness_area == "Major Cities of Australia"], na.rm = TRUE)
ref_post_r <- mean(lga_pre_post_remote$scripts_post[
  lga_pre_post_remote$remoteness_area == "Major Cities of Australia"], na.rm = TRUE)
ref_ratio_r <- ref_post_r / ref_pre_r

remote_eir <- lga_pre_post_remote %>%
  group_by(remoteness_area) %>%
  summarise(
    pre_mean  = mean(scripts_pre, na.rm = TRUE),
    post_mean = mean(scripts_post, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    group_ratio = post_mean / pre_mean,
    eir = group_ratio / ref_ratio_r,
    comparison = paste0(remoteness_area, " vs Major Cities"),
    dimension = "Remoteness"
  )

# Bootstrap CIs for remoteness
boot_remote_eir <- function() {
  lga_sample <- sample(unique(lga_pre_post_remote$lga_code), replace = TRUE)
  boot_data <- lga_pre_post_remote[lga_pre_post_remote$lga_code %in% lga_sample, ]

  ref_pre_b  <- mean(boot_data$scripts_pre[
    boot_data$remoteness_area == "Major Cities of Australia"], na.rm = TRUE)
  ref_post_b <- mean(boot_data$scripts_post[
    boot_data$remoteness_area == "Major Cities of Australia"], na.rm = TRUE)
  ref_ratio_b <- ref_post_b / ref_pre_b

  boot_data %>%
    group_by(remoteness_area) %>%
    summarise(pre = mean(scripts_pre, na.rm = TRUE),
              post = mean(scripts_post, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(eir = (post / pre) / ref_ratio_b) %>%
    select(remoteness_area, eir)
}

message("  Bootstrapping remoteness EIRs (", n_boot, " iterations)...")
boot_results_remote <- map_dfr(1:n_boot, ~ boot_remote_eir(), .id = "iter")

remote_eir_ci <- boot_results_remote %>%
  group_by(remoteness_area) %>%
  summarise(
    eir_lower = quantile(eir, 0.025, na.rm = TRUE),
    eir_upper = quantile(eir, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

remote_eir <- remote_eir %>%
  left_join(remote_eir_ci, by = "remoteness_area")

# Combine and save
all_eir <- bind_rows(
  seifa_eir %>% select(dimension, comparison, pre_mean, post_mean,
                        group_ratio, eir, eir_lower, eir_upper),
  remote_eir %>% select(dimension, comparison, pre_mean, post_mean,
                          group_ratio, eir, eir_lower, eir_upper)
)

write_csv(all_eir, file.path(tab_dir, "equity_impact_ratios.csv"))
message("  Saved equity_impact_ratios.csv")

# Print EIRs
message("\n  Equity Impact Ratios:")
for (i in seq_len(nrow(all_eir))) {
  row <- all_eir[i, ]
  message("    ", sprintf("%-45s", row$comparison),
          " EIR = ", sprintf("%.4f", row$eir),
          " [", sprintf("%.4f", row$eir_lower), ", ",
          sprintf("%.4f", row$eir_upper), "]")
}

# --- Fig 23: EIR bar chart ---
eir_plot_data <- all_eir %>%
  filter(comparison != "Q5 vs Q5",
         comparison != "Major Cities of Australia vs Major Cities") %>%
  mutate(comparison = fct_inorder(comparison))

fig23 <- ggplot(eir_plot_data, aes(x = eir, y = fct_rev(comparison))) +
  geom_vline(xintercept = 1.0, linewidth = 0.5, colour = "grey40",
             linetype = "dashed") +
  geom_errorbar(aes(xmin = eir_lower, xmax = eir_upper),
                width = 0.2, linewidth = 0.6) +
  geom_point(aes(colour = dimension), size = 3) +
  scale_colour_manual(values = c("SEIFA" = "steelblue",
                                  "Remoteness" = "darkorange"),
                      name = NULL) +
  labs(title = "Equity Impact Ratios: Prescription Rate Changes Relative to Reference",
       subtitle = "EIR > 1 = relatively greater increase (or smaller decrease) than reference group | 95% bootstrap CIs",
       x = "Equity Impact Ratio",
       y = NULL) +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(fig_dir, "fig23_equity_impact_ratios.png"),
       fig23, width = 10, height = 6, dpi = 300)
message("  Saved fig23_equity_impact_ratios.png")

# ==============================================================================
# 9. CHOROPLETH MAPS
# ==============================================================================
message("\n--- Choropleth maps ---")

# --- 9a. LGA-level pre/post % change (19-64) ---
lga_effects <- panel_1964 %>%
  mutate(period = if_else(post_reform == 1, "post", "pre")) %>%
  group_by(lga_code, period) %>%
  summarise(mean_scripts = mean(scripts_per_erp, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_scripts,
              names_prefix = "scripts_") %>%
  mutate(pct_change = (scripts_post - scripts_pre) / scripts_pre * 100)

write_csv(lga_effects, file.path(tab_dir, "lga_reform_effects.csv"))
message("  Saved lga_reform_effects.csv (", nrow(lga_effects), " LGAs)")

# Merge with shapefile
lga_map <- lga_sf %>%
  left_join(lga_effects, by = "lga_code") %>%
  left_join(lga_chars %>% select(lga_code, seifa_quintile), by = "lga_code")

# Crop to mainland + Tasmania (exclude external territories)
bbox_aus <- st_bbox(c(xmin = 112, ymin = -45, xmax = 155, ymax = -10),
                    crs = st_crs(lga_map))
lga_map_crop <- st_crop(lga_map, bbox_aus)

# --- Fig 24: Choropleth of reform effects ---
# Cap extreme values for colour scale
lga_map_crop <- lga_map_crop %>%
  mutate(pct_change_cap = pmin(pmax(pct_change, -30), 30))

fig24 <- ggplot(lga_map_crop) +
  geom_sf(aes(fill = pct_change_cap), colour = NA) +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0,
    limits = c(-30, 30),
    name = "% change",
    labels = function(x) paste0(x, "%"),
    oob = squish
  ) +
  labs(title = "LGA-Level Change in Prescription Rates After Copayment Reform (Ages 19-64)",
       subtitle = "% change in mean scripts per ERP: pre-reform (2018-2022) vs post-reform (2023-2025)") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

ggsave(file.path(fig_dir, "fig24_choropleth_reform_effects.png"),
       fig24, width = 10, height = 8, dpi = 300, bg = "white")
message("  Saved fig24_choropleth_reform_effects.png")

# --- Fig 25: SEIFA quintile context map ---
fig25 <- ggplot(lga_map_crop) +
  geom_sf(aes(fill = factor(seifa_quintile)), colour = NA) +
  scale_fill_brewer(
    palette = "RdYlGn", direction = -1,
    name = "SEIFA\nQuintile",
    labels = c("1" = "Q1 (Most\ndisadvantaged)",
               "2" = "Q2", "3" = "Q3", "4" = "Q4",
               "5" = "Q5 (Least\ndisadvantaged)"),
    na.value = "grey90"
  ) +
  labs(title = "SEIFA Index of Relative Socioeconomic Disadvantage by LGA",
       subtitle = "Quintile 1 = most disadvantaged, Quintile 5 = least disadvantaged") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

ggsave(file.path(fig_dir, "fig25_choropleth_seifa.png"),
       fig25, width = 10, height = 8, dpi = 300, bg = "white")
message("  Saved fig25_choropleth_seifa.png")

# ==============================================================================
# 10. SPATIAL AUTOCORRELATION
# ==============================================================================
message("\n--- Spatial autocorrelation ---")

# Need a clean spatial object with no missing pct_change and valid geometry
lga_for_moran <- lga_map_crop %>%
  filter(!is.na(pct_change), !st_is_empty(geometry))

message("  LGAs for Moran's I: ", nrow(lga_for_moran))

# Queen contiguity weights
nb <- poly2nb(lga_for_moran, queen = TRUE)

# Remove isolates (islands etc.) — use zero.policy
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Global Moran's I
moran_result <- moran.test(lga_for_moran$pct_change, lw, zero.policy = TRUE)

morans_i_df <- tibble(
  statistic     = "Global Moran's I",
  estimate      = moran_result$estimate["Moran I statistic"],
  expectation   = moran_result$estimate["Expectation"],
  variance      = moran_result$estimate["Variance"],
  z_score       = (moran_result$estimate["Moran I statistic"] -
                     moran_result$estimate["Expectation"]) /
                    sqrt(moran_result$estimate["Variance"]),
  p_value       = moran_result$p.value,
  n_lgas        = nrow(lga_for_moran),
  n_neighbours  = sum(card(nb) > 0)
)

write_csv(morans_i_df, file.path(tab_dir, "morans_i_results.csv"))
message("  Saved morans_i_results.csv")

message("  Global Moran's I = ", sprintf("%.4f", morans_i_df$estimate),
        " (z = ", sprintf("%.2f", morans_i_df$z_score),
        ", p = ", sprintf("%.4f", morans_i_df$p_value), ")")

# ==============================================================================
# 11. TABLE OUTPUTS SUMMARY
# ==============================================================================
message("\n--- All table outputs ---")
message("  outputs/tables/descriptive_seifa_pre_post.csv")
message("  outputs/tables/descriptive_remoteness_pre_post.csv")
message("  outputs/tables/did_seifa_results.csv")
message("  outputs/tables/did_remoteness_results.csv")
message("  outputs/tables/did_triple_diff_results.csv")
message("  outputs/tables/did_seifa_remoteness_interaction.csv")
message("  outputs/tables/equity_impact_ratios.csv")
message("  outputs/tables/lga_reform_effects.csv")
message("  outputs/tables/morans_i_results.csv")

# ==============================================================================
# 12. SUMMARY
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("PHASE 4: EQUITY ANALYSIS COMPLETE")
message(paste(rep("=", 70), collapse = ""))

message("\nTables saved (9):")
message("  descriptive_seifa_pre_post.csv")
message("  descriptive_remoteness_pre_post.csv")
message("  did_seifa_results.csv")
message("  did_remoteness_results.csv")
message("  did_triple_diff_results.csv")
message("  did_seifa_remoteness_interaction.csv")
message("  equity_impact_ratios.csv")
message("  lga_reform_effects.csv")
message("  morans_i_results.csv")

message("\nFigures saved (8):")
message("  fig18_seifa_time_series.png")
message("  fig19_remoteness_time_series.png")
message("  fig20_did_seifa_coefplot.png")
message("  fig21_did_remoteness_coefplot.png")
message("  fig22_triple_diff_heatmap.png")
message("  fig23_equity_impact_ratios.png")
message("  fig24_choropleth_reform_effects.png")
message("  fig25_choropleth_seifa.png")

message("\nKey findings:")
message("\n  DiD SEIFA interaction terms (post_reform x quintile, vs Q5):")
for (i in seq_len(nrow(seifa_interactions))) {
  row <- seifa_interactions[i, ]
  message("    ", sprintf("%-30s", row$term),
          " est = ", sprintf("%+.6f", row$estimate),
          " (p = ", sprintf("%.4f", row$p.value), ")")
}

message("\n  Triple-diff key coefficient:")
triple_q1 <- did_triple_tidy %>%
  filter(str_detect(term, "post_reform:is_general:seifa_quintile1"))
if (nrow(triple_q1) > 0) {
  message("    post_reform:is_general:seifa_quintile1 = ",
          sprintf("%+.6f", triple_q1$estimate),
          " [", sprintf("%.6f", triple_q1$conf.low), ", ",
          sprintf("%.6f", triple_q1$conf.high), "]",
          " (p = ", sprintf("%.4f", triple_q1$p.value), ")")
}

message("\n  Equity impact ratios:")
for (i in seq_len(nrow(all_eir))) {
  row <- all_eir[i, ]
  message("    ", sprintf("%-45s", row$comparison),
          " EIR = ", sprintf("%.4f", row$eir),
          " [", sprintf("%.4f", row$eir_lower), ", ",
          sprintf("%.4f", row$eir_upper), "]")
}

message("\n  Spatial autocorrelation:")
message("    Global Moran's I = ", sprintf("%.4f", morans_i_df$estimate),
        " (z = ", sprintf("%.2f", morans_i_df$z_score),
        ", p = ", sprintf("%.4f", morans_i_df$p_value), ")")
if (morans_i_df$p_value < 0.05) {
  message("    -> Significant spatial clustering of reform effects")
} else {
  message("    -> No significant spatial clustering detected")
}

message("\nInterpretation:")
message("  EIR > 1: relatively greater increase (or smaller decrease) vs reference")
message("  EIR < 1: relatively greater decrease vs reference")
message("  EIR = 1: no differential change")
message("  Positive DiD coefficient: group gained more (or lost less) than reference")

message("\nNew package dependencies (if not installed):")
message("  install.packages(c('data.table', 'fixest', 'spdep'))")

message("\nNext: Run 07_sensitivity.R (Phase 5)")
