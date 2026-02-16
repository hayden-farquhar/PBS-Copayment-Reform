# ==============================================================================
# 01b_descriptive_analysis.R
# Descriptive statistics and exploratory time series plots (Phase 1)
#
# Inputs:  data/processed/pbs_national_monthly.csv
#          data/processed/pbs_state_monthly.csv
#          data/processed/lga_characteristics.csv
#          data/processed/erp_lga.csv
#
# Outputs: outputs/tables/  — descriptive stats CSVs
#          outputs/figures/  — exploratory time series plots
# ==============================================================================

library(tidyverse)
library(scales)
library(lubridate)

# -- Paths --------------------------------------------------------------------
proj_dir <- here::here()
proc_dir <- file.path(proj_dir, "data", "processed")
fig_dir  <- file.path(proj_dir, "outputs", "figures")
tab_dir  <- file.path(proj_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# -- Theme --------------------------------------------------------------------
theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold"),
                  legend.position = "bottom"))

intervention_date <- as.Date("2023-01-01")
covid_start       <- as.Date("2020-03-01")
covid_end         <- as.Date("2021-12-31")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
message("Loading data...")

national <- read_csv(file.path(proc_dir, "pbs_national_monthly.csv"),
                     show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

state_df <- read_csv(file.path(proc_dir, "pbs_state_monthly.csv"),
                     show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

lga_chars <- read_csv(file.path(proc_dir, "lga_characteristics.csv"),
                      show_col_types = FALSE)

message("  National: ", format(nrow(national), big.mark = ","), " rows, ",
        n_distinct(national$atc_name), " ATC classes")
message("  State:    ", format(nrow(state_df), big.mark = ","), " rows")
message("  LGAs:     ", nrow(lga_chars), " with characteristics")

# ==============================================================================
# 2. DESCRIPTIVE STATISTICS TABLES
# ==============================================================================
message("\nGenerating descriptive statistics...")

# --- 2a. ATC class summary (national, Total age group) ---
atc_summary <- national %>%
  filter(age_group == "Total") %>%
  group_by(atc_name) %>%
  summarise(
    months_observed   = n(),
    mean_scripts_erp  = mean(scripts_per_erp, na.rm = TRUE),
    sd_scripts_erp    = sd(scripts_per_erp, na.rm = TRUE),
    mean_benefits_erp = mean(benefits_per_erp, na.rm = TRUE),
    sd_benefits_erp   = sd(benefits_per_erp, na.rm = TRUE),
    # Pre vs post reform
    mean_scripts_pre  = mean(scripts_per_erp[date < intervention_date], na.rm = TRUE),
    mean_scripts_post = mean(scripts_per_erp[date >= intervention_date], na.rm = TRUE),
    pct_change_scripts = (mean_scripts_post - mean_scripts_pre) / mean_scripts_pre * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_scripts_erp))

write_csv(atc_summary, file.path(tab_dir, "atc_class_summary.csv"))
message("  Saved atc_class_summary.csv (", nrow(atc_summary), " ATC classes)")

# Print top 15 by volume
message("\n  Top 15 ATC classes by mean scripts per ERP:")
atc_summary %>%
  head(15) %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  print(n = 15, width = 120)

# --- 2b. Pre/post reform comparison (national, Total age group) ---
pre_post <- national %>%
  filter(age_group == "Total") %>%
  mutate(period = if_else(date < intervention_date, "Pre-reform", "Post-reform")) %>%
  group_by(period) %>%
  summarise(
    n_months          = n_distinct(date),
    mean_scripts_erp  = mean(scripts_per_erp, na.rm = TRUE),
    sd_scripts_erp    = sd(scripts_per_erp, na.rm = TRUE),
    mean_benefits_erp = mean(benefits_per_erp, na.rm = TRUE),
    sd_benefits_erp   = sd(benefits_per_erp, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(pre_post, file.path(tab_dir, "pre_post_reform_summary.csv"))
message("  Saved pre_post_reform_summary.csv")

# --- 2c. Age group comparison ---
age_summary <- national %>%
  filter(age_group != "Total") %>%
  mutate(period = if_else(date < intervention_date, "Pre-reform", "Post-reform")) %>%
  group_by(age_group, period) %>%
  summarise(
    mean_scripts_erp  = mean(scripts_per_erp, na.rm = TRUE),
    mean_benefits_erp = mean(benefits_per_erp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = period,
              values_from = c(mean_scripts_erp, mean_benefits_erp)) %>%
  mutate(pct_change_scripts = (`mean_scripts_erp_Post-reform` -
                                 `mean_scripts_erp_Pre-reform`) /
           `mean_scripts_erp_Pre-reform` * 100)

write_csv(age_summary, file.path(tab_dir, "age_group_summary.csv"))
message("  Saved age_group_summary.csv")

# --- 2d. State comparison ---
state_summary <- state_df %>%
  filter(age_group == "Total") %>%
  mutate(period = if_else(date < intervention_date, "Pre-reform", "Post-reform")) %>%
  group_by(state, period) %>%
  summarise(
    mean_scripts_erp  = mean(scripts_per_erp, na.rm = TRUE),
    mean_benefits_erp = mean(benefits_per_erp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = period,
              values_from = c(mean_scripts_erp, mean_benefits_erp)) %>%
  mutate(pct_change_scripts = (`mean_scripts_erp_Post-reform` -
                                 `mean_scripts_erp_Pre-reform`) /
           `mean_scripts_erp_Pre-reform` * 100)

write_csv(state_summary, file.path(tab_dir, "state_summary.csv"))
message("  Saved state_summary.csv")

# --- 2e. LGA characteristics distribution ---
lga_distribution <- lga_chars %>%
  count(seifa_quintile, remoteness_area) %>%
  arrange(seifa_quintile, remoteness_area)

write_csv(lga_distribution, file.path(tab_dir, "lga_seifa_remoteness_distribution.csv"))
message("  Saved lga_seifa_remoteness_distribution.csv")

# ==============================================================================
# 3. TIME SERIES PLOTS
# ==============================================================================
message("\nGenerating plots...")

# Helper: add intervention and COVID shading to a ggplot
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

# --- 3a. Overall national prescription rate ---
p1_data <- national %>%
  filter(age_group == "Total") %>%
  group_by(date) %>%
  summarise(total_scripts_erp = sum(scripts_per_erp, na.rm = TRUE),
            .groups = "drop")

p1 <- ggplot(p1_data, aes(date, total_scripts_erp)) +
  geom_line(colour = "steelblue", linewidth = 0.4) +
  labs(title = "Total PBS Prescriptions per Capita (All ATC2 Classes)",
       subtitle = "National monthly rate, May 2002 – Dec 2025",
       x = NULL, y = "Scripts per ERP") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y")

p1 <- add_intervention_lines(p1)
ggsave(file.path(fig_dir, "fig01_national_total_scripts.png"),
       p1, width = 10, height = 5, dpi = 300)
message("  Saved fig01_national_total_scripts.png")

# --- 3b. By age group ---
p2_data <- national %>%
  filter(age_group != "Total") %>%
  group_by(date, age_group) %>%
  summarise(total_scripts_erp = sum(scripts_per_erp, na.rm = TRUE),
            .groups = "drop")

p2 <- ggplot(p2_data, aes(date, total_scripts_erp, colour = age_group)) +
  geom_line(linewidth = 0.4) +
  labs(title = "PBS Prescriptions per Capita by Age Group",
       subtitle = "National monthly rate, May 2002 – Dec 2025",
       x = NULL, y = "Scripts per ERP", colour = "Age group") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_colour_brewer(palette = "Set1")

p2 <- add_intervention_lines(p2)
ggsave(file.path(fig_dir, "fig02_national_by_age.png"),
       p2, width = 10, height = 5, dpi = 300)
message("  Saved fig02_national_by_age.png")

# --- 3c. Top 10 ATC classes over time ---
top10_atc <- atc_summary %>% head(10) %>% pull(atc_name)

p3_data <- national %>%
  filter(age_group == "Total", atc_name %in% top10_atc)

p3 <- ggplot(p3_data, aes(date, scripts_per_erp, colour = atc_name)) +
  geom_line(linewidth = 0.3) +
  labs(title = "PBS Prescriptions per Capita — Top 10 ATC2 Classes",
       subtitle = "National monthly rate, May 2002 – Dec 2025",
       x = NULL, y = "Scripts per ERP", colour = NULL) +
  scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
  guides(colour = guide_legend(ncol = 2))

p3 <- add_intervention_lines(p3)
ggsave(file.path(fig_dir, "fig03_top10_atc_classes.png"),
       p3, width = 12, height = 6, dpi = 300)
message("  Saved fig03_top10_atc_classes.png")

# --- 3d. Faceted top 6 ATC classes (clearer view) ---
top6_atc <- atc_summary %>% head(6) %>% pull(atc_name)

p4_data <- national %>%
  filter(age_group == "Total", atc_name %in% top6_atc) %>%
  mutate(atc_name = factor(atc_name, levels = top6_atc))

p4 <- ggplot(p4_data, aes(date, scripts_per_erp)) +
  geom_line(colour = "steelblue", linewidth = 0.3) +
  facet_wrap(~ atc_name, scales = "free_y", ncol = 2) +
  geom_vline(xintercept = intervention_date,
             linetype = "dashed", colour = "red", linewidth = 0.5) +
  annotate("rect", xmin = covid_start, xmax = covid_end,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  labs(title = "PBS Prescriptions per Capita — Top 6 ATC2 Classes",
       subtitle = "Red dashed line = January 2023 copayment reduction",
       x = NULL, y = "Scripts per ERP") +
  scale_x_date(date_breaks = "4 years", date_labels = "%Y")

ggsave(file.path(fig_dir, "fig04_top6_atc_faceted.png"),
       p4, width = 10, height = 8, dpi = 300)
message("  Saved fig04_top6_atc_faceted.png")

# --- 3e. Pre/post reform bar chart by ATC class (top 15) ---
top15_atc <- atc_summary %>% head(15) %>% pull(atc_name)

p5_data <- national %>%
  filter(age_group == "Total", atc_name %in% top15_atc) %>%
  mutate(period = if_else(date < intervention_date, "Pre-reform", "Post-reform")) %>%
  group_by(atc_name, period) %>%
  summarise(mean_scripts = mean(scripts_per_erp, na.rm = TRUE), .groups = "drop") %>%
  mutate(atc_name = factor(atc_name, levels = rev(top15_atc)))

p5 <- ggplot(p5_data, aes(mean_scripts, atc_name, fill = period)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(title = "Mean PBS Prescription Rate: Pre vs Post Reform",
       subtitle = "Top 15 ATC Level 2 classes by volume",
       x = "Mean scripts per ERP (monthly)", y = NULL, fill = NULL) +
  scale_fill_manual(values = c("Pre-reform" = "grey60", "Post-reform" = "steelblue"))

ggsave(file.path(fig_dir, "fig05_pre_post_bar.png"),
       p5, width = 10, height = 7, dpi = 300)
message("  Saved fig05_pre_post_bar.png")

# --- 3f. Percentage change by ATC class ---
p6_data <- atc_summary %>%
  head(20) %>%
  mutate(atc_name = fct_reorder(atc_name, pct_change_scripts),
         direction = if_else(pct_change_scripts >= 0, "Increase", "Decrease"))

p6 <- ggplot(p6_data, aes(pct_change_scripts, atc_name, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(title = "Percentage Change in Prescription Rate After Copayment Reform",
       subtitle = "Top 20 ATC2 classes, comparing pre vs post January 2023",
       x = "% change in mean scripts per ERP", y = NULL) +
  scale_fill_manual(values = c("Increase" = "steelblue", "Decrease" = "tomato"),
                    guide = "none")

ggsave(file.path(fig_dir, "fig06_pct_change_atc.png"),
       p6, width = 10, height = 7, dpi = 300)
message("  Saved fig06_pct_change_atc.png")

# --- 3g. By state ---
p7_data <- state_df %>%
  filter(age_group == "Total") %>%
  group_by(date, state) %>%
  summarise(total_scripts = sum(scripts_per_erp, na.rm = TRUE), .groups = "drop")

p7 <- ggplot(p7_data, aes(date, total_scripts)) +
  geom_line(colour = "steelblue", linewidth = 0.3) +
  facet_wrap(~ state, scales = "free_y", ncol = 3) +
  geom_vline(xintercept = intervention_date,
             linetype = "dashed", colour = "red", linewidth = 0.5) +
  annotate("rect", xmin = covid_start, xmax = covid_end,
           ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.3) +
  labs(title = "PBS Prescriptions per Capita by State/Territory",
       subtitle = "Red dashed line = January 2023 copayment reduction",
       x = NULL, y = "Scripts per ERP") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y")

ggsave(file.path(fig_dir, "fig07_state_faceted.png"),
       p7, width = 12, height = 8, dpi = 300)
message("  Saved fig07_state_faceted.png")

# --- 3h. Age group comparison: potential general vs concessional proxy ---
p8_data <- national %>%
  filter(age_group %in% c("19-64", "65+")) %>%
  group_by(date, age_group) %>%
  summarise(total_scripts = sum(scripts_per_erp, na.rm = TRUE), .groups = "drop")

p8 <- ggplot(p8_data, aes(date, total_scripts, colour = age_group)) +
  geom_line(linewidth = 0.4) +
  labs(title = "PBS Prescriptions: 19–64 (General Proxy) vs 65+ (Concessional Proxy)",
       subtitle = "Basis for CausalImpact control series",
       x = NULL, y = "Scripts per ERP", colour = "Age group") +
  scale_colour_manual(values = c("19-64" = "steelblue", "65+" = "darkorange"))

p8 <- add_intervention_lines(p8)
ggsave(file.path(fig_dir, "fig08_age_proxy_comparison.png"),
       p8, width = 10, height = 5, dpi = 300)
message("  Saved fig08_age_proxy_comparison.png")

# --- 3i. Seasonality check (monthly pattern) ---
p9_data <- national %>%
  filter(age_group == "Total",
         date >= "2015-01-01", date < "2023-01-01") %>%
  group_by(date) %>%
  summarise(total_scripts = sum(scripts_per_erp, na.rm = TRUE), .groups = "drop") %>%
  mutate(month = month(date, label = TRUE))

p9 <- ggplot(p9_data, aes(month, total_scripts)) +
  geom_boxplot(fill = "steelblue", alpha = 0.5, outlier.size = 1) +
  labs(title = "Monthly Seasonality in PBS Prescription Rates (Pre-Reform)",
       subtitle = "Distribution of total monthly scripts per ERP, 2015–2022",
       x = NULL, y = "Total scripts per ERP")

ggsave(file.path(fig_dir, "fig09_seasonality_boxplot.png"),
       p9, width = 8, height = 5, dpi = 300)
message("  Saved fig09_seasonality_boxplot.png")

# ==============================================================================
# 4. SUMMARY
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("PHASE 1 DESCRIPTIVE ANALYSIS COMPLETE")
message(paste(rep("=", 60), collapse = ""))

message("\nTables in: ", tab_dir)
for (f in list.files(tab_dir, pattern = "\\.csv$")) {
  message("  ", f)
}

message("\nFigures in: ", fig_dir)
for (f in list.files(fig_dir, pattern = "\\.png$")) {
  message("  ", f)
}

message("\nKey observations to check:")
message("  1. Is there a visible level shift at January 2023 in fig01?")
message("  2. Do the age groups diverge after reform in fig08?")
message("  3. Which ATC classes show the largest changes in fig06?")
message("  4. Is there strong seasonality that needs Fourier terms (fig09)?")
message("  5. How much COVID disruption is visible in the pre-period?")
