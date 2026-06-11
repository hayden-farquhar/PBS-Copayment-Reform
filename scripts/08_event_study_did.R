# ==============================================================================
# 08_event_study_did.R
# Event-study / triple-difference identification of the 2023 PBS copayment cut
#
# Four threats to identification motivate this analysis. The script builds the
# estimators that address each directly:
#
#   (1) Sensitivity to the choice of pre-intervention period / parallel-trends
#   (2) Need for an event-study (dynamic DiD) with a testable pre-trend
#   (3) COVID/Omicron confounding: a binary pandemic dummy ending in 2021 leaves
#       the Omicron-contaminated immediate pre-intervention window (Dec 2021 -
#       2022) uncontrolled and imposes a uniform national shock with no
#       geographic heterogeneity
#   (4) Contemporaneous GP-access confounders, esp. the Nov-2023 tripling of
#       bulk-billing incentives that targeted CONCESSION-CARD HOLDERS
#       (= the 65+ concessional control group, plus working-age Health Care
#        Card holders inside the 19-64 "general" arm)
#
# Design logic:
#   - Identification is the general (19-64, treated) vs concessional (65+,
#     control) difference-in-differences around the 1 Jan 2023 reform.
#   - Event study via fixest::i(event_time, is_general, ref = -1) with
#       unit FE  = lga_code^age_group   (absorbs all time-invariant level gaps)
#       time FE  = date                 (absorbs ALL national time shocks,
#                                        incl. Omicron, NON-PARAMETRICALLY)
#     -> the binary covid_period dummy is RETIRED in favour of month FE.
#   - Threat 3 (geographic heterogeneity): add state^date (and remoteness^date)
#     FE so pandemic disruption is allowed to differ by jurisdiction/remoteness.
#   - Threat 4: headline ATT on the CLEAN pre-Nov-2023 post-window, plus an
#     explicit concessional-specific bulk-billing shock term in the full window.
#   - Threat 1: re-estimate event study + static ATT across four pre-period
#     anchors (2015/2018/2020/2021) to show stability (or characterise drift).
#
# IMPORTANT caveat carried through (Hayden, 2026-06-09): the 65+/19-64 age proxy
# for concessional/general is LEAKY. Concession-card holders include Health Care
# Card holders (income-tested, any age), disability/carer/parenting pensioners,
# DVA, and Commonwealth Seniors Health Card holders - so a non-trivial share of
# the 19-64 "general" arm is actually concessional, and the Nov-2023 bulk-billing
# change touched BOTH arms (more in the 65+ arm). This attenuates the DiD and
# blurs the cleanliness of the control. Treated as a first-class limitation, not
# a footnote.
#
# Inputs:  data/processed/pbs_national_monthly.csv   (national, by ATC x age)
#          data/processed/pbs_lga_monthly.csv         (2.9 GB LGA panel)
#          data/processed/lga_characteristics.csv
#
# Outputs: outputs/tables/event_study_national.csv
#          outputs/tables/event_study_lga.csv
#          outputs/tables/event_study_pretrend_tests.csv
#          outputs/tables/did_att_specifications.csv
#          outputs/tables/did_preperiod_anchor_grid.csv
#          outputs/tables/event_study_lga_deseasonalised.csv
#          outputs/tables/did_att_deseasonalised.csv
#          outputs/tables/did_preperiod_anchor_grid_deseasonalised.csv
#          outputs/tables/first_stage_benefit_per_script.csv
#          outputs/tables/dose_response_classes.csv
#          outputs/tables/ddd_results.csv
#          outputs/tables/honestdid_sensitivity.csv
#          outputs/tables/placebo_in_time.csv
#          outputs/figures/fig29_event_study_national.png
#          outputs/figures/fig30_event_study_lga.png
#          outputs/figures/fig31_att_specification_curve.png
#          outputs/figures/fig32_preperiod_anchor_grid.png
#          outputs/figures/fig33_event_study_lga_deseasonalised.png
#          outputs/figures/fig34_preperiod_anchor_grid_deseasonalised.png
#          outputs/figures/fig35_first_stage_benefit_per_script.png
#          outputs/figures/fig36_dose_response_classes.png
#          outputs/figures/fig37_honestdid_sensitivity.png
#          outputs/figures/fig38_placebo_in_time.png
#
# Run AFTER 05_equity_analysis.R. Working directory = project root.
# ==============================================================================

# ==============================================================================
# 0. SETUP
# ==============================================================================
library(tidyverse)
library(lubridate)
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

# --- Key dates ---------------------------------------------------------------
intervention_date <- as.Date("2023-01-01")   # PBS general copayment cut
bulk_billing_date <- as.Date("2023-11-01")   # Tripling of bulk-billing incentive
omicron_start     <- as.Date("2021-12-01")   # Omicron wave hits Australia
omicron_end       <- as.Date("2022-12-31")   # (illustrative band, for plotting)

# Helper: signed event time in MONTHS relative to the reform (Jan 2023 = 0,
# Dec 2022 = -1, Feb 2023 = +1). Robust to day-of-month via a month index.
to_event_time <- function(d) {
  (year(d) - year(intervention_date)) * 12 + (month(d) - month(intervention_date))
}

# Reference event-month for the event study (the eve of the reform).
EVENT_REF <- -1L

# Clean post-window upper bound (event months 0..CLEAN_POST_MAX inclusive) that
# precede the Nov-2023 bulk-billing shock. Nov 2023 is event month +10.
CLEAN_POST_MAX <- 9L   # Jan 2023 .. Oct 2023

# Primary event-study display window (months pre/post). Pre-period anchors are
# varied separately in Section 5.
ES_LEAD_MAX <- 48L     # show up to 48 months of leads
ES_LAG_MAX  <- 35L     # post-period extends to Dec 2025

# ==============================================================================
# 1. NATIONAL EVENT STUDY (transparent visual; differenced contrast)
# ==============================================================================
# At the national level there are only two series (general vs concessional), so
# instead of a 2-cluster FE model we form the DiD contrast directly:
#     d_t = scripts_general_t - scripts_concessional_t
# and regress d_t on event-time dummies (ref = -1). This is the visual analogue
# of the LGA event study and matches the paper's existing national framing.
message("\n=== Section 1: National event study (differenced contrast) ===")

nat <- fread(file.path(proc_dir, "pbs_national_monthly.csv"),
             select = c("atc_name", "age_group", "date",
                        "scripts_per_erp"))
nat <- nat[age_group %in% c("19-64", "65+")]

# Overall (all-ATC) rate per age group per month
nat_overall <- nat[, .(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE)),
                   by = .(date, age_group)]
nat_overall[, date := as.Date(date)]

nat_wide <- dcast(nat_overall, date ~ age_group, value.var = "scripts_per_erp")
setnames(nat_wide, c("19-64", "65+"), c("general", "concessional"))
nat_wide[, diff_gc := general - concessional]
nat_wide[, event_time := to_event_time(date)]

# Restrict to the display window
nat_es <- nat_wide[event_time >= -ES_LEAD_MAX & event_time <= ES_LAG_MAX]
nat_es[, et_f := relevel(factor(event_time), ref = as.character(EVENT_REF))]

# Regress the differenced contrast on event-time dummies; HAC (Newey-West) SEs
# for autocorrelation, consistent with the paper's ITS approach.
fit_nat <- feols(diff_gc ~ i(event_time, ref = EVENT_REF), data = nat_es,
                 vcov = NW() ~ event_time)

es_nat <- broom::tidy(fit_nat, conf.int = TRUE) %>%
  filter(str_detect(term, "event_time::")) %>%
  mutate(event_time = as.integer(str_extract(term, "-?\\d+"))) %>%
  arrange(event_time)

# Add the reference row (0 by construction) so the plot has a point at -1
es_nat <- bind_rows(
  es_nat,
  tibble(term = "ref", estimate = 0, std.error = NA,
         statistic = NA, p.value = NA, conf.low = 0, conf.high = 0,
         event_time = EVENT_REF)
) %>% arrange(event_time)

write_csv(es_nat, file.path(tab_dir, "event_study_national.csv"))
message("  Saved event_study_national.csv (", nrow(es_nat), " event months)")

fig29 <- ggplot(es_nat, aes(x = event_time, y = estimate)) +
  annotate("rect", xmin = to_event_time(omicron_start),
           xmax = to_event_time(omicron_end), ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "red",
             linewidth = 0.5) +
  geom_vline(xintercept = to_event_time(bulk_billing_date) - 0.5,
             linetype = "dotted", colour = "purple", linewidth = 0.6) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "steelblue",
              alpha = 0.2) +
  geom_line(colour = "steelblue", linewidth = 0.5) +
  geom_point(colour = "steelblue", size = 1) +
  annotate("text", x = -0.5, y = Inf, label = "Reform\n(Jan 2023)",
           hjust = 1.05, vjust = 1.3, colour = "red", size = 3,
           fontface = "italic") +
  annotate("text", x = to_event_time(bulk_billing_date), y = Inf,
           label = "Bulk-billing\ntripling (Nov 2023)", hjust = -0.03,
           vjust = 1.3, colour = "purple", size = 3, fontface = "italic") +
  annotate("text", x = to_event_time(omicron_start) + 6, y = -Inf,
           label = "Omicron", vjust = -0.8, colour = "grey45", size = 3,
           fontface = "italic") +
  labs(title = "National Event Study: General-minus-Concessional Prescription Gap",
       subtitle = paste0("Differenced contrast (19-64 minus 65+), event time relative to reform; ref = -1 (Dec 2022). ",
                         "Newey-West 95% CI."),
       x = "Months relative to copayment reduction (0 = Jan 2023)",
       y = "Differential in scripts per ERP (general - concessional)")

ggsave(file.path(fig_dir, "fig29_event_study_national.png"),
       fig29, width = 12, height = 6, dpi = 300)
message("  Saved fig29_event_study_national.png")

# ==============================================================================
# 2. LGA EVENT STUDY (primary inference; high-powered, many clusters)
# ==============================================================================
message("\n=== Section 2: LGA event study (primary) ===")
message("  Loading LGA panel (2.9 GB) - subset of columns only...")

lga_raw <- fread(
  file.path(proc_dir, "pbs_lga_monthly.csv"),
  select = c("state", "lga_code", "date", "age_group", "scripts_per_erp",
             "remoteness_area")
)
lga_raw <- lga_raw[age_group %in% c("19-64", "65+")]

# Aggregate across ATC classes -> LGA x month x age_group overall rate
lga_agg <- lga_raw[
  , .(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE)),
  by = .(state, lga_code, date, age_group, remoteness_area)
]
rm(lga_raw); gc()

lga_agg[, date := as.Date(date)]
lga_agg[, lga_code := as.integer(lga_code)]
lga_agg[, is_general := as.integer(age_group == "19-64")]
lga_agg[, event_time := to_event_time(date)]

# Build the event-study estimation frame (full window kept; Section 5 varies it)
# Primary spec uses leads back to ES_LEAD_MAX so the figure is legible.
es_lga <- lga_agg[event_time >= -ES_LEAD_MAX & event_time <= ES_LAG_MAX]
message("  Event-study frame: ", format(nrow(es_lga), big.mark = ","), " rows, ",
        uniqueN(es_lga$lga_code), " LGAs")

# --- Primary event study -----------------------------------------------------
# i(event_time, is_general, ref = -1): general-vs-concessional differential at
# each event month relative to the eve of the reform.
# Unit FE lga_code^age_group absorbs the is_general main effect and all
# time-invariant level gaps; date FE absorbs ALL national time shocks (Omicron
# included). SEs clustered by LGA.
fit_lga <- feols(
  scripts_per_erp ~ i(event_time, is_general, ref = EVENT_REF)
  | lga_code^age_group + date,
  data = es_lga, cluster = ~lga_code
)

es_lga_tidy <- broom::tidy(fit_lga, conf.int = TRUE) %>%
  filter(str_detect(term, "event_time::")) %>%
  mutate(event_time = as.integer(str_extract(term, "-?\\d+"))) %>%
  arrange(event_time)

es_lga_tidy <- bind_rows(
  es_lga_tidy,
  tibble(term = "ref", estimate = 0, std.error = NA, statistic = NA,
         p.value = NA, conf.low = 0, conf.high = 0, event_time = EVENT_REF)
) %>% arrange(event_time)

write_csv(es_lga_tidy, file.path(tab_dir, "event_study_lga.csv"))
message("  Saved event_study_lga.csv")

# --- Joint pre-trend test ----------------------------------------------------
# H0: all pre-period interaction coefficients (event_time < -1) = 0.
pre_terms <- es_lga_tidy %>%
  filter(event_time < EVENT_REF, term != "ref") %>%
  pull(term)
pretrend_wald <- wald(fit_lga, keep = pre_terms)

# Post-period (treatment-active) joint test for completeness
post_terms <- es_lga_tidy %>% filter(event_time >= 0) %>% pull(term)
post_wald <- wald(fit_lga, keep = post_terms)

pretrend_tbl <- tibble(
  test = c("Pre-trend (event_time < -1 == 0)", "Post (event_time >= 0 == 0)"),
  stat = c(pretrend_wald$stat, post_wald$stat),
  df1  = c(pretrend_wald$df1, post_wald$df1),
  df2  = c(pretrend_wald$df2, post_wald$df2),
  p_value = c(pretrend_wald$p, post_wald$p)
)
write_csv(pretrend_tbl, file.path(tab_dir, "event_study_pretrend_tests.csv"))
message("  Pre-trend joint Wald p = ", sprintf("%.4f", pretrend_wald$p),
        " (p>0.05 supports parallel pre-trends)")
message("  Saved event_study_pretrend_tests.csv")

# --- Fig 30: LGA event study -------------------------------------------------
fig30 <- ggplot(es_lga_tidy, aes(x = event_time, y = estimate)) +
  annotate("rect", xmin = to_event_time(omicron_start),
           xmax = to_event_time(omicron_end), ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "red",
             linewidth = 0.5) +
  geom_vline(xintercept = to_event_time(bulk_billing_date) - 0.5,
             linetype = "dotted", colour = "purple", linewidth = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0,
                colour = "steelblue", alpha = 0.5) +
  geom_point(aes(colour = event_time >= 0), size = 1.4) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue"),
                      labels = c("FALSE" = "Pre-reform (placebo leads)",
                                 "TRUE" = "Post-reform (dynamic effects)"),
                      name = NULL) +
  annotate("text", x = -0.5, y = Inf, label = "Reform", hjust = 1.1,
           vjust = 1.4, colour = "red", size = 3, fontface = "italic") +
  annotate("text", x = to_event_time(bulk_billing_date), y = Inf,
           label = "Bulk-billing\ntripling", hjust = -0.05, vjust = 1.4,
           colour = "purple", size = 3, fontface = "italic") +
  labs(title = "LGA Event Study: Dynamic General-vs-Concessional DiD",
       subtitle = paste0("i(event_time, is_general, ref=-1) | LGA x age-group FE + month FE; ",
                         "clustered (LGA) 95% CI. Pre-trend Wald p = ",
                         sprintf("%.3f", pretrend_wald$p), "."),
       x = "Months relative to copayment reduction (0 = Jan 2023)",
       y = "Differential effect on scripts per ERP (general - concessional)")

ggsave(file.path(fig_dir, "fig30_event_study_lga.png"),
       fig30, width = 12, height = 6, dpi = 300)
message("  Saved fig30_event_study_lga.png")

# ==============================================================================
# 3. STATIC ATT SPECIFICATIONS
#    (threat 3: geography/Omicron; threat 4: bulk-billing confounder)
# ==============================================================================
message("\n=== Section 3: Static ATT under alternative specifications ===")

# Estimation frame: use a symmetric, defensible window (2019-2025 -> leads to
# -48). 'post' = reform active; 'clean_post' restricts post to pre-Nov-2023.
att_frame <- lga_agg[event_time >= -ES_LEAD_MAX & event_time <= ES_LAG_MAX]
att_frame[, post := as.integer(event_time >= 0)]
att_frame[, post_general := post * is_general]
# Clean (pre-bulk-billing) post indicator and its DiD interaction
att_frame[, clean_post := as.integer(event_time >= 0 & event_time <= CLEAN_POST_MAX)]
att_frame[, clean_post_general := clean_post * is_general]
# Concessional-specific bulk-billing shock: active from Nov 2023, for the
# control arm. Modelled as bb_post (touches 65+ via the base term) plus its
# interaction with is_general (differential for the 19-64 arm). See caveat:
# working-age HCC holders mean is_general is NOT bulk-billing-free.
att_frame[, bb_post := as.integer(date >= bulk_billing_date)]
att_frame[, bb_post_general := bb_post * is_general]

# Spec A: baseline TWFE DiD (month FE absorbs national COVID/Omicron)
specA <- feols(scripts_per_erp ~ post_general | lga_code^age_group + date,
               data = att_frame, cluster = ~lga_code)

# Spec B: + geographic heterogeneity in time shocks (state-by-month FE)
specB <- feols(scripts_per_erp ~ post_general | lga_code^age_group + date + state^date,
               data = att_frame, cluster = ~lga_code)

# Spec C: + remoteness-by-month FE (pandemic disruption varied by remoteness)
specC <- feols(scripts_per_erp ~ post_general
               | lga_code^age_group + date + remoteness_area^date,
               data = att_frame, cluster = ~lga_code)

# Spec D: full window but explicitly controlling the concessional bulk-billing
# shock (bb_post on the control arm + its general interaction)
specD <- feols(scripts_per_erp ~ post_general + bb_post + bb_post_general
               | lga_code^age_group + date,
               data = att_frame, cluster = ~lga_code)

# Spec E: CLEAN post-window (Jan-Oct 2023, pre-bulk-billing) - the headline
# uncontaminated ATT
specE <- feols(scripts_per_erp ~ clean_post_general | lga_code^age_group + date,
               data = att_frame[event_time <= CLEAN_POST_MAX | event_time < 0],
               cluster = ~lga_code)

extract_att <- function(fit, label, coef_name) {
  ti <- broom::tidy(fit, conf.int = TRUE)
  row <- ti[ti$term == coef_name, ]
  tibble(spec = label, term = coef_name, estimate = row$estimate,
         conf.low = row$conf.low, conf.high = row$conf.high,
         p.value = row$p.value, n = nobs(fit))
}

att_tbl <- bind_rows(
  extract_att(specA, "A: TWFE (month FE absorbs COVID)", "post_general"),
  extract_att(specB, "B: + state x month FE (geo heterogeneity)", "post_general"),
  extract_att(specC, "C: + remoteness x month FE", "post_general"),
  extract_att(specD, "D: full window, bulk-billing controlled", "post_general"),
  extract_att(specE, "E: clean pre-Nov-2023 post-window", "clean_post_general")
)
write_csv(att_tbl, file.path(tab_dir, "did_att_specifications.csv"))
message("  Saved did_att_specifications.csv")
print(att_tbl)

# --- Fig 31: ATT specification curve ----------------------------------------
att_plot <- att_tbl %>% mutate(spec = fct_inorder(spec))
fig31 <- ggplot(att_plot, aes(x = estimate, y = fct_rev(spec))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2,
                 colour = "steelblue") +
  geom_point(aes(colour = p.value < 0.05), size = 3) +
  scale_colour_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey60"),
                      labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                      name = NULL) +
  labs(title = "ATT Robustness: General-vs-Concessional DiD Across Specifications",
       subtitle = "Each row adds a more stringent control for time/geography shocks; point = ATT on scripts per ERP, 95% CI.",
       x = "Average treatment effect on the treated (scripts per ERP)", y = NULL)
ggsave(file.path(fig_dir, "fig31_att_specification_curve.png"),
       fig31, width = 10, height = 5, dpi = 300)
message("  Saved fig31_att_specification_curve.png")

# ==============================================================================
# 4. PRE-PERIOD ANCHOR GRID (threat 1: baseline-period sensitivity)
# ==============================================================================
message("\n=== Section 4: Pre-period anchor sensitivity grid ===")

# Re-estimate the CLEAN-window ATT (Spec E design) and the post-period event-
# study average under four pre-period start anchors. Stable estimates across
# anchors directly rebut "results depend on the baseline period selected".
anchors <- as.Date(c("2015-01-01", "2018-01-01", "2020-01-01", "2021-01-01"))

anchor_rows <- map_dfr(anchors, function(a0) {
  d <- lga_agg[date >= a0 & event_time <= ES_LAG_MAX]
  d[, clean_post := as.integer(event_time >= 0 & event_time <= CLEAN_POST_MAX)]
  d[, clean_post_general := clean_post * is_general]
  d[, post_general := as.integer(event_time >= 0) * is_general]
  d[, bb_post := as.integer(date >= bulk_billing_date)]
  d[, bb_post_general := bb_post * is_general]

  # Clean-window ATT
  f_clean <- feols(scripts_per_erp ~ clean_post_general | lga_code^age_group + date,
                   data = d[event_time <= CLEAN_POST_MAX | event_time < 0],
                   cluster = ~lga_code)
  rc <- broom::tidy(f_clean, conf.int = TRUE)
  rc <- rc[rc$term == "clean_post_general", ]

  # Full-window ATT with bulk-billing controlled
  f_full <- feols(scripts_per_erp ~ post_general + bb_post + bb_post_general
                  | lga_code^age_group + date, data = d, cluster = ~lga_code)
  rf <- broom::tidy(f_full, conf.int = TRUE)
  rf <- rf[rf$term == "post_general", ]

  bind_rows(
    tibble(anchor = format(a0, "%Y"), spec = "Clean pre-Nov-2023 ATT",
           estimate = rc$estimate, conf.low = rc$conf.low,
           conf.high = rc$conf.high, p.value = rc$p.value),
    tibble(anchor = format(a0, "%Y"), spec = "Full-window ATT (BB controlled)",
           estimate = rf$estimate, conf.low = rf$conf.low,
           conf.high = rf$conf.high, p.value = rf$p.value)
  )
})

write_csv(anchor_rows, file.path(tab_dir, "did_preperiod_anchor_grid.csv"))
message("  Saved did_preperiod_anchor_grid.csv")
print(anchor_rows)

fig32 <- ggplot(anchor_rows, aes(x = anchor, y = estimate, colour = spec,
                                 group = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  scale_colour_manual(values = c("Clean pre-Nov-2023 ATT" = "steelblue",
                                 "Full-window ATT (BB controlled)" = "darkorange"),
                      name = NULL) +
  labs(title = "Pre-Period Anchor Sensitivity of the DiD ATT",
       subtitle = "If pre-period choice drives the estimate, it swings across anchors; if robust, estimates cluster.",
       x = "Pre-period start (anchor year)",
       y = "ATT on scripts per ERP (95% CI)")
ggsave(file.path(fig_dir, "fig32_preperiod_anchor_grid.png"),
       fig32, width = 10, height = 5.5, dpi = 300)
message("  Saved fig32_preperiod_anchor_grid.png")

# ==============================================================================
# 6. DE-SEASONALISED SPECIFICATIONS
# ==============================================================================
# The saturated event study (Section 2) revealed a strong ~12-month, GENERAL-
# SPECIFIC seasonal cycle: the general-minus-concessional gap peaks every January
# and troughs every December. This is the PBS safety-net / copayment-year reset,
# which lands on 1 January and lifts the general arm specifically - exactly the
# reform date. With ref = -1 (Dec 2022, a seasonal trough) the t=0 (Jan 2023)
# contrast is mechanically inflated, and the pre-period anchor sensitivity in
# Section 4 is a direct symptom (moving the baseline re-phases the season).
#
# Fix: remove each (lga x age_group) series' own recurring month-of-year level,
# estimated from PRE-period data ONLY (event_time < 0), so the reform effect (a
# post-period deviation from the seasonal norm) is preserved while the recurring
# seasonal swing is netted out. y_adj = scripts_per_erp - pre-period seasonal
# mean for that LGA x arm x calendar-month. Re-run the event study, the static
# ATT, and the anchor grid on y_adj. Interpretation:
#   - flat de-seasonalised LEADS (near 0)            => parallel trends holds
#   - a sustained post-period STEP in the lags       => genuine reform effect
#   - leads still swinging / no clean post step       => effect not separable
#     from general-specific seasonality (honest null / design limitation)
# NOTE: seasonal FE cannot be ADDED to the event study - i(event_time,is_general)
# already saturates time-by-group, so month-of-year FE would be collinear and
# dropped. De-seasonalising the OUTCOME is the only non-collinear route.
message("\n=== Section 6: De-seasonalised specifications ===")

lga_agg[, month_of_year := month(date)]

# Subtract pre-period (lga x age_group x calendar-month) means. Pre-period is
# defined WITHIN whatever window dt spans (event_time < 0), so the same helper
# serves both the main frame and each anchor sub-window.
deseasonalise <- function(dt) {
  seas <- dt[event_time < 0,
             .(seas_mean = mean(scripts_per_erp, na.rm = TRUE)),
             by = .(lga_code, age_group, month_of_year)]
  out <- merge(dt, seas, by = c("lga_code", "age_group", "month_of_year"),
               all.x = TRUE, sort = FALSE)
  out[, y_adj := scripts_per_erp - seas_mean]
  out[]
}

ds_frame <- lga_agg[event_time >= -ES_LEAD_MAX & event_time <= ES_LAG_MAX]
ds_frame <- deseasonalise(ds_frame)
ds_frame <- ds_frame[!is.na(y_adj)]
ds_frame[, post := as.integer(event_time >= 0)]
ds_frame[, post_general := post * is_general]
ds_frame[, clean_post := as.integer(event_time >= 0 & event_time <= CLEAN_POST_MAX)]
ds_frame[, clean_post_general := clean_post * is_general]
ds_frame[, bb_post := as.integer(date >= bulk_billing_date)]
ds_frame[, bb_post_general := bb_post * is_general]

# --- 6a. De-seasonalised event study ----------------------------------------
fit_lga_ds <- feols(
  y_adj ~ i(event_time, is_general, ref = EVENT_REF)
  | lga_code^age_group + date,
  data = ds_frame, cluster = ~lga_code
)

es_ds_tidy <- broom::tidy(fit_lga_ds, conf.int = TRUE) %>%
  filter(str_detect(term, "event_time::")) %>%
  mutate(event_time = as.integer(str_extract(term, "-?\\d+"))) %>%
  arrange(event_time)
es_ds_tidy <- bind_rows(
  es_ds_tidy,
  tibble(term = "ref", estimate = 0, std.error = NA, statistic = NA,
         p.value = NA, conf.low = 0, conf.high = 0, event_time = EVENT_REF)
) %>% arrange(event_time)
write_csv(es_ds_tidy, file.path(tab_dir, "event_study_lga_deseasonalised.csv"))
message("  Saved event_study_lga_deseasonalised.csv")

pre_terms_ds  <- es_ds_tidy %>% filter(event_time < EVENT_REF, term != "ref") %>% pull(term)
post_terms_ds <- es_ds_tidy %>% filter(event_time >= 0) %>% pull(term)
pretrend_wald_ds <- wald(fit_lga_ds, keep = pre_terms_ds)
post_wald_ds     <- wald(fit_lga_ds, keep = post_terms_ds)
message("  De-seasonalised pre-trend joint Wald p = ",
        sprintf("%.4f", pretrend_wald_ds$p),
        "  (read MAGNITUDES too: huge N makes any wobble significant)")

# Persist the de-seasonalised joint pre-trend / post Wald F so the magnitude
# (not just the p-value embedded in fig33) is available for eTable S9.
pretrend_tbl_ds <- tibble(
  test = c("De-seasonalised pre-trend (event_time < -1 == 0)",
           "De-seasonalised post (event_time >= 0 == 0)"),
  stat = c(pretrend_wald_ds$stat, post_wald_ds$stat),
  df1  = c(pretrend_wald_ds$df1, post_wald_ds$df1),
  df2  = c(pretrend_wald_ds$df2, post_wald_ds$df2),
  p_value = c(pretrend_wald_ds$p, post_wald_ds$p)
)
write_csv(pretrend_tbl_ds,
          file.path(tab_dir, "event_study_pretrend_tests_deseasonalised.csv"))
message("  Saved event_study_pretrend_tests_deseasonalised.csv (de-seas Wald F)")

# Mean post-period dynamic coefficient (rough effect size on de-seasonalised y)
mean_post_ds <- mean(es_ds_tidy$estimate[es_ds_tidy$event_time >= 0], na.rm = TRUE)
mean_pre_ds  <- mean(abs(es_ds_tidy$estimate[es_ds_tidy$event_time < EVENT_REF &
                                               es_ds_tidy$term != "ref"]), na.rm = TRUE)
message("  Mean post-period coef (de-seas) = ", sprintf("%.4f", mean_post_ds),
        " vs mean |pre-period coef| = ", sprintf("%.4f", mean_pre_ds))

fig33 <- ggplot(es_ds_tidy, aes(x = event_time, y = estimate)) +
  annotate("rect", xmin = to_event_time(omicron_start),
           xmax = to_event_time(omicron_end), ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.35) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "red",
             linewidth = 0.5) +
  geom_vline(xintercept = to_event_time(bulk_billing_date) - 0.5,
             linetype = "dotted", colour = "purple", linewidth = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0,
                colour = "steelblue", alpha = 0.5) +
  geom_point(aes(colour = event_time >= 0), size = 1.4) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "steelblue"),
                      labels = c("FALSE" = "Pre-reform (placebo leads)",
                                 "TRUE" = "Post-reform (dynamic effects)"),
                      name = NULL) +
  annotate("text", x = -0.5, y = Inf, label = "Reform", hjust = 1.1,
           vjust = 1.4, colour = "red", size = 3, fontface = "italic") +
  labs(title = "De-seasonalised LGA Event Study: General-vs-Concessional DiD",
       subtitle = paste0("Outcome = scripts per ERP minus pre-period month-of-year mean (LGA x arm). ",
                         "Pre-trend Wald p = ", sprintf("%.3f", pretrend_wald_ds$p), "."),
       x = "Months relative to copayment reduction (0 = Jan 2023)",
       y = "De-seasonalised differential effect (scripts per ERP)")
ggsave(file.path(fig_dir, "fig33_event_study_lga_deseasonalised.png"),
       fig33, width = 12, height = 6, dpi = 300)
message("  Saved fig33_event_study_lga_deseasonalised.png")

# --- 6b. De-seasonalised static ATT -----------------------------------------
specA_ds <- feols(y_adj ~ post_general | lga_code^age_group + date,
                  data = ds_frame, cluster = ~lga_code)
specD_ds <- feols(y_adj ~ post_general + bb_post + bb_post_general
                  | lga_code^age_group + date,
                  data = ds_frame, cluster = ~lga_code)
specE_ds <- feols(y_adj ~ clean_post_general | lga_code^age_group + date,
                  data = ds_frame[event_time <= CLEAN_POST_MAX | event_time < 0],
                  cluster = ~lga_code)
att_ds_tbl <- bind_rows(
  extract_att(specA_ds, "A-ds: TWFE de-seasonalised", "post_general"),
  extract_att(specD_ds, "D-ds: full window, BB controlled", "post_general"),
  extract_att(specE_ds, "E-ds: clean pre-Nov-2023 window", "clean_post_general")
)
write_csv(att_ds_tbl, file.path(tab_dir, "did_att_deseasonalised.csv"))
message("  Saved did_att_deseasonalised.csv")
print(att_ds_tbl)

# --- 6c. De-seasonalised pre-period anchor grid -----------------------------
anchor_rows_ds <- map_dfr(anchors, function(a0) {
  d <- lga_agg[date >= a0 & event_time <= ES_LAG_MAX]
  d <- deseasonalise(d)
  d <- d[!is.na(y_adj)]
  d[, clean_post := as.integer(event_time >= 0 & event_time <= CLEAN_POST_MAX)]
  d[, clean_post_general := clean_post * is_general]
  d[, post_general := as.integer(event_time >= 0) * is_general]
  d[, bb_post := as.integer(date >= bulk_billing_date)]
  d[, bb_post_general := bb_post * is_general]

  f_clean <- feols(y_adj ~ clean_post_general | lga_code^age_group + date,
                   data = d[event_time <= CLEAN_POST_MAX | event_time < 0],
                   cluster = ~lga_code)
  rc <- broom::tidy(f_clean, conf.int = TRUE)
  rc <- rc[rc$term == "clean_post_general", ]

  f_full <- feols(y_adj ~ post_general + bb_post + bb_post_general
                  | lga_code^age_group + date, data = d, cluster = ~lga_code)
  rf <- broom::tidy(f_full, conf.int = TRUE)
  rf <- rf[rf$term == "post_general", ]

  bind_rows(
    tibble(anchor = format(a0, "%Y"), spec = "Clean pre-Nov-2023 ATT (de-seas)",
           estimate = rc$estimate, conf.low = rc$conf.low,
           conf.high = rc$conf.high, p.value = rc$p.value),
    tibble(anchor = format(a0, "%Y"), spec = "Full-window ATT (de-seas, BB ctrl)",
           estimate = rf$estimate, conf.low = rf$conf.low,
           conf.high = rf$conf.high, p.value = rf$p.value)
  )
})
write_csv(anchor_rows_ds,
          file.path(tab_dir, "did_preperiod_anchor_grid_deseasonalised.csv"))
message("  Saved did_preperiod_anchor_grid_deseasonalised.csv")
print(anchor_rows_ds)

fig34 <- ggplot(anchor_rows_ds, aes(x = anchor, y = estimate, colour = spec,
                                    group = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1,
                position = position_dodge(width = 0.3)) +
  geom_point(size = 3, position = position_dodge(width = 0.3)) +
  scale_colour_manual(values = c("Clean pre-Nov-2023 ATT (de-seas)" = "steelblue",
                                 "Full-window ATT (de-seas, BB ctrl)" = "darkorange"),
                      name = NULL) +
  labs(title = "De-seasonalised Pre-Period Anchor Sensitivity of the DiD ATT",
       subtitle = "De-seasonalised, the residual stays baseline-dependent: the ATT slides and flips sign across anchors.",
       x = "Pre-period start (anchor year)",
       y = "De-seasonalised ATT on scripts per ERP (95% CI)")
ggsave(file.path(fig_dir, "fig34_preperiod_anchor_grid_deseasonalised.png"),
       fig34, width = 10, height = 5.5, dpi = 300)
message("  Saved fig34_preperiod_anchor_grid_deseasonalised.png")

# ==============================================================================
# 7. FIRST-STAGE / DOSE-DELIVERED CHECK (national)
# ==============================================================================
# Was the reform's price reduction actually delivered to general patients? The
# government benefit per script = benefits_per_erp / scripts_per_erp. When the
# general copayment fell $42.50 -> $30, the government pays up to $12.50 more per
# (above-$30) general script, so general benefit/script should STEP UP at Jan
# 2023; concessional (copay unchanged) should not. A strong first stage + flat
# scripts (Sections 3/6) = the gold-standard inelastic-demand null: the policy
# was delivered but utilisation did not respond.
message("\n=== Section 7: First-stage / dose-delivered check ===")

nat_b <- fread(file.path(proc_dir, "pbs_national_monthly.csv"),
               select = c("atc_name", "age_group", "date",
                          "benefits_per_erp", "scripts_per_erp"))
nat_b <- nat_b[age_group %in% c("19-64", "65+")]
nat_b[, date := as.Date(date)]

nat_bs <- nat_b[, .(benefits = sum(benefits_per_erp, na.rm = TRUE),
                    scripts  = sum(scripts_per_erp,  na.rm = TRUE)),
                by = .(date, age_group)]
nat_bs[, bps := benefits / scripts]            # avg government benefit per script ($)
nat_bw <- dcast(nat_bs, date ~ age_group, value.var = "bps")
setnames(nat_bw, c("19-64", "65+"), c("gen_bps", "con_bps"))
nat_bw[, diff_bps := gen_bps - con_bps]
nat_bw[, event_time := to_event_time(date)]
nat_bw[, moy := month(date)]

# De-seasonalise the differenced benefit/script series (pre-period moy means)
seas_b <- nat_bw[event_time < 0, .(sm = mean(diff_bps, na.rm = TRUE)), by = moy]
nat_bw <- merge(nat_bw, seas_b, by = "moy", all.x = TRUE, sort = FALSE)
nat_bw[, diff_bps_adj := diff_bps - sm]
setorder(nat_bw, event_time)

fs_pre  <- nat_bw[event_time < 0, mean(diff_bps_adj, na.rm = TRUE)]
fs_post <- nat_bw[event_time >= 0 & event_time <= CLEAN_POST_MAX,
                  mean(diff_bps_adj, na.rm = TRUE)]
fs_jump <- fs_post - fs_pre
# Also the raw general step (more interpretable: the $ extra subsidy per general script)
gen_step <- nat_bw[event_time >= 0 & event_time <= CLEAN_POST_MAX, mean(gen_bps, na.rm = TRUE)] -
            nat_bw[event_time >= -12 & event_time < 0, mean(gen_bps, na.rm = TRUE)]
con_step <- nat_bw[event_time >= 0 & event_time <= CLEAN_POST_MAX, mean(con_bps, na.rm = TRUE)] -
            nat_bw[event_time >= -12 & event_time < 0, mean(con_bps, na.rm = TRUE)]

first_stage_tbl <- tibble(
  measure = c("General benefit/script step (clean post - 12mo pre, $)",
              "Concessional benefit/script step ($)",
              "De-seasonalised general-minus-concessional step ($)"),
  value = c(gen_step, con_step, fs_jump)
)
write_csv(nat_bw[, .(date, event_time, gen_bps, con_bps, diff_bps, diff_bps_adj)],
          file.path(tab_dir, "first_stage_benefit_per_script.csv"))
message("  First stage (general benefit/script step): ", sprintf("$%.2f", gen_step))
message("  Concessional step (placebo, should be ~0): ", sprintf("$%.2f", con_step))
message("  De-seas general-minus-concessional step:  ", sprintf("$%.2f", fs_jump))
print(first_stage_tbl)

fig35 <- ggplot(nat_bw[event_time >= -ES_LEAD_MAX & event_time <= ES_LAG_MAX]) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_line(aes(event_time, gen_bps, colour = "General (19-64)"), linewidth = 0.5) +
  geom_line(aes(event_time, con_bps, colour = "Concessional (65+)"), linewidth = 0.5) +
  scale_colour_manual(values = c("General (19-64)" = "steelblue",
                                 "Concessional (65+)" = "darkorange"), name = NULL) +
  annotate("text", x = -0.5, y = Inf, label = "Reform", hjust = 1.1, vjust = 1.4,
           colour = "red", size = 3, fontface = "italic") +
  labs(title = "First Stage: Government Benefit per Script (dose delivered?)",
       subtitle = paste0("If the reform was delivered, the general line steps up at Jan 2023 (govt pays the patient's ",
                         "$12.50). General step = ", sprintf("$%.2f", gen_step), "."),
       x = "Months relative to copayment reduction (0 = Jan 2023)",
       y = "Government benefit per script ($)")
ggsave(file.path(fig_dir, "fig35_first_stage_benefit_per_script.png"),
       fig35, width = 12, height = 6, dpi = 300)
message("  Saved fig35_first_stage_benefit_per_script.png")

# ==============================================================================
# 8. DOSE-RESPONSE ACROSS ATC CLASSES (national)
# ==============================================================================
# Treatment intensity ("bite") varies by drug price: the copayment cut only
# changes the patient price for scripts priced above $30. We proxy each class's
# delivered dose by the mechanical jump in GENERAL benefit/script at the reform,
# then ask whether classes with a bigger dose saw a bigger script response.
# Near-zero slope = inelastic demand / null robust across the price distribution.
# CAVEAT: ATC Level-2 aggregation blurs the within-class price distribution, so
# 'dose' is measured with error; read the slope's sign/magnitude, not precision.
message("\n=== Section 8: Dose-response across ATC classes ===")

nat_cls <- nat_b[, .(benefits = sum(benefits_per_erp, na.rm = TRUE),
                     scripts  = sum(scripts_per_erp,  na.rm = TRUE)),
                 by = .(atc_name, age_group, date)]
nat_cls[, date := as.Date(date)]
nat_cls[, event_time := to_event_time(date)]
nat_cls[, moy := month(date)]
nat_cls[, bps := benefits / scripts]
nat_cls[, period := fifelse(event_time < 0, "pre",
                     fifelse(event_time <= CLEAN_POST_MAX, "post_clean", "post_late"))]

# Per-class GENERAL dose = general benefit/script (clean post - pre)
gdose <- nat_cls[age_group == "19-64" & period %in% c("pre", "post_clean"),
                 .(bps = mean(bps, na.rm = TRUE)), by = .(atc_name, period)]
gdose_w <- dcast(gdose, atc_name ~ period, value.var = "bps")
gdose_w[, dose_bps := post_clean - pre]

# Per-class script DiD response: de-seasonalise the gen-minus-con scripts gap
sc <- nat_cls[, .(scripts = scripts), by = .(atc_name, age_group, date, event_time, moy)]
sc_w <- dcast(sc, atc_name + date + event_time + moy ~ age_group, value.var = "scripts")
setnames(sc_w, c("19-64", "65+"), c("gen_sc", "con_sc"))
sc_w[, diff_sc := gen_sc - con_sc]
seas_sc <- sc_w[event_time < 0, .(sm = mean(diff_sc, na.rm = TRUE)), by = .(atc_name, moy)]
sc_w <- merge(sc_w, seas_sc, by = c("atc_name", "moy"), all.x = TRUE, sort = FALSE)
sc_w[, diff_sc_adj := diff_sc - sm]
resp <- sc_w[, .(
  resp = mean(diff_sc_adj[event_time >= 0 & event_time <= CLEAN_POST_MAX], na.rm = TRUE) -
         mean(diff_sc_adj[event_time < 0], na.rm = TRUE),
  vol  = mean(gen_sc[event_time >= 0 & event_time <= CLEAN_POST_MAX], na.rm = TRUE)
), by = atc_name]

dose_resp <- merge(gdose_w[, .(atc_name, dose_bps)], resp, by = "atc_name")
dose_resp <- dose_resp[is.finite(dose_bps) & is.finite(resp)]
# Volume-weighted dose-response slope
dr_fit <- lm(resp ~ dose_bps, data = dose_resp, weights = vol)
dr_slope <- coef(dr_fit)[["dose_bps"]]
dr_p <- summary(dr_fit)$coefficients["dose_bps", "Pr(>|t|)"]
write_csv(dose_resp, file.path(tab_dir, "dose_response_classes.csv"))
message("  Dose-response slope (script response per $1 dose) = ",
        sprintf("%.4f", dr_slope), "  (p = ", sprintf("%.3f", dr_p), ")")
message("  Slope ~ 0 => no more response where the reform bit harder (inelastic).")

fig36 <- ggplot(dose_resp, aes(dose_bps, resp)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(aes(size = vol), alpha = 0.5, colour = "steelblue") +
  geom_smooth(aes(weight = vol), method = "lm", colour = "red", se = TRUE) +
  scale_size_continuous(guide = "none") +
  labs(title = "Dose-Response: Script Response vs Reform Bite by ATC Class",
       subtitle = paste0("x = mechanical general benefit/script jump ($ dose); y = de-seasonalised script DiD. ",
                         "Slope = ", sprintf("%.4f", dr_slope), " (p = ", sprintf("%.3f", dr_p), ")."),
       x = "Reform dose: general benefit/script jump ($, clean post - pre)",
       y = "De-seasonalised script response (general - concessional)")
ggsave(file.path(fig_dir, "fig36_dose_response_classes.png"),
       fig36, width = 10, height = 6, dpi = 300)
message("  Saved fig36_dose_response_classes.png")

# ==============================================================================
# 9. TRIPLE-DIFFERENCE (DDD): general x concessional x high/low reform bite
# ==============================================================================
# The strongest design here. Split classes into high-bite vs low-bite (top vs
# bottom dose tercile) and estimate the triple interaction. With age x date FE
# (absorbs the general-specific January seasonal spike, common across classes)
# AND bite-group x date FE (absorbs class-specific shocks), only the triple
# (general x post x high-bite) survives. Its identifying assumption - that the
# general-minus-concessional gap would have trended similarly across high- and
# low-bite classes absent the reform - is far weaker than simple parallel trends.
message("\n=== Section 9: Triple-difference (DDD) ===")

ddd_tbl <- tryCatch({
  terc <- quantile(dose_resp$dose_bps, c(1/3, 2/3), na.rm = TRUE)
  gdose_w[, bite := fifelse(dose_bps >= terc[2], "high",
                    fifelse(dose_bps <= terc[1], "low", "mid"))]
  bite_map <- gdose_w[, .(atc_name, bite)]

  # Estimate at the CLASS level (do NOT pre-collapse to 2 bite groups). Collapsing
  # left only 4 series (high/low x gen/con), so a cluster-robust VCOV was rank-
  # deficient and eigen() failed. Keeping classes gives ~tens of clusters and lets
  # us cluster on atc_name - the level at which 'bite' (drug price) is assigned.
  ddd <- merge(nat_cls, bite_map, by = "atc_name")
  ddd <- ddd[bite %in% c("high", "low") & event_time >= -ES_LEAD_MAX & is.finite(scripts)]
  ddd[, scripts_per_erp := scripts]
  ddd[, is_general := as.integer(age_group == "19-64")]
  ddd[, is_high := as.integer(bite == "high")]
  ddd[, post := as.integer(event_time >= 0)]
  ddd[, clean := as.integer(event_time >= 0 & event_time <= CLEAN_POST_MAX)]
  ddd[, ddd_term := post * is_general * is_high]
  ddd[, ddd_clean := clean * is_general * is_high]

  # atc_name^date absorbs every class-specific time shock (incl. bite-group x time);
  # age_group^date absorbs the general-specific January seasonal spike; atc_name^
  # age_group absorbs class x age levels. Only the triple (general x post x high)
  # survives. Cluster by class (atc_name), the unit bite intensity is assigned at.
  ddd_full <- feols(scripts_per_erp ~ ddd_term
                    | atc_name^age_group + age_group^date + atc_name^date,
                    data = ddd, cluster = ~atc_name)
  ddd_cln  <- feols(scripts_per_erp ~ ddd_clean
                    | atc_name^age_group + age_group^date + atc_name^date,
                    data = ddd[event_time <= CLEAN_POST_MAX | event_time < 0],
                    cluster = ~atc_name)

  # Wild cluster bootstrap p-values (restricted / WCR; Cameron-Gelbach-Miller 2008,
  # Roodman et al. 2019). Hand-rolled in base R + fixest so NO external package is
  # required (fwildclusterboot is not built for every R version). The triple-
  # difference clusters on ATC class and the high-bite (treated) tercile holds only
  # ~28 of 86 classes, so few-treated-cluster inference is warranted. Procedure:
  # (1) fit the restricted model imposing H0: beta = 0 (drop the DDD term); (2) draw
  # Rademacher weights w_g in {-1,+1}, one per ATC-class cluster; (3) form
  # y* = fitted_restricted + w_g * resid_restricted; (4) refit the FULL model on y*
  # and collect the cluster-robust t-stat; (5) p = share of |t*| >= |t0| over
  # B = 9,999 replications (with the +1/+1 correction). This refits feols ~10k times
  # per window and may take several minutes; lower WCB_B to trade precision for speed.
  FE_DDD <- "atc_name^age_group + age_group^date + atc_name^date"
  WCB_B  <- 9999
  wcb_p_manual <- function(d, param, B = WCB_B, seed = 20230101) {
    d <- as.data.frame(d)
    if (length(unique(d[[param]])) < 2) return(NA_real_)
    fit_r <- feols(as.formula(paste0("scripts_per_erp ~ 1 | ", FE_DDD)),
                   data = d, cluster = ~atc_name)
    used  <- obs(fit_r)                      # rows fixest actually kept
    d     <- d[used, , drop = FALSE]
    res_r <- as.numeric(resid(fit_r))
    fit_v <- d$scripts_per_erp - res_r       # restricted fitted values
    fit_full <- feols(as.formula(paste0("scripts_per_erp ~ ", param, " | ", FE_DDD)),
                      data = d, cluster = ~atc_name)
    t0 <- as.numeric(coef(fit_full)[param] / sqrt(diag(vcov(fit_full)))[param])
    groups <- unique(d$atc_name); idx <- match(d$atc_name, groups); ng <- length(groups)
    fb_form <- as.formula(paste0("y_star ~ ", param, " | ", FE_DDD))
    set.seed(seed); tstar <- numeric(B)
    for (b in seq_len(B)) {
      w <- sample(c(-1, 1), ng, replace = TRUE)        # Rademacher, per cluster
      d$y_star <- fit_v + w[idx] * res_r
      fb <- feols(fb_form, data = d, cluster = ~atc_name)
      tstar[b] <- as.numeric(coef(fb)[param] / sqrt(diag(vcov(fb)))[param])
      if (b %% 1000 == 0) message(sprintf("    WCB %s: %d/%d", param, b, B))
    }
    (1 + sum(abs(tstar) >= abs(t0), na.rm = TRUE)) / (1 + B)
  }
  wcb_full <- tryCatch(wcb_p_manual(ddd, "ddd_term"),
                       error = function(e) { message("  WCB full failed: ", conditionMessage(e)); NA_real_ })
  wcb_cln  <- tryCatch(wcb_p_manual(ddd[event_time <= CLEAN_POST_MAX | event_time < 0], "ddd_clean"),
                       error = function(e) { message("  WCB clean failed: ", conditionMessage(e)); NA_real_ })
  message(sprintf("  Wild cluster bootstrap p (restricted, B=%d): full = %s, clean = %s",
                  WCB_B, format(wcb_full), format(wcb_cln)))

  tbl <- bind_rows(
    broom::tidy(ddd_full, conf.int = TRUE) %>% filter(term == "ddd_term") %>%
      mutate(spec = "DDD full window (post)"),
    broom::tidy(ddd_cln, conf.int = TRUE) %>% filter(term == "ddd_clean") %>%
      mutate(spec = "DDD clean pre-Nov-2023 window")
  ) %>% select(spec, estimate, std.error, statistic, p.value, conf.low, conf.high)
  tbl$wcb_p <- c(wcb_full, wcb_cln)
  write_csv(tbl, file.path(tab_dir, "ddd_results.csv"))
  message("  DDD (general x post x high-bite) triple coefficient:")
  print(tbl)
  message("  Near-zero / n.s. triple => no reform effect even where bite was largest,")
  message("  net of seasonality AND class shocks (the cleanest available null).")
  tbl
}, error = function(e) {
  message("  Section 9 DDD failed: ", conditionMessage(e), " (skipped, continuing)")
  NULL
})

# ==============================================================================
# 10. HONESTDID SENSITIVITY (Rambachan & Roth) on de-seasonalised event study
# ==============================================================================
# Parallel trends is violated, so the DiD point estimate is not trustworthy as a
# point. HonestDiD asks: GIVEN the pre-period violations we actually observe, how
# large could the true post-period effect be under bounded extrapolations of the
# differential trend? If the robust confidence set brackets 0 (and excludes large
# positive effects), 'unidentified' becomes 'bounded near zero'. We use a compact
# clean window (-12..+9) so the relative-magnitudes reference is well-posed.
message("\n=== Section 10: HonestDiD sensitivity bounds ===")

honestdid_ok <- requireNamespace("HonestDiD", quietly = TRUE)
if (!honestdid_ok) {
  message("  HonestDiD not installed; attempting install...")
  try(install.packages("HonestDiD", repos = "https://cloud.r-project.org"),
      silent = TRUE)
  honestdid_ok <- requireNamespace("HonestDiD", quietly = TRUE)
}

if (honestdid_ok) {
  tryCatch({
    es_hd <- ds_frame[event_time >= -12 & event_time <= 9]
    fit_hd <- feols(y_adj ~ i(event_time, is_general, ref = EVENT_REF)
                    | lga_code^age_group + date, data = es_hd, cluster = ~lga_code)
    b <- coef(fit_hd); V <- vcov(fit_hd)
    idx <- grep("event_time::", names(b))
    et  <- as.integer(str_extract(names(b)[idx], "-?\\d+"))
    o <- order(et); idx <- idx[o]; et <- et[o]
    betahat <- as.numeric(b[idx]); sigma <- as.matrix(V[idx, idx])
    numPre  <- sum(et < EVENT_REF)
    numPost <- sum(et >= 0)
    l_vec <- rep(1 / numPost, numPost)

    orig <- HonestDiD::constructOriginalCS(betahat, sigma, numPre, numPost,
                                           l_vec = l_vec)
    sens <- HonestDiD::createSensitivityResults_relativeMagnitudes(
      betahat, sigma, numPre, numPost, l_vec = l_vec,
      Mbarvec = c(0, 0.5, 1, 1.5, 2))

    # Coerce to plain numerics: constructOriginalCS returns lb/ub as 1x1 matrix
    # cells, which become matrix-typed tibble columns and break printing/write_csv.
    hd_tbl <- bind_rows(
      tibble(Mbar = NA_real_,
             lb = as.numeric(orig$lb), ub = as.numeric(orig$ub),
             method = "Original (parallel-trends CI)"),
      tibble(Mbar = as.numeric(sens$Mbar),
             lb = as.numeric(sens$lb), ub = as.numeric(sens$ub),
             method = "Relative magnitudes")
    )
    write_csv(hd_tbl, file.path(tab_dir, "honestdid_sensitivity.csv"))
    message("  Average post-period effect (de-seasonalised), robust CIs:")
    print(as.data.frame(hd_tbl))
    excl0 <- with(hd_tbl, (lb > 0) | (ub < 0))
    if (any(excl0, na.rm = TRUE)) {
      message("  Robust CI EXCLUDES 0 at Mbar = ",
              paste(sprintf("%.1f", hd_tbl$Mbar[which(excl0)]), collapse = ", "),
              " (an effect survives that much pre-trend extrapolation)")
    } else {
      message("  Robust CI INCLUDES 0 at EVERY Mbar (0..2): the null is robust even ",
              "when the post differential trend is allowed to be up to 2x the ",
              "largest pre-period violation.")
    }
    hd_fig <- HonestDiD::createSensitivityPlot_relativeMagnitudes(sens, orig)
    ggsave(file.path(fig_dir, "fig37_honestdid_sensitivity.png"),
           hd_fig, width = 9, height = 5.5, dpi = 300)
    message("  Saved fig37_honestdid_sensitivity.png")
  }, error = function(e) {
    message("  Section 10 HonestDiD failed: ", conditionMessage(e), " (skipped, continuing)")
  })
} else {
  message("  HonestDiD unavailable; skipped (install manually to enable).")
}

# ==============================================================================
# 11. PLACEBO-IN-TIME (fake January reform dates)
# ==============================================================================
# The most legible proof of the seasonal artifact: every January shows a Dec->Jan
# jump in the general-minus-concessional gap (the safety-net reset), not just Jan
# 2023. If the true reform year is unremarkable among the placebo years, the raw
# 'effect' is seasonality, not the policy.
message("\n=== Section 11: Placebo-in-time (fake January reform dates) ===")

# Reuse the national differenced scripts gap (diff_gc) built in Section 1.
plac <- nat_wide[, .(date, diff_gc)]
plac[, yr := year(date)]
plac[, mo := month(date)]
jan_delta <- map_dfr(2016:2025, function(y) {
  jan <- plac[yr == y & mo == 1, diff_gc]
  dec <- plac[yr == (y - 1) & mo == 12, diff_gc]
  if (length(jan) == 1 && length(dec) == 1)
    tibble(year = y, dec_to_jan_jump = jan - dec) else NULL
})
jan_delta <- jan_delta %>% mutate(is_reform = year == 2023)
write_csv(jan_delta, file.path(tab_dir, "placebo_in_time.csv"))
message("  Dec->Jan jump in general-minus-concessional gap, by year:")
print(jan_delta)
reform_rank <- mean(jan_delta$dec_to_jan_jump <= jan_delta$dec_to_jan_jump[jan_delta$is_reform])
message("  2023 jump is at the ", sprintf("%.0f%%", 100 * reform_rank),
        " percentile of all Dec->Jan jumps (50% = perfectly unremarkable).")

fig38 <- ggplot(jan_delta, aes(factor(year), dec_to_jan_jump, fill = is_reform)) +
  geom_col() +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "red"),
                    labels = c("FALSE" = "Placebo Januaries",
                               "TRUE" = "Reform (Jan 2023)"), name = NULL) +
  labs(title = "Placebo-in-Time: Every January Shows the Safety-Net-Reset Jump",
       subtitle = "If Jan 2023 is unremarkable among placebo Januaries, the raw 'effect' is seasonal.",
       x = "Year", y = "Dec -> Jan jump in general-minus-concessional gap")
ggsave(file.path(fig_dir, "fig38_placebo_in_time.png"),
       fig38, width = 10, height = 5.5, dpi = 300)
message("  Saved fig38_placebo_in_time.png")

# ==============================================================================
# 12. SUMMARY
# ==============================================================================
message("\n", paste(rep("=", 70), collapse = ""))
message("EVENT-STUDY / DiD ANALYSIS COMPLETE")
message(paste(rep("=", 70), collapse = ""))
message("\nIdentification-threat coverage:")
message("  (1) Pre-period sensitivity -> Section 4 anchor grid (fig32)")
message("  (2) Event study           -> Sections 1-2 (fig29 national, fig30 LGA)")
message("  (3) COVID/Omicron + geo   -> month FE retires the binary dummy;")
message("                               Specs B/C add state^date & remoteness^date FE")
message("  (4) Bulk-billing confound -> Spec E clean window + Spec D explicit control")
message("  (S6) Seasonality artefact -> Section 6 de-seasonalises the outcome,")
message("                               re-runs event study + ATT + anchor grid")
message("  (S7) First stage / dose   -> Section 7 benefit-per-script step (fig35):")
message("                               did the reform actually change the price?")
message("  (S8) Dose-response        -> Section 8 per-class bite vs response (fig36)")
message("  (S9) Triple-difference    -> Section 9 DDD nets out seasonal + class shocks")
message("  (S10) HonestDiD bounds    -> Section 10 R-R sensitivity (fig37)")
message("  (S11) Placebo-in-time     -> Section 11 fake-January reforms (fig38)")
message("\nKey numbers to read off:")
message("  - RAW pre-trend Wald p (fig30): ", sprintf("%.4f", pretrend_wald$p),
        "  (rejected; leads dominated by January seasonality)")
message("  - DE-SEAS pre-trend Wald p (fig33): ", sprintf("%.4f", pretrend_wald_ds$p))
message("  - DE-SEAS mean post coef ", sprintf("%.4f", mean_post_ds),
        " vs mean |pre coef| ", sprintf("%.4f", mean_pre_ds),
        " (effect must clear the pre-period band to be real)")
message("  - Headline de-seasonalised clean ATT: see did_att_deseasonalised.csv")
message("  - FIRST STAGE: general benefit/script step (fig35) - if large & positive,")
message("    the price mechanism fired; flat scripts then => inelastic-demand null.")
message("  - DOSE-RESPONSE slope (fig36): ~0 across classes => no price gradient.")
message("  - DDD triple coefficient (ddd_results.csv): cleanest identified estimate.")
message("  - HonestDiD breakdown Mbar (fig37): how far parallel-trends must bend")
message("    before the CI includes 0 (small Mbar => fragile).")
message("  - PLACEBO: 2023 Dec->Jan jump percentile (Section 11) - if not extreme,")
message("    the 'reform' jump is indistinguishable from any other January.")
message("  - INTERPRETATION: a GENUINE null requires (a) a firing first stage, (b) a")
message("    flat dose-response, (c) a DDD near zero, AND (d) HonestDiD bounds that")
message("    stay near zero. Non-identification looks different: large first stage but")
message("    pre-trends/placebo that swamp any post step.")
message("\nCAVEAT to carry into the manuscript Methods/Limitations:")
message("  The 65+/19-64 age proxy is leaky - concession-card holders include")
message("  working-age Health Care Card / pensioner / DVA / CSHC holders, so the")
message("  control arm is impure and the Nov-2023 bulk-billing change touched both")
message("  arms. The DiD is therefore an attenuated lower bound on the contrast.")
message("\nTables: event_study_national.csv, event_study_lga.csv,")
message("        event_study_pretrend_tests.csv, did_att_specifications.csv,")
message("        did_preperiod_anchor_grid.csv, event_study_lga_deseasonalised.csv,")
message("        did_att_deseasonalised.csv, did_preperiod_anchor_grid_deseasonalised.csv,")
message("        first_stage_benefit_per_script.csv, dose_response_classes.csv,")
message("        ddd_results.csv, honestdid_sensitivity.csv, placebo_in_time.csv")
message("Figures: fig29-fig38")
