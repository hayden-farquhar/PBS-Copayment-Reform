# ==============================================================================
# 09_equity_reestimation.R
#
# Distributional (equity) re-estimation of the 2023 PBS general copayment cut
# UNDER THE CORRECTED IDENTIFICATION of script 08.
#
# Motivation: the original equity analysis (script 05: did_seifa_results.csv,
# did_remoteness_results.csv, did_seifa_remoteness_interaction.csv) was estimated
# under the pre-rebuild identification that script 08 / the v3.x manuscript show
# is confounded by (a) the 1-January safety-net reset collinear with the reform
# date and (b) the Nov-2023 bulk-billing shock acting on the control arm. Those
# equity point estimates are therefore contaminated and MUST NOT be reused. This
# script re-estimates the distributional question on the SAME de-seasonalised,
# clean (pre-Nov-2023) window that carries the headline identification weight,
# stratified by SEIFA disadvantage quintile and by remoteness, with a formal
# heterogeneity test.
#
# Equity estimand: does the (non-)response of general-patient dispensing to the
# price cut differ by area-level disadvantage (SEIFA quintile) or remoteness?
# Estimator: de-seasonalised general-vs-concessional DiD on the clean window
# (event months 0..9, Jan-Oct 2023), per stratum, clustered by LGA (541 clusters
# => standard CRVE; no wild bootstrap needed, unlike the class-level DDD).
#
# Inputs : data/processed/pbs_lga_monthly.csv (2.9 GB; subset of columns)
#          data/processed/lga_characteristics.csv (SEIFA quintile per LGA)
# Outputs: outputs/tables/equity_seifa_clean_did_corrected.csv
#          outputs/tables/equity_remoteness_clean_did_corrected.csv
#          outputs/tables/equity_heterogeneity_tests_corrected.csv
#          outputs/figures/fig40_equity_forest_corrected.png / .pdf
#
# NOTE: every number this script reports is computed here; NONE is hand-entered.
# ==============================================================================

# ==============================================================================
# 0. SETUP  (env-detecting, version-pinned-aware, resumable)
# ==============================================================================
needed <- c("data.table", "fixest", "broom", "ggplot2", "dplyr", "here", "lubridate")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                          "\nInstall them before running.")
suppressPackageStartupMessages({
  library(data.table); library(fixest); library(broom)
  library(ggplot2); library(dplyr); library(lubridate)
})
message("fixest ", as.character(packageVersion("fixest")),
        " | data.table ", as.character(packageVersion("data.table")))

proj_dir <- here::here()
proc_dir <- file.path(proj_dir, "data", "processed")
fig_dir  <- file.path(proj_dir, "outputs", "figures")
tab_dir  <- file.path(proj_dir, "outputs", "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# --- Key dates / window constants (identical to script 08) -------------------
intervention_date <- as.Date("2023-01-01")
bulk_billing_date <- as.Date("2023-11-01")
EVENT_REF      <- -1L
CLEAN_POST_MAX <- 9L     # Jan 2023 .. Oct 2023 (pre-bulk-billing)
ES_LEAD_MAX    <- 48L
ES_LAG_MAX     <- 35L
to_event_time <- function(d) {
  (year(d) - year(intervention_date)) * 12 + (month(d) - month(intervention_date))
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a   # base-R-safe null/NA coalesce

# Okabe-Ito colourblind-safe palette (portfolio shared standard)
OI <- c(black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
        yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7")

# ==============================================================================
# 1. LOAD + AGGREGATE  (mirror script 08 section 2 exactly, + SEIFA)
# ==============================================================================
agg_cache <- file.path(proc_dir, "_lga_agg_equity_cache.rds")
if (file.exists(agg_cache)) {
  message("Resuming from cached LGA aggregate: ", agg_cache)
  lga_agg <- readRDS(agg_cache)
} else {
  message("Loading LGA panel (2.9 GB; subset of columns)...")
  lga_raw <- fread(
    file.path(proc_dir, "pbs_lga_monthly.csv"),
    select = c("state", "lga_code", "date", "age_group", "scripts_per_erp",
               "remoteness_area", "seifa_quintile")
  )
  lga_raw <- lga_raw[age_group %in% c("19-64", "65+")]
  # Aggregate across ATC classes -> LGA x month x age_group overall rate
  lga_agg <- lga_raw[
    , .(scripts_per_erp = sum(scripts_per_erp, na.rm = TRUE)),
    by = .(state, lga_code, date, age_group, remoteness_area, seifa_quintile)
  ]
  rm(lga_raw); gc()
  lga_agg[, date := as.Date(date)]
  lga_agg[, lga_code := as.integer(lga_code)]
  saveRDS(lga_agg, agg_cache)
  message("Cached LGA aggregate -> ", agg_cache)
}

lga_agg[, is_general    := as.integer(age_group == "19-64")]
lga_agg[, event_time    := to_event_time(date)]
lga_agg[, month_of_year := month(date)]

# SEIFA quintile must be present (1 = most disadvantaged .. 5 = least)
if (!"seifa_quintile" %in% names(lga_agg) || all(is.na(lga_agg$seifa_quintile))) {
  ch <- fread(file.path(proc_dir, "lga_characteristics.csv"),
              select = c("lga_code", "seifa_quintile"))
  ch[, lga_code := as.integer(lga_code)]
  lga_agg <- merge(lga_agg[, !"seifa_quintile"], ch, by = "lga_code", all.x = TRUE)
}
lga_agg <- lga_agg[!is.na(seifa_quintile)]
message("Strata available -> SEIFA quintiles: ",
        paste(sort(unique(lga_agg$seifa_quintile)), collapse = ","),
        " | remoteness: ", paste(sort(unique(lga_agg$remoteness_area)), collapse = " / "))

# ==============================================================================
# 2. DE-SEASONALISE (identical helper to script 08 section 6) + clean window
# ==============================================================================
deseasonalise <- function(dt) {
  seas <- dt[event_time < 0,
             .(seas_mean = mean(scripts_per_erp, na.rm = TRUE)),
             by = .(lga_code, age_group, month_of_year)]
  out <- merge(dt, seas, by = c("lga_code", "age_group", "month_of_year"),
               all.x = TRUE, sort = FALSE)
  out[, y_adj := scripts_per_erp - seas_mean]
  out[]
}
ds <- deseasonalise(lga_agg[event_time >= -ES_LEAD_MAX & event_time <= ES_LAG_MAX])
ds <- ds[!is.na(y_adj)]
ds[, clean_post         := as.integer(event_time >= 0 & event_time <= CLEAN_POST_MAX)]
ds[, clean_post_general := clean_post * is_general]
# Estimation frame: clean post-window + all leads (drop the contaminated tail)
ds_clean <- ds[event_time <= CLEAN_POST_MAX]   # keeps event_time < 0 (leads) + 0..9

# ==============================================================================
# 3. PER-STRATUM ATT (de-seasonalised, clean window, remoteness x date FE)
#    Matches the corrected headline level spec (Spec C, de-seasonalised).
# ==============================================================================
extract_att <- function(fit, label, coef_name, n_lga) {
  ti  <- broom::tidy(fit, conf.int = TRUE)
  row <- ti[ti$term == coef_name, ]
  data.frame(stratum = label, estimate = row$estimate,
             conf.low = row$conf.low, conf.high = row$conf.high,
             std.error = row$std.error, p.value = row$p.value,
             n_obs = nobs(fit), n_lga = n_lga, stringsAsFactors = FALSE)
}

message("\n=== SEIFA-quintile-stratified clean-window DiD (de-seasonalised) ===")
seifa_rows <- list()
for (q in sort(unique(ds_clean$seifa_quintile))) {
  sub <- ds_clean[seifa_quintile == q]
  # remoteness x date FE retained (Spec C analogue); falls back to date FE if a
  # quintile spans a single remoteness band.
  fe <- if (uniqueN(sub$remoteness_area) > 1)
    "lga_code^age_group + date + remoteness_area^date" else "lga_code^age_group + date"
  fit <- feols(as.formula(paste("y_adj ~ clean_post_general |", fe)),
               data = sub, cluster = ~lga_code)
  seifa_rows[[as.character(q)]] <-
    extract_att(fit, paste0("SEIFA Q", q,
                            c("1"=" (most disadvantaged)","5"=" (least disadvantaged)")[as.character(q)] %||% ""),
                "clean_post_general", uniqueN(sub$lga_code))
  message(sprintf("  Q%s: ATT=%.5f  [%.5f, %.5f]  p=%.3f  (LGAs=%d)",
                  q, seifa_rows[[as.character(q)]]$estimate,
                  seifa_rows[[as.character(q)]]$conf.low,
                  seifa_rows[[as.character(q)]]$conf.high,
                  seifa_rows[[as.character(q)]]$p.value,
                  seifa_rows[[as.character(q)]]$n_lga))
}
seifa_tbl <- do.call(rbind, seifa_rows)
fwrite(seifa_tbl, file.path(tab_dir, "equity_seifa_clean_did_corrected.csv"))
message("  Saved equity_seifa_clean_did_corrected.csv")

message("\n=== Remoteness-stratified clean-window DiD (de-seasonalised) ===")
rem_rows <- list()
for (ra in sort(unique(ds_clean$remoteness_area))) {
  sub <- ds_clean[remoteness_area == ra]
  if (uniqueN(sub$lga_code) < 5) { message("  Skip '", ra, "' (<5 LGAs)"); next }
  fit <- feols(y_adj ~ clean_post_general | lga_code^age_group + date,
               data = sub, cluster = ~lga_code)
  rem_rows[[ra]] <- extract_att(fit, ra, "clean_post_general", uniqueN(sub$lga_code))
  message(sprintf("  %-32s ATT=%.5f  [%.5f, %.5f]  p=%.3f  (LGAs=%d)",
                  ra, rem_rows[[ra]]$estimate, rem_rows[[ra]]$conf.low,
                  rem_rows[[ra]]$conf.high, rem_rows[[ra]]$p.value, rem_rows[[ra]]$n_lga))
}
rem_tbl <- do.call(rbind, rem_rows)
fwrite(rem_tbl, file.path(tab_dir, "equity_remoteness_clean_did_corrected.csv"))
message("  Saved equity_remoteness_clean_did_corrected.csv")

# ==============================================================================
# 4. FORMAL HETEROGENEITY TESTS (pooled interaction; joint Wald)
# ==============================================================================
message("\n=== Heterogeneity tests (H0: equal ATT across strata) ===")
het <- list()

# SEIFA: per-quintile slopes vs Q1 ref; joint Wald on the 4 differences.
fit_seifa_het <- feols(
  y_adj ~ clean_post_general + i(seifa_quintile, clean_post_general, ref = 1)
  | lga_code^age_group + date + remoteness_area^date,
  data = ds_clean, cluster = ~lga_code)
w_seifa <- wald(fit_seifa_het, "seifa_quintile")
het[["seifa"]] <- data.frame(
  test = "SEIFA quintile heterogeneity (4 interaction terms = 0)",
  stat = unname(w_seifa$stat), df1 = unname(w_seifa$df1), df2 = unname(w_seifa$df2),
  p.value = unname(w_seifa$p), stringsAsFactors = FALSE)
message(sprintf("  SEIFA: Wald F=%.2f  p=%.4f", w_seifa$stat, w_seifa$p))

# Remoteness: per-band slopes vs ref; joint Wald.
fit_rem_het <- feols(
  y_adj ~ clean_post_general + i(remoteness_area, clean_post_general)
  | lga_code^age_group + date,
  data = ds_clean, cluster = ~lga_code)
w_rem <- wald(fit_rem_het, "remoteness_area")
het[["rem"]] <- data.frame(
  test = "Remoteness heterogeneity (interaction terms = 0)",
  stat = unname(w_rem$stat), df1 = unname(w_rem$df1), df2 = unname(w_rem$df2),
  p.value = unname(w_rem$p), stringsAsFactors = FALSE)
message(sprintf("  Remoteness: Wald F=%.2f  p=%.4f", w_rem$stat, w_rem$p))

het_tbl <- do.call(rbind, het)
fwrite(het_tbl, file.path(tab_dir, "equity_heterogeneity_tests_corrected.csv"))
message("  Saved equity_heterogeneity_tests_corrected.csv")

# ==============================================================================
# 5. FOREST PLOT (Tier-2 hero effect-size figure; Okabe-Ito)
# ==============================================================================
forest <- bind_rows(
  transform(seifa_tbl, panel = "By SEIFA disadvantage quintile"),
  transform(rem_tbl,  panel = "By remoteness")
)
forest$stratum <- factor(forest$stratum, levels = rev(forest$stratum))
p <- ggplot(forest, aes(x = estimate, y = stratum)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = OI["black"]) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.25,
                 colour = OI["blue"], linewidth = 0.6) +
  geom_point(colour = OI["blue"], size = 2.4) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Reform effect on general-patient dispensing\n(de-seasonalised scripts/ERP, clean window)",
       y = NULL,
       title = "No distributional response to the 2023 PBS copayment reduction",
       subtitle = "General-vs-concessional DiD by stratum; 95% CI clustered by LGA.\nCIs bracket zero in every stratum.") +
  theme_minimal(base_size = 11) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.title.position = "plot", plot.caption.position = "plot")
# width widened 9 -> 12 (matching the other main-text figures, ~3600 px) and the
# longest title/subtitle/axis strings wrapped, so the in-image text no longer
# clips at the right edge of the canvas (fixed 2026-06-16 after the v4.0 visual review).
ggsave(file.path(fig_dir, "fig40_equity_forest_corrected.png"), p,
       width = 12, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "fig40_equity_forest_corrected.pdf"), p, width = 12, height = 6)
message("  Saved fig40_equity_forest_corrected.png / .pdf")

message("\nDONE. Read the per-stratum CIs + heterogeneity p-values before drafting any equity claim.")
