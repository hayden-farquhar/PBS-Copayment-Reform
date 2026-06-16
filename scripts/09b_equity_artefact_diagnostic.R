# ==============================================================================
# 09b_equity_artefact_diagnostic.R
#
# Is the stratified "distributional response" (script 09) a real equity signal
# or the same seasonality/proxy artefact the national analysis already exposed?
#
# Decisive test (mirrors what killed the NATIONAL positive in script 08): vary
# the pre-period ANCHOR and check whether each stratum's clean-window ATT is
# sign-stable and magnitude-stable, on de-seasonalised vs raw outcomes. A real
# effect is anchor-robust; an artefact reverses/attenuates across anchors.
#
# Uses the cached LGA aggregate from script 09 (no 2.9 GB reload).
# Output: outputs/tables/equity_artefact_anchor_grid.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(fixest); library(broom); library(lubridate)
})
proj_dir <- here::here()
proc_dir <- file.path(proj_dir, "data", "processed")
tab_dir  <- file.path(proj_dir, "outputs", "tables")

intervention_date <- as.Date("2023-01-01"); CLEAN_POST_MAX <- 9L
to_event_time <- function(d) (year(d)-year(intervention_date))*12 + (month(d)-month(intervention_date))

cache <- file.path(proc_dir, "_lga_agg_equity_cache.rds")
stopifnot(file.exists(cache))
lga <- readRDS(cache)
lga[, is_general := as.integer(age_group == "19-64")]
lga[, event_time := to_event_time(date)]
lga[, month_of_year := month(date)]
if (!"seifa_quintile" %in% names(lga)) stop("seifa_quintile missing from cache")
lga <- lga[!is.na(seifa_quintile)]

# De-seasonalise within a given pre-period start (anchor): pre-period means use
# event_time in [anchor, -1]; clean window = event_time 0..9.
fit_stratum <- function(dt, anchor) {
  d <- dt[event_time >= anchor & event_time <= CLEAN_POST_MAX]
  seas <- d[event_time < 0,
            .(seas_mean = mean(scripts_per_erp, na.rm = TRUE)),
            by = .(lga_code, age_group, month_of_year)]
  d <- merge(d, seas, by = c("lga_code","age_group","month_of_year"), all.x = TRUE, sort = FALSE)
  d[, y_adj := scripts_per_erp - seas_mean]
  d <- d[!is.na(y_adj)]
  d[, clean_post := as.integer(event_time >= 0)]
  d[, clean_post_general := clean_post * is_general]
  out <- list()
  for (yv in c("y_adj","scripts_per_erp")) {
    f <- feols(as.formula(paste(yv, "~ clean_post_general | lga_code^age_group + date")),
               data = d, cluster = ~lga_code)
    r <- broom::tidy(f, conf.int = TRUE); r <- r[r$term == "clean_post_general", ]
    out[[yv]] <- data.frame(outcome = ifelse(yv=="y_adj","de-seasonalised","raw"),
                            estimate = r$estimate, conf.low = r$conf.low,
                            conf.high = r$conf.high, p.value = r$p.value)
  }
  do.call(rbind, out)
}

anchors <- c("-48"=-48L, "-36"=-36L, "-24"=-24L, "-12"=-12L, "-6"=-6L)  # pre-period starts
res <- list()
for (q in sort(unique(lga$seifa_quintile))) {
  for (a in names(anchors)) {
    rr <- fit_stratum(lga[seifa_quintile == q], anchors[[a]])
    rr$stratum <- paste0("SEIFA Q", q); rr$anchor_months <- a
    res[[paste0("s",q,a)]] <- rr
  }
}
for (ra in sort(unique(lga$remoteness_area))) {
  sub <- lga[remoteness_area == ra]
  if (uniqueN(sub$lga_code) < 5) next
  for (a in names(anchors)) {
    rr <- fit_stratum(sub, anchors[[a]])
    rr$stratum <- ra; rr$anchor_months <- a
    res[[paste0("r",substr(ra,1,4),a)]] <- rr
  }
}
grid <- do.call(rbind, res)
grid <- grid[, c("stratum","anchor_months","outcome","estimate","conf.low","conf.high","p.value")]
fwrite(grid, file.path(tab_dir, "equity_artefact_anchor_grid.csv"))

# Sign-stability summary on the de-seasonalised series across anchors
ds <- grid[grid$outcome == "de-seasonalised", ]
cat("\n=== De-seasonalised clean-window ATT across pre-period anchors ===\n")
for (s in unique(ds$stratum)) {
  sub <- ds[ds$stratum == s, ]
  signs <- sign(sub$estimate); flips <- length(unique(signs[signs != 0])) > 1
  sig <- sum(sub$p.value < 0.05, na.rm = TRUE)
  cat(sprintf("%-30s est range [% .4f, % .4f]  sign-flips=%s  sig@5%%=%d/%d\n",
              s, min(sub$estimate), max(sub$estimate), flips, sig, nrow(sub)))
}
cat("\nRule: a real effect is sign-stable + persistently significant across anchors;\n")
cat("an artefact flips sign or loses significance as the anchor moves.\n")
