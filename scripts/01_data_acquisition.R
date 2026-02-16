# ==============================================================================
# 01_data_acquisition.R
# Download and import AIHW PBS monthly data and ABS linkage files
#
# Data sources:
#   1. AIHW PBS Monthly Data (supplementary Excel tables)
#   2. ABS SEIFA 2021 (IRSD) at LGA level
#   3. ABS Remoteness Areas 2021 (SA1-to-RA allocation)
#   4. ABS LGA boundary shapefiles (2021 ASGS Edition 3)
#   5. ABS Estimated Resident Population by LGA
#
# Outputs: Raw files saved to data/raw/ and data/spatial/
# ==============================================================================

library(httr)
library(readxl)
library(utils)

# -- Paths --------------------------------------------------------------------
proj_dir  <- here::here()
raw_dir   <- file.path(proj_dir, "data", "raw")
spatial_dir <- file.path(proj_dir, "data", "spatial")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

# -- Helper: download with status reporting -----------------------------------
safe_download <- function(url, destfile, description = "") {
  if (file.exists(destfile)) {
    message("  [SKIP] Already exists: ", basename(destfile))
    return(invisible(TRUE))
  }
  message("  [DOWN] ", description, " -> ", basename(destfile))
  tryCatch({
    resp <- GET(url, write_disk(destfile, overwrite = TRUE),
                progress(), timeout(300))
    if (http_error(resp)) {
      warning("  [FAIL] HTTP ", status_code(resp), " for ", url)
      unlink(destfile)
      return(invisible(FALSE))
    }
    message("  [OK]   ", round(file.size(destfile) / 1e6, 1), " MB")
    invisible(TRUE)
  }, error = function(e) {
    warning("  [FAIL] ", e$message)
    unlink(destfile)
    invisible(FALSE)
  })
}

# ==============================================================================
# 1. AIHW PBS MONTHLY DATA
# ==============================================================================
# The AIHW provides PBS prescription data through:
#   - Interactive dashboard (Jan 2015 onwards)
#   - Downloadable supplementary Excel tables (historical + current)
#
# The dashboard is at:
#   https://www.aihw.gov.au/reports/medicines/pbs-monthly-data/contents/dashboard
#
# Data tables are at:
#   https://www.aihw.gov.au/reports/medicines/pbs-monthly-data/data
#
# Available files (from the /data page):
#   - Historical data tables by ATC Level 1 (May 2002 - Jan 2015): ~111 MB XLSX
#   - Historical data tables by ATC Level 2 (May 2002 - Jan 2015): ~408 MB XLSX
#   - Current supplementary data tables (Jan 2015 onwards): XLSX
#
# IMPORTANT: The AIHW data page uses JavaScript-rendered download links that
# cannot be reliably scraped. Manual download is the most reliable approach.
# ==============================================================================

message("\n=== 1. AIHW PBS Monthly Data ===")
message("
MANUAL DOWNLOAD REQUIRED for AIHW PBS data.

Steps:
  1. Go to: https://www.aihw.gov.au/reports/medicines/pbs-monthly-data/data
  2. Download ALL available supplementary data tables (Excel workbooks)
     - Look for 'Data tables' or 'Supplementary tables' download links
     - Historical tables (May 2002 - Jan 2015)
     - Current tables (Jan 2015 onwards)
  3. Save the downloaded .xlsx files to:
     ", raw_dir, "
  4. Also check the dashboard for any additional export options:
     https://www.aihw.gov.au/reports/medicines/pbs-monthly-data/contents/dashboard

Key data elements needed:
  - Monthly prescription counts by ATC Level 2 therapeutic class
  - Breakdown by patient category (General vs Concessional)
  - Geographic breakdown by LGA (or state/territory at minimum)
  - Prescription rates per capita (if available)
  - Time period: as far back as available through most recent month
")

# Check if any PBS data files exist
pbs_files <- list.files(raw_dir, pattern = "(?i)pbs|pharm|prescri",
                        full.names = TRUE)
if (length(pbs_files) > 0) {
  message("Found ", length(pbs_files), " PBS-related file(s) in data/raw/:")
  for (f in pbs_files) {
    message("  ", basename(f), " (", round(file.size(f) / 1e6, 1), " MB)")
  }
} else {
  message("No PBS data files found yet in data/raw/")
  message("Please download the files as described above before proceeding.")
}

# ==============================================================================
# 2. ABS SEIFA 2021 — IRSD at LGA Level
# ==============================================================================
message("\n=== 2. ABS SEIFA 2021 (IRSD by LGA) ===")

seifa_url <- paste0(
  "https://www.abs.gov.au/statistics/people/people-and-communities/",
  "socio-economic-indexes-areas-seifa-australia/2021/",
  "Local%20Government%20Area%2C%20Indexes%2C%20SEIFA%202021.xlsx"
)
seifa_file <- file.path(raw_dir, "SEIFA_2021_LGA.xlsx")

safe_download(seifa_url, seifa_file, "ABS SEIFA 2021 LGA Indexes")

# If automatic download fails, try alternative URL pattern
if (!file.exists(seifa_file)) {
  message("  Trying alternative URL...")
  seifa_url_alt <- paste0(
    "https://www.abs.gov.au/statistics/people/people-and-communities/",
    "socio-economic-indexes-areas-seifa-australia/latest-release/",
    "Local%20Government%20Area%2C%20Indexes%2C%20SEIFA%202021.xlsx"
  )
  safe_download(seifa_url_alt, seifa_file, "ABS SEIFA 2021 LGA (alt URL)")
}

if (!file.exists(seifa_file)) {
  message("
  MANUAL DOWNLOAD FALLBACK:
    1. Go to: https://www.abs.gov.au/statistics/people/people-and-communities/socio-economic-indexes-areas-seifa-australia/latest-release
    2. Under 'Data downloads', find 'Local Government Area, Indexes, SEIFA 2021'
    3. Download and save to: ", seifa_file)
}

# ==============================================================================
# 3. ABS REMOTENESS AREAS 2021 (SA1-to-RA Allocation)
# ==============================================================================
message("\n=== 3. ABS Remoteness Areas 2021 (SA1-to-RA Allocation) ===")

ra_url <- paste0(
  "https://www.abs.gov.au/statistics/standards/",
  "australian-statistical-geography-standard-asgs-edition-3/",
  "jul2021-jun2026/access-and-downloads/allocation-files/",
  "RA_2021_AUST.xlsx"
)
ra_file <- file.path(raw_dir, "RA_2021_AUST.xlsx")

safe_download(ra_url, ra_file, "ABS Remoteness Areas 2021 allocation")

if (!file.exists(ra_file)) {
  message("
  MANUAL DOWNLOAD FALLBACK:
    1. Go to: https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/allocation-files
    2. Download 'Remoteness Areas 2021' allocation file (RA_2021_AUST.xlsx)
    3. Save to: ", ra_file)
}

# Also download MB-to-LGA allocation for building the LGA-RA crosswalk
message("\n  Downloading MB-to-LGA allocation...")

mb_lga_url <- paste0(
  "https://www.abs.gov.au/statistics/standards/",
  "australian-statistical-geography-standard-asgs-edition-3/",
  "jul2021-jun2026/access-and-downloads/allocation-files/",
  "LGA_2021_AUST.xlsx"
)
mb_lga_file <- file.path(raw_dir, "LGA_2021_AUST_allocation.xlsx")

safe_download(mb_lga_url, mb_lga_file, "ABS MB-to-LGA 2021 allocation")

# Download MB main structure allocation (MB→SA1 mapping, needed for RA crosswalk)
message("\n  Downloading MB main structure allocation (for MB→SA1 mapping)...")

mb_sa1_url <- paste0(
  "https://www.abs.gov.au/statistics/standards/",
  "australian-statistical-geography-standard-asgs-edition-3/",
  "jul2021-jun2026/access-and-downloads/allocation-files/",
  "MB_2021_AUST.xlsx"
)
mb_sa1_file <- file.path(raw_dir, "MB_2021_AUST.xlsx")

safe_download(mb_sa1_url, mb_sa1_file, "ABS MB 2021 main structure allocation")

if (!file.exists(sa1_lga_file)) {
  message("
  MANUAL DOWNLOAD FALLBACK:
    1. Go to: https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/allocation-files
    2. Download 'Local Government Areas 2021' allocation file
    3. Save to: ", sa1_lga_file)
}

# ==============================================================================
# 4. ABS LGA BOUNDARY SHAPEFILES (2021 ASGS Edition 3)
# ==============================================================================
message("\n=== 4. ABS LGA Boundary Shapefiles 2021 ===")

lga_shp_url <- paste0(
  "https://www.abs.gov.au/statistics/standards/",
  "australian-statistical-geography-standard-asgs-edition-3/",
  "jul2021-jun2026/access-and-downloads/digital-boundary-files/",
  "LGA_2021_AUST_GDA2020_SHP.zip"
)
lga_shp_zip <- file.path(spatial_dir, "LGA_2021_AUST_GDA2020_SHP.zip")

safe_download(lga_shp_url, lga_shp_zip, "ABS LGA 2021 shapefiles")

# Extract shapefiles
if (file.exists(lga_shp_zip)) {
  lga_shp_dir <- file.path(spatial_dir, "LGA_2021")
  if (!dir.exists(lga_shp_dir)) {
    message("  Extracting shapefiles...")
    dir.create(lga_shp_dir, showWarnings = FALSE)
    unzip(lga_shp_zip, exdir = lga_shp_dir)
    message("  [OK] Extracted to ", lga_shp_dir)
  } else {
    message("  [SKIP] Already extracted")
  }
} else {
  message("
  MANUAL DOWNLOAD FALLBACK:
    1. Go to: https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files
    2. Download 'LGA_2021_AUST_GDA2020_SHP.zip'
    3. Save to: ", lga_shp_zip, "
    4. Extract to: ", file.path(spatial_dir, "LGA_2021"))
}

# ==============================================================================
# 5. ABS ESTIMATED RESIDENT POPULATION BY LGA
# ==============================================================================
message("\n=== 5. ABS Estimated Resident Population by LGA ===")

# ERP by LGA — ABS Regional Population publication
# 32180DS0004 = ERP by LGA, SUA, RA, CED, SED (2001-2024). Table 1 has LGA ERP.
erp_url <- paste0(
  "https://www.abs.gov.au/statistics/people/population/regional-population/",
  "2023-24/32180DS0004_2001-24.xlsx"
)
erp_file <- file.path(raw_dir, "ERP_LGA_2001_2024.xlsx")

safe_download(erp_url, erp_file, "ABS ERP by LGA 2001-2024")

if (!file.exists(erp_file)) {
  message("
  MANUAL DOWNLOAD FALLBACK:
    1. Go to: https://www.abs.gov.au/statistics/people/population/regional-population/latest-release
    2. Under 'Data downloads', find 'Population estimates by LGA' or similar
    3. Download the Excel data cube with LGA-level ERP from 2001 to latest
    4. Save to: ", erp_file)
}

# ==============================================================================
# SUMMARY
# ==============================================================================
message("\n", paste(rep("=", 60), collapse = ""))
message("DATA ACQUISITION SUMMARY")
message(paste(rep("=", 60), collapse = ""))

check_file <- function(path, label) {
  if (file.exists(path)) {
    sz <- round(file.size(path) / 1e6, 1)
    message("  [OK]   ", label, " (", sz, " MB)")
  } else {
    message("  [MISS] ", label, " -- manual download needed")
  }
}

check_file(seifa_file, "SEIFA 2021 LGA")
check_file(ra_file, "Remoteness Areas allocation")
check_file(mb_lga_file, "MB-to-LGA allocation")
check_file(mb_sa1_file, "MB main structure (MB→SA1)")
check_file(lga_shp_zip, "LGA shapefiles")
check_file(erp_file, "ERP by LGA")

pbs_files <- list.files(raw_dir, pattern = "(?i)pbs|pharm|prescri")
if (length(pbs_files) > 0) {
  message("  [OK]   AIHW PBS data (", length(pbs_files), " files)")
} else {
  message("  [MISS] AIHW PBS data -- manual download required (see above)")
}

message("\nNext step: Run 02_data_cleaning.py to wrangle data into tidy panels.")
