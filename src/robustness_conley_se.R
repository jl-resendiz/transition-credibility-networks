# ──────────────────────────────────────────────────────────────────────
# Conley (1999) Standard Errors for Channel Decomposition
#
# Replicates Table 2 of the paper using fixest::conley() for
# spatially-robust inference at 250, 500, and 1000 km cutoffs.
# ──────────────────────────────────────────────────────────────────────
library(fixest)
library(data.table)

# Resolve root: works with Rscript and source()
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("--file=", .args, value = TRUE)
if (length(.file_arg) > 0) {
  .script_dir <- dirname(normalizePath(sub("--file=", "", .file_arg)))
} else {
  .script_dir <- getwd()
}
root <- normalizePath(file.path(.script_dir, ".."))
derived <- file.path(root, "data", "derived")
raw     <- file.path(root, "data", "raw")
results <- file.path(root, "results", "metrics")

cat("Loading data...\n")

# Monthly returns
ret <- fread(file.path(derived, "returns", "monthly_returns.csv"))
ret[, ym := substr(datadate, 1, 7)]

# Fama-French monthly factors → vwretd
ff_lines <- readLines(file.path(raw, "factors", "F-F_Research_Data_Factors.csv"))
ff_rows <- list()
for (line in ff_lines) {
  parts <- trimws(strsplit(line, ",")[[1]])
  if (length(parts) >= 5 && nchar(parts[1]) == 6 && grepl("^[0-9]+$", parts[1])) {
    ym <- paste0(substr(parts[1], 1, 4), "-", substr(parts[1], 5, 6))
    mktrf <- as.numeric(parts[2])
    rf    <- as.numeric(parts[5])
    if (!is.na(mktrf) && !is.na(rf)) {
      ff_rows[[length(ff_rows) + 1]] <- data.table(ym = ym, vwretd = (mktrf + rf) / 100)
    }
  }
}
ff <- rbindlist(ff_rows)

# Merge to get AR
ret <- merge(ret, ff, by = "ym", all.x = TRUE)
ret[, ar := ret_monthly - vwretd]

# Events (first-mover retirements)
events <- fread(file.path(derived, "events", "coal_retirement_events.csv"))
fm_col <- intersect(c("first_mover", "is_first_mover"), names(events))
if (length(fm_col) > 0) events <- events[tolower(as.character(get(fm_col[1]))) %in% c("1", "true")]
# Require matched gvkeys
events <- events[!is.na(matched_gvkeys) & nchar(as.character(matched_gvkeys)) > 0]
# Build event_ym: prefer announcement_date > event_date > ret_year-06
for (col in c("event_date", "announcement_date")) {
  if (col %in% names(events)) {
    if (is.integer(events[[col]]) || is.numeric(events[[col]])) {
      events[, (col) := as.character(as.Date(get(col), origin = "1970-01-01"))]
    } else {
      events[, (col) := as.character(get(col))]
    }
  }
}
events[, event_ym := NA_character_]
if ("announcement_date" %in% names(events)) {
  events[nchar(announcement_date) >= 7, event_ym := substr(announcement_date, 1, 7)]
}
if ("event_date" %in% names(events)) {
  events[is.na(event_ym) & nchar(event_date) >= 7, event_ym := substr(event_date, 1, 7)]
}
if ("ret_year" %in% names(events)) {
  events[is.na(event_ym) & !is.na(ret_year), event_ym := paste0(ret_year, "-06")]
}
events <- events[!is.na(event_ym)]
# Create unique event_id
events[, event_id := paste0("ret_", seq_len(.N))]
# Parse matched gvkeys (semicolon-separated)
events[, gvkey := as.character(matched_gvkeys)]
cat(sprintf("  Events: %d first-mover retirements with dates\n", nrow(events)))

# Weight matrices
load_W <- function(path) {
  w <- fread(path)
  setnames(w, c("gvkey_i", "gvkey_j", "w_ij"))
  w[, w_ij := as.numeric(w_ij)]
  w
}
W_geo  <- load_W(file.path(derived, "networks", "weight_matrix_W_geo.csv"))
W_fuel <- load_W(file.path(derived, "networks", "weight_matrix_W_fuel.csv"))
W_reg  <- load_W(file.path(derived, "networks", "weight_matrix_W_regulatory.csv"))

# Firm centroids (lat/lon for Conley)
centroids <- fread(file.path(derived, "networks", "firm_centroids.csv"))
setnames(centroids, tolower(names(centroids)))
# Normalize column names
if ("centroid_lat" %in% names(centroids)) setnames(centroids, "centroid_lat", "lat")
if ("centroid_lon" %in% names(centroids)) setnames(centroids, "centroid_lon", "lon")

# Alpha panel
alpha <- fread(file.path(derived, "fundamentals", "firm_alpha_panel.csv"))

# Fundamentals (for SIC sector matching)
fundm <- fread(file.path(derived, "fundamentals", "firm_fundamentals.csv"))
# Get most recent SIC per gvkey
sic_col <- intersect(c("sic", "sich", "sic4"), names(fundm))
if (length(sic_col) > 0) {
  sic_map <- fundm[!is.na(get(sic_col[1])), .(sic = get(sic_col[1])[.N]), by = gvkey]
} else {
  sic_map <- data.table(gvkey = character(), sic = character())
}

# ── Build event-firm panel ───────────────────────────────────────────

add_months <- function(ym, delta) {
  y <- as.integer(substr(ym, 1, 4))
  m <- as.integer(substr(ym, 6, 7)) + delta
  while (m > 12) { y <- y + 1; m <- m - 12 }
  while (m < 1)  { y <- y - 1; m <- m + 12 }
  sprintf("%04d-%02d", y, m)
}

cat("Building event-firm pairs for [-1, +3] window...\n")
all_gvkeys <- unique(ret$gvkey)
obs_list <- list()

for (i in seq_len(nrow(events))) {
  ev <- events[i]
  ev_ym    <- ev$event_ym
  ev_id    <- ev$event_id
  # Event firm gvkeys (semicolon-separated matched_gvkeys)
  ev_gvkeys <- trimws(unlist(strsplit(as.character(ev$matched_gvkeys), ";")))
  ev_gvkeys <- ev_gvkeys[nchar(ev_gvkeys) > 0]
  if (length(ev_gvkeys) == 0) next
  ev_gvkey <- ev_gvkeys[1]  # primary event firm for W lookup

  # CAR window: [-1, +3] months
  months_needed <- sapply(-1:3, function(d) add_months(ev_ym, d))

  for (gk in all_gvkeys) {
    if (gk %in% ev_gvkeys) next
    # Get ARs in window
    firm_ret <- ret[gvkey == gk & ym %in% months_needed]
    if (nrow(firm_ret) < 3) next
    car <- sum(firm_ret$ar, na.rm = TRUE)

    # Spatial weights
    w_g <- W_geo[gvkey_i == ev_gvkey & gvkey_j == gk, w_ij]
    w_f <- W_fuel[gvkey_i == ev_gvkey & gvkey_j == gk, w_ij]
    w_r <- W_reg[gvkey_i == ev_gvkey & gvkey_j == gk, w_ij]
    if (length(w_g) == 0) w_g <- 0
    if (length(w_f) == 0) w_f <- 0
    if (length(w_r) == 0) w_r <- 0

    # Same sector
    ev_sic <- sic_map[gvkey == ev_gvkey, sic]
    gk_sic <- sic_map[gvkey == gk, sic]
    same_sector <- 0
    if (length(ev_sic) > 0 && length(gk_sic) > 0) {
      if (!is.na(ev_sic[1]) && !is.na(gk_sic[1]) &&
          nchar(ev_sic[1]) >= 2 && nchar(gk_sic[1]) >= 2) {
        same_sector <- as.integer(substr(ev_sic[1], 1, 2) == substr(gk_sic[1], 1, 2))
      }
    }

    # Centroid for this firm
    clat <- centroids[gvkey == gk, lat]
    clon <- centroids[gvkey == gk, lon]
    if (length(clat) == 0) { clat <- NA; clon <- NA }

    obs_list[[length(obs_list) + 1]] <- data.table(
      event_id = ev_id, gvkey = gk, car = car,
      w_geo = w_g[1], w_fuel = w_f[1], w_reg = w_r[1],
      same_sector = same_sector,
      lat = clat[1], lon = clon[1]
    )
  }
  if (i %% 20 == 0) cat(sprintf("  Processed %d/%d events\n", i, nrow(events)))
}

panel <- rbindlist(obs_list)
cat(sprintf("Panel: %d event-firm pairs\n", nrow(panel)))

# Drop rows with missing coordinates for Conley
panel_geo <- panel[!is.na(lat) & !is.na(lon)]
cat(sprintf("Panel with coordinates: %d rows\n", nrow(panel_geo)))

# ── Estimate ─────────────────────────────────────────────────────────

cat("\nEstimating channel decomposition...\n\n")

est <- feols(car ~ w_geo + w_fuel + w_reg + same_sector,
             data = panel_geo, lean = FALSE)

# Event-clustered (baseline)
se_event <- summary(est, cluster = ~event_id)

# Conley SEs at three cutoffs
cat("Computing Conley SEs (this may take a moment)...\n")
se_250  <- summary(est, vcov = conley(cutoff = 250,  distance = "spherical") ~ lat + lon)
se_500  <- summary(est, vcov = conley(cutoff = 500,  distance = "spherical") ~ lat + lon)
se_1000 <- summary(est, vcov = conley(cutoff = 1000, distance = "spherical") ~ lat + lon)

# ── Format results ───────────────────────────────────────────────────

vars <- c("w_geo", "w_fuel", "w_reg", "same_sector")

fmt <- function(sm, var) {
  ct <- coeftable(sm)
  idx <- which(rownames(ct) == var)
  if (length(idx) == 0) return(c(NA, NA, NA))
  c(ct[idx, 1], ct[idx, 2], ct[idx, 3])
}

lines <- c(
  "# Conley Standard Errors: Channel Decomposition",
  "",
  sprintf("N = %d event-firm pairs", nrow(panel_geo)),
  "",
  "| Variable | Coeff | Event-clustered SE (t) | Conley 250km SE (t) | Conley 500km SE (t) | Conley 1000km SE (t) |",
  "|---|---|---|---|---|---|"
)

for (v in vars) {
  ev <- fmt(se_event, v)
  c2 <- fmt(se_250, v)
  c5 <- fmt(se_500, v)
  c10 <- fmt(se_1000, v)
  lines <- c(lines, sprintf("| %s | %.3f | %.3f (%.2f) | %.3f (%.2f) | %.3f (%.2f) | %.3f (%.2f) |",
    v, ev[1],
    ev[2], ev[3],
    c2[2], c2[3],
    c5[2], c5[3],
    c10[2], c10[3]))
}

lines <- c(lines, "",
  "Notes: Conley (1999) SEs use Bartlett kernel with spherical (Haversine) distances.",
  "Firm coordinates are capacity-weighted plant centroids from GEM/GeoAsset.",
  "Dependent variable: cumulative abnormal return over [-1, +3] months.",
  "Market adjustment: vwretd from Fama-French factors.")

dir.create(results, recursive = TRUE, showWarnings = FALSE)
writeLines(lines, file.path(results, "conley_se.md"))

cat("\n=== Results ===\n")
cat(paste(lines, collapse = "\n"), "\n")
cat(sprintf("\nWritten to: %s\n", file.path(results, "conley_se.md")))
