# ──────────────────────────────────────────────────────────────────────
# Wild Cluster Bootstrap for Phase-Out DiD
#
# With only G=14 Tier-1 binding events, standard cluster-robust SEs
# are unreliable. This script implements:
#   1. CR2 bias-reduced SEs via clubSandwich (Bell-McCaffrey)
#   2. Manual wild cluster bootstrap with Rademacher weights
#   3. fixest two-way clustering as baseline comparison
#
# References:
#   Cameron, Gelbach & Miller (2008) J. Business & Economic Statistics
#   Webb (2023) Econometrics Journal
#   Pustejovsky & Tipton (2018) JBES [clubSandwich]
# ──────────────────────────────────────────────────────────────────────
library(fixest)
library(clubSandwich)
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

# Fama-French monthly → vwretd
ff_lines <- readLines(file.path(raw, "factors", "F-F_Research_Data_Factors.csv"))
ff_rows <- list()
for (line in ff_lines) {
  parts <- trimws(strsplit(line, ",")[[1]])
  if (length(parts) >= 5 && nchar(parts[1]) == 6 && grepl("^[0-9]+$", parts[1])) {
    ym <- paste0(substr(parts[1], 1, 4), "-", substr(parts[1], 5, 6))
    mktrf <- as.numeric(parts[2]); rf <- as.numeric(parts[5])
    if (!is.na(mktrf) && !is.na(rf))
      ff_rows[[length(ff_rows) + 1]] <- data.table(ym = ym, vwretd = (mktrf + rf) / 100)
  }
}
ff <- rbindlist(ff_rows)
ret <- merge(ret, ff, by = "ym", all.x = TRUE)
ret[, ar := ret_monthly - vwretd]

# Phase-out events (Tier-1 binding only)
events <- fread(file.path(derived, "events", "coal_phaseout_shocks_events.csv"))
# Detect tier column name
tier_col <- intersect(c("tier", "exogeneity_tier"), names(events))
if (length(tier_col) > 0) events <- events[get(tier_col[1]) == 1 | get(tier_col[1]) == "1"]
# Detect binding column
if ("binding" %in% names(events)) events <- events[tolower(as.character(binding)) %in% c("1", "true", "yes")]
# Convert event_date: fread may read dates as integer (days since epoch)
if (is.integer(events$event_date) || is.numeric(events$event_date)) {
  events[, event_date := as.character(as.Date(event_date, origin = "1970-01-01"))]
} else {
  events[, event_date := as.character(event_date)]
}
events <- events[!is.na(event_date) & nchar(event_date) >= 7]
events[, event_ym := substr(event_date, 1, 7)]
cat(sprintf("  Tier-1 binding events: %d\n", nrow(events)))

# Alpha panel → coal_share for treatment intensity
alpha <- fread(file.path(derived, "fundamentals", "firm_alpha_panel.csv"))

# ── Build panel ──────────────────────────────────────────────────────

add_months <- function(ym, delta) {
  y <- as.integer(substr(ym, 1, 4))
  m <- as.integer(substr(ym, 6, 7)) + delta
  while (m > 12) { y <- y + 1; m <- m - 12 }
  while (m < 1)  { y <- y - 1; m <- m + 12 }
  sprintf("%04d-%02d", y, m)
}

cat("Building event-firm-month panel...\n")
all_gvkeys <- unique(ret$gvkey)
obs <- list()

for (i in seq_len(nrow(events))) {
  ev <- events[i]
  ev_ym <- ev$event_ym
  ev_id <- if ("event_id" %in% names(ev)) ev$event_id else paste0("po_", i)
  ev_country <- if ("country" %in% names(ev)) ev$country else ""

  months_window <- sapply(-6:12, function(d) add_months(ev_ym, d))

  for (gk in all_gvkeys) {
    firm_ret <- ret[gvkey == gk & ym %in% months_window]
    if (nrow(firm_ret) < 6) next

    # Coal share as treatment intensity
    cs <- alpha[gvkey == gk, alpha]
    if (length(cs) == 0) next
    coal_share <- cs[1]
    if (is.na(coal_share)) next

    for (j in seq_len(nrow(firm_ret))) {
      tau <- which(months_window == firm_ret$ym[j]) - 7  # -6 maps to tau=-6
      post <- as.integer(tau >= 0)
      obs[[length(obs) + 1]] <- data.table(
        event_id = ev_id, gvkey = gk,
        ym = firm_ret$ym[j], tau = tau,
        ar = firm_ret$ar[j],
        coal_share = coal_share,
        post = post,
        exp_post = coal_share * post
      )
    }
  }
  if (i %% 5 == 0) cat(sprintf("  Processed %d/%d events\n", i, nrow(events)))
}

panel <- rbindlist(obs)
cat(sprintf("Panel: %d obs, %d events, %d firms\n",
            nrow(panel), length(unique(panel$event_id)), length(unique(panel$gvkey))))

# ── Estimation ───────────────────────────────────────────────────────

cat("\nEstimating DiD: AR ~ coal_share x Post | firm + month FE\n\n")

# fixest with firm + month FE
est <- feols(ar ~ exp_post | gvkey + ym, data = panel)

# (a) Event-clustered
se_event <- summary(est, cluster = ~event_id)

# (b) Two-way clustered (event + firm)
se_twoway <- summary(est, cluster = ~event_id + gvkey)

# (c) CR2 bias-reduced (clubSandwich) — designed for few clusters
cat("Computing CR2 (Bell-McCaffrey) SEs via clubSandwich...\n")
# clubSandwich works with lm objects; extract the demeaned model
panel[, ar_dm := ar - mean(ar), by = gvkey]
panel[, ar_dm := ar_dm - mean(ar_dm), by = ym]
panel[, ep_dm := exp_post - mean(exp_post), by = gvkey]
panel[, ep_dm := ep_dm - mean(ep_dm), by = ym]
lm_dm <- lm(ar_dm ~ ep_dm - 1, data = panel)
cr2 <- coef_test(lm_dm, vcov = "CR2", cluster = panel$event_id, test = "Satterthwaite")

# (d) Manual wild cluster bootstrap (Rademacher weights)
cat("Running wild cluster bootstrap (B=999, Rademacher weights)...\n")

beta_obs <- coef(est)["exp_post"]
cluster_ids <- unique(panel$event_id)
G <- length(cluster_ids)

# Restricted model (H0: beta=0) → just FE
est_r <- feols(ar ~ 1 | gvkey + ym, data = panel)
fitted_r <- fitted(est_r)
resid_u  <- residuals(est)

# Assign cluster index
panel[, cl_idx := match(event_id, cluster_ids)]

B <- 999
t_obs <- coeftable(se_event)["exp_post", "t value"]
t_boot <- numeric(B)

set.seed(42)
for (b in seq_len(B)) {
  # Rademacher weights: +1 or -1 per cluster
  w <- sample(c(-1, 1), G, replace = TRUE)
  # Bootstrap y*
  panel[, y_star := fitted_r + w[cl_idx] * resid_u]
  # Re-estimate
  est_b <- tryCatch(
    feols(y_star ~ exp_post | gvkey + ym, data = panel, warn = FALSE, notes = FALSE),
    error = function(e) NULL
  )
  if (is.null(est_b)) { t_boot[b] <- 0; next }
  se_b <- tryCatch({
    sm <- summary(est_b, cluster = ~event_id)
    coeftable(sm)["exp_post", "t value"]
  }, error = function(e) 0)
  t_boot[b] <- se_b
  if (b %% 100 == 0) cat(sprintf("  Bootstrap: %d/%d\n", b, B))
}

# Symmetric p-value
p_boot <- (1 + sum(abs(t_boot) >= abs(t_obs))) / (1 + B)

# ── Report ───────────────────────────────────────────────────────────

ct_ev <- coeftable(se_event)
ct_tw <- coeftable(se_twoway)

lines <- c(
  "# Wild Cluster Bootstrap: Phase-Out DiD",
  "",
  sprintf("G = %d Tier-1 binding events, N = %d, B = %d bootstrap replications", G, nrow(panel), B),
  "",
  "## exp_post (coal_share x Post) coefficient",
  "",
  "| Inference Method | SE | t-stat | p-value | Clusters |",
  "|---|---|---|---|---|",
  sprintf("| Event-clustered (CR1) | %.4f | %.2f | %.3f | %d |",
    ct_ev["exp_post", 2], ct_ev["exp_post", 3], ct_ev["exp_post", 4], G),
  sprintf("| Two-way (event + firm) | %.4f | %.2f | %.3f | %d + %d |",
    ct_tw["exp_post", 2], ct_tw["exp_post", 3], ct_tw["exp_post", 4],
    G, length(unique(panel$gvkey))),
  sprintf("| CR2 Bell-McCaffrey | %.4f | %.2f | %.3f | %d |",
    cr2$SE, cr2$tstat, cr2$p_Satt, G),
  sprintf("| Wild bootstrap (Rademacher) | — | %.2f (obs) | %.3f | %d |",
    t_obs, p_boot, G),
  "",
  sprintf("Coefficient: %.4f", beta_obs),
  "",
  "Notes: CR2 uses Satterthwaite degrees of freedom (clubSandwich).",
  "Wild bootstrap uses Rademacher weights (+1/-1) under H0: beta=0.",
  "Firm + month FE absorbed via fixest. Restricted model imposes beta=0."
)

dir.create(results, recursive = TRUE, showWarnings = FALSE)
writeLines(lines, file.path(results, "strategy3_phaseout_wild_bootstrap.md"))

cat("\n=== Results ===\n")
cat(paste(lines, collapse = "\n"), "\n")
cat(sprintf("\nWritten to: %s\n", file.path(results, "strategy3_phaseout_wild_bootstrap.md")))
