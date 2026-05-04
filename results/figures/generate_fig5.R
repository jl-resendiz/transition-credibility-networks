#!/usr/bin/env Rscript
# Figure 5: Cumulative gamma_fuel(tau) over the daily event-time window.
#
# Visualizes the four phases identified in §5 daily event-study:
#   (i)   Pre-event drift (taus -21 to -2)
#   (ii)  Announcement reaction (-1 to +1)
#   (iii) Post-announcement recovery (+2 to +10)
#   (iv)  Long-horizon resumption (+11 to +21)
#
# Input:  results/metrics/daily_event_time_path.md (parsed via lookup table)
# Output: results/figures/fig5_daily_path.pdf
#         manuscript/figures/fig5_daily_path.pdf

library(ggplot2)
library(sysfonts)
library(showtext)

font_add_google("Roboto", "roboto")
font_add_google("Roboto Condensed", "roboto_cond")
showtext_auto()

# Read the markdown table from daily_event_time_path.md
# (expects a markdown table; we parse manually)
path <- "results/metrics/daily_event_time_path.md"
if (!file.exists(path)) stop(paste("Missing:", path))

lines <- readLines(path)
# Pull table rows: lines starting with "|" and containing a leading numeric tau
data_rows <- grep("^\\| *[+-]?[0-9]+ *\\|", lines, value = TRUE)

# Parse: | tau | T | gamma_fuel | t_fuel | gamma_geo | t_geo |
parse_row <- function(line) {
  parts <- strsplit(gsub("\\|", "", line), ",|\\s+", perl = TRUE)[[1]]
  parts <- parts[parts != ""]
  list(
    tau = as.integer(parts[1]),
    T   = as.integer(parts[2]),
    g_fuel = as.numeric(gsub("\\+", "", parts[3])),
    t_fuel = as.numeric(gsub("\\+", "", parts[4])),
    g_geo  = as.numeric(gsub("\\+", "", parts[5])),
    t_geo  = as.numeric(gsub("\\+", "", parts[6]))
  )
}

df <- do.call(rbind.data.frame, lapply(data_rows, parse_row))
df <- df[!is.na(df$tau) & !is.na(df$g_fuel), ]
df <- df[order(df$tau), ]

# Compute cumulative gamma_fuel
df$cum_g_fuel <- cumsum(df$g_fuel)
df$cum_g_geo  <- cumsum(df$g_geo)

# Standard-error envelope (95% CI from individual-day NW SE; conservatively
# treat days as independent for the cumulative — slightly understates SE
# because of within-event correlation). For visualization only.
df$se_fuel <- abs(df$g_fuel / df$t_fuel)  # back out SE from t-stat
df$cum_se_fuel <- sqrt(cumsum(df$se_fuel^2))
df$cum_lo <- df$cum_g_fuel - 1.96 * df$cum_se_fuel
df$cum_hi <- df$cum_g_fuel + 1.96 * df$cum_se_fuel

# Phase shading
phases <- data.frame(
  xmin = c(-21.5, -1.5, +1.5, +10.5),
  xmax = c( -1.5, +1.5, +10.5, +21.5),
  phase = c("Pre-event drift", "Announcement", "Recovery", "Long-horizon resumption"),
  fill = c("#e8e8e8", "#d4a76a", "#cce8d4", "#e8d8d8")
)

p <- ggplot(df, aes(x = tau, y = cum_g_fuel)) +
  geom_rect(data = phases,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            inherit.aes = FALSE, alpha = 0.4) +
  scale_fill_identity() +
  geom_ribbon(aes(ymin = cum_lo, ymax = cum_hi), fill = "gray75", alpha = 0.4) +
  geom_hline(yintercept = 0, color = "gray60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.3, linetype = "dashed") +
  geom_line(color = "gray20", linewidth = 0.6) +
  geom_point(color = "gray20", size = 1.2) +
  scale_x_continuous(breaks = seq(-21, 21, 7),
                     limits = c(-21.5, 21.5)) +
  labs(x = "Trading-day offset from announcement (τ)",
       y = "Cumulative γ_fuel(τ)",
       title = NULL) +
  annotate("text", x = -11, y = max(df$cum_hi) * 0.95,
           label = "Pre-event drift", size = 2.8,
           family = "roboto_cond", color = "gray30") +
  annotate("text", x = 0, y = max(df$cum_hi) * 0.95,
           label = "Announce.", size = 2.8,
           family = "roboto_cond", color = "gray30") +
  annotate("text", x = 6, y = max(df$cum_hi) * 0.95,
           label = "Recovery", size = 2.8,
           family = "roboto_cond", color = "gray30") +
  annotate("text", x = 16, y = max(df$cum_hi) * 0.95,
           label = "Long horizon", size = 2.8,
           family = "roboto_cond", color = "gray30") +
  theme_minimal(base_family = "roboto", base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray93", linewidth = 0.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 9, color = "gray20"),
    axis.text  = element_text(size = 8, color = "gray40"),
    plot.margin = margin(8, 14, 8, 8)
  )

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("results/figures/fig5_daily_path.pdf", p,
       width = 5.5, height = 3.2, device = cairo_pdf)
ggsave("manuscript/figures/fig5_daily_path.pdf", p,
       width = 5.5, height = 3.2, device = cairo_pdf)
cat("Saved fig5_daily_path.pdf to results/figures/ and manuscript/figures/\n")
