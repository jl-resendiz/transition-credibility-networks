library(ggplot2)
library(sysfonts)
library(showtext)

font_add_google("Roboto", "roboto")
font_add_google("Roboto Condensed", "roboto_cond")
showtext_auto()

df <- read.csv("results/summaries/event_level_betas.csv")
df$below <- df$beta_fuel < df$beta_geo
n_below <- sum(df$below)
n_total <- nrow(df)

# Trim outlier for cleaner axes; note in caption
xlim <- c(-15, 12)
ylim <- c(-22, 12)

p <- ggplot(df, aes(x = beta_geo, y = beta_fuel)) +
  geom_hline(yintercept = 0, color = "gray85", linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "gray85", linewidth = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "21", color = "gray60", linewidth = 0.35) +
  geom_point(aes(fill = below), shape = 21, size = 1.6, stroke = 0.15,
             color = "white", alpha = 0.85) +
  scale_fill_manual(values = c("TRUE" = "gray25", "FALSE" = "gray72"), guide = "none") +
  # Annotations positioned safely inside the plot area
  annotate("text", x = -4, y = -20.5,
           label = paste0(n_below, " of ", n_total, " events (", round(100*n_below/n_total), "%)"),
           size = 2.8, family = "roboto_cond", color = "gray25", fontface = "italic") +
  annotate("text", x = 8, y = 10,
           label = paste0(n_total - n_below, " of ", n_total),
           size = 2.6, family = "roboto_cond", color = "gray55", fontface = "italic") +
  annotate("text", x = 9.5, y = 7.5, label = "fuel = geo",
           size = 2.1, family = "roboto_cond", color = "gray50", angle = 38) +
  scale_x_continuous(limits = xlim, breaks = seq(-15, 10, 5)) +
  scale_y_continuous(limits = ylim, breaks = seq(-20, 10, 5)) +
  labs(x = "Geographic coefficient (event-level)",
       y = "Fuel coefficient (event-level)") +
  theme_minimal(base_family = "roboto", base_size = 9) +
  theme(
    panel.grid.major = element_line(color = "gray93", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 8.5, color = "gray20"),
    axis.text = element_text(size = 7.5, color = "gray40"),
    plot.margin = margin(8, 14, 8, 8)
  )

ggsave("results/figures/fig1_fuel_vs_geo.pdf", p, width = 3.5, height = 3.5,
       device = cairo_pdf)
cat("Saved: results/figures/fig1_fuel_vs_geo.pdf\n")
