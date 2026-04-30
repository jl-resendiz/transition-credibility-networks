library(ggplot2)
library(sysfonts)
library(showtext)

font_add_google("Roboto", "roboto")
font_add_google("Roboto Condensed", "roboto_cond")
showtext_auto()

# Load event-time path data
df <- read.csv("results/summaries/event_time_betas.csv")
df$ci_lo <- df$beta_fuel - 1.96 * df$se_fuel
df$ci_hi <- df$beta_fuel + 1.96 * df$se_fuel
df$period <- ifelse(df$tau <= -2, "Pre",
              ifelse(df$tau <= 3, "Headline window", "Post"))
df$period <- factor(df$period,
                    levels = c("Pre", "Headline window", "Post"))

# Headline window shaded region
hl_band <- data.frame(xmin = -1.5, xmax = 3.5)

p <- ggplot(df, aes(x = tau, y = beta_fuel)) +
  geom_rect(data = hl_band, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax,
                ymin = -Inf, ymax = Inf),
            fill = "gray92", alpha = 0.5) +
  geom_hline(yintercept = 0, color = "gray60",
             linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, color = period),
                width = 0.25, linewidth = 0.4) +
  geom_point(aes(color = period, fill = period),
             shape = 21, size = 2, stroke = 0.4) +
  scale_color_manual(values = c("Pre" = "gray45",
                                "Headline window" = "black",
                                "Post" = "gray45")) +
  scale_fill_manual(values = c("Pre" = "gray85",
                               "Headline window" = "gray35",
                               "Post" = "gray85")) +
  scale_x_continuous(breaks = seq(-12, 6, 2),
                     limits = c(-12.5, 6.5)) +
  labs(x = "Event time (months)",
       y = expression("Fama-MacBeth " * beta[fuel](tau)),
       title = NULL) +
  theme_minimal(base_family = "roboto", base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "gray93",
                                    linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 8.5, color = "gray20"),
    axis.text = element_text(size = 7.5, color = "gray40"),
    plot.margin = margin(8, 14, 8, 8)
  )

ggsave("results/figures/fig4_event_time_path.pdf", p,
       width = 5.5, height = 3.2, device = cairo_pdf)
ggsave("manuscript/figures/fig4_event_time_path.pdf", p,
       width = 5.5, height = 3.2, device = cairo_pdf)
cat("Saved fig4_event_time_path.pdf to results/figures/ and manuscript/figures/\n")
