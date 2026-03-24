library(ggplot2)
library(sysfonts)
library(showtext)

font_add_google("Roboto", "roboto")
font_add_google("Roboto Condensed", "roboto_cond")
showtext_auto()

df <- data.frame(
  tercile = factor(
    c("2011\u20132013", "2013\u20132014", "2014\u20132022"),
    levels = rev(c("2011\u20132013", "2013\u20132014", "2014\u20132022"))
  ),
  beta = c(1.717, -3.578, -5.518),
  se   = c(2.292, 1.185, 1.541),
  n    = c(58, 58, 59),
  sig  = c(FALSE, TRUE, TRUE)
)
df$ci_lo <- df$beta - 1.96 * df$se
df$ci_hi <- df$beta + 1.96 * df$se

p <- ggplot(df, aes(x = beta, y = tercile)) +
  geom_vline(xintercept = 0, color = "gray75", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.12,
                 color = "gray50", linewidth = 0.35) +
  geom_point(aes(fill = sig), shape = 21, size = 3.2, stroke = 0.2, color = "white") +
  scale_fill_manual(values = c("TRUE" = "gray20", "FALSE" = "gray70"), guide = "none") +
  geom_text(aes(label = paste0(n, " events")),
            hjust = -0.25, size = 2.6, nudge_y = 0.18,
            family = "roboto_cond", color = "gray40") +
  labs(x = "Mean fuel coefficient (Fama-MacBeth)", y = NULL) +
  scale_x_continuous(breaks = seq(-10, 6, 2)) +
  theme_minimal(base_family = "roboto", base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray93", linewidth = 0.2),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 8, color = "gray20"),
    axis.text = element_text(size = 7.5, color = "gray40"),
    plot.margin = margin(8, 18, 8, 8)
  )

ggsave("results/figures/fig2_calendar_time.pdf", p, width = 3.5, height = 2.2,
       device = cairo_pdf)
cat("Saved: results/figures/fig2_calendar_time.pdf\n")
