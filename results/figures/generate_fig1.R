library(ggplot2)
library(tidyr)
library(sysfonts)
library(showtext)

font_add_google("Roboto", "roboto")
font_add_google("Roboto Condensed", "roboto_cond")
showtext_auto()

df <- read.csv("results/summaries/event_level_betas.csv")

# Reshape to long format
long <- pivot_longer(df, cols = c(beta_fuel, beta_geo),
                     names_to = "channel", values_to = "coefficient")
long$channel <- factor(long$channel,
                       levels = c("beta_fuel", "beta_geo"),
                       labels = c("Fuel-mix similarity",
                                  "Geographic proximity"))

# Summary statistics
mu_fuel <- mean(df$beta_fuel)
mu_geo  <- mean(df$beta_geo)

means <- data.frame(
  channel = factor(c("Fuel-mix similarity", "Geographic proximity"),
                   levels = c("Fuel-mix similarity",
                              "Geographic proximity")),
  mu = c(mu_fuel, mu_geo)
)

p <- ggplot(long, aes(x = coefficient)) +
  geom_histogram(binwidth = 2, fill = "gray35", color = "white",
                 linewidth = 0.15, alpha = 0.8) +
  geom_vline(xintercept = 0, color = "gray75", linewidth = 0.3) +
  geom_vline(data = means, aes(xintercept = mu),
             linetype = "dashed", color = "gray25", linewidth = 0.4) +
  geom_text(data = means,
            aes(x = mu - 0.8, y = Inf,
                label = paste0("mean = ", round(mu, 1))),
            vjust = 2, hjust = 1, size = 2.5, family = "roboto_cond",
            color = "gray25", fontface = "italic") +
  facet_wrap(~ channel, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = seq(-20, 10, 5),
                     limits = c(-22, 8)) +
  labs(x = "Event-level coefficient", y = "Count") +
  theme_minimal(base_family = "roboto", base_size = 9) +
  theme(
    panel.grid.major = element_line(color = "gray93",
                                    linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_text(size = 8.5, color = "gray20"),
    axis.text = element_text(size = 7.5, color = "gray40"),
    strip.text = element_text(size = 8.5, family = "roboto_cond",
                              color = "gray20", hjust = 0),
    plot.margin = margin(8, 14, 8, 8)
  )

ggsave("results/figures/fig1_fuel_vs_geo.pdf", p,
       width = 3.5, height = 3.5, device = cairo_pdf)
ggsave("manuscript/figures/fig1_fuel_vs_geo.pdf", p,
       width = 3.5, height = 3.5, device = cairo_pdf)
cat("Saved to results/figures/ and manuscript/figures/\n")
