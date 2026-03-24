# Figure Specifications

Three figures for the manuscript "When Coal Retires: The Propagation of Stranding Risk."
Design follows Schwabish (JEP 2014): show the data, reduce clutter, integrate text and graph.

---

## Figure 1: Event-Level Fuel vs. Geographic Coefficients

**Purpose:** Visualise the paper's central finding — fuel dominates geography in 82/117 events.

**Placement:** Section 4.1 (Fuel-Mix Similarity Dominates), after Table 1.

### Data source
- File: `results/summaries/event_level_betas.csv`
- Columns: `beta_fuel` (y-axis), `beta_geo` (x-axis)
- N = 117 events (rows in the CSV)

### Design

**Type:** Scatterplot with 45-degree reference line.

**Axes:**
- X-axis: β_geo (geographic coefficient), range approximately [-5, +10]
- Y-axis: β_fuel (fuel coefficient), range approximately [-25, +15]
- Both axes labeled with variable names and "(event-level Fama-MacBeth coefficient)"
- Zero lines in light gray for both axes

**Reference line:** 45-degree line (y = x) in medium gray, dashed. Points BELOW this line have β_fuel < β_geo (fuel more negative than geography). Label the line at the upper-right edge: "β_fuel = β_geo"

**Points:**
- All 117 events as circles, medium size
- Points below the 45-degree line: dark fill (these are the 82/117 = 70%)
- Points above the 45-degree line: light/hollow fill (these are the 35/117 = 30%)
- No individual point labels (too many events)

**Annotations:**
- Text annotation in the lower-left quadrant: "82 of 117 events (70%)" with an arrow or bracket pointing to the cluster below the line
- Text annotation in the upper region: "35 of 117 events (30%)"

**Title:** "Fuel-mix similarity produces more negative returns than geographic proximity in 70 percent of events"

**Notes:** "Each point represents one Fama-MacBeth cross-sectional regression. Points below the 45-degree line indicate events where the fuel coefficient is more negative than the geographic coefficient."

**Color:** Grayscale-safe. Dark gray fill for below-line points, light gray for above-line points. No color needed.

**Size:** Single column width (fits a journal column). Approximately 3.5 x 3.5 inches.

### Production (R code sketch)
```r
library(ggplot2)
df <- read.csv("results/summaries/event_level_betas.csv")
df$below_line <- df$beta_fuel < df$beta_geo

ggplot(df, aes(x = beta_geo, y = beta_fuel)) +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(fill = below_line), shape = 21, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("TRUE" = "gray30", "FALSE" = "gray80"), guide = "none") +
  annotate("text", x = -2, y = -20, label = "82 of 117 events (70%)",
           size = 3, hjust = 0, fontface = "italic") +
  annotate("text", x = 5, y = 10, label = "35 of 117 events (30%)",
           size = 3, hjust = 0, fontface = "italic") +
  labs(x = expression(hat(gamma)[geo] ~ "(event-level)"),
       y = expression(hat(gamma)[fuel] ~ "(event-level)")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

ggsave("results/figures/fig1_fuel_vs_geo.pdf", width = 3.5, height = 3.5)
```

---

## Figure 2: Calendar-Time Strengthening of the Fuel Coefficient

**Purpose:** Visualise that the fuel signal is absent early and strengthens monotonically.

**Placement:** Section 4.5 (The Fuel Signal Strengthens Over Time).

### Data (hardcoded from pipeline results)

| Tercile | Period | N events | Mean β_fuel | NW SE | t |
|---|---|---|---|---|---|
| T1 | 2011–2013 | 58 | +1.717 | 2.292 | 0.75 |
| T2 | 2013–2014 | 58 | -3.578 | 1.185 | -3.02 |
| T3 | 2014–2022 | 59 | -5.518 | 1.541 | -3.58 |

### Design

**Type:** Dot plot with 95% confidence intervals (horizontal bars). NOT a bar chart (Schwabish: dot plots use less ink and emphasise point estimates more clearly than bars).

**Axes:**
- X-axis: Mean fuel coefficient (β_fuel), range [-10, +6]
- Y-axis: Three categories, labeled "2011–2013 (T1)", "2013–2014 (T2)", "2014–2022 (T3)", ordered top to bottom (earliest at top)
- Vertical zero reference line in light gray

**Points:**
- Three dots (one per tercile), large filled circles
- 95% CI bars: beta ± 1.96 * NW_SE
  - T1: [+1.717 - 4.493, +1.717 + 4.493] = [-2.776, +6.210]
  - T2: [-3.578 - 2.323, -3.578 + 2.323] = [-5.901, -1.255]
  - T3: [-5.518 - 3.021, -5.518 + 3.021] = [-8.539, -2.497]
- T1 in light gray (insignificant — CI crosses zero)
- T2 and T3 in dark gray (significant — CIs do not cross zero)

**Annotations:**
- Next to each dot: N events in parentheses, e.g., "(58 events)"
- Significance stars next to T2 and T3 dots: "***"

**Title:** "The fuel coefficient strengthens over calendar time"

**Notes:** "Mean Fama-MacBeth fuel coefficient by event-year tercile. Horizontal bars show 95% confidence intervals using Newey-West standard errors. The fuel signal is absent in 2011–2013 and significant at the 1% level thereafter."

**Color:** Grayscale-safe. Dark circles for significant, light for insignificant.

**Size:** Single column width, approximately 3.5 x 2.5 inches.

### Production (R code sketch)
```r
library(ggplot2)
df <- data.frame(
  tercile = factor(c("2011–2013", "2013–2014", "2014–2022"),
                   levels = rev(c("2011–2013", "2013–2014", "2014–2022"))),
  beta = c(1.717, -3.578, -5.518),
  se = c(2.292, 1.185, 1.541),
  n = c(58, 58, 59),
  sig = c(FALSE, TRUE, TRUE)
)
df$ci_lo <- df$beta - 1.96 * df$se
df$ci_hi <- df$beta + 1.96 * df$se

ggplot(df, aes(x = beta, y = tercile)) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, color = "gray40") +
  geom_point(aes(fill = sig), shape = 21, size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "gray20", "FALSE" = "gray75"), guide = "none") +
  geom_text(aes(label = paste0("(", n, " events)")), hjust = -0.3, size = 2.8, nudge_y = 0.15) +
  labs(x = expression("Mean " * hat(gamma)[fuel]),
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave("results/figures/fig2_calendar_time.pdf", width = 3.5, height = 2.5)
```

---

## Figure 3: Global Map of Retirement Events and Sample Firms

**Purpose:** Visualise the paper's global coverage and the "distance does not insulate" framing.

**Placement:** Section 1 (Introduction) or Section 3.1 (Sample).

### Data sources
- Retirement events: `data/derived/events/coal_retirement_events.csv`
  - Filter: `is_first_mover == True` (179 events)
  - Columns: `lat`, `lon`, `capacity_mw`, `country`
- Firm centroids: `data/derived/networks/firm_centroids.csv`
  - Columns: `centroid_lat`, `centroid_lon`, `total_mw`
  - N = 414 firms with GPS coordinates

### Design

**Type:** World map (Robinson or Natural Earth projection) with two layers of points.

**Layer 1 — Retirement events:**
- Red/dark triangles (or diamonds), sized by `capacity_mw`
- Concentrated in Europe, North America, India, China, Australia
- Only first-mover events (179 points)

**Layer 2 — Sample firms:**
- Blue/light circles, sized by `total_mw`
- Scattered globally (80 countries)
- Smaller and lighter than retirement events to avoid visual dominance

**Background:**
- Light gray country boundaries
- No ocean fill (white background)
- No graticules (latitude/longitude grid lines — adds clutter)

**Annotations:**
- A text label in each major region: "35 retirements" (Europe), "135 retirements" (USA), "15 retirements" (India), etc. Only label the top 4-5 regions.
- No individual event or firm labels (too many points)

**Title:** "Coal retirements concentrate in a few regions; exposed utilities span the globe"

**Notes:** "Triangles mark 179 first-mover coal retirement events (2011–2022), sized by retired MW. Circles mark 414 sample utilities, sized by total installed capacity. Geographic distance between retirements and firms does not predict stock returns; fuel-mix similarity does."

**Color:** Two-color scheme: dark red/orange for retirements, steel blue for firms. Must be distinguishable in grayscale (use different shapes as backup). Colorblind-safe: avoid red-green pairs.

**Size:** Full page width, approximately 6.5 x 3.5 inches (landscape aspect for a world map).

### Production (R code sketch)
```r
library(ggplot2)
library(sf)
library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")

events <- read.csv("data/derived/events/coal_retirement_events.csv")
events <- events[events$is_first_mover == "True" & !is.na(events$lat), ]
firms <- read.csv("data/derived/networks/firm_centroids.csv")

ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray80", linewidth = 0.2) +
  geom_point(data = firms,
             aes(x = centroid_lon, y = centroid_lat, size = total_mw / 1000),
             color = "steelblue", alpha = 0.4, shape = 16) +
  geom_point(data = events,
             aes(x = lon, y = lat, size = capacity_mw / 100),
             color = "firebrick", alpha = 0.6, shape = 17) +
  scale_size_continuous(range = c(0.5, 4), guide = "none") +
  coord_sf(crs = "+proj=robin") +
  annotate("text", x = -100, y = -40, label = "▲ Retirement events\n● Sample utilities",
           size = 2.5, hjust = 0, color = "gray30") +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

ggsave("results/figures/fig3_world_map.pdf", width = 6.5, height = 3.5)
```

---

## LaTeX Integration

Add to `main.tex` preamble (already exists):
```latex
\usepackage{graphicx}
```

### Figure 1 placement (after Table 1 in Section 4.1):
```latex
\begin{figure}[H]
\centering
\includegraphics[width=\columnwidth]{../results/figures/fig1_fuel_vs_geo.pdf}
\caption{Fuel-mix similarity produces more negative returns than geographic proximity in 70 percent of events}
\label{fig:fuel_vs_geo}
\begin{flushleft}\footnotesize
Notes: Each point represents one Fama-MacBeth cross-sectional regression
(117 events with $\geq 20$ firms). Points below the 45-degree line indicate
events where $\hat{\gamma}_{\text{fuel}} < \hat{\gamma}_{\text{geo}}$.
Dark points: fuel more negative (82 events, 70\%). Light points: geography
more negative (35 events, 30\%).
\end{flushleft}
\end{figure}
```

### Figure 2 placement (in Section 4.5):
```latex
\begin{figure}[H]
\centering
\includegraphics[width=\columnwidth]{../results/figures/fig2_calendar_time.pdf}
\caption{The fuel coefficient strengthens over calendar time}
\label{fig:calendar_time}
\begin{flushleft}\footnotesize
Notes: Mean Fama-MacBeth fuel coefficient by event-year tercile. Horizontal
bars show 95\% confidence intervals using Newey-West standard errors. The
fuel signal is absent in 2011--2013 and significant at the 1 percent level
thereafter.
\end{flushleft}
\end{figure}
```

### Figure 3 placement (in Introduction or Section 3.1):
```latex
\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{../results/figures/fig3_world_map.pdf}
\caption{Coal retirements concentrate in a few regions; exposed utilities span the globe}
\label{fig:world_map}
\begin{flushleft}\footnotesize
Notes: Triangles mark 179 first-mover coal retirement events (2011--2022),
sized by retired MW. Circles mark 414 sample utilities, sized by total
installed capacity. Geographic distance does not predict stock returns
around retirements; fuel-mix similarity does.
\end{flushleft}
\end{figure}
```
