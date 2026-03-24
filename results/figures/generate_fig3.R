library(ggplot2)
library(sf)
library(rnaturalearth)
library(sysfonts)
library(showtext)

font_add_google("Roboto", "roboto")
font_add_google("Roboto Condensed", "roboto_cond")
showtext_auto()

robin <- "+proj=robin"
world <- ne_countries(scale = "medium", returnclass = "sf")
world_robin <- st_transform(world, robin)

events <- read.csv("data/derived/events/coal_retirement_events.csv",
                    stringsAsFactors = FALSE)
events <- events[events$is_first_mover == "True", ]
events <- events[!is.na(events$lat) & !is.na(events$lon), ]
events <- events[events$lat != 0 & events$lon != 0, ]
events_sf <- st_as_sf(events, coords = c("lon", "lat"), crs = 4326)
events_sf <- st_transform(events_sf, robin)
events_xy <- as.data.frame(st_coordinates(events_sf))

firms <- read.csv("data/derived/networks/firm_centroids.csv")
firms <- firms[!is.na(firms$centroid_lat) & !is.na(firms$centroid_lon), ]
firms <- firms[firms$centroid_lat != 0 & firms$centroid_lon != 0, ]
firms_sf <- st_as_sf(firms, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
firms_sf <- st_transform(firms_sf, robin)
firms_xy <- as.data.frame(st_coordinates(firms_sf))

# Two-color scheme: steel blue for firms, dark red for events
col_firm  <- "#6B8DAE"   # muted steel blue
col_event <- "#8B2500"   # dark brick red

p <- ggplot() +
  geom_sf(data = world_robin, fill = "gray97", color = "gray82", linewidth = 0.08) +
  geom_point(data = firms_xy, aes(x = X, y = Y),
             color = col_firm, size = 0.7, alpha = 0.45, shape = 16) +
  geom_point(data = events_xy, aes(x = X, y = Y),
             color = col_event, size = 1.2, alpha = 0.6, shape = 17) +
  # Manual legend with colored symbols
  annotate("point", x = -10500000, y = -4200000,
           shape = 17, size = 2.2, color = col_event) +
  annotate("text", x = -9600000, y = -4200000,
           label = paste0("Retirement events (", nrow(events_xy), ")"),
           size = 2.2, family = "roboto_cond", color = "gray25",
           hjust = 0, vjust = 0.5) +
  annotate("point", x = -10500000, y = -5100000,
           shape = 16, size = 1.5, color = col_firm) +
  annotate("text", x = -9600000, y = -5100000,
           label = paste0("Sample utilities (", nrow(firms_xy), ")"),
           size = 2.2, family = "roboto_cond", color = "gray25",
           hjust = 0, vjust = 0.5) +
  coord_sf(xlim = c(-12000000, 16000000), ylim = c(-5800000, 8500000),
           expand = FALSE) +
  theme_void(base_family = "roboto", base_size = 8) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave("results/figures/fig3_world_map.pdf", p, width = 6.5, height = 3.3,
       device = cairo_pdf)
cat("Saved: results/figures/fig3_world_map.pdf\n")
