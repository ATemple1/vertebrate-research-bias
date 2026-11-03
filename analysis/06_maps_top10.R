## ============================================================================
## 06_maps_top10.R
##
## Where is research effort concentrated, geographically and among species?
## Maps (A) species richness, (B) total publications and (C) mean
## publications per species on a global grid, and (D) shows the ten
## most-studied species in each of the six major classes.
##
## Produces: Figure — Map.pdf
##
## Input data (in data/):
##   global_richness.tif    — number of species whose range overlaps each cell
##   global_papersum.tif    — summed papers of the species in each cell
##   global_mean_papers.tif — mean papers per species in each cell
##   species_top10.csv      — ten most-studied species per major class
##
## Note: silhouettes are downloaded from PhyloPic (rphylopic); needs internet.
## ============================================================================

library(tidyverse)
library(rphylopic)
library(grid)
library(patchwork)
library(cowplot)
library(terra)
library(tidyterra)

seed <- 2009    # project-wide seed; identical in every script
set.seed(seed)
options(scipen = 999)
dir.create("output", showWarnings = FALSE)


# --- Maps (panels A-C) --------------------------------------------------------

## Cells with value 0 (no species in A; no publications in B and C) or Inf
## are set to NA so they plot as transparent
r1 <- rast("data/global_richness.tif")
r1[r1 == 0 | is.infinite(r1)] <- NA

r2 <- rast("data/global_papersum.tif")
r2[r2 == 0 | is.infinite(r2)] <- NA

r3 <- rast("data/global_mean_papers.tif")
r3[r3 == 0 | is.infinite(r3)] <- NA

## Panel A: species richness
p01 <- ggplot() +
  geom_spatraster(data = r1) +
  scale_fill_viridis_c(
    na.value = "transparent",
    trans = "log1p",
    breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000),
    limits = c(1,3000),
    guide = guide_colourbar(
      barheight = unit(2.5, "mm"),
      label.theme = element_text(margin = margin(t = 1)))
    ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5),
    legend.box.margin = margin(b = -17),
    legend.margin = margin(0, 0, 0, 0),
    text = element_text(size = 6)
  ) +
  labs( fill = "Species richness")
p01

## Panel B: total publications across the species in each cell
p02 <- ggplot() +
  geom_spatraster(data = r2) +
  scale_fill_viridis_c(
    na.value = "transparent",
    trans = "log1p",
    breaks = c(1, 10, 300, 10000, 300000),
    guide = guide_colourbar(
      barheight = unit(2.5, "mm"),
      label.theme = element_text(margin = margin(t = 1)))
    ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5),
    legend.box.margin = margin(b = -17),
    legend.margin = margin(0, 0, 0, 0),
    text = element_text(size = 6)
  ) +
  labs(fill = "Sum of publications across species")
p02

## Panel C: mean publications per species
p03 <- ggplot() +
  geom_spatraster(data = r3) +
  scale_fill_viridis_c(
    na.value = "transparent",
    trans = "log1p",
    breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000),
    limits = c(1,3000),
    guide = guide_colourbar(
      barheight = unit(2.5, "mm"),
      label.theme = element_text(margin = margin(t = 1)))
    ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title.position = "top",
    legend.title = element_text(hjust = 0.5),
    legend.box.margin = margin(b = -17),
    legend.margin = margin(0, 0, 0, 0),
    text = element_text(size = 6)
  ) +
  labs( fill = "Mean papers per species")
p03

## Stack panels A-C
p04 <- p01 / p02 / p03
p04


# --- Ten most-studied species per class (panel D) -----------------------------

d <- read.csv("data/species_top10.csv")

## Facet order: tetrapods, then fishes
class_order <- c("Mammalia", "Aves", "Reptilia",
                 "Amphibia", "Actinopterygii", "Chondrichthyes")

d <- d %>%
  mutate(class = factor(class, levels = class_order))

## Bars start at 300 papers; leading spaces in common_name offset the labels
## from the points
p05 <-
ggplot(d, aes(x = papers, y = reorder(rank, papers), colour = class)) +
  geom_point() +
  geom_text(aes(label = common_name), hjust = 0, size = 5/.pt) +
  geom_segment(aes(x = 300, xend = papers,
                 yend = reorder(rank, papers), colour = class),
             linewidth = 1) +

  facet_wrap(~ class, ncol = 1, scales = "free_y") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_colour_manual(values = c("#D08020", "#F0B840", "#785018", "#508070", "#60C8C0", "#28A0B0")) +
  scale_x_continuous(trans = "log1p", breaks = c(100, 300, 1000, 3000, 10000, 30000), expand = c(0, 0), limits = c(NA, 500000)) +
  labs(x = "Publications", y = "Ranking") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    strip.text = element_text(size = 6*1.2, face = "bold", hjust = 0, vjust = 1),
    strip.background = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
    )
p05

## PhyloPic silhouettes recoloured to the class palette (bird and mammal
## mirrored)
p_aves <- recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Buteo lagopus"), format = "raster"), fill = "#F0B840")
p_aves <- rasterGrob(p_aves[, ncol(p_aves):1, ]) ## flip bird

p_rept <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Cnemaspis niyomwanae"), format = "raster"), fill = "#785018"))

p_amph <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Fritziana fissilis"), format = "raster"), fill = "#508070"))

p_mamm <- recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Sylvilagus audubonii"), format = "raster"), fill = "#D08020")
p_mamm <- rasterGrob(p_mamm[, ncol(p_mamm):1, ]) ## flip mammal

p_acti <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Hemigrammus rhodostomus"), format = "raster"), fill = "#60C8C0"))

p_chon <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Carcharhinus obscurus"), format = "raster"), fill = "#28A0B0"))

## Add silhouettes to panel D
pT10 <- ggdraw(p05) +
  draw_grob(p_mamm, x = .75, y = .92, width = .12, height = .12) +
  draw_grob(p_aves, x = .75, y = .76, width = .12, height = .12) +
  draw_grob(p_rept, x = .75, y = .585, width = .15, height = .15) +
  draw_grob(p_amph, x = .75, y = .44, width = .12, height = .12) +
  draw_grob(p_acti, x = .75, y = .265, width = .15, height = .15) +
  draw_grob(p_chon, x = .73, y = .105, width = .18, height = .15)
pT10


# --- Figure assembly ----------------------------------------------------------

pMAP <- ggdraw() +
  draw_plot(p04,
            x = 0, y = 0, width = 2/3, height = 1) +
  draw_plot(pT10,
            x = 2/3, y = 0, width = 1/3, height = 1) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 10,
                  x = c(.01, .01, .01, 2/3+.01), y = c(.99, 2/3-.01, 1/3-.01, .99))
pMAP

ggsave("output/Map.pdf", pMAP, width = 18.4, height = 19.9, units = "cm", dpi = 1200, bg = "white")
