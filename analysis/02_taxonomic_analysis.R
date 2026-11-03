## ============================================================================
## 02_taxonomic_analysis.R
##
## Is research effort taxonomically structured? An intercept-only negative
## binomial model of per-species paper counts with nested taxonomic random
## intercepts (class/order/family/genus) gives the expected number of papers
## per species at each taxonomic level (posterior median, 50% and 95% CrIs).
## Figure: (A) class/order/family tree with rings of median-paper percentiles,
## (B) share of species vs share of publications by class, (C) expected
## publications per species by order against class-level estimates.
##
## Produces: Figure — Taxa_Bias.pdf
##           Table  — class_effects.csv, order_effects.csv,
##                    family_effects.csv, genus_effects.csv
##           Object — m_taxa.rds (used by 05_counterfactual_deficits.R)
##
## Input data (in data/):
##   species_data_final.csv — species-level research effort, traits, geography
##                            and IUCN status (64,707 species)
##
## Run time:   ~<N> h for the brms fit                     [placeholder]
## Note: silhouettes are downloaded from PhyloPic (rphylopic); needs internet.
## ============================================================================

library(tidyverse)
library(tidybayes)
library(dagitty)
library(brms)
library(cmdstanr)
library(ggforce)
library(cowplot)
library(data.tree)
library(phangorn)
library(ggtree)
library(grid)
library(rphylopic)

seed <- 2009    # project-wide seed; identical in every script
set.seed(seed)
options(scipen = 999)
dir.create("output", showWarnings = FALSE)


# --- Load data ----------------------------------------------------------------

d <- read.csv("data/species_data_final.csv")
sapply(d, class)

## DAG as specified in 01_dag_development.R, kept for reference; the
## taxonomic model below has no exposure and so needs no adjustment set.
dag <- dagitty(
"dag {
abundance [latent]
biome []
charisma [latent]
gdp []
papers [outcome]
range []
remoteness []
size []
status []
year []
biogeography []

abundance -> status
biome -> papers
biome -> range
biome -> remoteness
biome -> size
biome -> status
biome -> year
biome -> gdp
charisma -> papers
gdp -> papers
gdp -> status
range -> abundance
range -> gdp
range -> papers
range -> status
range -> year
remoteness -> abundance
remoteness -> papers
remoteness -> range
remoteness -> size
remoteness -> status
remoteness -> year
size -> abundance
size -> charisma
size -> papers
size -> range
size -> status
size -> year
status -> papers
year -> papers
year -> status
biogeography -> gdp
biogeography -> remoteness
biogeography -> papers
biogeography -> year
biogeography -> size
biogeography -> range
}
"
)

plot(graphLayout(dag))


# --- Priors -------------------------------------------------------------------

## Weakly informative priors, shared by every model in the analysis
## (03_causal_models_global.R and 04_causal_models_by_class.R add normal(0, 1)
## on slopes). They exclude implausible values but leave ~65,000 species to
## dominate the posterior (Gelman 2006, 2008; McElreath 2020). The intercept
## is centred on log1p(median papers) rather than on the mean (34 papers),
## which a few extremely well-studied species inflate. The half-t on the
## taxonomic SDs neither shrinks the few-group class level towards zero nor
## behaves erratically like vaguer priors; a scale of 1 keeps the summed
## variation across four taxonomic levels realistic.
priors_re <- c(
  prior(normal(0.69, 2.5),  class = "Intercept"), # log(median(papers) + 1) = log(2); typical species in a typical genus.
                                                  # 95%: ~0.02-270 papers; covers all class medians (1-6)
  prior(student_t(3, 0, 1), class = "sd"),        # half-t (brms truncates at 0); all taxonomic levels.
                                                  # Typical ~2-fold between-group differences, >20-fold allowed; heavy tail
                                                  # avoids shrinking clade effects (Gelman 2006, few-group class level)
  prior(exponential(1),     class = "shape")      # NB phi; mean 1, central 95%: 0.03-3.7
)


# --- Taxonomic model ----------------------------------------------------------

## One model gives class-, order-, family- and genus-level effects. Class is a
## random effect (not a fixed factor) for consistency with the causal models,
## so no reference class is needed.
m_taxa <- brm(
  papers ~
    1 +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors_re,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_taxa)

## Convergence
plot(m_taxa, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_taxa, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution

saveRDS(m_taxa, "output/m_taxa.rds")


# --- Expected papers per species by taxonomic level ---------------------------

# Use posterior_epred() with re_formula to obtain conditional fitted values
# at each taxonomic level. re_formula controls which REs are included;
# omitted REs are set to zero, i.e. predictions for an average lower-level
# group (conditional, not marginalised over lower levels). This is
# equivalent to the manual draw-by-draw stacking (e.g. exp(intercept +
# r_class + r_order)) but avoids fragile joins on level names that brms
# may encode differently from the raw data (e.g. spaces → dots).

# Helper: summarise a posterior_epred matrix to quantiles
summarise_epred <- function(epred_matrix, nd) {
  nd %>% mutate(
    mean   = colMeans(epred_matrix),
    median = apply(epred_matrix, 2, median),
    Q2.5   = apply(epred_matrix, 2, quantile, probs = 0.025),
    Q25    = apply(epred_matrix, 2, quantile, probs = 0.25),
    Q75    = apply(epred_matrix, 2, quantile, probs = 0.75),
    Q97.5  = apply(epred_matrix, 2, quantile, probs = 0.975)
  )
}

# Class: intercept + class RE
nd_class <- d %>% distinct(class)
class_ref <- summarise_epred(
  posterior_epred(m_taxa, newdata = nd_class,
                  re_formula = ~ (1 | class)),
  nd_class
)

# Order: intercept + class + order REs
nd_order <- d %>% distinct(class, order)
order_ref <- summarise_epred(
  posterior_epred(m_taxa, newdata = nd_order,
                  re_formula = ~ (1 | class) + (1 | order)),
  nd_order
)

# Family: intercept + class + order + family REs
nd_family <- d %>% distinct(class, order, family)
family_ref <- summarise_epred(
  posterior_epred(m_taxa, newdata = nd_family,
                  re_formula = ~ (1 | class) + (1 | order) + (1 | family)),
  nd_family
)

# Genus: all REs (re_formula = NULL uses full model formula)
nd_genus <- d %>% distinct(class, order, family, genus)
genus_ref <- summarise_epred(
  posterior_epred(m_taxa, newdata = nd_genus, re_formula = NULL),
  nd_genus
)

## Supplementary tables
write_csv(class_ref,  "output/class_effects.csv")
write_csv(order_ref,  "output/order_effects.csv")
write_csv(family_ref, "output/family_effects.csv")
write_csv(genus_ref,  "output/genus_effects.csv")


# --- Expected publications per species by order (panel C) ---------------------

# Prepare order plot data
order_plot_data <- order_ref %>%
  mutate(
    class = factor(class, levels = c(
      "AMPHIBIA", "REPTILIA", "AVES", "MAMMALIA",
      "ACTINOPTERYGII", "CHONDRICHTHYES", "MYXINI",
      "PETROMYZONTI", "SARCOPTERYGII"
    )),
    order = fct_reorder2(order, class, median)
  ) %>%
  arrange(class, median) %>%
  mutate(
    order = fct_inorder(order),
    row = factor(
      if_else(class %in% c("AMPHIBIA", "REPTILIA", "AVES", "MAMMALIA"),
              "Tetrapods", "Fishes"),
      levels = c("Tetrapods", "Fishes")
    )
  )

# Class-level reference estimates
class_ref_df <- class_ref %>%
  rename(class_est  = median,
         class_lo95 = Q2.5,  class_lo50 = Q25,
         class_hi50 = Q75,   class_hi95 = Q97.5)

# Y-axis positions computed WITHIN each facet
class_yrange <- order_plot_data %>%
  group_by(row) %>%
  mutate(y_num = as.numeric(factor(order,
    levels = levels(order)[levels(order) %in% order]))) %>%
  ungroup() %>%
  group_by(row, class) %>%
  summarise(ymin = min(y_num) - 0.5, ymax = max(y_num) + 0.5, .groups = "drop") %>%
  left_join(class_ref_df, by = "class")

# Boundary lines
vlines <- class_yrange %>%
  group_by(row) %>%
  arrange(ymin) %>%
  mutate(boundary = ymax) %>%
  slice(-n()) %>%
  ungroup()

# Plot
p01 <- ggplot(order_plot_data, aes(x = order, y = median)) +

  # Class-level reference ribbon — 95% CI (lighter)
  geom_rect(data = class_yrange,
            aes(xmin = ymin, xmax = ymax, ymin = class_lo95, ymax = class_hi95),
            fill = "grey85", alpha = 0.4, inherit.aes = FALSE) +

  # Class-level reference ribbon — 50% CI (darker)
  geom_rect(data = class_yrange,
            aes(xmin = ymin, xmax = ymax, ymin = class_lo50, ymax = class_hi50),
            fill = "grey75", alpha = 0.4, inherit.aes = FALSE) +

  # Class-level reference line (posterior median)
  geom_segment(data = class_yrange,
               aes(x = ymin, xend = ymax, y = class_est, yend = class_est),
               colour = "grey50", linewidth = 0.3, linetype = "solid",
               inherit.aes = FALSE) +

  # Class boundary lines
  geom_vline(data = vlines, aes(xintercept = boundary),
             colour = "grey80", linewidth = 0.3) +

  # 95% CI — thin line
  geom_linerange(aes(ymin = Q2.5, ymax = Q97.5, colour = class),
                 linewidth = 0.5, alpha = 0.8) +

  # 50% CI — thick line
  geom_linerange(aes(ymin = Q25, ymax = Q75, colour = class),
                 linewidth = 1) +

  # Posterior median point
  geom_point(colour = "black", size = 1) +

  # Scales
  scale_color_manual(values = c(
    "AMPHIBIA"       = "#508070",
    "REPTILIA"       = "#785018",
    "AVES"           = "#F0B840",
    "MAMMALIA"       = "#D08020",
    "ACTINOPTERYGII" = "#60C8C0",
    "CHONDRICHTHYES" = "#28A0B0",
    "MYXINI"         = "#1878A0",
    "PETROMYZONTI"   = "#145070",
    "SARCOPTERYGII"  = "#082840"
  )) +
  scale_x_discrete(labels = str_to_sentence) +
  scale_y_continuous(transform = "log1p",
                     breaks = c(0, 1, 3, 10, 30, 100, 300, 1000),
                     expand = expansion(mult = c(0.05, 0.2))) +
  facet_wrap(~ row, nrow = 1, scales = "free_y") +
  labs(x = "Order", y = "Expected publications per species") +
  coord_flip() +
  theme_classic(base_size = 6) +
  theme(
    text = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none"
  )

p01


# --- Share of species vs share of publications by class (panel B) -------------

# Prepare proportional data
dumbbell_data <- d %>%
  group_by(class) %>%
  summarise(
    n_species    = n(),
    total_papers = sum(papers),
    .groups      = "drop"
  ) %>%
  mutate(
    pct_species = 100 * n_species / sum(n_species),
    pct_papers  = 100 * total_papers / sum(total_papers),
    # Bias direction: positive = over-studied relative to richness
    bias        = pct_papers - pct_species,
    # Label for annotation
    class_label = str_to_sentence(class)
  ) %>%
  # Sort by species richness (most species-rich at top)
  arrange(n_species) %>%
  mutate(class_label = fct_inorder(class_label))

# Axis transformation
# Log10 percentage axis to resolve Sarcopterygii/Myxini/Petromyzonti
# while keeping Actinopterygii/Aves readable.
# Custom breaks spanning the range of the data.

axis_breaks <- c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30)
axis_labels <- c("0.01", "0.03", "0.1", "0.3", "1", "3", "10", "30")

# Plot
p02 <- ggplot(dumbbell_data, aes(y = class_label)) +
  # Connecting segment (the "dumbbell bar")
  geom_segment(
    aes(
      x    = pct_species,
      xend = pct_papers,
      yend = class_label
    ),
    colour    = "grey60",
    linewidth = 0.4
  ) +

  # Species share point
  geom_point(
    aes(x = pct_species),
    shape = 21,
    fill  = "white",
    colour = "grey30",
    size  = 1,
    stroke = 0.4
  ) +

  # Paper share point
  geom_point(
    aes(x = pct_papers),
    shape = 16,
    colour = "grey30",
    size  = 1
  ) +

  # Log-scaled x-axis
  scale_x_log10(
    breaks = axis_breaks,
    labels = axis_labels,
    expand = expansion(mult = c(0.05, 0.10))  # right padding for annotations
  ) +

  scale_y_discrete(
    expand = expansion(mult = c(0.05, 0.15))
  ) +

  # Inline legend via annotation (top of panel)
  annotate("point", x = 0.036, y = 9.7, shape = 21, fill = "white",
           colour = "grey30", size = 5/.pt, stroke = 0.4) +
  annotate("text",  x = 0.05, y = 9.7, label = "% of species",
           size = 5/.pt, hjust = 0, vjust = 0.5) +
  annotate("point", x = 1.5,   y = 9.7, shape = 16,
           colour = "grey30", size = 5/.pt) +
  annotate("text",  x = 2, y = 9.7, label = "% of publications",
           size = 5/.pt, hjust = 0, vjust = 0.5) +
  geom_text(
     aes(x = max(pct_papers) * 2, label = scales::comma(n_species)),
     size = 5/.pt, hjust = 0, colour = "grey40", fontface = "italic"
   ) +
   annotate("text", x = max(dumbbell_data$pct_papers) * 2, y = 9.7,
            label = "Species", size = 5/.pt, hjust = 0, fontface = "italic") +
  labs(
    x = "Share of vertebrate total (%)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 6) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.2)
  )

p02


# --- Taxonomic tree with research-effort rings (panel A) ----------------------

## Class/order/family tree built from the taxonomy (topology only)
tax <- d %>%
  distinct(class, order, family) %>%
  arrange(class, order, family) %>%
  mutate(pathString = paste("Vertebrata", class, order, family, sep = "/"))

tax_tree <- as.Node(tax)
tree <- as.phylo(tax_tree)

## Rings, inner to outer: percentile rank of median papers per species at
## class, order and family level. Gaps separate classes around the fan.
# Tip metadata
tip_data <- d %>%
  group_by(class, order, family) %>%
  summarise(
    n_species     = n(),
    median_papers = median(papers),
    .groups = "drop"
  )

# Median papers at each taxonomic level
class_med <- d %>%
  group_by(class) %>%
  summarise(class_papers = median(papers), .groups = "drop")

order_med <- d %>%
  group_by(class, order) %>%
  summarise(order_papers = median(papers), .groups = "drop")

# Percentile ranks (within each level)
ring_data <- tip_data %>%
  left_join(order_med, by = c("class", "order")) %>%
  left_join(class_med, by = "class") %>%
  mutate(
    class_rank  = percent_rank(class_papers),
    order_rank  = percent_rank(order_papers),
    family_rank = percent_rank(median_papers)
  ) %>%
  select(family, class_rank, order_rank, family_rank) %>%
  column_to_rownames("family")

# Base tree
p_tree <- ggtree(tree, layout = "fan", size = 0.15) %<+% tip_data

# Inter-class gaps
gap_size    <- 5
pd          <- p_tree$data
n_tips      <- ape::Ntip(tree)
tip_spacing <- pd %>% filter(isTip) %>% pull(y) %>% sort() %>% diff() %>% median()

tips_sorted <- pd %>%
  filter(isTip) %>%
  arrange(y) %>%
  mutate(
    class        = setNames(tip_data$class, tip_data$family)[label],
    class_change = class != lag(class, default = first(class)),
    offset       = cumsum(class_change) * gap_size * tip_spacing
  )

tip_offsets <- setNames(tips_sorted$offset, tips_sorted$node)

desc <- Descendants(tree, seq_len(n_tips + tree$Nnode), type = "tips")
node_offsets <- vapply(desc, function(d) mean(tip_offsets[as.character(d)], na.rm = TRUE), numeric(1))

pd$y         <- pd$y + node_offsets[pd$node]
p_tree$data  <- pd
p_tree       <- p_tree + theme(legend.position = "none")

# Heatmap rings
p03 <- gheatmap(
  p_tree, ring_data,
  offset = 0.5, width = 0.3,
  colnames = FALSE, color = NA
) +
  # Research effort gradient (Violet)
  scale_fill_gradientn(
    colours = c("#f0eeef", "#b8b0c0", "#9890A0", "#605870", "#302838"),
    name    = "Publications\n(percentile)"
    ) +
  scale_y_continuous(expand = expansion(add = gap_size * tip_spacing / 2)) +
  scale_color_identity() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    text = element_text(size = 6),
    legend.title = element_text(hjust = 0.5),
    legend.text = element_text(hjust = 0.5),
    legend.title.position = "top",
    legend.key = element_rect(fill = NA)
  )

# Lighten the tree branches
p03$data$color <- "grey80"
p03 <- p03 + aes(color = color)

p03

# Cleanup
rm(tip_data, class_med, order_med, ring_data, p_tree, pd, n_tips,
   tip_spacing, tips_sorted, tip_offsets, desc, node_offsets, gap_size)


# --- Figure assembly ----------------------------------------------------------

## PhyloPic silhouettes recoloured to the class palette; the bird and mammal
## images are mirrored to face the same way as the others
p_acti <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Hemigrammus rhodostomus"), format = "raster"), fill = "#60C8C0"))
p_aves <- recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Buteo lagopus"), format = "raster"), fill = "#F0B840")
p_aves <- rasterGrob(p_aves[, ncol(p_aves):1, ]) ## flip bird
p_rept <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Cnemaspis niyomwanae"), format = "raster"), fill = "#785018"))
p_amph <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Fritziana fissilis"), format = "raster"), fill = "#508070"))
p_mamm <- recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Sylvilagus audubonii"), format = "raster"), fill = "#D08020")
  p_mamm <- rasterGrob(p_mamm[, ncol(p_mamm):1, ]) ## flip mammal
p_chon <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Carcharhinus obscurus"), format = "raster"), fill = "#28A0B0"))
p_myxi <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Myxine australis"), format = "raster"), fill = "#1878A0"))
p_petr <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Mesomyzon mengae"), format = "raster"), fill = "#145070"))
p_sarc <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Latimeria chalumnae"), format = "raster"), fill = "#082840"))

pTB <-
  ggdraw() +
  draw_plot(p03,
            x = 0, y = 2/3, width = 6/13, height = 1/3) +
  draw_plot(p02,
            x = 6/13, y = 2/3, width = 6/13, height = 1/3) +
  draw_plot(p01,
            x = 0, y = 0, width = 1, height = 2/3) +
  draw_plot_label(label = c("A", "B", "C"), size = 10,
                  x = c(0, 6/15, 0), y = c(1, 1, 2/3)) +
  ## Actin
  draw_grob(p_acti, x = 0.32, y = 0.7, width = 0.05, height = 0.05) +
  draw_grob(p_acti, x = 0.94, y = 0.934, width = 0.05, height = 0.05) +
  draw_grob(p_acti, x = 0.94, y = 0.24, width = 0.05, height = 0.05) +
  ## Aves
  draw_grob(p_aves, x = 0.01, y = 0.817, width = 0.04, height = 0.04) +
  draw_grob(p_aves, x = 0.942, y = 0.912, width = 0.04, height = 0.04) +
  draw_grob(p_aves, x = 0.46, y = 0.244, width = 0.04, height = 0.04) +
  ## Reptilia
  draw_grob(p_rept, x = 0.285, y = 0.94, width = 0.05, height = 0.05) +
  draw_grob(p_rept, x = 0.94, y = 0.875, width = 0.05, height = 0.05) +
  draw_grob(p_rept, x = 0.46, y = 0.055, width = 0.05, height = 0.05) +
  ## Amphibia
  draw_grob(p_amph, x = 0.36, y = 0.91, width = 0.04, height = 0.04) +
  draw_grob(p_amph, x = 0.94, y = 0.846, width = 0.04, height = 0.04) +
  draw_grob(p_amph, x = 0.46, y = 0.027, width = 0.04, height = 0.04) +
  ## Mammalia
  draw_grob(p_mamm, x = 0.135, y = 0.945, width = 0.04, height = 0.04) +
  draw_grob(p_mamm, x = 0.94, y = 0.817, width = 0.04, height = 0.04) +
  draw_grob(p_mamm, x = 0.46, y = 0.525, width = 0.04, height = 0.04) +
  ## Chondricthyes
  draw_grob(p_chon, x = 0.40, y = 0.85, width = 0.06, height = 0.05) +
  draw_grob(p_chon, x = 0.935, y = 0.78, width = 0.06, height = 0.05) +
  draw_grob(p_chon, x = 0.935, y = 0.525, width = 0.06, height = 0.05) +
  ## Myxini
  draw_grob(p_myxi, x = 0.41, y = 0.805, width = 0.045, height = 0.05) +
  draw_grob(p_myxi, x = 0.94, y = 0.748, width = 0.045, height = 0.05) +
  draw_grob(p_myxi, x = 0.94, y = 0.59, width = 0.045, height = 0.05) +
  ## Petromyzonti
  draw_grob(p_petr, x = 0.43, y = 0.815, width = 0.05, height = 0.05) +
  draw_grob(p_petr, x = 0.94, y = 0.717, width = 0.05, height = 0.05) +
  draw_grob(p_petr, x = 0.94, y = 0.602, width = 0.05, height = 0.05) +
  ## Sarcopterygii
  draw_grob(p_sarc, x = 0.41, y = 0.83, width = 0.05, height = 0.05) +
  draw_grob(p_sarc, x = 0.94, y = 0.689, width = 0.05, height = 0.05) +
  draw_grob(p_sarc, x = 0.94, y = 0.623, width = 0.05, height = 0.05)
pTB

ggsave("output/Taxa_Bias.pdf", pTB, width = 12.1, height = 18, units = "cm", dpi = 1200, bg = "white")
