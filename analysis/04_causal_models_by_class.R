## ============================================================================
## 04_causal_models_by_class.R
##
## Do the drivers of research effort differ among vertebrate classes? Refits
## every causal model from 03_causal_models_global.R separately within each of
## the six major classes (Actinopterygii, Chondrichthyes, Amphibia, Reptilia,
## Aves, Mammalia). Predictors are standardised within class and the class
## random effect is dropped. Biome and IUCN status levels with fewer than
## min_n species in a class are left out of that class's biome or status
## model. Figure: (A) global IRRs, (B) class IRRs against the global estimate,
## (C) heat map of IRRs by class. E-values as in 03.
##
## Produces: Figure — Causal_Effects.pdf
##           Table  — class_model_effectsize.csv,
##                    class_model_pairwise_biome.csv,
##                    class_model_pairwise_status.csv, class_model_evalues.csv
##           Object — m_<driver>_<class>.rds, up to 48 fitted models, e.g.
##                    m_status_aves.rds (used by 05_counterfactual_deficits.R)
##
## Input data (in data/):
##   species_data_final.csv — species-level research effort, traits, geography
##                            and IUCN status (64,707 species)
##
## Depends on: output/global_model_effectsize.csv from 03_causal_models_global.R
## Run time:   ~2 days for the brms fits
## Note: silhouettes are downloaded from PhyloPic (rphylopic); needs internet.
## ============================================================================

library(tidyverse)
library(dagitty)
library(brms)
library(cmdstanr)
library(patchwork)
library(rphylopic)
library(grid)
library(cowplot)
library(rstan)
library(EValue)

seed <- 2009    # project-wide seed; identical in every script
set.seed(seed)
options(scipen = 999)
dir.create("output", showWarnings = FALSE)


# --- Load data ----------------------------------------------------------------

d <- read.csv("data/species_data_final.csv")
sapply(d, class)


# --- Classes and factor levels ------------------------------------------------

## Six classes with enough species for class-level models. Rare biome or
## status levels (< min_n species in a class) are dropped only from that
## class's biome or status model; those species stay in every other model.
## EX and EW are therefore kept wherever a class has enough of them.
classes <- c("ACTINOPTERYGII", "CHONDRICHTHYES", "AMPHIBIA",
             "REPTILIA", "AVES", "MAMMALIA")

min_n <- 3  ## minimum species per factor level to retain in biome/status models

## global level orderings (subsetted per class to levels meeting min_n)
## biome: last surviving level is the STZ dropped level
global_biome_order  <- c("terrestrial", "freshwater", "marine",
                         "terrestrial+freshwater", "terrestrial+marine",
                         "freshwater+marine")
## status: LC-CR, then DD, EW, EX; the last retained level in a class (EX
## wherever it meets min_n) is the sum-to-zero dropped level
global_status_order <- c("LC", "NT", "VU", "EN", "CR", "DD", "EW", "EX")


# --- DAG ----------------------------------------------------------------------

## As specified in 01_dag_development.R; adjustment sets are the same as in 03
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

## Same priors as 03_causal_models_global.R. Within a class the class random
## effect is redundant, so only order/family/genus intercepts remain.
priors_cls <- c(
  prior(normal(0.69, 2.5),  class = "Intercept"), # same prior in every class; class medians within 0.5 log units of centre
  prior(normal(0, 1),       class = "b"),         # predictors standardised within class: effect per 1 within-class SD
  prior(student_t(3, 0, 1), class = "sd"),        # half-t; order, family, genus
  prior(exponential(1),     class = "shape")      # NB phi
)


# --- Effect-size scale --------------------------------------------------------

## As in 03_causal_models_global.R: IRR = exp(beta); continuous effects per
## 2 SD (coefficients doubled at extraction), here 2 SD *within the class*;
## biome and IUCN status as sum-to-zero contrasts, so each level's IRR is
## relative to the grand mean of that class's retained levels, with the last
## retained level recovered as minus the sum of the others. Effect sizes are
## posterior medians with 50% and 95% CrIs and pd; pairwise contrasts are
## medians with 95% intervals and pd.

## One row per effect (or contrast) and class, accumulated across the models
std_eff <- data.frame()
pairwise_biome_all  <- data.frame()
pairwise_status_all <- data.frame()


# --- Geographic range ---------------------------------------------------------

## Range = range_km2, area of the IUCN/BirdLife range polygons. Each driver
## section loops over classes: subset, fit, save, store the model in the
## session (m_<driver>_<class>, used by the E-value section) and extract IRRs.

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "range", outcome = "papers")

for (cls in classes) {

  cat("\n--- Range:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## model
  m <- brm(
    papers ~
      scale(log1p(range_km2)) +
      scale(median_lat) +
      biome +
      scale(log1p(dist_min)) +
      scale(log1p(dist_extent)) +
      scale(log1p(body_mass_g)) +
      (1|order) + (1|family) + (1|genus),
    data = d_cls,
    family = negbinomial(),
    prior = priors_cls,
    chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = seed,
    backend = "cmdstanr"
  )

  assign(paste0("m_range_", cls_lbl), m)
  saveRDS(m, paste0("output/m_range_", cls_lbl, ".rds"))
  assign(paste0("r2_range_", cls_lbl), bayes_R2(m))
  cat("Saved: m_range_", cls_lbl, "\n")

  ## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
  std_eff <-
    rbind(std_eff,
          as_draws_df(m) %>%
            select(b = b_scalelog1prange_km2) %>%
            mutate(b2 = b * 2, irr = exp(b2)) %>%
            summarise(
              Parameter = "range_km2",
              class = cls,
              trans = "per 2 SD",
              median = median(b2),
              l95 = quantile(b2, 0.025), u95 = quantile(b2, 0.975),
              IRR = median(irr),
              IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
              IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
              pd = max(mean(b > 0), mean(b < 0)),
            ) %>%
            as_tibble()
    )

  rm(m)
}
gc()


# --- Remoteness ---------------------------------------------------------------

## Remoteness = minimum elevation/depth (dist_min) and its extent (dist_extent)

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "remoteness", outcome = "papers")

for (cls in classes) {

  cat("\n--- Remoteness:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## model
  m <- brm(
    papers ~
      scale(log1p(dist_min)) +
      scale(log1p(dist_extent)) +
      scale(median_lat) +
      biome +
      (1|order) + (1|family) + (1|genus),
    data = d_cls,
    family = negbinomial(),
    prior = priors_cls,
    chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = seed,
    backend = "cmdstanr"
  )

  assign(paste0("m_remoteness_", cls_lbl), m)
  saveRDS(m, paste0("output/m_remoteness_", cls_lbl, ".rds"))
  assign(paste0("r2_remoteness_", cls_lbl), bayes_R2(m))
  cat("Saved: m_remoteness_", cls_lbl, "\n")

  ## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd — dist_min
  std_eff <-
    rbind(std_eff,
          as_draws_df(m) %>%
            select(b = b_scalelog1pdist_min) %>%
            mutate(b2 = b * 2, irr = exp(b2)) %>%
            summarise(
              Parameter = "dist_min",
              class = cls,
              trans = "per 2 SD",
              median = median(b2),
              l95 = quantile(b2, 0.025), u95 = quantile(b2, 0.975),
              IRR = median(irr),
              IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
              IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
              pd = max(mean(b > 0), mean(b < 0)),
            ) %>%
            as_tibble()
    )

  ## extent
  std_eff <-
    rbind(std_eff,
          as_draws_df(m) %>%
            select(b = b_scalelog1pdist_extent) %>%
            mutate(b2 = b * 2, irr = exp(b2)) %>%
            summarise(
              Parameter = "dist_extent",
              class = cls,
              trans = "per 2 SD",
              median = median(b2),
              l95 = quantile(b2, 0.025), u95 = quantile(b2, 0.975),
              IRR = median(irr),
              IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
              IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
              pd = max(mean(b > 0), mean(b < 0)),
            ) %>%
            as_tibble()
    )

  rm(m)
}
gc()


# --- GDP ----------------------------------------------------------------------

## GDP = median_gdp_bil, summed GDP (median of the last 10 years) of the
## countries whose land or EEZ overlaps the species' range

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "gdp", outcome = "papers")

for (cls in classes) {

  cat("\n--- GDP:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## model
  m <- brm(
    papers ~
      scale(log1p(median_gdp_bil)) +
      scale(median_lat) +
      biome +
      scale(log1p(range_km2)) +
      (1|order) + (1|family) + (1|genus),
    data = d_cls,
    family = negbinomial(),
    prior = priors_cls,
    chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = seed,
    backend = "cmdstanr"
  )

  assign(paste0("m_gdp_", cls_lbl), m)
  saveRDS(m, paste0("output/m_gdp_", cls_lbl, ".rds"))
  assign(paste0("r2_gdp_", cls_lbl), bayes_R2(m))
  cat("Saved: m_gdp_", cls_lbl, "\n")

  ## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
  std_eff <-
    rbind(std_eff,
          as_draws_df(m) %>%
            select(b = b_scalelog1pmedian_gdp_bil) %>%
            mutate(b2 = b * 2, irr = exp(b2)) %>%
            summarise(
              Parameter = "median_gdp_bil",
              class = cls,
              trans = "per 2 SD",
              median = median(b2),
              l95 = quantile(b2, 0.025), u95 = quantile(b2, 0.975),
              IRR = median(irr),
              IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
              IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
              pd = max(mean(b > 0), mean(b < 0)),
            ) %>%
            as_tibble()
    )

  rm(m)
}
gc()


# --- Body size ----------------------------------------------------------------

## Size = body_mass_g, from species trait databases

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "size", outcome = "papers")

for (cls in classes) {

  cat("\n--- Size:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## model
  m <- brm(
    papers ~
      scale(log1p(body_mass_g)) +
      scale(median_lat) +
      biome +
      scale(log1p(dist_min)) +
      scale(log1p(dist_extent)) +
      (1|order) + (1|family) + (1|genus),
    data = d_cls,
    family = negbinomial(),
    prior = priors_cls,
    chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = seed,
    backend = "cmdstanr"
  )

  assign(paste0("m_size_", cls_lbl), m)
  saveRDS(m, paste0("output/m_size_", cls_lbl, ".rds"))
  assign(paste0("r2_size_", cls_lbl), bayes_R2(m))
  cat("Saved: m_size_", cls_lbl, "\n")

  ## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
  std_eff <-
    rbind(std_eff,
          as_draws_df(m) %>%
            select(b = b_scalelog1pbody_mass_g) %>%
            mutate(b2 = b * 2, irr = exp(b2)) %>%
            summarise(
              Parameter = "body_mass_g",
              class = cls,
              trans = "per 2 SD",
              median = median(b2),
              l95 = quantile(b2, 0.025), u95 = quantile(b2, 0.975),
              IRR = median(irr),
              IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
              IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
              pd = max(mean(b > 0), mean(b < 0)),
            ) %>%
            as_tibble()
    )

  rm(m)
}
gc()


# --- Year of description ------------------------------------------------------

## Year = year of description, from the species' taxonomic authority

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "year", outcome = "papers")

for (cls in classes) {

  cat("\n--- Year:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## model
  m <- brm(
    papers ~
      scale(year) +
      scale(median_lat) +
      biome +
      scale(log1p(range_km2)) +
      scale(log1p(dist_min)) +
      scale(log1p(dist_extent)) +
      scale(log1p(body_mass_g)) +
      (1|order) + (1|family) + (1|genus),
    data = d_cls,
    family = negbinomial(),
    prior = priors_cls,
    chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = seed,
    backend = "cmdstanr"
  )

  assign(paste0("m_year_", cls_lbl), m)
  saveRDS(m, paste0("output/m_year_", cls_lbl, ".rds"))
  assign(paste0("r2_year_", cls_lbl), bayes_R2(m))
  cat("Saved: m_year_", cls_lbl, "\n")

  ## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
  std_eff <-
    rbind(std_eff,
          as_draws_df(m) %>%
            select(b = b_scaleyear) %>%
            mutate(b2 = b * 2, irr = exp(b2)) %>%
            summarise(
              Parameter = "year",
              class = cls,
              trans = "per 2 SD",
              median = median(b2),
              l95 = quantile(b2, 0.025), u95 = quantile(b2, 0.975),
              IRR = median(irr),
              IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
              IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
              pd = max(mean(b > 0), mean(b < 0)),
            ) %>%
            as_tibble()
    )

  rm(m)
}
gc()


# --- Biogeography -------------------------------------------------------------

## Biogeography = median_lat (signed), capturing the Global North-South axis

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "biogeography", outcome = "papers")

for (cls in classes) {

  cat("\n--- Biogeography:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## model
  m <- brm(
    papers ~
      scale(median_lat) +
      (1|order) + (1|family) + (1|genus),
    data = d_cls,
    family = negbinomial(),
    prior = priors_cls,
    chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = seed,
    backend = "cmdstanr"
  )

  assign(paste0("m_biogeography_", cls_lbl), m)
  saveRDS(m, paste0("output/m_biogeography_", cls_lbl, ".rds"))
  assign(paste0("r2_biogeography_", cls_lbl), bayes_R2(m))
  cat("Saved: m_biogeography_", cls_lbl, "\n")

  ## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
  std_eff <-
    rbind(std_eff,
          as_draws_df(m) %>%
            select(b = b_scalemedian_lat) %>%
            mutate(b2 = b * 2, irr = exp(b2)) %>%
            summarise(
              Parameter = "median_lat",
              class = cls,
              trans = "per 2 SD",
              median = median(b2),
              l95 = quantile(b2, 0.025), u95 = quantile(b2, 0.975),
              IRR = median(irr),
              IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
              IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
              pd = max(mean(b > 0), mean(b < 0)),
            ) %>%
            as_tibble()
    )

  rm(m)
}
gc()


# --- Biome --------------------------------------------------------------------

## Multi-biome species form their own levels (e.g. freshwater+marine), as in
## 03. A class needs at least two retained levels for a biome model.

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "biome", outcome = "papers")

for (cls in classes) {

  cat("\n--- Biome:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## 1. Identify well-represented biome levels in this class
  biome_counts <- table(d_cls$biome)
  cat("Biome counts:\n")
  print(biome_counts)

  cls_biome_levels <- names(biome_counts[biome_counts >= min_n])
  cls_biome_levels <- intersect(global_biome_order, cls_biome_levels)  # preserve global ordering
  n_biome <- length(cls_biome_levels)

  cat("Retained biome levels (n >=", min_n, "):", paste(cls_biome_levels, collapse = ", "), "\n")
  if (n_biome >= 2) cat("Dropped biome level:", cls_biome_levels[n_biome], "\n")

  if (n_biome >= 2) {

    ## 2. Create biome-model data (species in rare biome levels excluded here only)
    d_biome <- d_cls %>%
      filter(biome %in% cls_biome_levels) %>%
      droplevels()
    d_biome$biome <- factor(d_biome$biome, levels = cls_biome_levels)
    contrasts(d_biome$biome) <- contr.sum(n_biome)
    contrasts(d_biome$biome)

    ## 3. Fit model
    m <- brm(
      papers ~
        biome +
        (1|order) + (1|family) + (1|genus),
      data = d_biome,
      family = negbinomial(),
      prior = priors_cls,
      chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
      control = list(adapt_delta = 0.99, max_treedepth = 12),
      seed = seed,
      backend = "cmdstanr"
    )

    assign(paste0("m_biome_", cls_lbl), m)
    saveRDS(m, paste0("output/m_biome_", cls_lbl, ".rds"))
    assign(paste0("r2_biome_", cls_lbl), bayes_R2(m))
    cat("Saved: m_biome_", cls_lbl, "\n")

    ## 4. Extract all biome effects including the dropped level
    draws_raw_biome <- as_draws_df(m) %>%
      as_tibble() %>%
      select(starts_with("b_biome"))

    draws_all_biome <- draws_raw_biome %>%
      mutate(
        b_biome_dropped = -rowSums(across(everything()))
      )

    estimated_names_biome <- paste0("b_biome", cls_biome_levels[1:(n_biome - 1)])
    all_names_biome <- c(estimated_names_biome, paste0("b_biome", cls_biome_levels[n_biome]))
    names(draws_all_biome) <- all_names_biome

    ## 5. IRR of every retained level vs the class grand mean
    std_eff <-
      rbind(
        std_eff,
        draws_all_biome %>%
          pivot_longer(everything(),
                       names_to = "Parameter", values_to = "b") %>%
          mutate(
            Parameter = sub("^b_biome", "", Parameter),
            Parameter = paste(Parameter, "vs grand mean"),
            irr = exp(b)
          ) %>%
          group_by(Parameter) %>%
          summarise(
            class = cls,
            trans = "none (factor, sum-to-zero)",
            median = median(b),
            l95 = quantile(b, 0.025), u95 = quantile(b, 0.975),
            IRR = median(irr),
            IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
            IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
            pd = max(mean(b > 0), mean(b < 0)),
            .groups = "drop"
          )
      )

    ## 6. Pairwise contrasts between all biome levels
    pairwise_biome_all <-
      rbind(
        pairwise_biome_all,
        combn(cls_biome_levels, 2, simplify = FALSE) %>%
          map_dfr(function(pair) {
            a <- pair[1]
            b <- pair[2]
            diff <- draws_all_biome[[paste0("b_biome", a)]] - draws_all_biome[[paste0("b_biome", b)]]
            irr_diff <- exp(diff)
            tibble(
              class = cls,
              contrast = paste(a, "vs", b),
              diff_median = median(diff),
              diff_l95 = quantile(diff, 0.025),
              diff_u95 = quantile(diff, 0.975),
              IRR_ratio = median(irr_diff),
              IRR_ratio_l95 = quantile(irr_diff, 0.025),
              IRR_ratio_u95 = quantile(irr_diff, 0.975),
              pd = max(mean(diff > 0), mean(diff < 0)),
              sig_95 = !between(0, quantile(diff, 0.025), quantile(diff, 0.975))
            )
          })
      )

    ## clean up
    rm(m, d_biome, draws_raw_biome, draws_all_biome, estimated_names_biome, all_names_biome)

  } else {
    cat("WARNING: <2 biome levels meet min_n — skipping biome model for", cls, "\n")
  }

  rm(biome_counts, cls_biome_levels, n_biome)
}
gc()


# --- IUCN status --------------------------------------------------------------

## Status = IUCN Red List category; levels retained per class as for biome

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "status", outcome = "papers")

for (cls in classes) {

  cat("\n--- Status:", cls, "---\n")
  cls_lbl <- tolower(cls)
  d_cls   <- d %>% filter(class == cls) %>% droplevels()
  cat("n =", nrow(d_cls), "\n")

  ## 1. Identify well-represented status levels in this class
  status_counts <- table(d_cls$status)
  cat("Status counts:\n")
  print(status_counts)

  cls_status_levels <- names(status_counts[status_counts >= min_n])
  cls_status_levels <- intersect(global_status_order, cls_status_levels)  # preserve global ordering
  n_status <- length(cls_status_levels)

  cat("Retained status levels (n >=", min_n, "):", paste(cls_status_levels, collapse = ", "), "\n")
  cat("Dropped status level:", cls_status_levels[n_status], "\n")

  ## 2. Create status-model data (species in rare status levels excluded here only)
  d_status <- d_cls %>%
    filter(status %in% cls_status_levels) %>%
    droplevels()
  d_status$status <- factor(d_status$status, levels = cls_status_levels)
  contrasts(d_status$status) <- contr.sum(n_status)

  ## 3. Check what R will estimate
  contrasts(d_status$status)

  ## 4. Fit model
  m <- brm(
    papers ~
      status +
      biome +
      scale(log1p(median_gdp_bil)) +
      scale(log1p(range_km2)) +
      scale(log1p(dist_min)) +
      scale(log1p(dist_extent)) +
      scale(log1p(body_mass_g)) +
      scale(year) +
      (1|order) + (1|family) + (1|genus),
    data = d_status,
    family = negbinomial(),
    prior = priors_cls,
    chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    seed = seed,
    backend = "cmdstanr"
  )

  assign(paste0("m_status_", cls_lbl), m)
  saveRDS(m, paste0("output/m_status_", cls_lbl, ".rds"))
  assign(paste0("r2_status_", cls_lbl), bayes_R2(m))
  cat("Saved: m_status_", cls_lbl, "\n")

  ## 5. Extract all status effects including the dropped level
  draws_raw_status <- as_draws_df(m) %>%
    as_tibble() %>%
    select(starts_with("b_status"))

  # brms names sum-to-zero columns as status1, status2, ..., status(k-1)
  # These correspond to levels 1:(k-1) in the factor ordering
  # The dropped level (last retained level) = negative sum of all others
  draws_all_status <- draws_raw_status %>%
    mutate(
      b_status_dropped = -rowSums(across(everything()))
    )

  # Rename columns to actual level names
  estimated_names_status <- paste0("b_status", cls_status_levels[1:(n_status - 1)])
  all_names_status <- c(estimated_names_status, paste0("b_status", cls_status_levels[n_status]))
  names(draws_all_status) <- all_names_status

  ## Store the retained levels for this class (kept for reference; not used below)
  assign(paste0("status_levels_", cls_lbl), cls_status_levels)

  ## 6. IRR of every retained level vs the class grand mean
  std_eff <-
    rbind(
      std_eff,
      draws_all_status %>%
        pivot_longer(everything(),
                     names_to = "Parameter", values_to = "b") %>%
        mutate(
          Parameter = sub("^b_status", "", Parameter),
          Parameter = paste(Parameter, "vs grand mean"),
          irr = exp(b)
        ) %>%
        group_by(Parameter) %>%
        summarise(
          class = cls,
          trans = "none (factor, sum-to-zero)",
          median = median(b),
          l95 = quantile(b, 0.025), u95 = quantile(b, 0.975),
          IRR = median(irr),
          IRR_lower95 = quantile(irr, 0.025), IRR_upper95 = quantile(irr, 0.975),
          IRR_lower50 = quantile(irr, 0.25),  IRR_upper50 = quantile(irr, 0.75),
          pd = max(mean(b > 0), mean(b < 0)),
          .groups = "drop"
        )
    )

  ## 7. Pairwise contrasts between all status levels
  pairwise_status_all <-
    rbind(
      pairwise_status_all,
      combn(cls_status_levels, 2, simplify = FALSE) %>%
        map_dfr(function(pair) {
          a <- pair[1]
          b <- pair[2]
          diff <- draws_all_status[[paste0("b_status", a)]] - draws_all_status[[paste0("b_status", b)]]
          irr_diff <- exp(diff)
          tibble(
            class = cls,
            contrast = paste(a, "vs", b),
            diff_median = median(diff),
            diff_l95 = quantile(diff, 0.025),
            diff_u95 = quantile(diff, 0.975),
            IRR_ratio = median(irr_diff),
            IRR_ratio_l95 = quantile(irr_diff, 0.025),
            IRR_ratio_u95 = quantile(irr_diff, 0.975),
            pd = max(mean(diff > 0), mean(diff < 0)),
            sig_95 = !between(0, quantile(diff, 0.025), quantile(diff, 0.975))
          )
        })
    )

  ## clean up
  rm(m, d_status, status_counts, cls_status_levels, n_status,
     draws_raw_status, draws_all_status, estimated_names_status, all_names_status)
}
gc()

## Accumulated effect sizes and pairwise contrasts
std_eff
pairwise_biome_all
pairwise_status_all


# --- Save effect sizes --------------------------------------------------------

write.csv(std_eff, "output/class_model_effectsize.csv", row.names = FALSE)
write.csv(pairwise_biome_all, "output/class_model_pairwise_biome.csv", row.names = FALSE)
write.csv(pairwise_status_all, "output/class_model_pairwise_status.csv", row.names = FALSE)


# --- Figure: global and class-level effects (Causal_Effects.pdf) --------------

## Complete every class x parameter combination (NA where a level was not
## modelled in a class) so all panels share the same rows
std_eff_all <- std_eff %>%
  expand(class, Parameter) %>%
  left_join(std_eff)

## Global estimates from 03_causal_models_global.R
std_eff_global <- read.csv("output/global_model_effectsize.csv")

## Panel A: global model. Rows: continuous drivers in a fixed order (low to
## high IRR in the reported fit; update lev_order after a refit), then biome,
## then status
lev_order <- c(
  "year", "dist_min", "median_lat", "dist_extent", "median_gdp_bil",  "range_km2", "body_mass_g", # continuous rows
  "terrestrial vs grand mean", "freshwater vs grand mean", "marine vs grand mean", "terrestrial+freshwater vs grand mean",
  "terrestrial+marine vs grand mean", "freshwater+marine vs grand mean", # biome
  "DD vs grand mean", "LC vs grand mean", "NT vs grand mean", "VU vs grand mean", "EN vs grand mean", "CR vs grand mean",
  "EW vs grand mean", "EX vs grand mean" # status
  )

p01 <-
  ggplot(std_eff_global, aes(y = Parameter, x = IRR)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "black") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "black") +
  geom_point(size = 1, colour = "black") +
  annotate("text", x = Inf, y = 14.5, label = "Biome \n(vs grand mean) ",
           hjust = 1, vjust = 1.1, colour = "black", size = 5/.pt) +
  annotate("text", x = Inf, y = 8.5, label = "IUCN Status \n(vs grand mean) ",
           hjust = 1, vjust = 1.1, colour = "black", size = 5/.pt) +
  scale_x_log10() +
  scale_y_discrete(limits = rev(lev_order),
                   labels = c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD",
                              "Freshwater-\nmarine", "Terrestrial-\nmarine",
                              "Terrestrial-\nfreshwater", "Marine", "Freshwater",
                              "Terrestrial", "Body mass", "Geographic\n range",
                              "Gross domestic\nproduct", "Elevation/depth\nextent",
                              "Median\nlatitude", "Elevation/depth\nminimum",
                              "Year of\ndescription" ))  + # rev() so lev_order reads top-to-bottom
  labs(x = "IRR", title = "Global model") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", colour = "black"))
p01

## Panel B: Mammalia (class estimate; | marks the global estimate)
p02 <-
  ggplot(std_eff_all %>% filter(class == "MAMMALIA"),
         aes(y = Parameter, x = IRR)) +
  ## annotation lines
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  # global model
  geom_point(data = std_eff_global, colour = "black", size = 2, shape = "|") +
  # class model
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "#D08020") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "#D08020") +
  geom_point(size = 1, shape = 21, colour = "#D08020", fill = "black") +
  scale_x_log10() +
  scale_y_discrete(limits = rev(lev_order),
                   labels = c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD",
                              "F-M", "T-M", "T-F",
                              "Mar.", "Frw.", "Ter.",
                              "Mas.", "Rng.", "GDP",
                              "Ext.", "Lat.",
                              "Min.", "Year"))  + # rev() so lev_order reads top-to-bottom
  labs(x = "IRR", title = "Mammalia") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title = element_blank(),
        plot.title = element_text(face = "bold"))
p02

## Panel B: Aves (class estimate; | marks the global estimate)
p03 <-
  ggplot(std_eff_all %>% filter(class == "AVES"),
         aes(y = Parameter, x = IRR)) +
  ## annotation lines
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  # global model
  geom_point(data = std_eff_global, colour = "black", size = 2, shape = "|") +
  # class model
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "#F0B840") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "#F0B840") +
  geom_point(size = 1, shape = 21, colour = "#F0B840", fill = "black") +
  scale_y_discrete(limits = rev(lev_order)) +
  scale_x_log10() +
  labs(x = "IRR", title = "Aves") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(face = "bold"))
p03

## Panel B: Reptilia (class estimate; | marks the global estimate)
p04 <-
  ggplot(std_eff_all %>% filter(class == "REPTILIA"),
         aes(y = Parameter, x = IRR)) +
  ## annotation lines
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  # global model
  geom_point(data = std_eff_global, colour = "black", size = 2, shape = "|") +
  # class model
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "#785018") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "#785018") +
  geom_point(size = 1, shape = 21, colour = "#785018", fill = "black") +
  scale_y_discrete(limits = rev(lev_order)) +
  scale_x_log10() +
  labs(x = "IRR", title = "Reptilia") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(face = "bold"))
p04

## Panel B: Amphibia (class estimate; | marks the global estimate)
p05 <-
  ggplot(std_eff_all %>% filter(class == "AMPHIBIA"),
         aes(y = Parameter, x = IRR)) +
  ## annotation lines
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  # global model
  geom_point(data = std_eff_global, colour = "black", size = 2, shape = "|") +
  # class model
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "#508070") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "#508070") +
  geom_point(size = 1, shape = 21, colour = "#508070", fill = "black") +
  scale_x_log10() +
  scale_y_discrete(limits = rev(lev_order),
                   labels = c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD",
                              "F-M", "T-M", "T-F",
                              "Mar.", "Frw.", "Ter.",
                              "Mas.", "Rng.", "GDP",
                              "Ext.", "Lat.",
                              "Min.", "Year"))  + # rev() so lev_order reads top-to-bottom
  labs(x = "IRR", title = "Amphibia") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold"))
p05

## Panel B: Actinopterygii (class estimate; | marks the global estimate)
p06 <-
  ggplot(std_eff_all %>% filter(class == "ACTINOPTERYGII"),
         aes(y = Parameter, x = IRR)) +
  ## annotation lines
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  # global model
  geom_point(data = std_eff_global, colour = "black", size = 2, shape = "|") +
  # class model
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "#60C8C0") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "#60C8C0") +
  geom_point(size = 1, shape = 21, colour = "#60C8C0", fill = "black") +
  scale_y_discrete(limits = rev(lev_order)) +
  scale_x_log10() +
  labs(x = "IRR", title = "Actinopterygii") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(face = "bold"))
p06

## Panel B: Chondrichthyes (class estimate; | marks the global estimate)
p07 <-
  ggplot(std_eff_all %>% filter(class == "CHONDRICHTHYES"),
         aes(y = Parameter, x = IRR)) +
  ## annotation lines
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  # global model
  geom_point(data = std_eff_global, colour = "black", size = 2, shape = "|") +
  # class model
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "#28A0B0") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "#28A0B0") +
  geom_point(size = 1, shape = 21, colour = "#28A0B0", fill = "black") +
  scale_y_discrete(limits = rev(lev_order)) +
  scale_x_log10() +
  labs(x = "IRR", title = "Chondrichthyes") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(face = "bold"))
p07

## Panel C: heat map of IRRs (columns: global, then classes); grey = not
## estimated in that class, or IRR outside the 0.05-20 colour limits
p08 <- full_join(std_eff_all, std_eff_global) %>%
  mutate(class  = case_when(is.na(class) ~ "Global",
                            .default = class)) %>%
  ggplot(aes(x = class, y = Parameter, fill = IRR)) +
  geom_tile() +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  scale_x_discrete(limits = c("Global", "MAMMALIA", "AVES", "REPTILIA", "AMPHIBIA", "ACTINOPTERYGII", "CHONDRICHTHYES")) +
  scale_y_discrete(limits = rev(lev_order),
                   labels = c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD",
                              "Freshwater-\nMarine", "Terrestrial-\nMarine",
                              "Terrestrial-\nFreshwater", "Marine", "Freshwater",
                              "Terrestrial", "Body Mass", "Geographic\n Range",
                              "Gross Domestic\nProduct", "Elevation/Depth\nExtent",
                              "Median\nLatitude", "Elevation/Depth\nMinimum",
                              "Year of\nDescription" ))  + # rev() so lev_order reads top-to-bottom
  scale_fill_gradient2(
    low = "#6B5FB5",
    mid = "#f7f7f7",
    high = "#A65C6A",
    na.value = "grey90",
    midpoint = 1.0, # Forces the center of the gradient to 1
    transform = "log10",
    limits = c(.05, 20),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  guides(fill = guide_colourbar(
    theme = theme(
  legend.key.width  = unit(4, "lines"),
  legend.key.height = unit(.5, "lines"),
  legend.text.position = "top"
))) +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.text = element_text(margin = margin(b = 2.5)),  # space between labels and bar
        legend.margin = margin(-36.5, 0, 0, 0),
        legend.position="bottom",
        legend.title = element_blank())
p08

## PhyloPic silhouettes recoloured to the class palette (bird and mammal
## mirrored); a black circle marks the global model
p_glob <- circleGrob(r = 0.5, gp = gpar(fill = "black", col = "black"))

p_aves <- recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Buteo lagopus"), format = "raster"), fill = "#F0B840")
p_aves <- rasterGrob(p_aves[, ncol(p_aves):1, ]) ## flip bird

p_rept <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Cnemaspis niyomwanae"), format = "raster"), fill = "#785018"))

p_amph <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Fritziana fissilis"), format = "raster"), fill = "#508070"))

p_mamm <- recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Sylvilagus audubonii"), format = "raster"), fill = "#D08020")
p_mamm <- rasterGrob(p_mamm[, ncol(p_mamm):1, ]) ## flip mammal

p_acti <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Hemigrammus rhodostomus"), format = "raster"), fill = "#60C8C0"))

p_chon <- rasterGrob(recolor_phylopic(get_phylopic(uuid = get_uuid(name = "Carcharhinus obscurus"), format = "raster"), fill = "#28A0B0"))

## Assemble: A = global, B = class panels, C = heat map
pCB <-
  ggdraw((p01 | (p02 | p03 | p04) / (p05 | p06 | p07) | p08) +
  plot_layout(ncol = 3, widths = c(1, 3, .5))) +
  ## markings for subplots of std_eff
  draw_grob(p_glob, x = .235, y = .935, width = .037, height = .037) +
  draw_grob(p_mamm, x = .44, y = .935, width = .05, height = .05) +
  draw_grob(p_aves, x = .63, y = .935, width = .05, height = .05) +
  draw_grob(p_rept, x = .82, y = .93, width = .05, height = .05) +
  draw_grob(p_amph, x = .44, y = .48, width = .04, height = .04) +
  draw_grob(p_acti, x = .63, y = .475, width = .05, height = .05) +
  draw_grob(p_chon, x = .82, y = .475, width = .06, height = .05) +
  ## labels for heatmap
  draw_grob(p_glob, x = .888, y = .94, width = .016, height = .016) +
  draw_grob(editGrob(p_mamm, vp = viewport(angle = 90)), x = .897, y = .94, width = .02, height = .02) +
  draw_grob(editGrob(p_aves, vp = viewport(angle = 90)), x = .91, y = .94, width = .02, height = .02) +
  draw_grob(editGrob(p_rept, vp = viewport(angle = 90)), x = .922, y = .94, width = .02, height = .02) +
  draw_grob(editGrob(p_amph, vp = viewport(angle = 90)), x = .937, y = .94, width = .017, height = .017) +
  draw_grob(editGrob(p_acti, vp = viewport(angle = 90)), x = .948, y = .94, width = .02, height = .02) +
  draw_grob(editGrob(p_chon, vp = viewport(angle = 90)), x = .96, y = .94, width = .02, height = .02) +
  ## plot labels
  draw_plot_label(label = c("A", "B", "C"), size = 10,
                  x = c(.01, .28, .87), y = c(.99, .99, .99)) +
  draw_label("IRR", x = 0.935, y = 0.075, size = 6)

pCB

ggsave("output/Causal_Effects.pdf", pCB, width = 18.4, height = 12, units = "cm", dpi = 1200, bg = "white")


# --- Sensitivity to unmeasured confounding (E-values) -------------------------

## As in 03_causal_models_global.R, per class. The class models are taken
## from the session objects created above (m_<driver>_<class>).
evalue_nb <- function(model, predictor) {
  k     <- if (startsWith(predictor, "scale")) 2 else 1   # 2 s.d. for continuous
  fe    <- fixef(model, summary = TRUE, robust = TRUE)
  beta  <- fe[predictor, "Estimate"] * k
  lower <- fe[predictor, "Q2.5"] * k
  upper <- fe[predictor, "Q97.5"] * k
  evalues.RR(
    est = exp(beta),
    lo  = exp(lower),
    hi  = exp(upper)
  )
}

evalue_nb_dropped <- function(model, factor_prefix) {
  draws <- as_draws_df(model) %>%
    as_tibble() %>%
    select(starts_with(paste0("b_", factor_prefix)))

  dropped_draws <- -rowSums(draws)

  beta  <- median(dropped_draws)
  lower <- quantile(dropped_draws, 0.025)
  upper <- quantile(dropped_draws, 0.975)

  evalues.RR(
    est = exp(beta),
    lo  = exp(lower),
    hi  = exp(upper)
  )
}

## Apply across models — per class
evalue_table_all <- data.frame()

for (cls in classes) {

  cls_lbl <- tolower(cls)

  ## retrieve class-level models
  m_size_cls          <- get(paste0("m_size_", cls_lbl))
  m_range_cls         <- get(paste0("m_range_", cls_lbl))
  m_remoteness_cls    <- get(paste0("m_remoteness_", cls_lbl))
  m_gdp_cls           <- get(paste0("m_gdp_", cls_lbl))
  m_biogeography_cls  <- get(paste0("m_biogeography_", cls_lbl))
  m_year_cls          <- get(paste0("m_year_", cls_lbl))
  m_status_cls        <- get(paste0("m_status_", cls_lbl))

  ## build model_specs — continuous predictors
  model_specs <- list(
    list(model = m_size_cls,         predictor = "scalelog1pbody_mass_g"),
    list(model = m_range_cls,        predictor = "scalelog1prange_km2"),
    list(model = m_remoteness_cls,   predictor = "scalelog1pdist_min"),
    list(model = m_remoteness_cls,   predictor = "scalelog1pdist_extent"),
    list(model = m_gdp_cls,          predictor = "scalelog1pmedian_gdp_bil"),
    list(model = m_biogeography_cls, predictor = "scalemedian_lat"),
    list(model = m_year_cls,         predictor = "scaleyear")
  )

  ## biome STZ coefficients (biome1, biome2, ..., biome(k-1))
  if (exists(paste0("m_biome_", cls_lbl))) {
    m_biome_cls <- get(paste0("m_biome_", cls_lbl))
    biome_coefs <- rownames(fixef(m_biome_cls))
    biome_coefs <- biome_coefs[grepl("^biome", biome_coefs)]
    for (bc in biome_coefs) {
      model_specs <- c(model_specs, list(list(model = m_biome_cls, predictor = bc)))
    }
  }

  ## status STZ coefficients (status1, status2, ..., status(k-1))
  status_coefs <- rownames(fixef(m_status_cls))
  status_coefs <- status_coefs[grepl("^status", status_coefs)]
  for (sc in status_coefs) {
    model_specs <- c(model_specs, list(list(model = m_status_cls, predictor = sc)))
  }

  ## compute e-values for estimated coefficients
  evalue_results <- lapply(model_specs, function(x) {
    ev <- evalue_nb(x$model, x$predictor)
    ci <- ev["E-values", c("lower", "upper")]
    data.frame(
      class        = cls,
      predictor    = x$predictor,
      rr           = ev["RR", "point"],
      evalue_point = ev["E-values", "point"],
      evalue_ci    = unname(ci[!is.na(ci)][1])
    )
  })
  evalue_results <- do.call(rbind, evalue_results)

  ## Add dropped levels (last level of each sum-to-zero factor)
  # Biome dropped level
  if (exists(paste0("m_biome_", cls_lbl))) {
    ev_biome_dropped <- evalue_nb_dropped(m_biome_cls, "biome")
    ci_biome <- ev_biome_dropped["E-values", c("lower", "upper")]
    evalue_results <- rbind(evalue_results, data.frame(
      class        = cls,
      predictor    = "biome_dropped",
      rr           = ev_biome_dropped["RR", "point"],
      evalue_point = ev_biome_dropped["E-values", "point"],
      evalue_ci    = unname(ci_biome[!is.na(ci_biome)][1])
    ))
  }

  # Status dropped level
  ev_status_dropped <- evalue_nb_dropped(m_status_cls, "status")
  ci_status <- ev_status_dropped["E-values", c("lower", "upper")]
  evalue_results <- rbind(evalue_results, data.frame(
    class        = cls,
    predictor    = "status_dropped",
    rr           = ev_status_dropped["RR", "point"],
    evalue_point = ev_status_dropped["E-values", "point"],
    evalue_ci    = unname(ci_status[!is.na(ci_status)][1])
  ))

  evalue_table_all <- rbind(evalue_table_all, evalue_results)
}

evalue_table_all
write.csv(evalue_table_all, "output/class_model_evalues.csv")

## Clean up
rm(model_specs, evalue_results, m_size_cls, m_range_cls, m_remoteness_cls, m_gdp_cls,
   m_biogeography_cls, m_year_cls, m_status_cls)
if (exists("m_biome_cls")) rm(m_biome_cls)
if (exists("ev_biome_dropped")) rm(ev_biome_dropped, ci_biome)
rm(ev_status_dropped, ci_status)
gc()
