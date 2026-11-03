## ============================================================================
## 03_causal_models_global.R
##
## What drives differences in research effort among species? Each driver in
## the DAG (01_dag_development.R) gets its own negative binomial model of
## per-species paper counts across all vertebrates. Each model adjusts only for
## that driver's backdoor adjustment set and includes taxonomic random
## intercepts. Continuous effects are incidence rate ratios (IRRs) per 2 SD;
## biome and IUCN status are sum-to-zero contrasts against the grand mean,
## with all pairwise contrasts. E-values give robustness to unmeasured
## confounding.
##
## Produces: Table  — global_model_effectsize.csv (used by
##                    04_causal_models_by_class.R), pairwise_biome.csv,
##                    pairwise_status.csv, global_model_evalues.csv
##           Object — m_range.rds, m_remoteness.rds, m_gdp.rds, m_size.rds,
##                    m_year.rds, m_biogeography.rds, m_biome.rds,
##                    m_status.rds (fitted models; not read by later scripts)
##
## Input data (in data/):
##   species_data_final.csv — species-level research effort, traits, geography
##                            and IUCN status (64,707 species)
##
## Run time:   ~<N> h for the eight brms fits              [placeholder]
## ============================================================================

library(tidyverse)
library(dagitty)
library(brms)
library(cmdstanr)
library(patchwork)
library(rstan)
library(EValue)

seed <- 2009    # project-wide seed; identical in every script
set.seed(seed)
options(scipen = 999)
dir.create("output", showWarnings = FALSE)


# --- Load data ----------------------------------------------------------------

d <- read.csv("data/species_data_final.csv")
sapply(d, class)


# --- DAG ----------------------------------------------------------------------

## As specified in 01_dag_development.R. adjustmentSets() below prints the
## minimal backdoor adjustment set for each driver; each model's covariates
## are that set.
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

## As in 02_taxonomic_analysis.R, plus normal(0, 1) on slopes.
priors <- c(
  prior(normal(0.69, 2.5),  class = "Intercept"), # as in 02_taxonomic_analysis.R
  prior(normal(0, 1),       class = "b"),         # log scale; IRR ~1/7-7 per 1 SD (1/50-50 per 2 SD);
                                                  # same scale for sum-to-zero IUCN and biome contrasts
  prior(student_t(3, 0, 1), class = "sd"),        # half-t; same prior for class, order, family, genus
                                                  # (no belief about relative importance encoded)
  prior(exponential(1),     class = "shape")      # NB phi
)


# --- Effect-size scale --------------------------------------------------------

## Effects are reported as incidence rate ratios (IRR = exp(beta)), the
## multiplicative change in expected paper count. Continuous predictors are
## standardised (1 SD) inside each model formula, and their coefficients are
## doubled at extraction (exp(2 * beta)), so each IRR is the effect of a 2-SD
## increase (Gelman 2008). Biome and IUCN status use sum-to-zero contrasts, so
## each level's IRR is relative to the grand mean across levels (the
## unweighted average of the level effects), not to a reference level; the
## last level, which the contrasts drop, is recovered as minus the sum of the
## others. Pairwise level-vs-level contrasts are derived from the same
## posterior draws. Effect sizes are summarised as posterior medians with 50%
## and 95% credible intervals and the probability of direction (pd); pairwise
## contrasts as medians with 95% intervals and pd. All IRRs are read
## against 1; continuous and categorical effects share the multiplicative
## scale, but their magnitudes are only approximately comparable.

## One row per effect, accumulated across the models below
std_eff <- data.frame()


# --- Geographic range ---------------------------------------------------------

## Range = range_km2, area of the IUCN/BirdLife range polygons

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "range", outcome = "papers")

## model
m_range <- brm(
  papers ~
    scale(log1p(range_km2)) +
    scale(median_lat) +
    biome +
    scale(log1p(dist_min)) +
    scale(log1p(dist_extent)) +
    scale(log1p(body_mass_g)) +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_range)
## Convergence
plot(m_range, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_range, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_range, "output/m_range.rds")

## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
std_eff <-
  rbind(std_eff,
        as_draws_df(m_range) %>%
          select(b = b_scalelog1prange_km2) %>%
          mutate(b2 = b * 2, irr = exp(b2)) %>%
          summarise(
            Parameter = "range_km2",
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


# --- Remoteness ---------------------------------------------------------------

## Remoteness = minimum elevation/depth (dist_min) and the extent between
## minimum and maximum (dist_extent), from IUCN listings. A dist_min x
## dist_extent interaction was credible but negligible and not favoured by
## LOO, so the additive model is used.

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "remoteness", outcome = "papers")

## model
m_remoteness <- brm(
  papers ~
    scale(log1p(dist_min)) +
    scale(log1p(dist_extent)) +
    scale(median_lat) +
    biome +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_remoteness)
## Convergence
plot(m_remoteness, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_remoteness, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_remoteness, "output/m_remoteness.rds")

## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
## min
std_eff <-
  rbind(std_eff,
        as_draws_df(m_remoteness) %>%
          select(b = b_scalelog1pdist_min) %>%
          mutate(b2 = b * 2, irr = exp(b2)) %>%
          summarise(
            Parameter = "dist_min",
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
        as_draws_df(m_remoteness) %>%
          select(b = b_scalelog1pdist_extent) %>%
          mutate(b2 = b * 2, irr = exp(b2)) %>%
          summarise(
            Parameter = "dist_extent",
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


# --- GDP ----------------------------------------------------------------------

## GDP = median_gdp_bil, summed GDP (median of the last 10 years) of the
## countries whose land or EEZ overlaps the species' range

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "gdp", outcome = "papers")

## model
m_gdp <- brm(
  papers ~
    scale(log1p(median_gdp_bil)) +
    scale(median_lat) +
    biome +
    scale(log1p(range_km2)) +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_gdp)
## Convergence
plot(m_gdp, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_gdp, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_gdp, "output/m_gdp.rds")

## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
std_eff <-
  rbind(std_eff,
        as_draws_df(m_gdp) %>%
          select(b = b_scalelog1pmedian_gdp_bil) %>%
          mutate(b2 = b * 2, irr = exp(b2)) %>%
          summarise(
            Parameter = "median_gdp_bil",
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


# --- Body size ----------------------------------------------------------------

## Size = body_mass_g, from species trait databases

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "size", outcome = "papers")

## model
m_size <- brm(
  papers ~
    scale(log1p(body_mass_g)) +
    scale(median_lat) +
    biome +
    scale(log1p(dist_min)) +
    scale(log1p(dist_extent)) +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_size)
## Convergence
plot(m_size, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_size, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_size, "output/m_size.rds")

## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
std_eff <-
  rbind(std_eff,
        as_draws_df(m_size) %>%
          select(b = b_scalelog1pbody_mass_g) %>%
          mutate(b2 = b * 2, irr = exp(b2)) %>%
          summarise(
            Parameter = "body_mass_g",
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


# --- Year of description ------------------------------------------------------

## Year = year of description, from the species' taxonomic authority

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "year", outcome = "papers")

## model
m_year <- brm(
  papers ~
    scale(year) +
    scale(median_lat) +
    biome +
    scale(log1p(range_km2)) +
    scale(log1p(dist_min)) +
    scale(log1p(dist_extent)) +
    scale(log1p(body_mass_g)) +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_year)
## Convergence
plot(m_year, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_year, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_year, "output/m_year.rds")

## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
std_eff <-
  rbind(std_eff,
        as_draws_df(m_year) %>%
          select(b = b_scaleyear) %>%
          mutate(b2 = b * 2, irr = exp(b2)) %>%
          summarise(
            Parameter = "year",
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


# --- Biogeography -------------------------------------------------------------

## Biogeography = median_lat (signed), capturing the Global North-South axis

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "biogeography", outcome = "papers")

## model
m_biogeography <- brm(
  papers ~
    scale(median_lat) +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_biogeography)
## Convergence
plot(m_biogeography, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_biogeography, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_biogeography, "output/m_biogeography.rds")

## IRR per 2 SD: posterior median, 50% and 95% CrIs, pd
std_eff <-
  rbind(std_eff,
        as_draws_df(m_biogeography) %>%
          select(b = b_scalemedian_lat) %>%
          mutate(b2 = b * 2, irr = exp(b2)) %>%
          summarise(
            Parameter = "median_lat",
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


# --- Biome --------------------------------------------------------------------

## Biomes are not mutually exclusive, so species in more than one biome form
## their own level (e.g. freshwater+marine) of a single factor, rather than
## being coded as overlapping binary indicators

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "biome", outcome = "papers")

## 1. Set up contrasts
# Biome: sum-to-zero (exposure of interest)
# freshwater+marine is last, so it is the level recovered from the draws
d$biome <- factor(d$biome, levels = c("terrestrial", "freshwater", "marine",
                                       "terrestrial+freshwater", "terrestrial+marine",
                                       "freshwater+marine"))
contrasts(d$biome) <- contr.sum(nlevels(d$biome))

## 2. Check what R will estimate
contrasts(d$biome)
biome_levels <- levels(d$biome)
n_biome <- length(biome_levels)

## 3. Fit model
m_biome <- brm(
  papers ~
    biome +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_biome)
## Convergence
plot(m_biome, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_biome, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_biome, "output/m_biome.rds")

## 4. Extract all biome effects including the dropped level
draws_raw_biome <- as_draws_df(m_biome) %>%
  as_tibble() %>%
  select(starts_with("b_biome"))

draws_all_biome <- draws_raw_biome %>%
  mutate(
    b_biome_dropped = -rowSums(across(everything()))
  )

estimated_names_biome <- paste0("b_biome", biome_levels[1:(n_biome - 1)])
all_names_biome <- c(estimated_names_biome, paste0("b_biome", biome_levels[n_biome]))
names(draws_all_biome) <- all_names_biome

## 5. IRR of every level vs the grand mean
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
pairwise_biome <- combn(biome_levels, 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    a <- pair[1]
    b <- pair[2]
    diff <- draws_all_biome[[paste0("b_biome", a)]] - draws_all_biome[[paste0("b_biome", b)]]
    irr_diff <- exp(diff)
    tibble(
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

## clean up
rm(biome_levels, n_biome, draws_raw_biome, draws_all_biome, estimated_names_biome, all_names_biome)
gc()


# --- IUCN status --------------------------------------------------------------

## Status = IUCN Red List category (DD, LC, NT, VU, EN, CR, EW, EX). Biome,
## a covariate here, keeps the sum-to-zero contrasts set in the biome section
## above, so run the script top to bottom.

## Minimal backdoor adjustment set (the model's covariates)
adjustmentSets(dag, exposure = "status", outcome = "papers")

## 1. Set up contrasts
# Status: sum-to-zero (exposure of interest)
# Put EX last — least central to the narrative,
# and it gets recovered from posterior draws anyway
d$status <- factor(d$status, levels = c("DD", "LC", "NT", "VU", "EN", "CR", "EW", "EX"))
contrasts(d$status) <- contr.sum(nlevels(d$status))

## 2. Check what R will estimate
# Columns 1–7 get coefficients; EX (row 8) is the negative sum
contrasts(d$status)

# Store level names for later labelling
status_levels <- levels(d$status)
n_status <- length(status_levels)

## 3. Fit model
m_status <- brm(
  papers ~
    status +
    biome +
    scale(log1p(median_gdp_bil)) +
    scale(log1p(range_km2)) +
    scale(log1p(dist_min)) +
    scale(log1p(dist_extent)) +
    scale(log1p(body_mass_g)) +
    scale(year) +
    (1|class) + (1|order) + (1|family) + (1|genus),
  data = d,
  family = negbinomial(),
  prior = priors,
  chains = 4, cores = 4, iter = 4000, warmup = 2000, threads = threading(3),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed = seed,
  backend = "cmdstanr"
)

summary(m_status)
## Convergence
plot(m_status, variable = c("^b_", "^sd_", "^shape$"), regex = TRUE)  # trace plots
## Posterior predictive checks
pp_check(m_status, type = "dens_overlay", ndraws = 50) + scale_x_log10()   # whole distribution
saveRDS(m_status, "output/m_status.rds")

## 4. Extract all status effects including the dropped level
draws_raw_status <- as_draws_df(m_status) %>%
  as_tibble() %>%
  select(starts_with("b_status"))

# brms names sum-to-zero columns as status1, status2, ..., status(k-1)
# These correspond to levels 1:(k-1) in the factor ordering
# The dropped level (EX, level k) = negative sum of all others
draws_all_status <- draws_raw_status %>%
  mutate(
    b_status_dropped = -rowSums(across(everything()))
  )

# Rename columns to actual level names
estimated_names_status <- paste0("b_status", status_levels[1:(n_status - 1)])
all_names_status <- c(estimated_names_status, paste0("b_status", status_levels[n_status]))
names(draws_all_status) <- all_names_status

## 5. IRR of every level vs the grand mean
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

## 6. Pairwise contrasts between all status levels
pairwise_status <- combn(status_levels, 2, simplify = FALSE) %>%
  map_dfr(function(pair) {
    a <- pair[1]
    b <- pair[2]
    diff <- draws_all_status[[paste0("b_status", a)]] - draws_all_status[[paste0("b_status", b)]]
    irr_diff <- exp(diff)
    tibble(
      contrast = paste(a, "vs", b),
      diff_meidan = median(diff),
      diff_l95 = quantile(diff, 0.025),
      diff_u95 = quantile(diff, 0.975),
      IRR_ratio = median(irr_diff),
      IRR_ratio_l95 = quantile(irr_diff, 0.025),
      IRR_ratio_u95 = quantile(irr_diff, 0.975),
      pd = max(mean(diff > 0), mean(diff < 0)),
      sig_95 = !between(0, quantile(diff, 0.025), quantile(diff, 0.975))
    )
  })

## clean up
rm(status_levels, n_status, draws_raw_status, draws_all_status, estimated_names_status, all_names_status)
gc()


# --- Global effect sizes plot (diagnostic; not saved) -------------------------

## Rows: continuous drivers in a fixed order (low to high IRR in the reported
## fit; update lev_order after a refit), then biome, then status
lev_order <- c(
  "year", "dist_min", "median_lat", "dist_extent", "median_gdp_bil",  "range_km2", "body_mass_g", # continuous rows
  "terrestrial vs grand mean", "freshwater vs grand mean", "marine vs grand mean", "terrestrial+freshwater vs grand mean", "terrestrial+marine vs grand mean", "freshwater+marine vs grand mean", # biome
  "DD vs grand mean", "LC vs grand mean", "NT vs grand mean", "VU vs grand mean", "EN vs grand mean", "CR vs grand mean", "EW vs grand mean", "EX vs grand mean" # status
  )

p01 <-
  ggplot(std_eff, aes(y = Parameter, x = IRR)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = .5, colour = "grey75") +
  geom_linerange(aes(xmin = IRR_lower95, xmax = IRR_upper95),
                 linewidth = 0.5, alpha = 0.8, colour = "#9890A0") +
  geom_linerange(aes(xmin = IRR_lower50, xmax = IRR_upper50),
                 linewidth = 1, colour = "#9890A0") +
  geom_point(size = 1, colour = "black") +
  geom_hline(yintercept = 14.5, linetype = 3, linewidth = .5, colour = "grey75") +
  geom_hline(yintercept = 8.5, linetype = 3, linewidth = .5, colour = "grey75") +
  annotate("text", x = Inf, y = 21.5, label = "Numeric ",
           hjust = 1, vjust = 1, colour = "black", size = 5/.pt) +
  annotate("text", x = Inf, y = 14.5, label = "Biome \n(vs grand mean) ",
           hjust = 1, vjust = 2, colour = "black", size = 5/.pt) +
  annotate("text", x = Inf, y = 8.5, label = "IUCN Status \n(vs grand mean) ",
           hjust = 1, vjust = 2, colour = "black", size = 5/.pt) +
  scale_x_log10() +
  scale_y_discrete(limits = rev(lev_order),
                   labels = c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD",
                              "Freshwater-\nMarine", "Terrestrial-\nMarine", "Terrestrial-\nFreshwater", "Marine", "Freshwater", "Terrestrial",
                              "Body Mass", "Geographic\n Range", "Gross Domestic\nProduct", "Elevation/Depth\nExtent", "Median\nLatitude",
                              "Elevation/Depth\nMinimum", "Year of\nDescription" ))  + # rev() so lev_order reads top-to-bottom
  labs(x = "IRR") +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.title.y = element_blank())
p01


# --- Save effect sizes --------------------------------------------------------

write.csv(std_eff, "output/global_model_effectsize.csv")
write.csv(pairwise_biome, "output/pairwise_biome.csv")
write.csv(pairwise_status, "output/pairwise_status.csv")


# --- Sensitivity to unmeasured confounding (E-values) -------------------------

## The E-value is the minimum strength of association (risk-ratio scale) an
## unmeasured confounder would need with both the driver and research effort
## to fully explain away an effect (VanderWeele & Ding 2017); larger = more
## robust. IRRs are treated as risk ratios, continuous effects per 2 SD, from
## posterior medians and 95% CrIs. evalue_ci is the E-value for the interval
## limit nearest the null. For each sum-to-zero factor, the dropped level is
## recovered from the draws (evalue_nb_dropped).
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

## Apply across models
model_specs <- list(
  list(model = m_size, predictor = "scalelog1pbody_mass_g"),
  list(model = m_range, predictor = "scalelog1prange_km2"),
  list(model = m_remoteness, predictor = "scalelog1pdist_min"),
  list(model = m_remoteness, predictor = "scalelog1pdist_extent"),
  list(model = m_gdp, predictor = "scalelog1pmedian_gdp_bil"),
  list(model = m_biogeography, predictor = "scalemedian_lat"),
  list(model = m_year, predictor = "scaleyear"),
  ## biome
  list(model = m_biome, predictor = "biome1"), # Terrestrial
  list(model = m_biome, predictor = "biome2"), # Freshwater
  list(model = m_biome, predictor = "biome3"), # Marine
  list(model = m_biome, predictor = "biome4"), # Terrestrial-Freshwater
  list(model = m_biome, predictor = "biome5"), # Terrestrial-Marine
  ## status
  list(model = m_status, predictor = "status1"), # DD
  list(model = m_status, predictor = "status2"), # LC
  list(model = m_status, predictor = "status3"), # NT
  list(model = m_status, predictor = "status4"), # VU
  list(model = m_status, predictor = "status5"), # EN
  list(model = m_status, predictor = "status6"), # CR
  list(model = m_status, predictor = "status7") # EW
)

evalue_results <- lapply(model_specs, function(x) {
  ev <- evalue_nb(x$model, x$predictor)
  ci <- ev["E-values", c("lower", "upper")]
  data.frame(
    predictor    = x$predictor,
    rr           = ev["RR", "point"],
    evalue_point = ev["E-values", "point"],
    evalue_ci    = unname(ci[!is.na(ci)][1])
  )
})
evalue_results <- do.call(rbind, evalue_results)

## Add dropped levels (last level of each sum-to-zero factor)
# Biome dropped level
ev_biome_dropped <- evalue_nb_dropped(m_biome, "biome")
ci_biome <- ev_biome_dropped["E-values", c("lower", "upper")]
evalue_results <- rbind(evalue_results, data.frame(
  predictor    = "biome_dropped",
  rr           = ev_biome_dropped["RR", "point"],
  evalue_point = ev_biome_dropped["E-values", "point"],
  evalue_ci    = unname(ci_biome[!is.na(ci_biome)][1])
))

# Status dropped level
ev_status_dropped <- evalue_nb_dropped(m_status, "status")
ci_status <- ev_status_dropped["E-values", c("lower", "upper")]
evalue_results <- rbind(evalue_results, data.frame(
  predictor    = "status_dropped",
  rr           = ev_status_dropped["RR", "point"],
  evalue_point = ev_status_dropped["E-values", "point"],
  evalue_ci    = unname(ci_status[!is.na(ci_status)][1])
))

evalue_results
write.csv(evalue_results, "output/global_model_evalues.csv")

## Clean up
rm(model_specs, ev_biome_dropped, ci_biome, ev_status_dropped, ci_status)
gc()
