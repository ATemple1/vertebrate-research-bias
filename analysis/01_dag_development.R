## ============================================================================
## 01_dag_development.R
##
## Specifies the hypothesised causal structure linking species traits,
## geography, socioeconomics and IUCN status to research effort (a directed
## acyclic graph, DAG), and checks it against the data by testing each
## conditional independence the DAG implies. Each test compares nested linear
## mixed models with and without the focal variable (likelihood-ratio test and
## marginal partial R2), with taxonomic random intercepts.
##
## Produces: printed output only (implied independencies, likelihood-ratio
##           tests, partial R2). Histograms and the DAG plot are written to
##           Rplots.pdf under Rscript.
##
## Input data (in data/):
##   species_data_final.csv — species-level research effort, traits, geography
##                            and IUCN status (64,707 species)
## ============================================================================

library(tidyverse)
library(dagitty)
library(mgcv)
library(lme4)
library(performance)

seed <- 2009    # project-wide seed; identical in every script
set.seed(seed)
options(scipen = 999)
dir.create("output", showWarnings = FALSE)


# --- Load data ----------------------------------------------------------------

d <- read.csv("data/species_data_final.csv")
sapply(d, class)

d <- d %>%
  mutate(across(where(is.character), as.factor))

## Distribution checks. Predictors marked log1p are right-skewed and are
## log1p-transformed wherever they enter a model (dist_max is not modelled
## directly; it enters through dist_extent).
hist(d$papers) # right-skewed count; negative binomial response, untransformed
hist(d$year) # untransformed
hist(d$median_gdp_bil) # log1p
hist(d$range_km2) # log1p
hist(d$body_mass_g) # log1p
hist(d$dist_min) # log1p
hist(d$dist_max) # log1p
hist(d$dist_extent) # log1p


# --- DAG specification --------------------------------------------------------

## Abundance and charisma are latent: both plausibly matter but neither can be
## measured across ~65,000 species. DAG nodes map to data columns as follows:
##   biogeography = median_lat        gdp   = median_gdp_bil
##   remoteness   = dist_min, dist_extent
##   range        = range_km2         size  = body_mass_g
##   status, year, biome = columns of the same name; papers = outcome
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


# --- Conditional independence tests -------------------------------------------

## The DAG implies five conditional independencies (dagitty abbreviations:
## bggr = biogeography, biom = biome, stts = status, rang = range,
## rmtn = remoteness). If the data contradict one, the DAG is mis-specified.
iCI <- impliedConditionalIndependencies(dag)
print(impliedConditionalIndependencies(dag))

## Each test fits the response on its conditioning set without (null) and with
## the focal variable, compares the two by likelihood-ratio test (ML fits,
## REML = FALSE), and computes the focal variable's partial R2 from the
## marginal (fixed-effect) R2: (R2_full - R2_null) / (1 - R2_null).
## A priori failure criterion (Methods), applied when reading the printed
## results: p < 0.05 AND partial R2 > 0.05, because with
## ~65,000 species trivially small associations are statistically detectable.
## Taxonomic random intercepts (class/order/family/genus) stand in for
## phylogeny, which is not resolved for all vertebrates and would be
## computationally prohibitive at this scale; random slopes are not needed
## because only overall patterns are of interest.

## Test 1: bggr _||_ biom (median latitude independent of biome; no
##   conditioning set, so the null model has random intercepts only)
## Null model
m_CI1_lm_null <- lmer(
  median_lat ~
    1 +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d %>% filter(!is.na(biome)),
  REML = FALSE
)

## Full model: + focal variable
m_CI1_lm <- lmer(
  median_lat ~
    biome +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d %>% filter(!is.na(biome)),
  REML = FALSE
)

anova(m_CI1_lm_null, m_CI1_lm)

CT1_pR2 <- (r2_nakagawa(m_CI1_lm)$R2_marginal - r2_nakagawa(m_CI1_lm_null)$R2_marginal) / (1 - r2_nakagawa(m_CI1_lm_null)$R2_marginal)
CT1_pR2

## Test 2: bggr _||_ stts | biom, gdp, rang, rmtn, size, year
##   (median latitude independent of IUCN status given the other variables)
## Null model: conditioning set only
m_CI2_lm_null <- lmer(
  median_lat ~
    as.factor(biome) +
    scale(log1p(median_gdp_bil)) +
    scale(log1p(range_km2)) +
    scale(log1p(dist_min)) +
    scale(log1p(dist_extent)) +
    scale(log1p(body_mass_g)) +
    scale(year) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

## Full model: + focal variable
m_CI2_lm <- lmer(
  median_lat ~
    as.factor(status) +
    as.factor(biome) +
    scale(log1p(median_gdp_bil)) +
    scale(log1p(range_km2)) +
    scale(log1p(dist_min)) +
    scale(log1p(dist_extent)) +
    scale(log1p(body_mass_g)) +
    scale(year) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

anova(m_CI2_lm_null, m_CI2_lm)

CT2_pR2 <- (r2_nakagawa(m_CI2_lm)$R2_marginal - r2_nakagawa(m_CI2_lm_null)$R2_marginal) / (1 - r2_nakagawa(m_CI2_lm_null)$R2_marginal)
CT2_pR2

## Test 3: gdp _||_ rmtn | bggr, biom, rang
##   (GDP independent of remoteness; tested separately for minimum
##   depth/elevation and depth/elevation extent against a shared null)
## Null model: conditioning set only (shared by both remoteness tests)
m_CI3_lm_null <- lmer(
  log1p(median_gdp_bil) ~
    scale(median_lat) +
    as.factor(biome) +
    scale(log1p(range_km2)) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

## Full model: + minimum depth/elevation
m_CI3_lm_min <- lmer(
  log1p(median_gdp_bil) ~
    scale(log1p(dist_min)) +
    scale(median_lat) +
    as.factor(biome) +
    scale(log1p(range_km2)) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

## Full model: + depth/elevation extent
m_CI3_lm_ext <- lmer(
  log1p(median_gdp_bil) ~
    scale(log1p(dist_extent)) +
    scale(median_lat) +
    as.factor(biome) +
    scale(log1p(range_km2)) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

anova(m_CI3_lm_null, m_CI3_lm_min)
anova(m_CI3_lm_null, m_CI3_lm_ext)

CT3_pR2_min <- (r2_nakagawa(m_CI3_lm_min)$R2_marginal - r2_nakagawa(m_CI3_lm_null)$R2_marginal) / (1 - r2_nakagawa(m_CI3_lm_null)$R2_marginal)
CT3_pR2_min
CT3_pR2_ext <- (r2_nakagawa(m_CI3_lm_ext)$R2_marginal - r2_nakagawa(m_CI3_lm_null)$R2_marginal) / (1 - r2_nakagawa(m_CI3_lm_null)$R2_marginal)
CT3_pR2_ext

## Test 4: gdp _||_ size | bggr, biom, rang (GDP independent of body mass)
## Null model: conditioning set only
m_CI4_lm_null <- lmer(
  log1p(median_gdp_bil) ~
    scale(median_lat) +
    as.factor(biome) +
    scale(log1p(range_km2)) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

## Full model: + focal variable
m_CI4_lm <- lmer(
  log1p(median_gdp_bil) ~
    scale(log1p(body_mass_g)) +
    scale(median_lat) +
    as.factor(biome) +
    scale(log1p(range_km2)) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

anova(m_CI4_lm_null, m_CI4_lm)

CT2_pR4 <- (r2_nakagawa(m_CI4_lm)$R2_marginal - r2_nakagawa(m_CI4_lm_null)$R2_marginal) / (1 - r2_nakagawa(m_CI4_lm_null)$R2_marginal)
CT2_pR4

## Test 5: gdp _||_ year | bggr, biom, rang (GDP independent of year of
##   description)
## Null model: conditioning set only
m_CI5_lm_null <- lmer(
  log1p(median_gdp_bil) ~
    scale(median_lat) +
    as.factor(biome) +
    scale(log1p(range_km2)) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

## Full model: + focal variable
m_CI5_lm <- lmer(
  log1p(median_gdp_bil) ~
    scale(year) +
    scale(median_lat) +
    as.factor(biome) +
    scale(log1p(range_km2)) +
    (1 | class) +
    (1 | order) +
    (1 | family) +
    (1 | genus),
  data = d,
  REML = FALSE
)

anova(m_CI5_lm_null, m_CI5_lm)

CT2_pR5 <- (r2_nakagawa(m_CI5_lm)$R2_marginal - r2_nakagawa(m_CI5_lm_null)$R2_marginal) / (1 - r2_nakagawa(m_CI5_lm_null)$R2_marginal)
CT2_pR5
