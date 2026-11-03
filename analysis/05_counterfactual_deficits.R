## ============================================================================
## 05_counterfactual_deficits.R
##
## How much research does each bias cost? For each causal driver and class,
## posterior expected papers per species under current conditions (X) are
## compared with a counterfactual in which the species sits at the advantaged
## end of that driver (Y = X * adjustment ratio from the class model's
## coefficients). Deficit = Y - X, carried through 1,000 posterior draws and
## summed per class and across classes. Counterfactuals:
##   IUCN status — research rises from the class's own LC level by a per-step
##                 multiplier (lambda) fitted to the Aves IUCN gradient; DD is
##                 brought to LC; Aves is the reference and has no deficit
##   Biome       — every species in its class's most-researched pure biome
##   Continuous  — every species at its class's 95th (or 5th) percentile
##   Taxonomic   — every species given the Mammalia class effect (m_taxa)
## Each driver is a separate counterfactual, so deficits are not additive
## across drivers.
##
## Produces: Figure — Deficit.pdf
##           Table  — reference_levels_summary.csv, summary_class_deficits.csv,
##                    summary_all_deficits.csv, class_deficit_summary.csv,
##                    species_deficit_summary.csv,
##                    species_taxa_deficit_summary.csv
##
## Input data (in data/):
##   species_data_final.csv — species-level research effort, traits, geography
##                            and IUCN status (64,707 species)
##
## Depends on: output/m_<driver>_<class>.rds from 04_causal_models_by_class.R
##             output/m_taxa.rds from 02_taxonomic_analysis.R
## Run time:   ~<N> h; memory-heavy (posterior_epred for every species)
##                                                          [placeholder]
## Note: silhouettes are downloaded from PhyloPic (rphylopic); needs internet.
## ============================================================================

library(brms)
library(posterior)
library(tidyverse)
library(ggridges)
library(rphylopic)
library(grid)
library(patchwork)
library(cowplot)

seed <- 2009    # project-wide seed; identical in every script
set.seed(seed)
options(scipen = 999)
dir.create("output", showWarnings = FALSE)


# --- Parameters ---------------------------------------------------------------

classes <- c("ACTINOPTERYGII", "CHONDRICHTHYES", "AMPHIBIA",
             "REPTILIA", "AVES", "MAMMALIA")

# Level orderings and min_n as in 04_causal_models_by_class.R, so each model's
# retained levels and its dropped (last) sum-to-zero level can be rebuilt
global_status_order <- c("LC", "NT", "VU", "EN", "CR", "DD", "EW", "EX")
global_biome_order  <- c("terrestrial", "freshwater", "marine",
                         "terrestrial+freshwater", "terrestrial+marine",
                         "freshwater+marine")
min_n <- 3

# Subsample posterior draws to reduce memory (8000 full -> 1000 sufficient)
n_draws_use <- 1000


# --- Load data ----------------------------------------------------------------

d <- read.csv("data/species_data_final.csv")


# --- IUCN status deficit ------------------------------------------------------

## Reference curve: research should rise by a constant multiplier (lambda)
## per IUCN step above LC (LC = 0, NT = 1, VU = 2, EN = 3, CR = 4). log(lambda)
## is the slope of a log-linear fit through the Aves LC-CR IRRs, fitted within
## each posterior draw so that its uncertainty propagates. Each class is
## anchored at its own LC level, because papers already published cannot be
## removed: the deficit is how far research falls short of scaling up from
## LC. DD species are brought to LC (with threat unknown they should receive
## at least the baseline attention given to LC species); EW and EX are
## excluded. Per species and draw: X = posterior_epred() under current
## conditions, Y = X * adjustment ratio, deficit = Y - X. Negative deficits
## (a class above the reference) are kept.
## Aves defines lambda and so has no IUCN deficit: estimating one against its
## own fitted slope would be circular, and because the Aves curve is mildly
## concave the LC-anchored line overshoots its NT-EN IRRs, which would
## wrongly imply an Aves shortfall. Anchoring Aves at NT or VU instead was not
## adopted: it would anchor Aves differently from the other classes, the
## choice between NT and VU is arbitrary, and the resulting deficit is small.

## Status models from 04_causal_models_by_class.R
status_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_status_", cls_lbl, ".rds")
  cat("Loading:", fp, "\n")
  status_models[[cls]] <- readRDS(fp)
}

# Steps above LC for each status (DD handled separately)
status_steps <- tibble(
  status = c("DD", "LC", "NT", "VU", "EN", "CR"),
  steps  = c(NA, 0:4)
)

## Reference curve: per-step multiplier (lambda) from Aves. Recover all Aves
## IUCN log(IRR) draws, rebuilding the retained levels and the dropped
## sum-to-zero level exactly as in 04_causal_models_by_class.R.
d_aves <- d %>% filter(class == "AVES") %>% droplevels()
aves_status_counts <- table(d_aves$status)
aves_model_levels  <- names(aves_status_counts[aves_status_counts >= min_n])
aves_model_levels  <- intersect(global_status_order, aves_model_levels)
n_aves_levels      <- length(aves_model_levels)

cat("Aves model status levels:", paste(aves_model_levels, collapse = ", "), "\n")
cat("Dropped level:", aves_model_levels[n_aves_levels], "\n")

draws_aves <- as_draws_df(status_models[["AVES"]]) %>%
  as_tibble() %>%
  select(starts_with("b_status"))

# Recover dropped level
draws_aves <- draws_aves %>%
  mutate(b_status_dropped = -rowSums(across(everything())))
names(draws_aves) <- paste0("b_status", aves_model_levels)

## log(lambda) per posterior draw: least-squares slope of log(IRR) on rank
## (LC-CR = 1-5), computed in closed form for all draws at once
# Extract LC-CR columns only (in rank order)
lc_cr_cols <- paste0("b_status", c("LC", "NT", "VU", "EN", "CR"))
draws_lc_cr <- draws_aves[, lc_cr_cols]

# Fit log(IRR) ~ rank per draw
# rank = 1:5 (LC to CR), design matrix is constant across draws
rank_vec  <- 1:5
rank_mean <- mean(rank_vec)
rank_ss   <- sum((rank_vec - rank_mean)^2)  # = 10

# Slope per draw: sum((rank - mean) * y) / sum((rank - mean)^2)
# This is vectorised across draws (rows of the matrix)
rank_dev <- rank_vec - rank_mean  # c(-2, -1, 0, 1, 2)

log_lambda_draws_full <- as.numeric(
  (as.matrix(draws_lc_cr) %*% rank_dev) / rank_ss
)

# Subsample to n_draws_use
lambda_idx <- sort(sample(length(log_lambda_draws_full), n_draws_use))
log_lambda_draws <- log_lambda_draws_full[lambda_idx]

# Summary (from full posterior)
lambda_draws <- exp(log_lambda_draws_full)

cat(sprintf("Lambda posterior summary:\n"))
cat(sprintf("  Median: %.4f\n", median(lambda_draws)))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
            quantile(lambda_draws, 0.025), quantile(lambda_draws, 0.975)))

# Also report median IRRs for reference
aves_irrs <- draws_lc_cr %>%
  pivot_longer(everything(), names_to = "col", values_to = "log_irr") %>%
  mutate(status = sub("^b_status", "", col)) %>%
  group_by(status) %>%
  summarise(irr_aves = median(exp(log_irr)), .groups = "drop") %>%
  mutate(
    status = factor(status, levels = c("LC", "NT", "VU", "EN", "CR")),
    rank   = as.numeric(status)
  ) %>%
  arrange(rank)

cat("\nAves median IRRs (from grand mean):\n")
print(aves_irrs)

cat("\nImplied multipliers relative to LC (at median lambda):\n")
print(status_steps %>% filter(!is.na(steps)) %>%
        mutate(multiplier = median(lambda_draws)^steps))

rm(draws_aves, draws_lc_cr, d_aves, aves_status_counts,
   aves_model_levels, n_aves_levels, rank_dev, rank_ss, rank_mean,
   log_lambda_draws_full, lambda_idx)

## Deficit for one class: rebuild all IUCN log(IRR) draws from the
## sum-to-zero coefficients; build the LC-anchored reference curve per draw
## from that draw's lambda; compute the adjustment ratio per species and draw;
## X = posterior_epred(), Y = X * ratio, deficit = Y - X; then per-draw totals
## by class and by class x status for the density plots.
compute_iucn_deficit <- function(model, species_data, class_name,
                                 log_lambda_draws, status_steps,
                                 global_status_order, min_n,
                                 n_draws_use = 1000) {

  cat("\n--- IUCN deficit:", class_name, "---\n")

  ## 1. Detect which status levels were in the model (as in 04)
  status_counts     <- table(species_data$status)
  cls_model_levels  <- names(status_counts[status_counts >= min_n])
  cls_model_levels  <- intersect(global_status_order, cls_model_levels)
  n_model_levels    <- length(cls_model_levels)

  cat("Model status levels:", paste(cls_model_levels, collapse = ", "), "\n")
  cat("Dropped level:", cls_model_levels[n_model_levels], "\n")

  ## 2. Extract posterior draws for all status coefficients (subsampled)
  draws <- as_draws_df(model) %>% as_tibble()
  n_draws_full <- nrow(draws)
  draws_idx <- sort(sample(n_draws_full, min(n_draws_use, n_draws_full)))
  n_draws <- length(draws_idx)

  draws_status <- draws[draws_idx, ] %>% select(starts_with("b_status"))

  # Recover dropped level (last in class ordering)
  draws_status <- draws_status %>%
    mutate(b_status_dropped = -rowSums(across(everything())))
  names(draws_status) <- paste0("b_status", cls_model_levels)

  ## 3. LC anchor draws (log scale)
  log_irr_LC <- draws_status[["b_statusLC"]]

  ## 4. Summarise class-specific vs reference IRRs for reporting (DD-CR only)
  class_irrs <- draws_status %>%
    pivot_longer(everything(), names_to = "col", values_to = "log_irr") %>%
    mutate(status = sub("^b_status", "", col)) %>%
    filter(status %in% c("LC", "NT", "VU", "EN", "CR", "DD")) %>%
    group_by(status) %>%
    summarise(
      irr_class    = median(exp(log_irr)),
      irr_class_lo = quantile(exp(log_irr), 0.025),
      irr_class_hi = quantile(exp(log_irr), 0.975),
      .groups = "drop"
    ) %>%
    left_join(status_steps, by = "status") %>%
    mutate(
      irr_ref_anchored = median(exp(log_irr_LC)) * median(exp(log_lambda_draws))^steps,
      class = class_name
    )

  cat("Class vs reference IRRs:\n")
  print(class_irrs)

  ## 5. Get model training data, exclude EW/EX, map status to steps
  d_model <- model$data
  d_model$.row <- seq_len(nrow(d_model))

  d_model <- d_model %>%
    filter(!status %in% c("EW", "EX")) %>%
    left_join(status_steps, by = "status")

  cat("Species for deficit calculation:", nrow(d_model), "\n")

  ## 6. Compute adjustment ratio per species per posterior draw
  ##    For LC: steps = 0, ref = LC -> ratio = 1 (exactly, every draw)
  ##    For DD: brought to LC (minimum baseline research attention)
  ##    No floor — negatives retained for plotting flexibility

  n_sp <- nrow(d_model)
  adj_matrix <- matrix(1, nrow = n_sp, ncol = n_draws)

  for (i in seq_len(n_sp)) {
    sp_status <- as.character(d_model$status[i])
    log_irr_class_s <- draws_status[[paste0("b_status", sp_status)]]

    if (sp_status == "DD") {
      ## DD -> LC: lacking information to assess threat, DD species should
      ## receive at minimum the baseline research attention given to LC
      adj_matrix[i, ] <- exp(log_irr_LC - log_irr_class_s)
    } else {
      ## log_lambda_draws varies per draw — reference curve has uncertainty
      log_irr_ref_s <- log_irr_LC + d_model$steps[i] * log_lambda_draws
      adj_matrix[i, ] <- exp(log_irr_ref_s - log_irr_class_s)
    }
  }

  ## 7. Posterior expected counts under current conditions
  cat("Computing posterior_epred (this may take a moment)...\n")
  X <- posterior_epred(model, subset = draws_idx)  # n_draws x n_model_species
  X <- X[, d_model$.row, drop = FALSE]              # subset to retained species

  ## 8. Counterfactual and deficit per draw
  ##    Negative deficit = class exceeds reference at that status
  X_t <- t(X)                                  # n_species x n_draws
  Y_t <- X_t * adj_matrix                      # counterfactual papers
  deficit_draws <- Y_t - X_t                   # negative = exceeds reference

  ## 9. Summarise per species
  d_model <- d_model %>%
    mutate(
      predicted_papers    = apply(X_t, 1, median),
      predicted_papers_lo = apply(X_t, 1, quantile, 0.025),
      predicted_papers_hi = apply(X_t, 1, quantile, 0.975),
      adj_median = apply(adj_matrix, 1, median),
      adj_lo     = apply(adj_matrix, 1, quantile, 0.025),
      adj_hi     = apply(adj_matrix, 1, quantile, 0.975),
      counterfactual_papers    = apply(Y_t, 1, median),
      counterfactual_papers_lo = apply(Y_t, 1, quantile, 0.025),
      counterfactual_papers_hi = apply(Y_t, 1, quantile, 0.975),
      iucn_deficit      = apply(deficit_draws, 1, median),
      iucn_deficit_lo   = apply(deficit_draws, 1, quantile, 0.025),
      iucn_deficit_hi   = apply(deficit_draws, 1, quantile, 0.975),
      iucn_deficit_mean = apply(deficit_draws, 1, mean),
      prob_deficit      = apply(deficit_draws, 1, function(x) mean(x > 0)),
      class             = class_name
    )

  ## 10. Aggregate posterior distributions for density plotting
  ##     Sum deficit across species per draw -> full posterior of total deficit
  status_vec <- as.character(d_model$status)

  # Per class: total deficit distribution (one value per draw)
  class_dist <- tibble(
    class = class_name,
    draw  = seq_len(n_draws),
    total_deficit = colSums(deficit_draws)
  )

  # Per class x status: total and mean deficit distributions
  class_status_dist <- map_dfr(
    unique(status_vec),
    function(s) {
      idx <- which(status_vec == s)
      tibble(
        class  = class_name,
        status = s,
        draw   = seq_len(n_draws),
        total_deficit = colSums(deficit_draws[idx, , drop = FALSE]),
        mean_deficit  = total_deficit / length(idx),
        n_species     = length(idx)
      )
    }
  )

  ## 11. Flatten deficit draws with draw index for downstream subsampling
  flat_draws <- tibble(
    class = class_name,
    draw  = rep(seq_len(n_draws), each = n_sp),
    deficit_value = as.vector(deficit_draws)
  )

  list(species = d_model, class_irrs = class_irrs,
       class_dist = class_dist, class_status_dist = class_status_dist,
       flat_draws = flat_draws)
}

## Run for every class (the Aves IUCN deficit is computed but dropped later,
## as Aves is the reference)
irr_comparison   <- data.frame()
class_dists      <- data.frame()
cls_status_dists <- data.frame()
iucn_flat_draws  <- list()

for (cls in classes) {

  d_cls <- d %>% filter(class == cls) %>% droplevels()

  res <- compute_iucn_deficit(
    model               = status_models[[cls]],
    species_data        = d_cls,
    class_name          = cls,
    log_lambda_draws    = log_lambda_draws,
    status_steps        = status_steps,
    global_status_order = global_status_order,
    min_n               = min_n,
    n_draws_use         = n_draws_use
  )

  irr_comparison       <- bind_rows(irr_comparison, res$class_irrs)
  class_dists          <- bind_rows(class_dists, res$class_dist)
  cls_status_dists     <- bind_rows(cls_status_dists, res$class_status_dist)
  iucn_flat_draws[[cls]] <- res$flat_draws

  gc()
}

iucn_flat_draws <- bind_rows(iucn_flat_draws)
rm(status_models)
gc()

## Diagnostics

## Class IRRs vs the LC-anchored reference curve
irr_comparison %>%
  filter(status %in% c("LC", "NT", "VU", "EN", "CR")) %>%
  mutate(status = factor(status, levels = c("LC", "NT", "VU", "EN", "CR"))) %>%
  ggplot(aes(x = status, y = irr_class, colour = class, group = class)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = irr_class_lo, ymax = irr_class_hi),
                width = 0.15, alpha = 0.4) +
  geom_line(aes(y = irr_ref_anchored), linetype = "dashed", alpha = 0.6) +
  scale_y_log10() +
  facet_wrap(~class, scales = "free_y") +
  labs(x = "IUCN Status", y = "IRR (from grand mean)",
       title = "Class-specific IUCN effects vs. LC-anchored reference",
       subtitle = "Dashed = reference curve (class LC x lambda^steps)") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")

## Posterior density of each class's total IUCN deficit
class_dists %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total IUCN research deficit",
       subtitle = "Sum across all species per posterior draw") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")


# --- Biome deficit ------------------------------------------------------------

## Counterfactual: every species in its class's most-researched pure biome
## (terrestrial, freshwater or marine, highest median IRR). Species already in
## that biome have zero deficit. Multi-biome levels are not used as the anchor
## (small n and inflated IRRs) and get a deficit only if less researched than
## the anchor (ratio floored at 1).

## Biome models from 04_causal_models_by_class.R (a class without one is skipped)
biome_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_biome_", cls_lbl, ".rds")
  if (file.exists(fp)) {
    cat("Loading:", fp, "\n")
    biome_models[[cls]] <- readRDS(fp)
  } else {
    cat("WARNING: no biome model for", cls, "\n")
  }
}

## Deficit function
compute_biome_deficit <- function(model, species_data, class_name,
                                   global_biome_order, min_n,
                                   n_draws_use = 1000) {

  cat("\n--- Biome deficit:", class_name, "---\n")

  ## 1. Detect biome levels (as in 04)
  biome_counts     <- table(species_data$biome)
  cls_biome_levels <- names(biome_counts[biome_counts >= min_n])
  cls_biome_levels <- intersect(global_biome_order, cls_biome_levels)
  n_biome_levels   <- length(cls_biome_levels)

  cat("Biome levels:", paste(cls_biome_levels, collapse = ", "), "\n")
  cat("Dropped level:", cls_biome_levels[n_biome_levels], "\n")

  ## 2. Extract and subsample coefficient draws
  draws <- as_draws_df(model) %>% as_tibble()
  n_draws_full <- nrow(draws)
  draws_idx <- sort(sample(n_draws_full, min(n_draws_use, n_draws_full)))
  n_draws <- length(draws_idx)

  draws_biome <- draws[draws_idx, ] %>% select(starts_with("b_biome"))
  draws_biome <- draws_biome %>%
    mutate(b_biome_dropped = -rowSums(across(everything())))
  names(draws_biome) <- paste0("b_biome", cls_biome_levels)

  ## 3. Find anchor biome — restrict to pure biomes only
  ##    Combined biomes (e.g. freshwater+marine) often have small n and
  ##    high IRR driven by sparse data. Restricting the anchor to
  ##    {terrestrial, freshwater, marine} avoids inflated counterfactuals.
  pure_biomes <- c("terrestrial", "freshwater", "marine")
  pure_in_model <- intersect(pure_biomes, cls_biome_levels)

  if (length(pure_in_model) == 0) {
    stop("No pure biomes (terrestrial/freshwater/marine) present in ", class_name)
  }

  pure_cols <- paste0("b_biome", pure_in_model)
  median_log_irrs_pure <- sapply(draws_biome[, pure_cols, drop = FALSE], median)
  best_idx   <- which.max(median_log_irrs_pure)
  best_biome <- pure_in_model[best_idx]
  log_irr_best <- draws_biome[[paste0("b_biome", best_biome)]]

  cat("Anchor biome:", best_biome,
      sprintf("(median IRR = %.3f) [restricted to pure biomes]\n",
              exp(median(log_irr_best))))

  ## 4. Summarise class-specific biome IRRs for reporting
  class_biome_irrs <- draws_biome %>%
    pivot_longer(everything(), names_to = "col", values_to = "log_irr") %>%
    mutate(biome = sub("^b_biome", "", col)) %>%
    group_by(biome) %>%
    summarise(
      irr_class    = median(exp(log_irr)),
      irr_class_lo = quantile(exp(log_irr), 0.025),
      irr_class_hi = quantile(exp(log_irr), 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      anchor_biome = best_biome,
      irr_anchor   = exp(median(log_irr_best)),
      class        = class_name
    )

  cat("Biome IRRs:\n")
  print(class_biome_irrs)

  ## 5. Get model data — all species in the model
  d_model <- model$data
  d_model$.row <- seq_len(nrow(d_model))

  n_sp <- nrow(d_model)
  cat("Species for deficit calculation:", n_sp, "\n")

  ## 6. Compute adjustment ratio per species per draw
  ##    Species in anchor biome -> adj = 1, deficit = 0
  ##    Combined-biome species: adj floored at 1 (like above-anchor
  ##    species in continuous variables — no negative deficit)
  adj_matrix <- matrix(1, nrow = n_sp, ncol = n_draws)

  for (i in seq_len(n_sp)) {
    sp_biome <- as.character(d_model$biome[i])
    if (sp_biome == best_biome) next
    log_irr_current <- draws_biome[[paste0("b_biome", sp_biome)]]

    if (sp_biome %in% pure_biomes) {
      ## Pure biome below anchor — standard adjustment
      adj_matrix[i, ] <- exp(log_irr_best - log_irr_current)
    } else {
      ## Combined biome — floor adj at 1: deficit only if less
      ## researched than the pure-biome anchor, otherwise no deficit
      adj_matrix[i, ] <- pmax(exp(log_irr_best - log_irr_current), 1)
    }
  }

  ## 7. Posterior expected counts under current conditions
  cat("Computing posterior_epred...\n")
  X <- posterior_epred(model, subset = draws_idx)
  X_t <- t(X)
  Y_t <- X_t * adj_matrix
  deficit_draws <- Y_t - X_t

  ## 8. Summarise per species
  d_model <- d_model %>%
    mutate(
      predicted_papers    = apply(X_t, 1, median),
      predicted_papers_lo = apply(X_t, 1, quantile, 0.025),
      predicted_papers_hi = apply(X_t, 1, quantile, 0.975),
      adj_median = apply(adj_matrix, 1, median),
      adj_lo     = apply(adj_matrix, 1, quantile, 0.025),
      adj_hi     = apply(adj_matrix, 1, quantile, 0.975),
      counterfactual_papers    = apply(Y_t, 1, median),
      counterfactual_papers_lo = apply(Y_t, 1, quantile, 0.025),
      counterfactual_papers_hi = apply(Y_t, 1, quantile, 0.975),
      biome_deficit      = apply(deficit_draws, 1, median),
      biome_deficit_lo   = apply(deficit_draws, 1, quantile, 0.025),
      biome_deficit_hi   = apply(deficit_draws, 1, quantile, 0.975),
      biome_deficit_mean = apply(deficit_draws, 1, mean),
      prob_deficit       = apply(deficit_draws, 1, function(x) mean(x > 0)),
      class              = class_name
    )

  ## 9. Aggregate posterior distributions for density plotting
  biome_vec <- as.character(d_model$biome)

  class_dist <- tibble(
    class = class_name,
    draw  = seq_len(n_draws),
    total_deficit = colSums(deficit_draws)
  )

  class_biome_dist <- map_dfr(
    unique(biome_vec),
    function(b) {
      idx <- which(biome_vec == b)
      tibble(
        class  = class_name,
        biome  = b,
        draw   = seq_len(n_draws),
        total_deficit = colSums(deficit_draws[idx, , drop = FALSE]),
        mean_deficit  = total_deficit / length(idx),
        n_species     = length(idx)
      )
    }
  )

  ## 10. Flatten deficit draws with draw index for downstream subsampling
  flat_draws <- tibble(
    class = class_name,
    draw  = rep(seq_len(n_draws), each = n_sp),
    deficit_value = as.vector(deficit_draws)
  )

  list(species = d_model, class_biome_irrs = class_biome_irrs,
       class_dist = class_dist, class_biome_dist = class_biome_dist,
       flat_draws = flat_draws)
}

## Run for every class
biome_irr_comparison   <- data.frame()
biome_class_dists      <- data.frame()
biome_cls_biome_dists  <- data.frame()
biome_flat_draws       <- list()

for (cls in names(biome_models)) {

  d_cls <- d %>% filter(class == cls) %>% droplevels()

  res <- compute_biome_deficit(
    model              = biome_models[[cls]],
    species_data       = d_cls,
    class_name         = cls,
    global_biome_order = global_biome_order,
    min_n              = min_n,
    n_draws_use        = n_draws_use
  )

  biome_irr_comparison       <- bind_rows(biome_irr_comparison, res$class_biome_irrs)
  biome_class_dists          <- bind_rows(biome_class_dists, res$class_dist)
  biome_cls_biome_dists      <- bind_rows(biome_cls_biome_dists, res$class_biome_dist)
  biome_flat_draws[[cls]]    <- res$flat_draws

  gc()
}

biome_flat_draws <- bind_rows(biome_flat_draws)
rm(biome_models)
gc()

## Diagnostics

## Biome IRRs by class
biome_irr_comparison %>%
  ggplot(aes(x = reorder(biome, irr_class), y = irr_class, colour = class)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = irr_class_lo, ymax = irr_class_hi),
                width = 0.2, alpha = 0.5) +
  coord_flip() +
  facet_wrap(~class, scales = "free_x") +
  labs(x = NULL, y = "IRR (from grand mean)",
       title = "Biome effects on research effort by class",
       subtitle = "Anchor biome (highest IRR) defines the reference") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")

## Posterior density of each class's total biome deficit
biome_class_dists %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total biome research deficit") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")


# --- Continuous drivers: shared deficit function ------------------------------

## Counterfactual: every species at the most-advantaged value of the driver in
## its class: the class's 95th percentile if more of the driver means more
## papers, or its 5th percentile if less does. Direction comes from the sign
## of the class's posterior median coefficient, so drivers acting in opposite
## directions in different classes are handled automatically. Species at or
## beyond the anchor keep a ratio of 1 (no deficit); for the others,
## ratio = exp(beta * (anchor - x)), with x on the model's standardised scale.
compute_continuous_deficit <- function(model, class_name,
                                        coef_name, var_name, var_label,
                                        transform_fn = identity,
                                        n_draws_use = 1000,
                                        anchor_percentile = 0.95) {

  cat("\n---", var_label, "deficit:", class_name, "---\n")

  ## 1. Extract and subsample coefficient draws
  draws <- as_draws_df(model) %>% as_tibble()
  n_draws_full <- nrow(draws)
  draws_idx <- sort(sample(n_draws_full, min(n_draws_use, n_draws_full)))
  n_draws <- length(draws_idx)

  beta_draws <- draws[[paste0("b_", coef_name)]][draws_idx]

  ## 2. Determine bias direction from median coefficient
  beta_median <- median(beta_draws)
  bias_positive <- beta_median > 0
  direction <- ifelse(bias_positive, "positive", "negative")

  cat(sprintf("Coefficient median: %.4f (%s bias)\n", beta_median, direction))

  ## 3. Compute scaled predictor values (replicating what brms does)
  d_model <- model$data
  d_model$.row <- seq_len(nrow(d_model))

  raw_vals    <- d_model[[var_name]]
  transformed <- transform_fn(raw_vals)
  t_mean      <- mean(transformed)
  t_sd        <- sd(transformed)
  scaled_vals <- (transformed - t_mean) / t_sd

  ## 4. Set anchor at class-specific percentile
  if (bias_positive) {
    anchor_val <- quantile(scaled_vals, anchor_percentile)
    disadvantaged <- scaled_vals < anchor_val
    cat(sprintf("Anchor: %.2fth percentile (scaled = %.3f)\n",
                anchor_percentile * 100, anchor_val))
  } else {
    anchor_val <- quantile(scaled_vals, 1 - anchor_percentile)
    disadvantaged <- scaled_vals > anchor_val
    cat(sprintf("Anchor: %.2fth percentile (scaled = %.3f)\n",
                (1 - anchor_percentile) * 100, anchor_val))
  }

  n_sp <- nrow(d_model)
  n_disadvantaged <- sum(disadvantaged)
  cat(sprintf("Species: %d total, %d disadvantaged\n", n_sp, n_disadvantaged))

  ## 5. Compute adjustment ratio per species per draw
  ##    adj = exp(beta * (anchor - x_i)) for disadvantaged species
  ##    Species at or past anchor keep adj = 1
  adj_matrix <- matrix(1, nrow = n_sp, ncol = n_draws)

  for (i in which(disadvantaged)) {
    delta <- anchor_val - scaled_vals[i]
    adj_matrix[i, ] <- exp(beta_draws * delta)
  }

  ## 6. Posterior expected counts
  cat("Computing posterior_epred...\n")
  X <- posterior_epred(model, subset = draws_idx)
  X_t <- t(X)
  Y_t <- X_t * adj_matrix
  deficit_draws <- Y_t - X_t

  ## 7. Summarise per species
  deficit_col      <- paste0(var_label, "_deficit")
  deficit_col_lo   <- paste0(var_label, "_deficit_lo")
  deficit_col_hi   <- paste0(var_label, "_deficit_hi")
  deficit_col_mean <- paste0(var_label, "_deficit_mean")

  d_model <- d_model %>%
    mutate(
      predicted_papers    = apply(X_t, 1, median),
      predicted_papers_lo = apply(X_t, 1, quantile, 0.025),
      predicted_papers_hi = apply(X_t, 1, quantile, 0.975),
      adj_median = apply(adj_matrix, 1, median),
      adj_lo     = apply(adj_matrix, 1, quantile, 0.025),
      adj_hi     = apply(adj_matrix, 1, quantile, 0.975),
      counterfactual_papers    = apply(Y_t, 1, median),
      counterfactual_papers_lo = apply(Y_t, 1, quantile, 0.025),
      counterfactual_papers_hi = apply(Y_t, 1, quantile, 0.975),
      deficit      = apply(deficit_draws, 1, median),
      deficit_lo   = apply(deficit_draws, 1, quantile, 0.025),
      deficit_hi   = apply(deficit_draws, 1, quantile, 0.975),
      deficit_mean = apply(deficit_draws, 1, mean),
      prob_deficit = apply(deficit_draws, 1, function(x) mean(x > 0)),
      class        = class_name,
      bias_direction = direction,
      anchor_percentile_used = ifelse(bias_positive,
                                       anchor_percentile,
                                       1 - anchor_percentile)
    )

  ## 8. Aggregate posterior distributions for density plotting
  class_dist <- tibble(
    class = class_name,
    draw  = seq_len(n_draws),
    total_deficit = colSums(deficit_draws)
  )

  ## 9. Summary of anchor and coefficient for reporting
  class_summary <- tibble(
    class     = class_name,
    variable  = var_label,
    coef_median = beta_median,
    coef_lo   = quantile(beta_draws, 0.025),
    coef_hi   = quantile(beta_draws, 0.975),
    direction = direction,
    anchor_percentile = ifelse(bias_positive,
                                anchor_percentile,
                                1 - anchor_percentile),
    anchor_raw = ifelse(bias_positive,
                         quantile(raw_vals, anchor_percentile),
                         quantile(raw_vals, 1 - anchor_percentile)),
    n_species = n_sp,
    n_disadvantaged = n_disadvantaged
  )

  ## 10. Flatten deficit draws with draw index for downstream subsampling
  flat_draws <- tibble(
    class = class_name,
    draw  = rep(seq_len(n_draws), each = n_sp),
    deficit_value = as.vector(deficit_draws)
  )

  list(species = d_model, class_dist = class_dist,
       class_summary = class_summary, flat_draws = flat_draws)
}


# --- Biogeography deficit -----------------------------------------------------

## Biogeography models from 04_causal_models_by_class.R
biogeography_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_biogeography_", cls_lbl, ".rds")
  cat("Loading:", fp, "\n")
  biogeography_models[[cls]] <- readRDS(fp)
}

## Run for every class
biogeo_class_dists  <- data.frame()
biogeo_summaries    <- data.frame()
biogeo_flat_draws   <- list()

for (cls in classes) {

  res <- compute_continuous_deficit(
    model              = biogeography_models[[cls]],
    class_name         = cls,
    coef_name          = "scalemedian_lat",
    var_name           = "median_lat",
    var_label          = "biogeography",
    transform_fn       = identity,
    n_draws_use        = n_draws_use,
    anchor_percentile  = 0.95
  )

  biogeo_class_dists         <- bind_rows(biogeo_class_dists, res$class_dist)
  biogeo_summaries           <- bind_rows(biogeo_summaries, res$class_summary)
  biogeo_flat_draws[[cls]]   <- res$flat_draws

  gc()
}

biogeo_flat_draws <- bind_rows(biogeo_flat_draws)
rm(biogeography_models)
gc()

## Diagnostics

## Coefficient direction and anchor value per class
biogeo_summaries %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

## Posterior density of each class's total biogeography deficit
biogeo_class_dists %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total biogeography research deficit") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")


# --- GDP deficit --------------------------------------------------------------

## GDP models from 04_causal_models_by_class.R
gdp_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_gdp_", cls_lbl, ".rds")
  cat("Loading:", fp, "\n")
  gdp_models[[cls]] <- readRDS(fp)
}

## Run for every class
gdp_class_dists  <- data.frame()
gdp_summaries    <- data.frame()
gdp_flat_draws   <- list()

for (cls in classes) {

  res <- compute_continuous_deficit(
    model              = gdp_models[[cls]],
    class_name         = cls,
    coef_name          = "scalelog1pmedian_gdp_bil",
    var_name           = "median_gdp_bil",
    var_label          = "gdp",
    transform_fn       = log1p,
    n_draws_use        = n_draws_use,
    anchor_percentile  = 0.95
  )

  gdp_class_dists         <- bind_rows(gdp_class_dists, res$class_dist)
  gdp_summaries           <- bind_rows(gdp_summaries, res$class_summary)
  gdp_flat_draws[[cls]]   <- res$flat_draws

  gc()
}

gdp_flat_draws <- bind_rows(gdp_flat_draws)
rm(gdp_models)
gc()

## Diagnostics

## Coefficient direction and anchor value per class
gdp_summaries %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

## Posterior density of each class's total GDP deficit
gdp_class_dists %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total GDP research deficit") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")


# --- Range deficit ------------------------------------------------------------

## Range models from 04_causal_models_by_class.R
range_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_range_", cls_lbl, ".rds")
  cat("Loading:", fp, "\n")
  range_models[[cls]] <- readRDS(fp)
}

## Run for every class
range_class_dists  <- data.frame()
range_summaries    <- data.frame()
range_flat_draws   <- list()

for (cls in classes) {

  res <- compute_continuous_deficit(
    model              = range_models[[cls]],
    class_name         = cls,
    coef_name          = "scalelog1prange_km2",
    var_name           = "range_km2",
    var_label          = "range",
    transform_fn       = log1p,
    n_draws_use        = n_draws_use,
    anchor_percentile  = 0.95
  )

  range_class_dists         <- bind_rows(range_class_dists, res$class_dist)
  range_summaries           <- bind_rows(range_summaries, res$class_summary)
  range_flat_draws[[cls]]   <- res$flat_draws

  gc()
}

range_flat_draws <- bind_rows(range_flat_draws)
rm(range_models)
gc()

## Diagnostics

## Coefficient direction and anchor value per class
range_summaries %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

## Posterior density of each class's total range deficit
range_class_dists %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total range research deficit") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")


# --- Remoteness deficit -------------------------------------------------------

## dist_min and dist_extent come from the same model; each deficit holds the
## other at its observed value

## Remoteness models from 04_causal_models_by_class.R
remoteness_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_remoteness_", cls_lbl, ".rds")
  cat("Loading:", fp, "\n")
  remoteness_models[[cls]] <- readRDS(fp)
}

## Run for every class: minimum depth/elevation (dist_min)
distmin_class_dists  <- data.frame()
distmin_summaries    <- data.frame()
distmin_flat_draws   <- list()

for (cls in classes) {

  res <- compute_continuous_deficit(
    model              = remoteness_models[[cls]],
    class_name         = cls,
    coef_name          = "scalelog1pdist_min",
    var_name           = "dist_min",
    var_label          = "dist_min",
    transform_fn       = log1p,
    n_draws_use        = n_draws_use,
    anchor_percentile  = 0.95
  )

  distmin_class_dists         <- bind_rows(distmin_class_dists, res$class_dist)
  distmin_summaries           <- bind_rows(distmin_summaries, res$class_summary)
  distmin_flat_draws[[cls]]   <- res$flat_draws

  gc()
}

distmin_flat_draws <- bind_rows(distmin_flat_draws)
gc()

## Run for every class: depth/elevation extent (dist_extent)
distext_class_dists  <- data.frame()
distext_summaries    <- data.frame()
distext_flat_draws   <- list()

for (cls in classes) {

  res <- compute_continuous_deficit(
    model              = remoteness_models[[cls]],
    class_name         = cls,
    coef_name          = "scalelog1pdist_extent",
    var_name           = "dist_extent",
    var_label          = "dist_extent",
    transform_fn       = log1p,
    n_draws_use        = n_draws_use,
    anchor_percentile  = 0.95
  )

  distext_class_dists         <- bind_rows(distext_class_dists, res$class_dist)
  distext_summaries           <- bind_rows(distext_summaries, res$class_summary)
  distext_flat_draws[[cls]]   <- res$flat_draws

  gc()
}

distext_flat_draws <- bind_rows(distext_flat_draws)
rm(remoteness_models)
gc()

## Diagnostics

## Coefficient direction and anchor value per class
cat("--- dist_min ---\n")
distmin_summaries %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

cat("\n--- dist_extent ---\n")
distext_summaries %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

## Posterior density of each class's total remoteness deficit
bind_rows(
  distmin_class_dists %>% mutate(variable = "dist_min"),
  distext_class_dists %>% mutate(variable = "dist_extent")
) %>%
  ggplot(aes(x = total_deficit, fill = variable, colour = variable)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of remoteness research deficit",
       subtitle = "dist_min and dist_extent computed independently") +
  theme_classic() +
  theme(text = element_text(size = 8))


# --- Size deficit -------------------------------------------------------------

## Size models from 04_causal_models_by_class.R
size_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_size_", cls_lbl, ".rds")
  cat("Loading:", fp, "\n")
  size_models[[cls]] <- readRDS(fp)
}

## Run for every class
size_class_dists  <- data.frame()
size_summaries    <- data.frame()
size_flat_draws   <- list()

for (cls in classes) {

  res <- compute_continuous_deficit(
    model              = size_models[[cls]],
    class_name         = cls,
    coef_name          = "scalelog1pbody_mass_g",
    var_name           = "body_mass_g",
    var_label          = "size",
    transform_fn       = log1p,
    n_draws_use        = n_draws_use,
    anchor_percentile  = 0.95
  )

  size_class_dists         <- bind_rows(size_class_dists, res$class_dist)
  size_summaries           <- bind_rows(size_summaries, res$class_summary)
  size_flat_draws[[cls]]   <- res$flat_draws

  gc()
}

size_flat_draws <- bind_rows(size_flat_draws)
rm(size_models)
gc()

## Diagnostics

## Coefficient direction and anchor value per class
size_summaries %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

## Posterior density of each class's total size deficit
size_class_dists %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total body size research deficit") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")


# --- Year deficit -------------------------------------------------------------

## Year models from 04_causal_models_by_class.R
year_models <- list()
for (cls in classes) {
  cls_lbl <- tolower(cls)
  fp <- paste0("output/m_year_", cls_lbl, ".rds")
  cat("Loading:", fp, "\n")
  year_models[[cls]] <- readRDS(fp)
}

## Run for every class
year_class_dists  <- data.frame()
year_summaries    <- data.frame()
year_flat_draws   <- list()

for (cls in classes) {

  res <- compute_continuous_deficit(
    model              = year_models[[cls]],
    class_name         = cls,
    coef_name          = "scaleyear",
    var_name           = "year",
    var_label          = "year",
    transform_fn       = identity,
    n_draws_use        = n_draws_use,
    anchor_percentile  = 0.95
  )

  year_class_dists         <- bind_rows(year_class_dists, res$class_dist)
  year_summaries           <- bind_rows(year_summaries, res$class_summary)
  year_flat_draws[[cls]]   <- res$flat_draws

  gc()
}

year_flat_draws <- bind_rows(year_flat_draws)
rm(year_models)
gc()

## Diagnostics

## Coefficient direction and anchor value per class
year_summaries %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()

## Posterior density of each class's total year deficit
year_class_dists %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total year-of-description research deficit") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")


# --- Reference (anchor) levels ------------------------------------------------

## Biome: anchor biome per class (from biome_irr_comparison)
biome_refs <- biome_irr_comparison %>%
  distinct(class, anchor_biome) %>%
  transmute(
    class     = class,
    variable  = "biome",
    direction = NA_character_,
    anchor_percentile = NA_real_,
    anchor_value      = anchor_biome,
    anchor_type       = "categorical (highest IRR)"
  )

## Continuous variables: anchor value per class (from individual summary tables)
continuous_refs <- bind_rows(
  biogeo_summaries, gdp_summaries, range_summaries,
  distmin_summaries, distext_summaries,
  size_summaries, year_summaries
) %>%
  transmute(
    class     = class,
    variable  = variable,
    direction = direction,
    anchor_percentile = anchor_percentile,
    anchor_value      = as.character(round(anchor_raw, 4)),
    anchor_type       = paste0(anchor_percentile * 100, "th percentile")
  )

reference_table <- bind_rows(biome_refs, continuous_refs) %>%
  arrange(variable, class)

write.csv(reference_table, "output/reference_levels_summary.csv", row.names = FALSE)


# --- Deficit totals by driver and class ---------------------------------------

## Totals are summed within each posterior draw, then summarised (median,
## mean, 50% and 95% CrIs), so class medians need not add up to the
## all-class median. Aves x IUCN status is excluded (reference class).

## Helper: collapse one flat_draws frame to draw-level totals per class
summarise_flat <- function(flat_draws, var_label) {
  flat_draws %>%
    group_by(class, draw) %>%
    summarise(draw_total = sum(deficit_value), .groups = "drop") %>%
    mutate(variable = var_label)
}

## Summarise each variable separately (~6000 rows each), then bind
draw_totals <- bind_rows(
  summarise_flat(iucn_flat_draws,    "IUCN Status"),
  summarise_flat(biome_flat_draws,   "Biome"),
  summarise_flat(biogeo_flat_draws,  "Biogeography"),
  summarise_flat(gdp_flat_draws,     "GDP"),
  summarise_flat(range_flat_draws,   "Range"),
  summarise_flat(distmin_flat_draws, "Remoteness (min)"),
  summarise_flat(distext_flat_draws, "Remoteness (extent)"),
  summarise_flat(size_flat_draws,    "Size"),
  summarise_flat(year_flat_draws,    "Year")
) %>%
  filter(!(class == "AVES" & variable == "IUCN Status"))

## Total deficit per class per bias
summary_class_variable <- draw_totals %>%
  group_by(variable, class) %>%
  summarise(
    total_deficit_median = median(draw_total),
    total_deficit_mean = mean(draw_total),
    total_deficit_lo = quantile(draw_total, 0.025),
    total_deficit_lo50 = quantile(draw_total, 0.25),
    total_deficit_hi50 = quantile(draw_total, 0.75),
    total_deficit_hi = quantile(draw_total, 0.975),
    .groups = "drop"
  ) %>%
  arrange(variable, class)

## Total deficit across all classes per bias
summary_by_variable <- draw_totals %>%
  group_by(variable, draw) %>%
  summarise(draw_total = sum(draw_total), .groups = "drop_last") %>%
  summarise(
    total_deficit_median = median(draw_total),
    total_deficit_mean = mean(draw_total),
    total_deficit_lo = quantile(draw_total, 0.025),
    total_deficit_lo50 = quantile(draw_total, 0.25),
    total_deficit_hi50 = quantile(draw_total, 0.75),
    total_deficit_hi = quantile(draw_total, 0.975),
    .groups = "drop"
  ) %>%
  arrange(desc(total_deficit_median))

rm(draw_totals)
gc()

## Save
write.csv(summary_class_variable, "output/summary_class_deficits.csv", row.names = FALSE)
write.csv(summary_by_variable, "output/summary_all_deficits.csv", row.names = FALSE)


# --- Taxonomic (class) deficit ------------------------------------------------

## Research effort differs between classes independently of species traits.
## Using the taxonomic model from 02 (papers ~ 1 + (1|class) + (1|order) +
## (1|family) + (1|genus)), the counterfactual gives every species the
## Mammalia class effect while keeping its own order, family and genus
## effects. The class effect is additive on the log scale, so the ratio
## exp(RE_class[Mammalia] - RE_class[c]) is constant within a class and the
## Mammalia deficit is zero by definition.

## Taxonomic model from 02_taxonomic_analysis.R
taxa_model <- readRDS("output/m_taxa.rds")

compute_class_deficit <- function(model, anchor_class = "MAMMALIA",
                                   n_draws_use = 1000) {

  cat("\n--- Class deficit (anchor:", anchor_class, ") ---\n")

  ## 1. Extract class-level RE draws only (not the order, family or genus REs)
  draws <- as_draws_df(model) %>% as_tibble()
  n_draws_full <- nrow(draws)
  draws_idx <- sort(sample(n_draws_full, min(n_draws_use, n_draws_full)))
  n_draws <- length(draws_idx)

  class_re_cols <- grep("^r_class\\[", names(draws), value = TRUE)
  class_re_draws <- draws[draws_idx, class_re_cols]

  # Parse class names from column names
  class_names <- gsub("r_class\\[(.+),Intercept\\]", "\\1", class_re_cols)
  names(class_re_draws) <- class_names

  cat("Classes found:", paste(class_names, collapse = ", "), "\n")

  ## 2. Anchor class
  median_res <- sapply(class_re_draws, median)
  cat("Class RE medians (IRR scale):\n")
  print(round(exp(sort(median_res, decreasing = TRUE)), 3))

  anchor_re <- class_re_draws[[anchor_class]]
  cat("Anchor:", anchor_class,
      sprintf("(median IRR = %.3f)\n", exp(median(anchor_re))))

  ## 3. Class RE summary for reporting
  class_re_summary <- class_re_draws %>%
    pivot_longer(everything(), names_to = "class", values_to = "re") %>%
    group_by(class) %>%
    summarise(
      irr       = median(exp(re)),
      irr_lo    = quantile(exp(re), 0.025),
      irr_hi    = quantile(exp(re), 0.975),
      .groups   = "drop"
    ) %>%
    mutate(anchor = anchor_class)

  cat("Class IRRs:\n")
  print(class_re_summary)

  ## 4. Model data and adjustment matrix
  d_model <- model$data
  d_model$.row <- seq_len(nrow(d_model))
  n_sp <- nrow(d_model)
  cat("Species for deficit calculation:", n_sp, "\n")

  adj_matrix <- matrix(1, nrow = n_sp, ncol = n_draws)

  for (cls in setdiff(class_names, anchor_class)) {
    idx <- which(d_model$class == cls)
    if (length(idx) == 0) next
    cls_re <- class_re_draws[[cls]]
    adj_vec <- exp(anchor_re - cls_re)  # length n_draws
    adj_matrix[idx, ] <- matrix(adj_vec, nrow = length(idx),
                                 ncol = n_draws, byrow = TRUE)
  }

  ## 5. Posterior expected counts
  cat("Computing posterior_epred (full model, this may take a moment)...\n")
  X <- posterior_epred(model, subset = draws_idx)  # n_draws x n_sp
  X_t <- t(X)                                      # n_sp x n_draws
  Y_t <- X_t * adj_matrix
  deficit_draws <- Y_t - X_t

  ## 6. Flatten for plotting
  flat_draws <- tibble(
    class = rep(as.character(d_model$class), times = n_draws),
    draw  = rep(seq_len(n_draws), each = n_sp),
    deficit_value = as.vector(deficit_draws)
  )

  ## 7. Per-class total deficit distribution (one value per draw per class)
  class_dist <- map_dfr(unique(as.character(d_model$class)), function(cls) {
    idx <- which(d_model$class == cls)
    tibble(
      class = cls,
      draw  = seq_len(n_draws),
      total_deficit = colSums(deficit_draws[idx, , drop = FALSE])
    )
  })

  rm(X, X_t, Y_t, adj_matrix, deficit_draws)
  gc()

  list(flat_draws = flat_draws, class_dist = class_dist,
       class_re_summary = class_re_summary)
}

## Run for all classes, then keep the six focal classes
taxa_res <- compute_class_deficit(
  model        = taxa_model,
  anchor_class = "MAMMALIA",
  n_draws_use  = n_draws_use
)

taxa_flat_draws    <- taxa_res$flat_draws
taxa_class_dist    <- taxa_res$class_dist
taxa_re_summary    <- taxa_res$class_re_summary

taxa_flat_draws <- taxa_res$flat_draws %>%
  filter(class %in% classes)

taxa_class_dist <- taxa_res$class_dist %>%
  filter(class %in% classes)

rm(taxa_model, taxa_res)
gc()

## Diagnostics

## Class random-effect IRRs
taxa_re_summary %>%
  mutate(class = fct_reorder(class, irr)) %>%
  ggplot(aes(x = class, y = irr)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = irr_lo, ymax = irr_hi), width = 0.2) +
  coord_flip() +
  labs(x = NULL, y = "IRR (from grand mean)",
       title = "Class-level random effects on research effort",
       subtitle = paste("Anchor:", unique(taxa_re_summary$anchor))) +
  theme_classic() +
  theme(text = element_text(size = 10))

## Posterior density of each class's total taxonomic deficit
class_order <- c("MAMMALIA", "AVES", "REPTILIA",
                 "AMPHIBIA", "ACTINOPTERYGII", "CHONDRICHTHYES")

taxa_class_dist %>%
  mutate(class = factor(class, levels = class_order)) %>%
  ggplot(aes(x = total_deficit, fill = class, colour = class)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  facet_wrap(~class, scales = "free") +
  labs(x = "Total deficit (papers)", y = "Posterior density",
       title = "Posterior distribution of total class-level research deficit",
       subtitle = "Counterfactual: all species at Mammalia-level class effect") +
  theme_classic() +
  theme(text = element_text(size = 8), legend.position = "none")

## Summary

## Per-class totals summarised across draws; "Total" sums the classes
## within each draw first
taxa_summary <- taxa_class_dist %>%
  group_by(class) %>%
  summarise(
    total_deficit_median = median(total_deficit),
    total_deficit_mean = mean(total_deficit),
    total_deficit_lo = quantile(total_deficit, 0.025),
    total_deficit_lo50 = quantile(total_deficit, 0.25),
    total_deficit_hi50 = quantile(total_deficit, 0.75),
    total_deficit_hi = quantile(total_deficit, 0.975),
    .groups = "drop"
  ) %>%
  arrange(desc(total_deficit_mean))

total_row <- taxa_class_dist %>%
  group_by(draw) %>%
  summarise(total_deficit = sum(total_deficit), .groups = "drop") %>%
  summarise(
    class = "Total",
    total_deficit_median = median(total_deficit),
    total_deficit_mean = mean(total_deficit),
    total_deficit_lo = quantile(total_deficit, 0.025),
    total_deficit_lo50 = quantile(total_deficit, 0.25),
    total_deficit_hi50 = quantile(total_deficit, 0.75),
    total_deficit_hi = quantile(total_deficit, 0.975)
  )

taxa_summary <- bind_rows(taxa_summary, total_row)

write.csv(taxa_summary, "output/class_deficit_summary.csv", row.names = FALSE)

rm(taxa_class_dist)
gc()


# --- Figure: per-species deficits (Deficit.pdf) -------------------------------

## 1. Subsample 100 of 1000 draws per species for plotting
set.seed(seed)
plot_draw_idx <- sort(sample(n_draws_use, 100))

all_flat_draws <- bind_rows(
  iucn_flat_draws    %>% mutate(variable = "IUCN Status"),
  biome_flat_draws   %>% mutate(variable = "Biome"),
  biogeo_flat_draws  %>% mutate(variable = "Biogeography"),
  gdp_flat_draws     %>% mutate(variable = "GDP"),
  range_flat_draws   %>% mutate(variable = "Range"),
  distmin_flat_draws %>% mutate(variable = "Remoteness (min)"),
  distext_flat_draws %>% mutate(variable = "Remoteness (extent)"),
  size_flat_draws    %>% mutate(variable = "Size"),
  year_flat_draws    %>% mutate(variable = "Year")
) %>%
  filter(draw %in% plot_draw_idx) %>%
  select(-draw)

## Free per-variable flat draws (summary data section already ran)
rm(iucn_flat_draws, biome_flat_draws, biogeo_flat_draws, gdp_flat_draws,
   range_flat_draws, distmin_flat_draws, distext_flat_draws,
   size_flat_draws, year_flat_draws)
gc()

## 2. Remove Aves x IUCN Status (reference group)
all_flat_draws <- all_flat_draws %>%
  filter(!(class == "AVES" & variable == "IUCN Status"))

## Display order of drivers and classes
var_order <- rev(c("Year", "Remoteness (extent)", "Biogeography", "GDP",
               "Remoteness (min)", "Range",  "Size",
               "Biome", "IUCN Status"))

class_order <- c("MAMMALIA", "AVES", "REPTILIA",
                 "AMPHIBIA", "ACTINOPTERYGII", "CHONDRICHTHYES")

all_flat_draws <- all_flat_draws %>%
  mutate(
    variable = factor(variable, levels = var_order),
    class    = factor(class, levels = class_order)
  )

## Panel A: each ridge pools species heterogeneity and posterior uncertainty.
## Only positive deficits are plotted (log1p axis).
p01 <-
  ggplot(all_flat_draws %>% filter(deficit_value > 0), aes(x = deficit_value, y = variable, fill = class)) +
  geom_density_ridges(
    scale = 2, rel_min_height = 0.01,
    alpha = 1, bandwidth = 0.5,
    colour = "black", linewidth = .5) +
  geom_text(data = data.frame(class = factor("AVES", levels = levels(all_flat_draws$class)),
                              variable = "IUCN Status", deficit_value = 0.5),
          aes(x = deficit_value, y = variable, label = "Reference Level"),
          inherit.aes = FALSE, size = 2, hjust = 0, vjust = -0.5, colour = "grey30") +
  scale_fill_manual(values = c("#D08020", "#F0B840", "#785018", "#508070", "#60C8C0", "#28A0B0")) +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 1, 3, 10, 30, 100, 300, 1000),
    limits = c(0, 2000),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = expansion(add = c(0, 2.1)),
                   labels = c("IUCN status", "Biome", "Body mass", "Geographic\nrange", "Elevation/depth\nminimum",
                              "Gross domestic\nproduct", "Median\nlatitude", "Elevation/depth\nextent", "Year of\ndescription")) +
  facet_wrap(~ class, nrow = 1, scales = "free_x", labeller = labeller(class = str_to_sentence)) +
  labs(x = "Publication deficit per species", y = NULL) +
  theme_classic() +
  theme(
    text = element_text(size = 6),
    strip.text  = element_text(size = 6*1.2, face = "bold", hjust = 0, vjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
    )
p01

## Per-species deficit by class and driver (median, mean, 50% and 95%
## intervals), from the 100 plotting draws and positive deficits only
species_deficit_summary <- all_flat_draws %>%
  filter(deficit_value > 0) %>%
  group_by(class, variable) %>%
  summarise(sp_deficit_med = median(deficit_value),
            sp_deficit_mean = mean(deficit_value),
            sp_deficit_lo = quantile(deficit_value, 0.025),
            sp_deficit_lo50 = quantile(deficit_value, 0.25),
            sp_deficit_hi50 = quantile(deficit_value, 0.75),
            sp_deficit_hi = quantile(deficit_value, 0.975))

write.csv(species_deficit_summary, "output/species_deficit_summary.csv", row.names = FALSE)

## Panel B: per-species taxonomic deficit relative to Mammalia, from 100
## posterior draws

## Subsample 100 draws for plotting
set.seed(seed)
plot_draw_idx_taxa <- sort(sample(n_draws_use, 100))

taxa_plot_draws <- taxa_flat_draws %>%
  filter(draw %in% plot_draw_idx_taxa) %>%
  select(-draw) %>%
  mutate(class = factor(class, levels = class_order))

## Plot
p02 <-
  taxa_plot_draws %>%
  filter(class != "MAMMALIA") %>%
  ggplot(aes(x = deficit_value, colour = class), fill = NA) +
  geom_density(bw = 0.5, linewidth = 2) +
    scale_colour_manual(values = c(
    "AVES"           = "#F0B840",
    "REPTILIA"       = "#785018",
    "AMPHIBIA"       = "#508070",
    "ACTINOPTERYGII" = "#60C8C0",
    "CHONDRICHTHYES" = "#28A0B0"
  )) +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 1, 3, 10, 30, 100, 300, 1000),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Publication deficit per species relative to Mammalia",
       y = "Relative density") +
  coord_cartesian(xlim = c(0, 2000)) +
  theme_classic() +
  theme(text = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")
p02

## As above for the taxonomic deficit (100 plotting draws, positive only)
species_taxa_deficit_summary <- taxa_plot_draws %>%
  filter(deficit_value > 0) %>%
  group_by(class) %>%
  summarise(sp_deficit_med = median(deficit_value),
            sp_deficit_mean = mean(deficit_value),
            sp_deficit_lo = quantile(deficit_value, 0.025),
            sp_deficit_lo50 = quantile(deficit_value, 0.25),
            sp_deficit_hi50 = quantile(deficit_value, 0.75),
            sp_deficit_hi = quantile(deficit_value, 0.975))

write.csv(species_taxa_deficit_summary, "output/species_taxa_deficit_summary.csv", row.names = FALSE)

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

## Assemble: A = per-driver deficits, B = taxonomic deficit
pCD <- ggdraw(free(p01, type = "label") / free(p02, type = "label")) +
  ## panel A silhouettes
  draw_grob(p_mamm, x = .20, y = .9, width = .05, height = .05) +
  draw_grob(p_aves, x = .345, y = .9, width = .05, height = .05) +
  draw_grob(p_rept, x = .49, y = .9, width = .05, height = .05) +
  draw_grob(p_amph, x = .635, y = .9, width = .04, height = .04) +
  draw_grob(p_acti, x = .78, y = .89, width = .05, height = .05) +
  draw_grob(p_chon, x = .925, y = .89, width = .06, height = .05) +
  ## panel B silhouettes
  draw_grob(p_rept, x = .12, y = .42, width = .05, height = .05) +
  draw_grob(p_aves, x = .29, y = .42, width = .05, height = .05) +
  draw_grob(p_amph, x = .38, y = .42, width = .04, height = .04) +
  draw_grob(p_chon, x = .55, y = .35, width = .06, height = .05) +
  draw_grob(p_acti, x = .83, y = .16, width = .05, height = .05) +
  draw_plot_label(label = c("A", "B"), size = 10,
                  x = c(.01, .01), y = c(.99, .49))
pCD

ggsave("output/Deficit.pdf", pCD, width = 18.4, height = 12, units = "cm", dpi = 1200, bg = "white")
