# ============================================================
# 11_seasonal_catch_allocation_diagnostics.R
# Seasonal catch-allocation diagnostics for fishery conditioning
#
# Stock / analysis:
#   European anchovy ane.27.9aS (Gulf of Cádiz)
#   Shortcut MSE in FLBEIA + FLasher
#
# Role in the workflow:
#   10_build_fleet_reference.R
#           ↓
#   11_seasonal_catch_allocation_diagnostics.R
#           ↓
#   reference fishery-conditioning documentation
#   and candidate seasonal-allocation sensitivity scenarios
#
# Purpose:
#   Diagnose the historical seasonal allocation of catches represented
#   by the 100 bootstrap-conditioned OMs, verify the reference seasonal
#   shares used for projection, describe historical catch-at-age by
#   quarter, and construct candidate seasonal-allocation scenarios for
#   subsequent sensitivity analyses.
#
# Inputs:
#   data/Rdata/ane_stock_conditioned_100iter.rds
#   data/Rdata/fleets_conditioned_1000iter.rds
#   data/Rdata/fleets_ctrl_reference_1000iter.rds
#   data/Rdata/iteration_map_1000.csv
#
# Outputs for reference fishery conditioning:
#   outputs/fishery_conditioning/tables/
#   outputs/fishery_conditioning/figures/
#
# Outputs for sensitivity diagnostics:
#   outputs/fishery_sensitivity/tables/
#   outputs/fishery_sensitivity/figures/
#
# Important:
#   - Historical diagnostics are based on the 100 unique conditioned OMs.
#   - The 10 replicated projection trajectories per OM are not treated as
#     independent historical OMs.
#   - The reference seasonal allocation is the 2015-2024 mean implemented
#     by script 10.
#   - PCA and clustering are diagnostics for candidate sensitivity
#     scenarios; they do not redefine the reference OM.
#   - k = 3 is retained from the original diagnostic script and is not
#     interpreted here as a validated methodological decision.
# ============================================================

rm(list = ls())

# ============================================================
# 0. PACKAGES
# ============================================================

library(FLCore)
library(FLBEIA)
library(FLFishery)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(here)
library(strucchange)
library(changepoint)

# ============================================================
# 1. SETTINGS AND PATHS
# ============================================================

rds_dir <- file.path(here(), "data", "Rdata")

ref_fig_dir <- file.path(here(), "outputs", "fishery_conditioning", "figures")
ref_tab_dir <- file.path(here(), "outputs", "fishery_conditioning", "tables")
sens_fig_dir <- file.path(here(), "outputs", "fishery_sensitivity", "figures")
sens_tab_dir <- file.path(here(), "outputs", "fishery_sensitivity", "tables")

dir.create(ref_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ref_tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sens_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sens_tab_dir, recursive = TRUE, showWarnings = FALSE)

first.yr <- 1989
ass.yr <- 2024
proj.yrs <- 2025:2054
hist_years <- first.yr:ass.yr
n_om <- 100

# ============================================================
# 2. LOAD CURRENT FISHERY-CONDITIONING OBJECTS
# ============================================================

ane.stock.100 <- readRDS(file.path(rds_dir, "ane_stock_conditioned_100iter.rds"))
fleets.1000 <- readRDS(file.path(rds_dir, "fleets_conditioned_1000iter.rds"))
fleets.ctrl.1000 <- readRDS(file.path(rds_dir, "fleets_ctrl_reference_1000iter.rds"))
iteration_map_1000 <- read.csv(file.path(rds_dir, "iteration_map_1000.csv"))

stopifnot(nrow(iteration_map_1000) == 1000)
stopifnot(n_distinct(iteration_map_1000$om) == n_om)
stopifnot(all(table(iteration_map_1000$om) == 10))

# ============================================================
# 3. IDENTIFY THE 100 UNIQUE HISTORICAL OMs
# ============================================================
# The 1000-iteration fleet contains ten replicated trajectories
# for each of the 100 historical OMs.
#
# Historical diagnostics use one representative iteration per OM.
# The FLFleetExt object is kept intact. Selection is performed
# after extracting FLQuant quantities from the fleet.

representative_iter <- iteration_map_1000 |>
  group_by(om) |>
  slice_min(iter1000, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(om)

source_iter_100 <- representative_iter$iter1000

stopifnot(nrow(representative_iter) == n_om)
stopifnot(length(source_iter_100) == n_om)
stopifnot(identical(representative_iter$om, seq_len(n_om)))

fleets <- fleets.1000
fleets.ctrl <- fleets.ctrl.1000

# ============================================================
# 4. EXTRACT HISTORICAL SEASONAL CATCH PROPORTIONS
# ============================================================
# Historical catch is first extracted from the complete
# 1000-iteration fleet object.
#
# One representative iteration is then selected for each of
# the 100 historical OMs using iteration_map_1000.

catch_hist_1000 <- catchWStock(fleets, stock = "ANE")
catch_season_1000 <- quantSums(catch_hist_1000)
catch_season_df <- as.data.frame(catch_season_1000) |>
  as_tibble() |>
  rename(catch = data) |>
  mutate(
    year = as.integer(as.character(year)),
    season = as.character(season),
    iter1000 = as.integer(as.character(iter)),
    catch = as.numeric(catch)) |>
  filter(year %in% hist_years) |>
  filter(iter1000 %in% source_iter_100)

# ============================================================
# MAP REPRESENTATIVE ITERATIONS TO HISTORICAL OMs
# ============================================================

iter_to_om <- representative_iter |> dplyr::select(iter1000, om)
catch_season_df <- catch_season_df |> left_join(iter_to_om, by = "iter1000")

stopifnot(n_distinct(catch_season_df$iter1000) == n_om)
stopifnot(n_distinct(catch_season_df$om) == n_om)
stopifnot(!any(is.na(catch_season_df$om)))

# ============================================================
# CALCULATE ANNUAL SEASONAL CATCH PROPORTIONS
# ============================================================

seasonal_prop_df <- catch_season_df |>
  group_by(om, year) |>
  mutate(
    annual_catch = sum(catch, na.rm = TRUE),
    prop = catch / annual_catch) |>
  ungroup() |>
  mutate(prop = ifelse(is.finite(prop), prop, NA_real_))

seasonal_prop_check <- seasonal_prop_df |> group_by(om, year) |>
  summarise(sum_prop = sum(prop, na.rm = TRUE),.groups = "drop")

range(seasonal_prop_check$sum_prop)
n_distinct(seasonal_prop_df$om)
n_distinct(seasonal_prop_df$year)
range(seasonal_prop_df$year)

# ============================================================
# 5. SUMMARISE HISTORICAL SEASONAL ALLOCATION
# ============================================================

seasonal_summary <- seasonal_prop_df |>
  group_by(season) |>
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    cv_prop = sd_prop / mean_prop,
    min_prop = min(prop, na.rm = TRUE),
    p05_prop = quantile(prop, 0.05, na.rm = TRUE),
    p50_prop = quantile(prop, 0.50, na.rm = TRUE),
    p95_prop = quantile(prop, 0.95, na.rm = TRUE),
    max_prop = max(prop, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(seasonal_summary,file.path(ref_tab_dir, "seasonal_catch_allocation_summary.csv"))

# ============================================================
# 5.1 FIGURE — HISTORICAL SEASONAL ALLOCATION
# ============================================================

p_ts <- ggplot(seasonal_prop_df,aes(x = year, y = prop, colour = season)) +
  geom_line() +
  geom_point(size = 1) +
  labs(x = "Year",y = "Proportion of annual catch",
    colour = "Quarter",title = "Historical seasonal catch allocation") +
  theme_bw()

ggsave(file.path(ref_fig_dir, "seasonal_catch_allocation_timeseries.png"),p_ts,width = 9,height = 5,dpi = 300)

# ============================================================
# 6. CHECK REFERENCE PROJECTION SEASONAL SHARES
# ============================================================
# Script 10 assigns the mean seasonal catch shares from 2015-2024 to
# every projection year. Recalculate those means from the historical
# fleet and compare them with fleets.ctrl.


recent10_years <- (ass.yr - 9):ass.yr

recent10_by_om <- seasonal_prop_df |>
  filter(year %in% recent10_years) |>
  group_by(om, season) |>
  summarise(
    historical_recent10_mean = mean(prop, na.rm = TRUE),
    .groups = "drop"
  )

seasonal_share_projection <- as.data.frame(fleets.ctrl$seasonal.share[[1]]) |>
  as_tibble() |>
  rename(reference_share = data) |>
  mutate(
    year = as.integer(as.character(year)),
    season = as.character(season),
    iter1000 = as.integer(as.character(iter))
  ) |>
  filter(year == min(proj.yrs)) |>
  filter(iter1000 %in% source_iter_100) |>
  left_join(iter_to_om, by = "iter1000") |>
  dplyr::select(om, season, reference_share)

stopifnot(n_distinct(seasonal_share_projection$om) == n_om)
stopifnot(!any(is.na(seasonal_share_projection$om)))

seasonal_share_check <- recent10_by_om |>
  left_join(
    seasonal_share_projection,
    by = c("om", "season")) |>
  mutate(abs_diff = abs(historical_recent10_mean - reference_share))

seasonal_share_max_diff <- max(seasonal_share_check$abs_diff, na.rm = TRUE)

stopifnot(seasonal_share_max_diff < 1e-10)

write_csv(seasonal_share_check,file.path(ref_tab_dir, "reference_seasonal_share_check_100OMs.csv"))

# ============================================================
# 7. YEAR-LEVEL SERIES FOR TEMPORAL DIAGNOSTICS

# ============================================================
# Trend, breakpoint, change-point, PCA and clustering analyses classify
# historical years. Therefore, a single seasonal allocation vector is
# required for each year. It is calculated here as the mean seasonal
# proportion across the 100 conditioned OMs.

seasonal_year_df <- seasonal_prop_df |>
  group_by(year, season) |>
  summarise(
    prop = mean(prop, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
# 8. DIAGNOSE TEMPORAL TRENDS BY SEASON
# ============================================================

trend_results <- seasonal_year_df |>
  group_by(season) |>
  group_modify(~ {
    
    fit <- lm(prop ~ year, data = .x)
    sfit <- summary(fit)
    
    tibble(
      intercept = coef(fit)[1],
      slope = coef(fit)[2],
      p_value_slope = coef(sfit)[2, 4],
      r_squared = sfit$r.squared)
  }) |>
  ungroup()

write_csv(trend_results,file.path(sens_tab_dir, "seasonal_catch_allocation_trend_results.csv"))

p_trend <- ggplot(seasonal_year_df,aes(x = year, y = prop)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ season, scales = "free_y") +
  labs(x = "Year",y = "Proportion of annual catch",
    title = "Linear trend in seasonal catch allocation") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_trends.png"),
  p_trend,width = 9,height = 6,dpi = 300)

# ============================================================
# 9. DIAGNOSE STRUCTURAL BREAKS BY SEASON
# ============================================================

break_results <- list()
break_fitted  <- list()

for (s in sort(unique(seasonal_year_df$season))) {
  
  dat_s <- seasonal_year_df |>
    filter(season == s) |>
    arrange(year) |>
    filter(!is.na(prop))
  
  bp <- breakpoints(prop ~ 1, data = dat_s)
  bp_opt <- breakpoints(prop ~ 1, data = dat_s, breaks = which.min(BIC(bp)))
  
  bp_years <- dat_s$year[bp_opt$breakpoints]
  bp_years <- bp_years[!is.na(bp_years)]
  
  break_results[[s]] <- tibble(
    season = s,
    n_breaks = length(bp_years),
    break_years = paste(bp_years, collapse = ", "),
    bic_min = min(BIC(bp), na.rm = TRUE))
  
  dat_s$regime <- breakfactor(bp_opt)
  dat_s$break_years <- paste(bp_years, collapse = ", ")
  
  break_fitted[[s]] <- dat_s
}

break_results_df <- bind_rows(break_results)
break_fitted_df  <- bind_rows(break_fitted)

write_csv(break_results_df,
  file.path(sens_tab_dir, "seasonal_catch_allocation_breakpoints.csv"))

write_csv(break_fitted_df,file.path(sens_tab_dir, "seasonal_catch_allocation_breakpoint_regimes.csv"))

p_breaks <- ggplot(break_fitted_df,aes(x = year, y = prop, colour = regime)) +
  geom_point(size = 1.8) +
  geom_line() +
  facet_wrap(~ season, scales = "free_y") +
  labs(x = "Year",y = "Proportion of annual catch",
    colour = "Regime",title = "Structural break analysis of seasonal catch allocation") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_breakpoints.png"),
  p_breaks,width = 9,height = 6,dpi = 300)

# ============================================================
# 10. DIAGNOSE CHANGE POINTS IN MEAN AND VARIANCE
# ============================================================

cpt_results <- list()

for (s in sort(unique(seasonal_year_df$season))) {
  
  dat_s <- seasonal_year_df |>
    filter(season == s) |>
    arrange(year) |>
    filter(!is.na(prop))
  
  x <- dat_s$prop
  
  cpt_mv <- cpt.meanvar(
    x,
    method = "PELT",
    penalty = "MBIC",
    class = TRUE)
  
  cpt_pos <- cpts(cpt_mv)
  cpt_years <- dat_s$year[cpt_pos]
  
  cpt_results[[s]] <- tibble(
    season = s,
    n_changepoints = length(cpt_years),
    changepoint_years = paste(cpt_years, collapse = ", "))
}

cpt_results_df <- bind_rows(cpt_results)

write_csv(cpt_results_df,file.path(sens_tab_dir, "seasonal_catch_allocation_changepoints_meanvar.csv"))

# ============================================================
# 11. BUILD YEAR × SEASON MATRIX
# ============================================================

seasonal_wide <- seasonal_year_df |>
  dplyr::select(year, season, prop) |>
  pivot_wider(
    names_from = season,
    values_from = prop,
    names_prefix = "Q") |>
  arrange(year)

write_csv(seasonal_wide,file.path(sens_tab_dir, "seasonal_catch_allocation_wide.csv"))

X <- seasonal_wide |>
  dplyr::select(starts_with("Q")) |>
  as.data.frame()

rownames(X) <- seasonal_wide$year

# ============================================================
# 12. PCA OF HISTORICAL SEASONAL ALLOCATION
# ============================================================

pca_fit <- prcomp(X, center = TRUE, scale. = TRUE)

pca_df <- as_tibble(pca_fit$x[, 1:2], rownames = "year") |>
  mutate(year = as.integer(year))

pca_var <- tibble(PC = paste0("PC", seq_along(pca_fit$sdev)),
                  variance_explained = pca_fit$sdev^2 / sum(pca_fit$sdev^2))

write_csv(pca_df,
  file.path(sens_tab_dir, "seasonal_catch_allocation_pca_scores.csv"))

write_csv(pca_var,
  file.path(sens_tab_dir, "seasonal_catch_allocation_pca_variance.csv"))

p_pca <- ggplot(pca_df,aes(x = PC1, y = PC2, label = year)) +
  geom_point(size = 2) +
  geom_text(vjust = -0.6, size = 3) +
  labs(x = "PC1",y = "PC2",title = "PCA of seasonal catch allocation") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_pca.png"),
  p_pca,width = 8,height = 6,dpi = 300)

# ============================================================
# 13. HIERARCHICAL CLUSTERING OF YEARS
# ============================================================

dist_mat <- dist(scale(X))
hc <- hclust(dist_mat, method = "ward.D2")

# Choose 3 clusters as a diagnostic default.
# This can be changed after inspecting the dendrogram.
k <- 3

cluster_df <- tibble(year = seasonal_wide$year,cluster = factor(cutree(hc, k = k)))

seasonal_cluster_df <- seasonal_year_df |>
  left_join(cluster_df, by = "year")

cluster_summary <- seasonal_cluster_df |>
  group_by(cluster, season) |>
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    n_years = n_distinct(year),
    years = paste(sort(unique(year)), collapse = ", "),
    .groups = "drop")

write_csv(cluster_df,
  file.path(sens_tab_dir, "seasonal_catch_allocation_year_clusters.csv"))

write_csv(cluster_summary,
  file.path(sens_tab_dir, "seasonal_catch_allocation_cluster_summary.csv"))

png(filename = file.path(sens_fig_dir, "seasonal_catch_allocation_dendrogram.png"),
  width = 1800,height = 1200,res = 200)

plot(hc,main = "Hierarchical clustering of seasonal catch allocation",xlab = "Year",sub = "")
rect.hclust(hc, k = k, border = 2:4)
dev.off()

p_cluster_ts <- ggplot(seasonal_cluster_df,aes(x = year,y = prop,colour = cluster)) +
  geom_point(size = 2) +
  facet_wrap(~ season,scales = "free_y") +
  labs(x = "Year",y = "Proportion of annual catch",
       colour = "Cluster",title = "Seasonal catch allocation clusters") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_clusters_timeseries.png"),
  p_cluster_ts,width = 9,height = 6,dpi = 300)

p_cluster_bar <- ggplot(cluster_summary,aes(x = season, y = mean_prop, fill = cluster)) +
  geom_col(position = "dodge") +
  labs(x = "Quarter",y = "Mean proportion of annual catch",
    fill = "Cluster",title = "Mean seasonal allocation by cluster") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_cluster_means.png"),p_cluster_bar,
width = 8,height = 5,dpi = 300)

#===============================================================================
# PCA coloured by hierarchical cluster
#===============================================================================

pca_plot_df <- pca_df |>
  left_join(cluster_df, by = "year")

p_pca_cluster <- ggplot(pca_plot_df,aes(x = PC1,y = PC2,colour = cluster,label = year)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6, size = 3) +
  labs(x = "PC1",y = "PC2",colour = "Cluster",title = "PCA of seasonal catch allocation by cluster") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_pca_by_cluster.png"),
  p_pca_cluster,width = 8,height = 6,dpi = 300)

p_pca_cluster_ellipse <- ggplot(pca_plot_df,aes(x = PC1, y = PC2, colour = cluster)) +
  stat_ellipse(aes(group = cluster), linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(aes(label = year), vjust = -0.6, size = 3) +
  labs(x = "PC1",y = "PC2",colour = "Cluster",title = "PCA of seasonal catch allocation by cluster") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_pca_by_cluster_ellipse.png"),
  p_pca_cluster_ellipse,width = 8,height = 6,dpi = 300)

pca_var |> mutate(variance_explained = round(variance_explained * 100, 1))


# ============================================================
# 14. COMPARE HISTORICAL, RECENT-10 AND RECENT-3 PERIODS
# ============================================================

period_df <- seasonal_year_df |>
  mutate(
    period = case_when(
      year %in% (ass.yr - 2):ass.yr ~ "recent3",
      year %in% (ass.yr - 9):ass.yr ~ "recent10",
      TRUE ~ "historical_previous"))

period_summary <- period_df |>
  group_by(period, season) |>
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    n_years = n_distinct(year),
    .groups = "drop")

write_csv(period_summary,
  file.path(sens_tab_dir, "seasonal_catch_allocation_period_summary.csv"))

p_period <- ggplot(
  period_summary,
  aes(x = season, y = mean_prop, fill = period)) +
  geom_col(position = "dodge") +
  labs(
    x = "Quarter",
    y = "Mean proportion of annual catch",
    fill = "Period",
    title = "Seasonal catch allocation by period") +
  theme_bw()

ggsave(file.path(sens_fig_dir, "seasonal_catch_allocation_period_comparison.png"),
  p_period,width = 8,height = 5,dpi = 300)

# ============================================================
# 15. BUILD CANDIDATE SEASONAL-ALLOCATION SCENARIOS
# ============================================================

scenario_historical <- seasonal_year_df |>
  group_by(season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = "historical")

scenario_recent10 <- seasonal_year_df |>
  filter(year %in% (ass.yr - 9):ass.yr) |>
  group_by(season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = "recent10")

scenario_recent3 <- seasonal_year_df |>
  filter(year %in% (ass.yr - 2):ass.yr) |>
  group_by(season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = "recent3")

scenario_clusters <- seasonal_cluster_df |>
  group_by(cluster, season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = paste0("cluster_", cluster)) |>
  dplyr::select(scenario, season, prop)

candidate_scenarios <- bind_rows(scenario_historical,
                                 scenario_recent10,
                                 scenario_recent3,
                                 scenario_clusters) |>
                      group_by(scenario) |>
                      mutate(prop = prop / sum(prop, na.rm = TRUE)) |>
                      ungroup()

write_csv(candidate_scenarios,file.path(sens_tab_dir, "candidate_seasonal_catch_allocation_scenarios.csv"))

p_candidates <- ggplot(candidate_scenarios,aes(x = season, y = prop, fill = scenario)) +
                geom_col(position = "dodge") +
                labs(x = "Quarter",y = "Catch allocation proportion",
                fill = "Scenario",title = "Candidate seasonal catch allocation scenarios") +
                theme_bw()

ggsave(file.path(sens_fig_dir, "candidate_seasonal_catch_allocation_scenarios.png"),
  p_candidates,width = 10,height = 6,dpi = 300)

# ============================================================
# 16. SUMMARISE SENSITIVITY DIAGNOSTICS
# ============================================================

diagnostic_summary <- seasonal_summary |>
  dplyr::select(season,mean_prop,cv_prop,p05_prop,p95_prop) |>
  left_join(trend_results,by = "season") |>
  left_join(break_results_df,by = "season") |>
  left_join(cpt_results_df,by = "season")

write_csv(diagnostic_summary,file.path(sens_tab_dir,"seasonal_catch_allocation_diagnostic_summary.csv"))

diagnostic_summary
candidate_scenarios
# ============================================================
# 17. HISTORICAL CATCH-AT-AGE COMPOSITION BY QUARTER
# ============================================================
# Catch-at-age is extracted from the complete 1000-iteration
# fleet object. One representative iteration per historical OM
# is then selected using the same mapping applied above.

catch_n_hist_1000 <- landings.n(fleets$SEINE@metiers$ALL@catches$ANE) +
                      discards.n(fleets$SEINE@metiers$ALL@catches$ANE)

catch_age_season_hist <- as.data.frame(catch_n_hist_1000) |>
                          as_tibble() |>
                          rename(catch_n = data) |>
                          mutate(
                          age = as.character(age),
                          year = as.integer(as.character(year)),
                          season = as.character(season),
                          iter1000 = as.integer(as.character(iter)),
                          catch_n = as.numeric(catch_n)) |>
                          filter(year %in% hist_years) |>
                          filter(iter1000 %in% source_iter_100) |>
                          left_join(iter_to_om, by = "iter1000")

stopifnot(n_distinct(catch_age_season_hist$iter1000) == n_om)
stopifnot(n_distinct(catch_age_season_hist$om) == n_om)
stopifnot(!any(is.na(catch_age_season_hist$om)))


catch_age_season_prop <- catch_age_season_hist |>
                          group_by(om, year, season) |>
                          mutate(
                          total_catch_n = sum(catch_n, na.rm = TRUE),
                          prop_age = catch_n / total_catch_n) |>
                          ungroup() |>
                          mutate(prop_age = ifelse(is.finite(prop_age), prop_age, NA_real_))


catch_age_prop_check <- catch_age_season_prop |>
                        group_by(om, year, season) |>
                        summarise(
                        sum_prop_age = sum(prop_age, na.rm = TRUE),
                        .groups = "drop")

range(catch_age_prop_check$sum_prop_age)
n_distinct(catch_age_season_prop$om)
range(catch_age_season_prop$year)

catch_age_season_summary <- catch_age_season_prop |>
                            group_by(season, age) |>
                            summarise(
                            mean_prop_age = mean(prop_age, na.rm = TRUE),
                            sd_prop_age = sd(prop_age, na.rm = TRUE),
                            .groups = "drop")

write_csv(catch_age_season_summary, file.path(ref_tab_dir, "historical_catch_at_age_composition_by_quarter.csv"))

p_catch_age_season <- ggplot(catch_age_season_summary,aes(x = season, y = mean_prop_age, fill = age)) +
                      geom_col(position = "stack") +
                      labs( x = "Quarter",y = "Mean catch-at-age proportion",
                      fill = "Age",title = "Historical catch-at-age composition by quarter") +
                      theme_bw()

ggsave(file.path(ref_fig_dir, "historical_catch_at_age_composition_by_quarter.png"),p_catch_age_season,width = 8,height = 5,dpi = 300)

catch_age_season_table <- catch_age_season_summary |>
                          dplyr::select(season, age, mean_prop_age) |>
                          pivot_wider(names_from = age,values_from = mean_prop_age,names_prefix = "Age_")

write_csv(catch_age_season_table,file.path(ref_tab_dir, "historical_catch_at_age_composition_by_quarter_table.csv"))

catch_age_season_table

# ============================================================
# 18. COMPLETION MESSAGE
# ============================================================

message("Reference fishery-conditioning outputs saved in: ", dirname(ref_fig_dir))
message("Sensitivity diagnostics saved in: ", dirname(sens_fig_dir))

