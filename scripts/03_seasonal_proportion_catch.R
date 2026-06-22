################################################################################
# Seasonal catch allocation diagnostics
# Anchovy 9aS OM
################################################################################

rm(list = ls())

library(FLCore)
library(FLBEIA)
library(FLFishery)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(here)

# Optional packages
# install.packages(c("strucchange", "changepoint", "factoextra"))
library(strucchange)
library(changepoint)

wd <- here()
setwd(wd)

out_dir <- "outputs/mse/diagnostics/seasonal_catch_allocation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#===============================================================================
# 1. Load MSE input objects
#===============================================================================

load("data/mse/test_100iter_10year.RData")
# This file should contain: fleets, fleets.ctrl, first.yr, ass.yr, proj.yrs, etc.

hist_years <- first.yr:ass.yr

#===============================================================================
# 2. Extract historical seasonal catch proportions
#===============================================================================

catch_hist <- catchWStock(fleets, stock = "ANE")

catch_season <- quantSums(catch_hist)

catch_season_df <- as.data.frame(iterMeans(catch_season)) |>
  as_tibble() |>
  rename(catch = data) |>
  mutate(
    year = as.integer(as.character(year)),
    season = as.character(season),
    catch = as.numeric(catch)
  ) |>
  filter(year %in% hist_years)

seasonal_prop_df <- catch_season_df |>
  group_by(year) |>
  mutate(
    annual_catch = sum(catch, na.rm = TRUE),
    prop = catch / annual_catch
  ) |>
  ungroup() |>
  mutate(
    prop = ifelse(is.finite(prop), prop, NA_real_)
  )

write_csv(
  seasonal_prop_df,
  file.path(out_dir, "historical_seasonal_catch_proportions.csv")
)

#===============================================================================
# 3. Summary statistics
#===============================================================================

seasonal_summary <- seasonal_prop_df |>
  group_by(season) |>
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop   = sd(prop, na.rm = TRUE),
    cv_prop   = sd_prop / mean_prop,
    min_prop  = min(prop, na.rm = TRUE),
    p05_prop  = quantile(prop, 0.05, na.rm = TRUE),
    p50_prop  = quantile(prop, 0.50, na.rm = TRUE),
    p95_prop  = quantile(prop, 0.95, na.rm = TRUE),
    max_prop  = max(prop, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  seasonal_summary,
  file.path(out_dir, "seasonal_catch_allocation_summary.csv")
)

#===============================================================================
# 4. Time series plot
#===============================================================================

p_ts <- ggplot(
  seasonal_prop_df,
  aes(x = year, y = prop, colour = season)
) +
  geom_line() +
  geom_point(size = 1.8) +
  labs(
    x = "Year",
    y = "Proportion of annual catch",
    colour = "Quarter",
    title = "Historical seasonal catch allocation"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_timeseries.png"),
  p_ts,
  width = 9,
  height = 5,
  dpi = 300
)

#===============================================================================
# 5. Trend analysis by season
#===============================================================================

trend_results <- seasonal_prop_df |>
  group_by(season) |>
  group_modify(~ {
    
    fit <- lm(prop ~ year, data = .x)
    sfit <- summary(fit)
    
    tibble(
      intercept = coef(fit)[1],
      slope = coef(fit)[2],
      p_value_slope = coef(sfit)[2, 4],
      r_squared = sfit$r.squared
    )
  }) |>
  ungroup()

write_csv(
  trend_results,
  file.path(out_dir, "seasonal_catch_allocation_trend_results.csv")
)

p_trend <- ggplot(
  seasonal_prop_df,
  aes(x = year, y = prop)
) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ season, scales = "free_y") +
  labs(
    x = "Year",
    y = "Proportion of annual catch",
    title = "Linear trend in seasonal catch allocation"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_trends.png"),
  p_trend,
  width = 9,
  height = 6,
  dpi = 300
)

#===============================================================================
# 6. Structural break analysis by season
#===============================================================================

break_results <- list()
break_fitted  <- list()

for (s in sort(unique(seasonal_prop_df$season))) {
  
  dat_s <- seasonal_prop_df |>
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
    bic_min = min(BIC(bp), na.rm = TRUE)
  )
  
  dat_s$regime <- breakfactor(bp_opt)
  dat_s$break_years <- paste(bp_years, collapse = ", ")
  
  break_fitted[[s]] <- dat_s
}

break_results_df <- bind_rows(break_results)
break_fitted_df  <- bind_rows(break_fitted)

write_csv(
  break_results_df,
  file.path(out_dir, "seasonal_catch_allocation_breakpoints.csv")
)

write_csv(
  break_fitted_df,
  file.path(out_dir, "seasonal_catch_allocation_breakpoint_regimes.csv")
)

p_breaks <- ggplot(
  break_fitted_df,
  aes(x = year, y = prop, colour = regime)
) +
  geom_point(size = 1.8) +
  geom_line() +
  facet_wrap(~ season, scales = "free_y") +
  labs(
    x = "Year",
    y = "Proportion of annual catch",
    colour = "Regime",
    title = "Structural break analysis of seasonal catch allocation"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_breakpoints.png"),
  p_breaks,
  width = 9,
  height = 6,
  dpi = 300
)

#===============================================================================
# 7. Change-point analysis in mean and variance
#===============================================================================

cpt_results <- list()

for (s in sort(unique(seasonal_prop_df$season))) {
  
  dat_s <- seasonal_prop_df |>
    filter(season == s) |>
    arrange(year) |>
    filter(!is.na(prop))
  
  x <- dat_s$prop
  
  cpt_mv <- cpt.meanvar(
    x,
    method = "PELT",
    penalty = "MBIC",
    class = TRUE
  )
  
  cpt_pos <- cpts(cpt_mv)
  cpt_years <- dat_s$year[cpt_pos]
  
  cpt_results[[s]] <- tibble(
    season = s,
    n_changepoints = length(cpt_years),
    changepoint_years = paste(cpt_years, collapse = ", ")
  )
}

cpt_results_df <- bind_rows(cpt_results)

write_csv(
  cpt_results_df,
  file.path(out_dir, "seasonal_catch_allocation_changepoints_meanvar.csv")
)

#===============================================================================
# 8. Wide matrix of seasonal proportions
#===============================================================================

seasonal_wide <- seasonal_prop_df |>
  select(year, season, prop) |>
  pivot_wider(
    names_from = season,
    values_from = prop,
    names_prefix = "Q"
  ) |>
  arrange(year)

write_csv(
  seasonal_wide,
  file.path(out_dir, "seasonal_catch_allocation_wide.csv")
)

X <- seasonal_wide |>
  select(starts_with("Q")) |>
  as.data.frame()

rownames(X) <- seasonal_wide$year

#===============================================================================
# 9. PCA
#===============================================================================

pca_fit <- prcomp(X, center = TRUE, scale. = TRUE)

pca_df <- as_tibble(pca_fit$x[, 1:2], rownames = "year") |>
  mutate(year = as.integer(year))

pca_var <- tibble(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  variance_explained = pca_fit$sdev^2 / sum(pca_fit$sdev^2)
)

write_csv(
  pca_df,
  file.path(out_dir, "seasonal_catch_allocation_pca_scores.csv")
)

write_csv(
  pca_var,
  file.path(out_dir, "seasonal_catch_allocation_pca_variance.csv")
)

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, label = year)
) +
  geom_point(size = 2) +
  geom_text(vjust = -0.6, size = 3) +
  labs(
    x = "PC1",
    y = "PC2",
    title = "PCA of seasonal catch allocation"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_pca.png"),
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

#===============================================================================
# 10. Hierarchical clustering of years
#===============================================================================

dist_mat <- dist(scale(X))
hc <- hclust(dist_mat, method = "ward.D2")

# Choose 3 clusters as a diagnostic default.
# This can be changed after inspecting the dendrogram.
k <- 3

cluster_df <- tibble(
  year = seasonal_wide$year,
  cluster = factor(cutree(hc, k = k))
)

seasonal_cluster_df <- seasonal_prop_df |>
  left_join(cluster_df, by = "year")

cluster_summary <- seasonal_cluster_df |>
  group_by(cluster, season) |>
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    n_years = n_distinct(year),
    years = paste(sort(unique(year)), collapse = ", "),
    .groups = "drop"
  )

write_csv(
  cluster_df,
  file.path(out_dir, "seasonal_catch_allocation_year_clusters.csv")
)

write_csv(
  cluster_summary,
  file.path(out_dir, "seasonal_catch_allocation_cluster_summary.csv")
)

png(
  filename = file.path(out_dir, "seasonal_catch_allocation_dendrogram.png"),
  width = 1800,
  height = 1200,
  res = 200
)
plot(
  hc,
  main = "Hierarchical clustering of seasonal catch allocation",
  xlab = "Year",
  sub = ""
)
rect.hclust(hc, k = k, border = 2:4)
dev.off()

p_cluster_ts <- ggplot(
  seasonal_cluster_df,
  aes(x = year, y = prop, colour = cluster)
) +
  geom_point(size = 2) +
  geom_line(aes(group = season)) +
  facet_wrap(~ season, scales = "free_y") +
  labs(
    x = "Year",
    y = "Proportion of annual catch",
    colour = "Cluster",
    title = "Seasonal catch allocation clusters"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_clusters_timeseries.png"),
  p_cluster_ts,
  width = 9,
  height = 6,
  dpi = 300
)

p_cluster_bar <- ggplot(
  cluster_summary,
  aes(x = season, y = mean_prop, fill = cluster)
) +
  geom_col(position = "dodge") +
  labs(
    x = "Quarter",
    y = "Mean proportion of annual catch",
    fill = "Cluster",
    title = "Mean seasonal allocation by cluster"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_cluster_means.png"),
  p_cluster_bar,
  width = 8,
  height = 5,
  dpi = 300
)

#===============================================================================
# PCA coloured by hierarchical cluster
#===============================================================================

pca_plot_df <- pca_df |>
  left_join(cluster_df, by = "year")

p_pca_cluster <- ggplot(
  pca_plot_df,
  aes(
    x = PC1,
    y = PC2,
    colour = cluster,
    label = year
  )
) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6, size = 3) +
  labs(
    x = "PC1",
    y = "PC2",
    colour = "Cluster",
    title = "PCA of seasonal catch allocation by cluster"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_pca_by_cluster.png"),
  p_pca_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

p_pca_cluster_ellipse <- ggplot(
  pca_plot_df,
  aes(
    x = PC1,
    y = PC2,
    colour = cluster,
    label = year
  )
) +
  stat_ellipse(aes(group = cluster), linewidth = 0.7) +
  geom_point(size = 3) +
  geom_text(vjust = -0.6, size = 3) +
  labs(
    x = "PC1",
    y = "PC2",
    colour = "Cluster",
    title = "PCA of seasonal catch allocation by cluster"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_pca_by_cluster_ellipse.png"),
  p_pca_cluster_ellipse,
  width = 8,
  height = 6,
  dpi = 300
)

pca_var |> 
  mutate(variance_explained = round(variance_explained * 100, 1))

ggplot(
  seasonal_cluster_df ,
  aes(x = year,y = prop,fill = cluster)) +
  geom_col() +
  facet_wrap(~ season, scales = "free_x") +
  labs(
    x = "Year",
    y = "Proportion of annual catch",
    fill = "Cluster",
    title = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
#===============================================================================
# 11. Period comparison: historical, recent10, recent3
#===============================================================================

period_df <- seasonal_prop_df |>
  mutate(
    period = case_when(
      year %in% (ass.yr - 2):ass.yr ~ "recent3",
      year %in% (ass.yr - 9):ass.yr ~ "recent10",
      TRUE ~ "historical_previous"
    )
  )

period_summary <- period_df |>
  group_by(period, season) |>
  summarise(
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop = sd(prop, na.rm = TRUE),
    n_years = n_distinct(year),
    .groups = "drop"
  )

write_csv(
  period_summary,
  file.path(out_dir, "seasonal_catch_allocation_period_summary.csv")
)

p_period <- ggplot(
  period_summary,
  aes(x = season, y = mean_prop, fill = period)
) +
  geom_col(position = "dodge") +
  labs(
    x = "Quarter",
    y = "Mean proportion of annual catch",
    fill = "Period",
    title = "Seasonal catch allocation by period"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_allocation_period_comparison.png"),
  p_period,
  width = 8,
  height = 5,
  dpi = 300
)

#===============================================================================
# 12. Candidate scenario values
#===============================================================================

scenario_historical <- seasonal_prop_df |>
  group_by(season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = "historical")

scenario_recent10 <- seasonal_prop_df |>
  filter(year %in% (ass.yr - 9):ass.yr) |>
  group_by(season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = "recent10")

scenario_recent3 <- seasonal_prop_df |>
  filter(year %in% (ass.yr - 2):ass.yr) |>
  group_by(season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = "recent3")

scenario_clusters <- seasonal_cluster_df |>
  group_by(cluster, season) |>
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  mutate(scenario = paste0("cluster_", cluster)) |>
  select(scenario, season, prop)

candidate_scenarios <- bind_rows(
  scenario_historical,
  scenario_recent10,
  scenario_recent3,
  scenario_clusters
) |>
  group_by(scenario) |>
  mutate(prop = prop / sum(prop, na.rm = TRUE)) |>
  ungroup()

write_csv(
  candidate_scenarios,
  file.path(out_dir, "candidate_seasonal_catch_allocation_scenarios.csv")
)

p_candidates <- ggplot(
  candidate_scenarios,
  aes(x = season, y = prop, fill = scenario)
) +
  geom_col(position = "dodge") +
  labs(
    x = "Quarter",
    y = "Catch allocation proportion",
    fill = "Scenario",
    title = "Candidate seasonal catch allocation scenarios"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "candidate_seasonal_catch_allocation_scenarios.png"),
  p_candidates,
  width = 10,
  height = 6,
  dpi = 300
)

#===============================================================================
# 13. Decision table
#===============================================================================

decision_table <- seasonal_summary |>
  select(season, mean_prop, cv_prop, p05_prop, p95_prop) |>
  left_join(trend_results, by = "season") |>
  left_join(break_results_df, by = "season") |>
  left_join(cpt_results_df, by = "season") |>
  mutate(
    high_variability = cv_prop > 0.30,
    trend_detected = p_value_slope < 0.05,
    breakpoint_detected = n_breaks > 0,
    changepoint_detected = n_changepoints > 0,
    include_sensitivity = high_variability |
      trend_detected |
      breakpoint_detected |
      changepoint_detected
  )

write_csv(
  decision_table,
  file.path(out_dir, "seasonal_catch_allocation_decision_table.csv")
)

decision_table
candidate_scenarios

message("Seasonal catch allocation diagnostics saved in: ", out_dir)

#==============================================================================
# Historical catch-at-age composition by quarter
#==============================================================================

catch_hist <- catchWStock(fleets, stock = "ANE")

catch_n_hist <- landings.n(fleets$SEINE@metiers$ALL@catches$ANE) +
  discards.n(fleets$SEINE@metiers$ALL@catches$ANE)

catch_age_season_hist <- as.data.frame(iterMeans(catch_n_hist)) |>
  as_tibble() |>
  rename(catch_n = data) |>
  mutate(
    age = as.character(age),
    year = as.integer(as.character(year)),
    season = as.character(season),
    catch_n = as.numeric(catch_n)
  ) |>
  filter(year %in% hist_years)

catch_age_season_prop <- catch_age_season_hist |>
  group_by(year, season) |>
  mutate(
    total_catch_n = sum(catch_n, na.rm = TRUE),
    prop_age = catch_n / total_catch_n
  ) |>
  ungroup() |>
  mutate(
    prop_age = ifelse(is.finite(prop_age), prop_age, NA_real_)
  )

catch_age_season_summary <- catch_age_season_prop |>
  group_by(season, age) |>
  summarise(
    mean_prop_age = mean(prop_age, na.rm = TRUE),
    sd_prop_age = sd(prop_age, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  catch_age_season_summary,
  file.path(out_dir, "historical_catch_at_age_composition_by_quarter.csv")
)

p_catch_age_season <- ggplot(
  catch_age_season_summary,
  aes(x = season, y = mean_prop_age, fill = age)
) +
  geom_col(position = "stack") +
  labs(
    x = "Quarter",
    y = "Mean catch-at-age proportion",
    fill = "Age",
    title = "Historical catch-at-age composition by quarter"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "historical_catch_at_age_composition_by_quarter.png"),
  p_catch_age_season,
  width = 8,
  height = 5,
  dpi = 300
)

catch_age_season_table <- catch_age_season_summary |>
  select(season, age, mean_prop_age) |>
  pivot_wider(
    names_from = age,
    values_from = mean_prop_age,
    names_prefix = "Age_"
  )

write_csv(
  catch_age_season_table,
  file.path(out_dir, "historical_catch_at_age_composition_by_quarter_table.csv")
)

catch_age_season_table



