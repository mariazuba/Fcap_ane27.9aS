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

wd <- here()
setwd(wd)

res_dir <- "outputs/mse/res"
out_dir <- "outputs/mse/diagnostics/fcap_biological_range"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hist_years <- 1989:2024
proj_years <- 2025:2044

flq_to_df <- function(x, value_name = "value") {
  as.data.frame(x) %>%
    as_tibble() %>%
    rename(!!value_name := data) %>%
    mutate(
      age = as.character(age),
      year = as.integer(as.character(year)),
      season = as.character(season),
      iter = as.integer(as.character(iter))
    )
}

extract_fcap <- function(file) {
  f <- basename(file)
  x <- str_extract(f, "Fcap_[0-9]+p[0-9]+")
  x <- str_remove(x, "Fcap_")
  as.numeric(str_replace(x, "p", "."))
}

extract_scenario <- function(file) {
  f <- basename(file)
  str_extract(f, "sc[0-9]+")
}

extract_mse_quantities <- function(file) {
  
  message("Reading: ", file)
  
  mse_res <- readRDS(file)
  
  biol  <- mse_res$biols$ANE
  fleet <- mse_res$fleets$SEINE
  metier <- fleet@metiers$ALL
  catch_obj <- metier@catches$ANE
  
  catch_n <- catch_obj@landings.n + catch_obj@discards.n
  stock_n <- biol@n
  
  exploitation_n <- catch_n / stock_n
  exploitation_n[!is.finite(exploitation_n)] <- NA
  
  catch_df <- flq_to_df(catch_n, "catch_n")
  expl_df  <- flq_to_df(exploitation_n, "u_age")
  stock_df <- flq_to_df(stock_n, "stock_n")
  
  stock_df <- stock_df %>%
    mutate(
      scenario = extract_scenario(file),
      Fcap = extract_fcap(file),
      cohort = year - as.integer(age)
    )
  
  catch_df <- catch_df %>%
    mutate(
      scenario = extract_scenario(file),
      Fcap = extract_fcap(file),
      cohort = year - as.integer(age)
    )
  
  expl_df <- expl_df %>%
    mutate(
      scenario = extract_scenario(file),
      Fcap = extract_fcap(file),
      cohort = year - as.integer(age)
    )
  
  list(
    catch = catch_df,
    exploitation = expl_df,
    stock = stock_df
  )
}

files <- list.files(
  res_dir,
  pattern = "sc[0-9]+_Fcap_.*\\.rds$",
  full.names = TRUE
)

files

all_data <- lapply(files, extract_mse_quantities)

catch_all <- bind_rows(lapply(all_data, `[[`, "catch"))
expl_all  <- bind_rows(lapply(all_data, `[[`, "exploitation"))
stock_all <- bind_rows(lapply(all_data, `[[`, "stock"))

# Catch-at-age anual
catch_age_annual <- catch_all %>%
  group_by(scenario, Fcap, year, age, iter) %>%
  summarise(catch_n = sum(catch_n, na.rm = TRUE), .groups = "drop")

# catch-at-age composition
prop_age <- catch_age_annual %>%
  group_by(scenario, Fcap, year, iter) %>%
  mutate(
    total_catch_n = sum(catch_n, na.rm = TRUE),
    prop_age = catch_n / total_catch_n
  ) %>%
  ungroup() %>%
  mutate(prop_age = ifelse(is.finite(prop_age), prop_age, NA))


# Referencia histórica
hist_prop_age <- prop_age %>%
  filter(year %in% hist_years) %>%
  group_by(age) %>%
  summarise(
    hist_mean = mean(prop_age, na.rm = TRUE),
    hist_p05  = quantile(prop_age, 0.05, na.rm = TRUE),
    hist_p95  = quantile(prop_age, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

hist_prop_age


# Proyección por Fcap
proj_prop_age <- prop_age %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, age) %>%
  summarise(
    proj_mean = mean(prop_age, na.rm = TRUE),
    proj_p05  = quantile(prop_age, 0.05, na.rm = TRUE),
    proj_p95  = quantile(prop_age, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

comparison_age <- proj_prop_age %>%
  left_join(hist_prop_age, by = "age") %>%
  mutate(
    outside_hist_range = proj_mean < hist_p05 | proj_mean > hist_p95,
    deviation_from_hist = abs(proj_mean - hist_mean)
  )

comparison_age

# Dependencia en cohortes jóvenes: edad 0-1
young_dep <- prop_age %>%
  mutate(age_group = case_when(
    age %in% c("0", "1") ~ "age0_1",
    age %in% c("2", "3") ~ "age2_3",
    TRUE ~ "other"
  )) %>%
  group_by(scenario, Fcap, year, iter, age_group) %>%
  summarise(prop = sum(prop_age, na.rm = TRUE), .groups = "drop")

hist_young <- young_dep %>%
  filter(year %in% hist_years, age_group == "age0_1") %>%
  summarise(
    hist_mean_age0_1 = mean(prop, na.rm = TRUE),
    hist_p95_age0_1  = quantile(prop, 0.95, na.rm = TRUE)
  )

proj_young <- young_dep %>%
  filter(year %in% proj_years, age_group == "age0_1") %>%
  group_by(Fcap) %>%
  summarise(
    proj_mean_age0_1 = mean(prop, na.rm = TRUE),
    proj_p95_age0_1  = quantile(prop, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hist_mean_age0_1 = hist_young$hist_mean_age0_1,
    hist_p95_age0_1  = hist_young$hist_p95_age0_1,
    excessive_young_dependence = proj_mean_age0_1 > hist_p95_age0_1
  )

proj_young

# Variabilidad de capturas por edad
cv_catch_age <- catch_age_annual %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, age, iter) %>%
  summarise(
    mean_catch = mean(catch_n, na.rm = TRUE),
    sd_catch = sd(catch_n, na.rm = TRUE),
    cv_catch = sd_catch / mean_catch,
    .groups = "drop"
  ) %>%
  mutate(cv_catch = ifelse(is.finite(cv_catch), cv_catch, NA)) %>%
  group_by(Fcap, age) %>%
  summarise(
    mean_cv = mean(cv_catch, na.rm = TRUE),
    median_cv = median(cv_catch, na.rm = TRUE),
    p95_cv = quantile(cv_catch, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

cv_catch_age

# F-at-age aproximado: explotación por edad
u_age_annual <- expl_all %>%
  group_by(scenario, Fcap, year, age, iter) %>%
  summarise(u_age = mean(u_age, na.rm = TRUE), .groups = "drop")

u_age_summary <- u_age_annual %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, age) %>%
  summarise(
    mean_u = mean(u_age, na.rm = TRUE),
    median_u = median(u_age, na.rm = TRUE),
    p95_u = quantile(u_age, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

u_age_summary

# Tabla final para decidir rango realista
fcap_bio_score <- comparison_age %>%
  group_by(Fcap) %>%
  summarise(
    n_ages_outside_hist = sum(outside_hist_range, na.rm = TRUE),
    mean_deviation_age_comp = mean(deviation_from_hist, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(proj_young, by = "Fcap") %>%
  mutate(
    biologically_realistic = case_when(
      excessive_young_dependence ~ FALSE,
      n_ages_outside_hist >= 2 ~ FALSE,
      TRUE ~ TRUE
    )
  )

fcap_bio_score

write_csv(fcap_bio_score, file.path(out_dir, "fcap_biological_realism_score.csv"))
write_csv(comparison_age, file.path(out_dir, "catch_age_composition_vs_historical.csv"))
write_csv(cv_catch_age, file.path(out_dir, "cv_catch_at_age_by_fcap.csv"))
write_csv(u_age_summary, file.path(out_dir, "exploitation_proxy_at_age_by_fcap.csv"))

# Figuras
p1 <- prop_age %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, age) %>%
  summarise(prop_age = mean(prop_age, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = factor(Fcap), y = prop_age, fill = age)) +
  geom_col(position = "stack") +
  labs(
    x = "Fcap",
    y = "Mean catch-at-age proportion",
    fill = "Age",
    title = "Projected catch-at-age composition by Fcap"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "catch_at_age_composition_by_fcap.png"),
       p1, width = 8, height = 5, dpi = 300)

p2 <- cv_catch_age %>%
  ggplot(aes(x = factor(Fcap), y = mean_cv, group = age, colour = age)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Fcap",
    y = "Mean CV of catch-at-age",
    colour = "Age",
    title = "Variability of catch-at-age by Fcap"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "cv_catch_at_age_by_fcap.png"),
       p2, width = 8, height = 5, dpi = 300)

p3 <- u_age_summary %>%
  ggplot(aes(x = factor(Fcap), y = mean_u, group = age, colour = age)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Fcap",
    y = "Mean exploitation proxy catch.n / stock.n",
    colour = "Age",
    title = "Approximate exploitation-at-age by Fcap"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "exploitation_proxy_at_age_by_fcap.png"),
       p3, width = 8, height = 5, dpi = 300)



# ============================================================
# 1. Historical vs projected catch-at-age composition
# ============================================================

hist_prop_age_ref <- prop_age %>%
  filter(year %in% hist_years) %>%
  group_by(age) %>%
  summarise(
    hist_mean = mean(prop_age, na.rm = TRUE),
    hist_p05  = quantile(prop_age, 0.05, na.rm = TRUE),
    hist_p95  = quantile(prop_age, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

proj_prop_age_ref <- prop_age %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, age) %>%
  summarise(
    proj_mean = mean(prop_age, na.rm = TRUE),
    proj_p05  = quantile(prop_age, 0.05, na.rm = TRUE),
    proj_p95  = quantile(prop_age, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

catch_age_hist_comparison <- proj_prop_age_ref %>%
  left_join(hist_prop_age_ref, by = "age") %>%
  mutate(
    outside_hist_p95 = proj_mean > hist_p95,
    outside_hist_p05 = proj_mean < hist_p05,
    outside_hist_range = outside_hist_p95 | outside_hist_p05,
    ratio_to_hist_mean = proj_mean / hist_mean,
    diff_from_hist_mean = proj_mean - hist_mean
  )

write_csv(
  catch_age_hist_comparison,
  file.path(out_dir, "historical_vs_projected_catch_age_composition.csv")
)

catch_age_hist_comparison

# ============================================================
# 2. Historical vs projected exploitation proxy at age
#    u_age = catch.n / stock.n
# ============================================================

u_age_annual <- expl_all %>%
  group_by(scenario, Fcap, year, age, iter) %>%
  summarise(
    u_age = mean(u_age, na.rm = TRUE),
    .groups = "drop"
  )

hist_u_age_ref <- u_age_annual %>%
  filter(year %in% hist_years) %>%
  group_by(age) %>%
  summarise(
    hist_mean_u = mean(u_age, na.rm = TRUE),
    hist_p05_u  = quantile(u_age, 0.05, na.rm = TRUE),
    hist_p95_u  = quantile(u_age, 0.95, na.rm = TRUE),
    hist_max_u  = max(u_age, na.rm = TRUE),
    .groups = "drop"
  )

proj_u_age_ref <- u_age_annual %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, age) %>%
  summarise(
    proj_mean_u = mean(u_age, na.rm = TRUE),
    proj_p05_u  = quantile(u_age, 0.05, na.rm = TRUE),
    proj_p95_u  = quantile(u_age, 0.95, na.rm = TRUE),
    proj_max_u  = max(u_age, na.rm = TRUE),
    .groups = "drop"
  )

u_age_hist_comparison <- proj_u_age_ref %>%
  left_join(hist_u_age_ref, by = "age") %>%
  mutate(
    mean_above_hist_p95 = proj_mean_u > hist_p95_u,
    p95_above_hist_p95  = proj_p95_u > hist_p95_u,
    mean_above_hist_max = proj_mean_u > hist_max_u,
    ratio_to_hist_mean_u = proj_mean_u / hist_mean_u,
    diff_from_hist_mean_u = proj_mean_u - hist_mean_u
  )

write_csv(
  u_age_hist_comparison,
  file.path(out_dir, "historical_vs_projected_exploitation_proxy_at_age.csv")
)

u_age_hist_comparison

# ============================================================
# 3. Biological realism score by Fcap
# ============================================================

fcap_hist_score <- u_age_hist_comparison %>%
  group_by(Fcap) %>%
  summarise(
    n_ages_mean_u_above_hist_p95 = sum(mean_above_hist_p95, na.rm = TRUE),
    n_ages_p95_u_above_hist_p95  = sum(p95_above_hist_p95, na.rm = TRUE),
    max_ratio_to_hist_mean_u = max(ratio_to_hist_mean_u, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    catch_age_hist_comparison %>%
      group_by(Fcap) %>%
      summarise(
        n_ages_catch_comp_outside_hist = sum(outside_hist_range, na.rm = TRUE),
        max_abs_diff_catch_comp = max(abs(diff_from_hist_mean), na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Fcap"
  ) %>%
  mutate(
    biological_flag = case_when(
      n_ages_mean_u_above_hist_p95 >= 2 ~ "outside historical exploitation range",
      n_ages_catch_comp_outside_hist >= 2 ~ "outside historical catch-at-age composition",
      TRUE ~ "within historical reference"
    )
  )

write_csv(
  fcap_hist_score,
  file.path(out_dir, "fcap_historical_biological_score.csv")
)

fcap_hist_score



p_hist_catch <- catch_age_hist_comparison %>%
  ggplot(aes(x = factor(Fcap), y = proj_mean, fill = age)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(ymin = hist_p05, ymax = hist_p95),
    position = position_dodge(width = 0.9),
    width = 0.25,
    linewidth = 0.6
  ) +
  facet_wrap(~ age, scales = "free_y") +
  labs(
    x = "Fcap",
    y = "Projected mean catch-at-age proportion",
    title = "Projected catch-at-age composition compared with historical range",
    subtitle = "Bars = projected mean; error bars = historical 5–95% range"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "historical_vs_projected_catch_age_composition.png"),
  p_hist_catch,
  width = 10,
  height = 6,
  dpi = 300
)

# Exploitation-at-age vs histórico
p_hist_u <- u_age_hist_comparison %>%
  ggplot(aes(x = factor(Fcap), y = proj_mean_u, fill = age)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(ymin = hist_p05_u, ymax = hist_p95_u),
    position = position_dodge(width = 0.9),
    width = 0.25,
    linewidth = 0.6
  ) +
  facet_wrap(~ age, scales = "free_y") +
  labs(
    x = "Fcap",
    y = "Projected mean exploitation proxy at age",
    title = "Projected exploitation-at-age compared with historical range",
    subtitle = "Bars = projected mean; error bars = historical 5–95% range"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "historical_vs_projected_exploitation_proxy_at_age.png"),
  p_hist_u,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# 4. Cohort survivorship diagnostics
#    N2/N1 and N3/N1
# ============================================================

stock_age_annual <- stock_all %>%
  group_by(scenario, Fcap, year, age, iter, cohort) %>%
  summarise(
    stock_n = mean(stock_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(age = as.integer(age))

# Wide format: N at age by year, iter and Fcap
stock_wide <- stock_age_annual %>%
  filter(age %in% c(1, 2, 3)) %>%
  select(scenario, Fcap, year, iter, age, stock_n) %>%
  pivot_wider(
    names_from = age,
    values_from = stock_n,
    names_prefix = "N_age"
  ) %>%
  mutate(
    N2_N1 = N_age2 / N_age1,
    N3_N1 = N_age3 / N_age1,
    N3_N2 = N_age3 / N_age2,
    N_old_N1 = (N_age2 + N_age3) / N_age1
  ) %>%
  mutate(
    across(
      c(N2_N1, N3_N1, N3_N2, N_old_N1),
      ~ ifelse(is.finite(.x), .x, NA_real_)
    )
  )

# Historical reference
hist_survivorship_ref <- stock_wide %>%
  filter(year %in% hist_years) %>%
  summarise(
    hist_mean_N2_N1 = mean(N2_N1, na.rm = TRUE),
    hist_p05_N2_N1  = quantile(N2_N1, 0.05, na.rm = TRUE),
    hist_p95_N2_N1  = quantile(N2_N1, 0.95, na.rm = TRUE),
    
    hist_mean_N3_N1 = mean(N3_N1, na.rm = TRUE),
    hist_p05_N3_N1  = quantile(N3_N1, 0.05, na.rm = TRUE),
    hist_p95_N3_N1  = quantile(N3_N1, 0.95, na.rm = TRUE),
    
    hist_mean_N3_N2 = mean(N3_N2, na.rm = TRUE),
    hist_p05_N3_N2  = quantile(N3_N2, 0.05, na.rm = TRUE),
    hist_p95_N3_N2  = quantile(N3_N2, 0.95, na.rm = TRUE),
    
    hist_mean_N_old_N1 = mean(N_old_N1, na.rm = TRUE),
    hist_p05_N_old_N1  = quantile(N_old_N1, 0.05, na.rm = TRUE),
    hist_p95_N_old_N1  = quantile(N_old_N1, 0.95, na.rm = TRUE)
  )

hist_survivorship_ref

# Comparación proyeccion vs históricos
proj_survivorship_ref <- stock_wide %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap) %>%
  summarise(
    proj_mean_N2_N1 = mean(N2_N1, na.rm = TRUE),
    proj_p05_N2_N1  = quantile(N2_N1, 0.05, na.rm = TRUE),
    proj_p95_N2_N1  = quantile(N2_N1, 0.95, na.rm = TRUE),
    
    proj_mean_N3_N1 = mean(N3_N1, na.rm = TRUE),
    proj_p05_N3_N1  = quantile(N3_N1, 0.05, na.rm = TRUE),
    proj_p95_N3_N1  = quantile(N3_N1, 0.95, na.rm = TRUE),
    
    proj_mean_N3_N2 = mean(N3_N2, na.rm = TRUE),
    proj_p05_N3_N2  = quantile(N3_N2, 0.05, na.rm = TRUE),
    proj_p95_N3_N2  = quantile(N3_N2, 0.95, na.rm = TRUE),
    
    proj_mean_N_old_N1 = mean(N_old_N1, na.rm = TRUE),
    proj_p05_N_old_N1  = quantile(N_old_N1, 0.05, na.rm = TRUE),
    proj_p95_N_old_N1  = quantile(N_old_N1, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hist_mean_N2_N1 = hist_survivorship_ref$hist_mean_N2_N1,
    hist_p05_N2_N1  = hist_survivorship_ref$hist_p05_N2_N1,
    hist_p95_N2_N1  = hist_survivorship_ref$hist_p95_N2_N1,
    
    hist_mean_N3_N1 = hist_survivorship_ref$hist_mean_N3_N1,
    hist_p05_N3_N1  = hist_survivorship_ref$hist_p05_N3_N1,
    hist_p95_N3_N1  = hist_survivorship_ref$hist_p95_N3_N1,
    
    hist_mean_N3_N2 = hist_survivorship_ref$hist_mean_N3_N2,
    hist_p05_N3_N2  = hist_survivorship_ref$hist_p05_N3_N2,
    hist_p95_N3_N2  = hist_survivorship_ref$hist_p95_N3_N2,
    
    hist_mean_N_old_N1 = hist_survivorship_ref$hist_mean_N_old_N1,
    hist_p05_N_old_N1  = hist_survivorship_ref$hist_p05_N_old_N1,
    hist_p95_N_old_N1  = hist_survivorship_ref$hist_p95_N_old_N1
  ) %>%
  mutate(
    N2_N1_below_hist_p05 = proj_mean_N2_N1 < hist_p05_N2_N1,
    N3_N1_below_hist_p05 = proj_mean_N3_N1 < hist_p05_N3_N1,
    N3_N2_below_hist_p05 = proj_mean_N3_N2 < hist_p05_N3_N2,
    N_old_N1_below_hist_p05 = proj_mean_N_old_N1 < hist_p05_N_old_N1
  )

proj_survivorship_ref

write_csv(
  proj_survivorship_ref,
  file.path(out_dir, "historical_vs_projected_cohort_survivorship.csv")
)

# Pasarlo a formato largo para figuras
surv_long <- stock_wide %>%
  select(scenario, Fcap, year, iter, N2_N1, N3_N1, N3_N2, N_old_N1) %>%
  pivot_longer(
    cols = c(N2_N1, N3_N1, N3_N2, N_old_N1),
    names_to = "indicator",
    values_to = "value"
  )

hist_surv_long <- surv_long %>%
  filter(year %in% hist_years) %>%
  group_by(indicator) %>%
  summarise(
    hist_mean = mean(value, na.rm = TRUE),
    hist_p05 = quantile(value, 0.05, na.rm = TRUE),
    hist_p95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

proj_surv_long <- surv_long %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, indicator) %>%
  summarise(
    proj_mean = mean(value, na.rm = TRUE),
    proj_p05 = quantile(value, 0.05, na.rm = TRUE),
    proj_p95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(hist_surv_long, by = "indicator") %>%
  mutate(
    below_hist_p05 = proj_mean < hist_p05,
    above_hist_p95 = proj_mean > hist_p95,
    outside_hist_range = below_hist_p05 | above_hist_p95
  )

write_csv(
  proj_surv_long,
  file.path(out_dir, "cohort_survivorship_long.csv")
)

# Figuras
p_surv <- proj_surv_long %>%
  ggplot(aes(x = factor(Fcap), y = proj_mean)) +
  geom_col() +
  geom_errorbar(
    aes(ymin = hist_p05, ymax = hist_p95),
    width = 0.25,
    linewidth = 0.6
  ) +
  facet_wrap(~ indicator, scales = "free_y") +
  labs(
    x = "Fcap",
    y = "Projected mean cohort survivorship ratio",
    title = "Projected cohort survivorship compared with historical range",
    subtitle = "Bars = projected mean; error bars = historical 5–95% range"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "historical_vs_projected_cohort_survivorship.png"),
  p_surv,
  width = 10,
  height = 6,
  dpi = 300
)

surv_score <- proj_surv_long %>%
  group_by(Fcap) %>%
  summarise(
    n_survivorship_indicators_below_hist_p05 =
      sum(below_hist_p05, na.rm = TRUE),
    n_survivorship_indicators_outside_hist =
      sum(outside_hist_range, na.rm = TRUE),
    .groups = "drop"
  )

fcap_hist_score_extended <- fcap_hist_score %>%
  left_join(surv_score, by = "Fcap") %>%
  mutate(
    biological_flag_extended = case_when(
      n_ages_mean_u_above_hist_p95 >= 2 &
        n_survivorship_indicators_below_hist_p05 >= 1 ~
        "outside historical exploitation and reduced cohort survivorship",
      
      n_ages_mean_u_above_hist_p95 >= 2 ~
        "outside historical exploitation range",
      
      n_survivorship_indicators_below_hist_p05 >= 1 ~
        "reduced cohort survivorship",
      
      TRUE ~
        "within historical biological envelope"
    )
  )

write_csv(
  fcap_hist_score_extended,
  file.path(out_dir, "fcap_historical_biological_score_extended.csv")
)


# ============================================================
# COHORT PERSISTENCE ANALYSIS
# ============================================================

stock_age_annual <- stock_all %>%
  group_by(scenario, Fcap, year, age, iter) %>%
  summarise(stock_n = mean(stock_n, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    age = as.integer(age),
    cohort = year - age
  )

# Follow same cohort across ages
cohort_wide <- stock_age_annual %>%
  filter(age %in% 0:3) %>%
  select(scenario, Fcap, cohort, age, iter, stock_n) %>%
  pivot_wider(
    names_from = age,
    values_from = stock_n,
    names_prefix = "N_age"
  ) %>%
  mutate(
    surv_1_to_2 = N_age2 / N_age1,
    surv_2_to_3 = N_age3 / N_age2,
    surv_1_to_3 = N_age3 / N_age1,
    old_fraction = (N_age2 + N_age3) / (N_age0 + N_age1 + N_age2 + N_age3)
  ) %>%
  mutate(
    across(
      c(surv_1_to_2, surv_2_to_3, surv_1_to_3, old_fraction),
      ~ ifelse(is.finite(.x), .x, NA_real_)
    )
  )

# Historical reference
hist_cohort_ref <- cohort_wide %>%
  filter(cohort %in% hist_years) %>%
  summarise(
    hist_mean_surv_1_to_2 = mean(surv_1_to_2, na.rm = TRUE),
    hist_p05_surv_1_to_2  = quantile(surv_1_to_2, 0.05, na.rm = TRUE),
    hist_p95_surv_1_to_2  = quantile(surv_1_to_2, 0.95, na.rm = TRUE),
    
    hist_mean_surv_2_to_3 = mean(surv_2_to_3, na.rm = TRUE),
    hist_p05_surv_2_to_3  = quantile(surv_2_to_3, 0.05, na.rm = TRUE),
    hist_p95_surv_2_to_3  = quantile(surv_2_to_3, 0.95, na.rm = TRUE),
    
    hist_mean_surv_1_to_3 = mean(surv_1_to_3, na.rm = TRUE),
    hist_p05_surv_1_to_3  = quantile(surv_1_to_3, 0.05, na.rm = TRUE),
    hist_p95_surv_1_to_3  = quantile(surv_1_to_3, 0.95, na.rm = TRUE),
    
    hist_mean_old_fraction = mean(old_fraction, na.rm = TRUE),
    hist_p05_old_fraction  = quantile(old_fraction, 0.05, na.rm = TRUE),
    hist_p95_old_fraction  = quantile(old_fraction, 0.95, na.rm = TRUE)
  )

# Projected cohort persistence
proj_cohort_ref <- cohort_wide %>%
  filter(cohort %in% proj_years) %>%
  group_by(Fcap) %>%
  summarise(
    proj_mean_surv_1_to_2 = mean(surv_1_to_2, na.rm = TRUE),
    proj_mean_surv_2_to_3 = mean(surv_2_to_3, na.rm = TRUE),
    proj_mean_surv_1_to_3 = mean(surv_1_to_3, na.rm = TRUE),
    proj_mean_old_fraction = mean(old_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hist_p05_surv_1_to_2 = hist_cohort_ref$hist_p05_surv_1_to_2,
    hist_p05_surv_2_to_3 = hist_cohort_ref$hist_p05_surv_2_to_3,
    hist_p05_surv_1_to_3 = hist_cohort_ref$hist_p05_surv_1_to_3,
    hist_p05_old_fraction = hist_cohort_ref$hist_p05_old_fraction,
    
    surv_1_to_2_below_hist = proj_mean_surv_1_to_2 < hist_p05_surv_1_to_2,
    surv_2_to_3_below_hist = proj_mean_surv_2_to_3 < hist_p05_surv_2_to_3,
    surv_1_to_3_below_hist = proj_mean_surv_1_to_3 < hist_p05_surv_1_to_3,
    old_fraction_below_hist = proj_mean_old_fraction < hist_p05_old_fraction
  )

write_csv(
  proj_cohort_ref,
  file.path(out_dir, "cohort_persistence_by_fcap.csv")
)

proj_cohort_ref

cohort_long <- cohort_wide %>%
  select(
    scenario, Fcap, cohort, iter,
    surv_1_to_2, surv_2_to_3, surv_1_to_3, old_fraction
  ) %>%
  pivot_longer(
    cols = c(surv_1_to_2, surv_2_to_3, surv_1_to_3, old_fraction),
    names_to = "indicator",
    values_to = "value"
  )

hist_cohort_long <- cohort_long %>%
  filter(cohort %in% hist_years) %>%
  group_by(indicator) %>%
  summarise(
    hist_p05 = quantile(value, 0.05, na.rm = TRUE),
    hist_p95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

proj_cohort_long <- cohort_long %>%
  filter(cohort %in% proj_years) %>%
  group_by(Fcap, indicator) %>%
  summarise(
    proj_mean = mean(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(hist_cohort_long, by = "indicator") %>%
  mutate(
    below_hist_p05 = proj_mean < hist_p05
  )

p_cohort <- ggplot(
  proj_cohort_long,
  aes(x = factor(Fcap), y = proj_mean)
) +
  geom_col() +
  geom_errorbar(
    aes(ymin = hist_p05, ymax = hist_p95),
    width = 0.25,
    linewidth = 0.6
  ) +
  facet_wrap(~ indicator, scales = "free_y") +
  labs(
    x = "Fcap",
    y = "Projected mean value",
    title = "Cohort persistence compared with historical range",
    subtitle = "Bars = projected mean; error bars = historical 5–95% range"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "cohort_persistence_by_fcap.png"),
  p_cohort,
  width = 10,
  height = 6,
  dpi = 300)

# ============================================================
# 7. Critical years and amplification of weak cohorts
# ============================================================

# 7.1 Summary of low recruitment years by Fcap and age
low_rec_summary <- low_rec_years %>%
  group_by(Fcap, age) %>%
  summarise(
    n_low_years = n_distinct(year),
    mean_low_n = mean(median_n, na.rm = TRUE),
    median_low_n = median(median_n, na.rm = TRUE),
    low_years = paste(sort(unique(year)), collapse = ", "),
    .groups = "drop"
  )

# Compare low recruitment abundance relative to Fcap = 1
low_rec_relative <- low_rec_summary %>%
  left_join(
    low_rec_summary %>%
      filter(Fcap == 1) %>%
      select(age, mean_low_n_Fcap1 = mean_low_n),
    by = "age"
  ) %>%
  mutate(
    ratio_to_Fcap1 = mean_low_n / mean_low_n_Fcap1,
    reduction_vs_Fcap1 = 1 - ratio_to_Fcap1
  )

write_csv(
  low_rec_relative,
  file.path(out_dir, "low_recruitment_years_relative_to_Fcap1.csv")
)

low_rec_relative

# Alta explotación durante años debiles
# ============================================================
# 8. Exploitation during weak recruitment years
# ============================================================

weak_years <- low_rec_years %>%
  distinct(Fcap, year)

u_weak_vs_normal <- u_age_annual %>%
  filter(year %in% proj_years) %>%
  left_join(
    weak_years %>% mutate(weak_year = TRUE),
    by = c("Fcap", "year")
  ) %>%
  mutate(
    weak_year = ifelse(is.na(weak_year), FALSE, weak_year)
  ) %>%
  group_by(Fcap, weak_year) %>%
  summarise(
    mean_u = mean(u_age, na.rm = TRUE),
    median_u = median(u_age, na.rm = TRUE),
    p95_u = quantile(u_age, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = weak_year,
    values_from = c(mean_u, median_u, p95_u),
    names_glue = "{.value}_{ifelse(weak_year, 'weak', 'normal')}"
  ) %>%
  mutate(
    amplification_mean_u = mean_u_weak / mean_u_normal,
    amplification_p95_u = p95_u_weak / p95_u_normal
  )

write_csv(
  u_weak_vs_normal,
  file.path(out_dir, "exploitation_amplification_in_weak_years.csv")
)

u_weak_vs_normal

# Persistencia de cohortes relativa a Fcap = 1
# ============================================================
# 9. Cohort persistence relative to Fcap = 1
# ============================================================

cohort_resilience <- cohort_wide %>%
  filter(cohort %in% proj_years) %>%
  group_by(Fcap) %>%
  summarise(
    mean_surv_1_to_2 = mean(surv_1_to_2, na.rm = TRUE),
    mean_surv_2_to_3 = mean(surv_2_to_3, na.rm = TRUE),
    mean_surv_1_to_3 = mean(surv_1_to_3, na.rm = TRUE),
    mean_old_fraction = mean(old_fraction, na.rm = TRUE),
    
    p05_surv_1_to_2 = quantile(surv_1_to_2, 0.05, na.rm = TRUE),
    p05_surv_2_to_3 = quantile(surv_2_to_3, 0.05, na.rm = TRUE),
    p05_surv_1_to_3 = quantile(surv_1_to_3, 0.05, na.rm = TRUE),
    p05_old_fraction = quantile(old_fraction, 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rel_surv_1_to_2 = mean_surv_1_to_2 / mean_surv_1_to_2[Fcap == 1],
    rel_surv_2_to_3 = mean_surv_2_to_3 / mean_surv_2_to_3[Fcap == 1],
    rel_surv_1_to_3 = mean_surv_1_to_3 / mean_surv_1_to_3[Fcap == 1],
    rel_old_fraction = mean_old_fraction / mean_old_fraction[Fcap == 1],
    
    reduction_surv_1_to_2 = 1 - rel_surv_1_to_2,
    reduction_surv_2_to_3 = 1 - rel_surv_2_to_3,
    reduction_surv_1_to_3 = 1 - rel_surv_1_to_3,
    reduction_old_fraction = 1 - rel_old_fraction
  )

write_csv(
  cohort_resilience,
  file.path(out_dir, "cohort_resilience_relative_to_Fcap1.csv")
)

cohort_resilience

# reducción de persistencia

cohort_resilience_long <- cohort_resilience %>%
  select(
    Fcap,
    reduction_surv_1_to_2,
    reduction_surv_2_to_3,
    reduction_surv_1_to_3,
    reduction_old_fraction
  ) %>%
  pivot_longer(
    cols = -Fcap,
    names_to = "indicator",
    values_to = "reduction"
  )

p_resilience <- ggplot(
  cohort_resilience_long,
  aes(x = factor(Fcap), y = reduction)
) +
  geom_col() +
  facet_wrap(~ indicator, scales = "free_y") +
  labs(
    x = "Fcap",
    y = "Reduction relative to Fcap = 1",
    title = "Reduction in cohort persistence relative to Fcap = 1"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "cohort_resilience_reduction_relative_to_Fcap1.png"),
  p_resilience,
  width = 10,
  height = 6,
  dpi = 300
)

# ============================================================
# 10. Final biological interpretation table
# ============================================================

final_bio_summary <- fcap_hist_score_extended %>%
  left_join(
    cohort_resilience %>%
      select(
        Fcap,
        rel_surv_1_to_2,
        rel_surv_2_to_3,
        rel_surv_1_to_3,
        rel_old_fraction
      ),
    by = "Fcap"
  ) %>%
  left_join(
    u_weak_vs_normal %>%
      select(Fcap, amplification_mean_u, amplification_p95_u),
    by = "Fcap"
  ) %>%
  mutate(
    biological_interpretation = case_when(
      Fcap <= 1 ~
        "within historical biological envelope",
      
      Fcap > 1 & Fcap < 1.5 ~
        "transition zone; slight increase in exploitation and reduced cohort persistence",
      
      Fcap >= 1.5 ~
        "outside historical exploitation envelope; reduced cohort persistence",
      
      TRUE ~ NA_character_
    )
  )

write_csv(
  final_bio_summary,
  file.path(out_dir, "final_biological_realism_summary_by_Fcap.csv")
)

final_bio_summary

as.data.frame(final_bio_summary)


# ============================================================
# Simple Figure A: ratio to historical P95
# ============================================================

u_ratio_p95 <- u_age_hist_comparison %>%
  mutate(
    age = factor(age, levels = c("0", "1", "2", "3")),
    ratio_to_hist_p95 = proj_mean_u / hist_p95_u
  )

p_u_ratio <- ggplot(
  u_ratio_p95,
  aes(x = factor(Fcap), y = ratio_to_hist_p95, group = age, colour = age)
) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  labs(
    x = "Fcap",
    y = "Projected mean U-at-age / historical P95",
    colour = "Age",
    title = "Exploitation-at-age relative to the historical upper range",
    subtitle = "Values above 1 exceed the historical 95th percentile"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "exploitation_ratio_to_historical_p95.png"),
  p_u_ratio,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# Simple Figure B: number of ages exceeding historical P95
# ============================================================

u_exceed_count <- u_age_hist_comparison %>%
  group_by(Fcap) %>%
  summarise(
    n_ages_above_hist_p95 = sum(mean_above_hist_p95, na.rm = TRUE),
    .groups = "drop"
  )

p_u_count <- ggplot(
  u_exceed_count,
  aes(x = factor(Fcap), y = n_ages_above_hist_p95)
) +
  geom_col() +
  scale_y_continuous(
    breaks = 0:4,
    limits = c(0, 4)
  ) +
  labs(
    x = "Fcap",
    y = "Number of age classes",
    title = "Age classes exceeding the historical exploitation range",
    subtitle = "Count of ages where projected mean U-at-age exceeds historical P95"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "number_ages_exceeding_historical_p95.png"),
  p_u_count,
  width = 7,
  height = 5,
  dpi = 300
)

#######################################################################3
u_age_season <- expl_all %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, season, age) %>%
  summarise(
    proj_mean_u = mean(u_age, na.rm = TRUE),
    .groups = "drop"
  )

hist_u_age_season <- expl_all %>%
  filter(year %in% hist_years) %>%
  group_by(season, age) %>%
  summarise(
    hist_p95_u = quantile(u_age, 0.95, na.rm = TRUE),
    hist_mean_u = mean(u_age, na.rm = TRUE),
    .groups = "drop"
  )

u_age_season_compare <- u_age_season %>%
  left_join(hist_u_age_season, by = c("season", "age")) %>%
  mutate(
    ratio_to_hist_p95 = proj_mean_u / hist_p95_u,
    above_hist_p95 = ratio_to_hist_p95 > 1
  )

p_season_ratio <- ggplot(
  u_age_season_compare,
  aes(x = factor(Fcap), y = ratio_to_hist_p95, group = age, colour = age)
) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line() +
  geom_point() +
  facet_wrap(~ season) +
  labs(
    x = "Fcap",
    y = "Projected U / historical seasonal P95",
    colour = "Age",
    title = "Seasonal exploitation-at-age relative to historical range"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_exploitation_ratio_to_hist_p95.png"),
  p_season_ratio,
  width = 10,
  height = 6,
  dpi = 300
)

################################################################################
season_exceed_count <- u_age_season_compare %>%
  group_by(Fcap, season) %>%
  summarise(
    n_age_classes_above_p95 = sum(above_hist_p95, na.rm = TRUE),
    .groups = "drop"
  )

p_season_count <- ggplot(
  season_exceed_count,
  aes(x = factor(Fcap), y = n_age_classes_above_p95)
) +
  geom_col() +
  facet_wrap(~ season) +
  scale_y_continuous(breaks = 0:4, limits = c(0, 4)) +
  labs(
    x = "Fcap",
    y = "Number of ages above historical P95",
    title = "Seasonal age classes exceeding historical exploitation range"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_number_ages_exceeding_hist_p95.png"),
  p_season_count,
  width = 10,
  height = 6,
  dpi = 300
)
################################################################################
season_catch_comp <- catch_all %>%
  filter(year %in% proj_years) %>%
  group_by(Fcap, season, age, iter) %>%
  summarise(catch_n = sum(catch_n, na.rm = TRUE), .groups = "drop") %>%
  group_by(Fcap, season, iter) %>%
  mutate(
    prop_age = catch_n / sum(catch_n, na.rm = TRUE)
  ) %>%
  ungroup()

season_catch_comp_summary <- season_catch_comp %>%
  group_by(Fcap, season, age) %>%
  summarise(
    mean_prop = mean(prop_age, na.rm = TRUE),
    .groups = "drop"
  )

p_season_catch_comp <- ggplot(
  season_catch_comp_summary,
  aes(x = factor(Fcap), y = mean_prop, fill = age)
) +
  geom_col(position = "fill") +
  facet_wrap(~ season) +
  labs(
    x = "Fcap",
    y = "Catch-at-age composition",
    fill = "Age",
    title = "Projected seasonal catch-at-age composition by Fcap"
  ) +
  theme_bw()

ggsave(
  file.path(out_dir, "seasonal_catch_at_age_composition_by_fcap.png"),
  p_season_catch_comp,
  width = 10,
  height = 6,
  dpi = 300
)
