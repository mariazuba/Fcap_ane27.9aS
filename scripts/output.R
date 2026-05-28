rm(list = ls())

library(FLBEIA)
library(FLCore)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(here)

# ============================================================
# 1. CONFIGURACIÓN
# ============================================================

wd <- here()
setwd(wd)

load("data/FLBEIA_inputs/FLBEIA_inputs_Fcap_ane9aS.RData")
load("outputs/Fcap_results/FLBEIAshiny_Fcap_inputs.RData")

results_list <- readRDS("outputs/Fcap_results/results_list_Fcap.rds")
summary_Fcap <- readRDS("outputs/Fcap_results/summary_Fcap.rds")

dir.create("outputs/Fcap_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Fcap_tables",  recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PLOTS PARA EXPLICAR FLEET / FISHERY DYNAMICS
# ============================================================

library(FLBEIA)
library(FLCore)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("report/figs", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Seasonal catch share used to distribute annual TAC
# ------------------------------------------------------------
seasonal_share_df <- as.data.frame(fleets.ctrl$seasonal.share[[1]])
names(seasonal_share_df)[names(seasonal_share_df) == "data"] <- "share"

p_share <- seasonal_share_df %>%
  mutate(
    year = as.numeric(as.character(year)),
    season = as.factor(season)) %>%
  filter(year %in% c(hist.yrs, proj.yrs)) %>%
  group_by(year, season) %>%
  summarise(share = mean(share, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = share, colour = season)) +
  annotate(
    "rect",
    xmin = proj.yr,
    xmax = max(proj.yrs),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = proj.yr, linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 11) +
  labs(
    title = "Seasonal catch share",
    subtitle = "",
    x = NULL,
    y = "Proportion of annual catch",
    colour = "Season")

ggsave("report/figs/fleet_seasonal_catch_share.png",p_share, width = 7,height = 4,dpi = 300)



# ------------------------------------------------------------
# 3. Landings weight-at-age used in the fleet
# ------------------------------------------------------------

lwa_df <- as.data.frame(landings.wt(fleets$SEINE@metiers$ALL@catches$ANE))
names(lwa_df)[names(lwa_df) == "data"] <- "landings_wt"

p_lwt <- lwa_df %>%
  mutate(year = as.numeric(as.character(year))) %>%
  group_by(year, age, season) %>%
  summarise(landings_wt = mean(landings_wt, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = year, y = landings_wt, colour = as.factor(age))) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ season, scales = "free_y") +
  geom_vline(xintercept = proj.yr, linetype = "dashed") +
  theme_bw(base_size = 11) +
  labs(
    title = "Fleet landings weight-at-age",
    subtitle = "",
    x = NULL,
    y = "Landings weight",
    colour = "Age")

ggsave("report/figs/fleet_landings_weight_age_season.png",
       p_lwt, width = 8, height = 5, dpi = 300)


# ============================================================
# Catch-at-age by season
# ============================================================

catch_n_df <- as.data.frame(
  landings.n(fleets$SEINE@metiers$ALL@catches$ANE)
)

names(catch_n_df)[names(catch_n_df) == "data"] <- "catch_n"

p_catch_age <- catch_n_df %>%
  mutate(
    year = as.numeric(as.character(year)),
    age = as.factor(age),
    season = as.factor(season)
  ) %>%
  filter(year %in% c(hist.yrs, proj.yrs)) %>%
  group_by(year, age, season) %>%
  summarise(
    catch_n = mean(catch_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = year, y = catch_n, colour = age)) +
  annotate(
    "rect",
    xmin = proj.yr,
    xmax = max(proj.yrs),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08
  ) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = proj.yr, linetype = "dashed") +
  facet_wrap(~ season, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Catch-at-age by season",
    subtitle = "",
    x = NULL,
    y = "Catch in numbers",
    colour = "Age"
  )

ggsave(
  "report/figs/fleet_catch_at_age_by_season.png",
  p_catch_age,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# PLOT: Advice component - TAC scenarios
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("report/figs", recursive = TRUE, showWarnings = FALSE)

# Historical ICES TAC used in the advice object
tac_hist <- data.frame(
  year = 2019:2024,
  TAC = c(5278, 8856, 9459, 4383, 1892, 1733),
  scenario = "Historical ICES TAC"
)

# Fixed TAC scenarios
tac_scenarios <- expand.grid(
  year = proj.yrs,
  scenario = c("TAC = 0", "TAC = 8000", "TAC = 10000")
) %>%
  mutate(
    TAC = case_when(
      scenario == "TAC = 0" ~ 0,
      scenario == "TAC = 8000" ~ 8000,
      scenario == "TAC = 10000" ~ 10000
    )
  )

tac_df <- bind_rows(tac_hist, tac_scenarios)

p_advice_tac <- ggplot(tac_df, aes(x = year, y = TAC, colour = scenario)) +
  annotate(
    "rect",
    xmin = proj.yr,
    xmax = max(proj.yrs),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.08
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = proj.yr, linetype = "dashed") +
  theme_bw(base_size = 11) +
  labs(
    title = "Historical ICES TAC and fixed TAC scenarios used in the MSE",
    subtitle = "",
    x = NULL,
    y = "TAC (t)",
    colour = "Scenario"
  )

ggsave(
  "report/figs/advice_TAC_scenarios.png",
  p_advice_tac,
  width = 7,
  height = 4,
  dpi = 300
)



# ============================================================
# 2. FUNCIÓN PARA EXTRAER RESULTADOS POR ESCENARIO
# ============================================================

extract_Fcap_outputs <- function(res_i, scenario, Blim, Bpa) {
  
  biol_i <- res_i$biols$ANE
  
  # SSB por season
  ssb_all <- quantSums(
    biol_i@n *
      biol_i@wt *
      predict(biol_i@mat)
  )
  
  ssb_s2 <- ssb_all[, , , "2", , ]
  ssb_df <- as.data.frame(ssb_s2)
  names(ssb_df)[names(ssb_df) == "data"] <- "SSB"
  
  ssb_df <- ssb_df %>%
    mutate(
      scenario = scenario,
      Blim = Blim,
      Bpa  = Bpa
    )
  
  # TAC
  tac_df <- as.data.frame(res_i$advice$TAC)
  names(tac_df)[names(tac_df) == "data"] <- "TAC"
  
  tac_df <- tac_df %>%
    mutate(scenario = scenario)
  
  # F advice si existe
  if (!is.null(res_i$advice$Fadv)) {
    fadv_df <- as.data.frame(res_i$advice$Fadv)
    names(fadv_df)[names(fadv_df) == "data"] <- "Fadv"
    fadv_df <- fadv_df %>%
      mutate(scenario = scenario)
  } else {
    fadv_df <- NULL
  }
  
  list(
    ssb  = ssb_df,
    tac  = tac_df,
    fadv = fadv_df
  )
}

# ============================================================
# 3. EXTRAER TODO
# ============================================================

out_list <- lapply(names(results_list), function(sc) {
  extract_Fcap_outputs(
    res_i    = results_list[[sc]],
    scenario = sc,
    Blim     = Blim,
    Bpa      = Bpa
  )
})

names(out_list) <- names(results_list)

ssb_all <- bind_rows(lapply(out_list, `[[`, "ssb"))
tac_all <- bind_rows(lapply(out_list, `[[`, "tac"))
fadv_all <- bind_rows(lapply(out_list, `[[`, "fadv"))

# ============================================================
# 4. TABLAS RESUMEN
# ============================================================

ssb_summary <- ssb_all %>%
  group_by(scenario, year) %>%
  summarise(
    ssb_q05 = quantile(SSB, 0.05, na.rm = TRUE),
    ssb_q50 = quantile(SSB, 0.50, na.rm = TRUE),
    ssb_q95 = quantile(SSB, 0.95, na.rm = TRUE),
    risk_Blim = mean(SSB < Blim, na.rm = TRUE),
    risk_Bpa  = mean(SSB < Bpa,  na.rm = TRUE),
    Blim = first(Blim),
    Bpa  = first(Bpa),
    .groups = "drop"
  )

tac_summary <- tac_all %>%
  group_by(scenario, year) %>%
  summarise(
    tac_q05 = quantile(TAC, 0.05, na.rm = TRUE),
    tac_q50 = quantile(TAC, 0.50, na.rm = TRUE),
    tac_q95 = quantile(TAC, 0.95, na.rm = TRUE),
    tac_mean = mean(TAC, na.rm = TRUE),
    .groups = "drop"
  )

if (!is.null(fadv_all)) {
  fadv_summary <- fadv_all %>%
    group_by(scenario, year) %>%
    summarise(
      fadv_q05 = quantile(Fadv, 0.05, na.rm = TRUE),
      fadv_q50 = quantile(Fadv, 0.50, na.rm = TRUE),
      fadv_q95 = quantile(Fadv, 0.95, na.rm = TRUE),
      fadv_mean = mean(Fadv, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  fadv_summary <- NULL
}

risk_summary <- risk_season2_proj %>%
  group_by(scenario) %>%
  summarise(
    risk_Blim_max  = max(pBlim, na.rm = TRUE),
    risk_Blim_mean = mean(pBlim, na.rm = TRUE),
    risk_Bpa_max   = max(pBpa, na.rm = TRUE),
    risk_Bpa_mean  = mean(pBpa, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 4b. F DESDE advice$Fadv
# ============================================================

if (!is.null(fadv_summary)) {
  
  f_summary <- fadv_summary %>%
    transmute(
      scenario = scenario,
      year = as.numeric(as.character(year)),
      f_q05 = fadv_q05,
      f_q50 = fadv_q50,
      f_q95 = fadv_q95,
      period = ifelse(year <= last.obs.yr, "Historical", "Projection")
    )
  
  write_csv(
    f_summary,
    "outputs/Fcap_tables/f_from_Fadv_summary.csv"
  )
  
} else {
  
  warning("No existe advice$Fadv. No se crea f_summary.")
  f_summary <- NULL
}


# ============================================================
# 4c. PANEL BIO DESDE results_list:
# rec season 3, SSB season 2, harvest age 3, TAC
# ============================================================

extract_panel_outputs <- function(res_i, scenario) {
  
  biol_i <- res_i$biols$ANE
  
  keep_cols <- function(x, variable_name) {
    df <- as.data.frame(x)
    names(df)[names(df) == "data"] <- "value"
    
    df %>%
      mutate(variable = variable_name) %>%
      select(year, iter, value, variable)
  }
  
  # Recruitment: age 0, season 3
  rec <- biol_i@n["0", , , "3", , ]
  rec_df <- keep_cols(rec, "rec")
  
  # SSB: season 2
  ssb_all <- quantSums(
    biol_i@n *
      biol_i@wt *
      predict(biol_i@mat)
  )
  ssb_s2 <- ssb_all[, , , "2", , ]
  ssb_df <- keep_cols(ssb_s2, "ssb")
  
  # F / harvest: age 3
  harv <- harvest(biol_i)
  fbar <- harv["3", , , , , ]
  f_df <- keep_cols(fbar, "f")
  
  # TAC / catch advice
  tac_df <- keep_cols(res_i$advice$TAC, "catch")
  
  bind_rows(
    rec_df,
    ssb_df,
    f_df,
    tac_df
  ) %>%
    mutate(
      scenario = scenario,
      year = as.numeric(as.character(year))
    )
}

panel_raw <- bind_rows(lapply(names(results_list), function(sc) {
  extract_panel_outputs(results_list[[sc]], sc)
}))

panel_all <- panel_raw %>%
  group_by(scenario, variable, year) %>%
  summarise(
    q05 = quantile(value, 0.05, na.rm = TRUE),
    q50 = quantile(value, 0.50, na.rm = TRUE),
    q95 = quantile(value, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    period = ifelse(year <= last.obs.yr, "Historical", "Projection"),
    variable = factor(
      variable,
      levels = c("rec", "ssb", "f", "catch"),
      labels = c("Recruitment season 3", "SSB season 2", "Harvest age 3", "TAC")
    )
  )

write_csv(
  panel_all,
  "outputs/Fcap_tables/panel_rec_ssbS2_harvest_TAC_summary.csv"
)
# ============================================================
# 5. GRÁFICAS
# ============================================================

# ---- SSB season 2 ----

p_ssb <- ggplot(
  ssb_summary,
  aes(
    x = as.numeric(as.character(year)),
    y = ssb_q50,
    colour = scenario,
    fill = scenario
  )
) +
  geom_ribbon(
    aes(ymin = ssb_q05, ymax = ssb_q95),
    alpha = 0.20,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  geom_hline(aes(yintercept = Blim), linetype = 2, colour = "red") +
  geom_hline(aes(yintercept = Bpa),  linetype = 3, colour = "blue") +
  labs(
    title = "SSB in season 2 by Fcap scenario",
    subtitle = "Median and 5–95% interval across iterations",
    x = "Year",
    y = "SSB season 2",
    colour = "Scenario",
    fill = "Scenario"
  ) +
  theme_bw()

ggsave(
  "outputs/Fcap_figures/SSB_season2_by_Fcap.png",
  p_ssb,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- TAC ----

p_tac <- ggplot(
  tac_summary,
  aes(
    x = as.numeric(as.character(year)),
    y = tac_q50,
    colour = scenario,
    fill = scenario
  )
) +
  geom_ribbon(
    aes(ymin = tac_q05, ymax = tac_q95),
    alpha = 0.20,
    colour = NA
  ) +
  geom_line(linewidth = 1) +
  labs(
    title = "TAC advice by Fcap scenario",
    subtitle = "Median and 5–95% interval across iterations",
    x = "Year",
    y = "TAC",
    colour = "Scenario",
    fill = "Scenario"
  ) +
  theme_bw()

ggsave(
  "outputs/Fcap_figures/TAC_by_Fcap.png",
  p_tac,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- Risk Blim ----

p_risk_blim <- ggplot(
  ssb_summary,
  aes(
    x = as.numeric(as.character(year)),
    y = risk_Blim,
    colour = scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "red") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Risk of SSB season 2 falling below Blim",
    subtitle = "Dashed line = 5% precautionary threshold",
    x = "Year",
    y = "P(SSB season 2 < Blim)",
    colour = "Scenario"
  ) +
  theme_bw()

ggsave(
  "outputs/Fcap_figures/Risk_Blim_by_Fcap.png",
  p_risk_blim,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- Risk Bpa ----

p_risk_bpa <- ggplot(
  ssb_summary,
  aes(
    x = as.numeric(as.character(year)),
    y = risk_Bpa,
    colour = scenario
  )
) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Risk of SSB season 2 falling below Bpa",
    x = "Year",
    y = "P(SSB season 2 < Bpa)",
    colour = "Scenario"
  ) +
  theme_bw()

ggsave(
  "outputs/Fcap_figures/Risk_Bpa_by_Fcap.png",
  p_risk_bpa,
  width = 10,
  height = 6,
  dpi = 300
)

# ---- Fadv si existe ----

if (!is.null(fadv_summary)) {
  
  p_fadv <- ggplot(
    fadv_summary,
    aes(
      x = as.numeric(as.character(year)),
      y = fadv_q50,
      colour = scenario,
      fill = scenario
    )
  ) +
    geom_ribbon(
      aes(ymin = fadv_q05, ymax = fadv_q95),
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(linewidth = 1) +
    labs(
      title = "F advice by Fcap scenario",
      subtitle = "Median and 5–95% interval across iterations",
      x = "Year",
      y = "F advice",
      colour = "Scenario",
      fill = "Scenario"
    ) +
    theme_bw()
  
  ggsave(
    "outputs/Fcap_figures/Fadv_by_Fcap.png",
    p_fadv,
    width = 10,
    height = 6,
    dpi = 300
  )
}

# ---- Relación Fcap vs riesgo máximo ----

risk_plot_df <- risk_summary %>%
  mutate(
    Fcap = as.numeric(gsub("Fcap_", "", scenario))
  )

p_fcap_risk <- ggplot(
  risk_plot_df,
  aes(x = Fcap, y = risk_Blim_max)
) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "red") +
  labs(
    title = "Maximum risk below Blim by Fcap",
    subtitle = "SSB evaluated in season 2",
    x = "Fcap",
    y = "Maximum P(SSB season 2 < Blim)"
  ) +
  theme_bw()

ggsave(
  "outputs/Fcap_figures/Fcap_vs_risk_Blim.png",
  p_fcap_risk,
  width = 8,
  height = 5,
  dpi = 300
)

# ---- Panel rec, SSB season 2, F, catch ----

p_panel <- ggplot() +
  geom_ribbon(
    data = panel_all %>% filter(period == "Projection"),
    aes(
      x = year,
      ymin = q05,
      ymax = q95,
      fill = scenario
    ),
    alpha = 0.20,
    colour = NA
  ) +
  geom_line(
    data = panel_all %>% filter(period == "Historical"),
    aes(
      x = year,
      y = q50,
      group = interaction(scenario, variable)
    ),
    colour = "grey40",
    alpha = 0.35,
    linewidth = 0.7
  ) +
  geom_line(
    data = panel_all %>% filter(period == "Projection"),
    aes(
      x = year,
      y = q50,
      colour = scenario
    ),
    linewidth = 1
  ) +
  geom_vline(
    xintercept = last.obs.yr,
    linetype = 2
  ) +
  facet_wrap(
    ~ variable,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Biological indicators by Fcap scenario",
    subtitle = "SSB evaluated in spawning season 2",
    x = "Year",
    y = NULL,
    colour = "Scenario",
    fill = "Scenario"
  ) +
  theme_bw()

ggsave(
  "outputs/Fcap_figures/Bio_panel_rec_ssbS2_f_catch_by_Fcap.png",
  p_panel,
  width = 11,
  height = 8,
  dpi = 300
)

# ============================================================
# 6. MOSTRAR RESULTADOS EN CONSOLA
# ============================================================

print(summary_Fcap)
print(risk_summary)

message("Done. Tables saved in outputs/Fcap_tables")
message("Figures saved in outputs/Fcap_figures")

