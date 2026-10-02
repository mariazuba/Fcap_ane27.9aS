# ============================================================
# 13_VALIDATE_REFERENCE_OM.R
# Diagnostics for open-loop validation of the reference OM
# Cleaned diagnostic version: 2026-09-25
#
# Purpose:
#   Evaluate the structural and dynamic behaviour of the
#   reference operating model before implementation of the
#   management procedure.
#
# Validation scenarios:
#   C0    - No fishing
#   C7000 - Constant annual catch advice of 7,000 t
#
# Diagnostics include:
#   - Historical-to-projection continuity
#   - Weight-at-age and natural mortality
#   - Catchability and seasonal fishery structure
#   - Recruitment and Beverton-Holt dynamics
#   - SSB, fishing mortality and realised catch
#   - Unfished biological behaviour
#   - Stock-fishery response under constant catch
#   - Extreme-trajectory diagnostics
#
# INPUTS:
#   data/Rdata/
#     - biols_conditioned_1000iter.rds
#     - fleets_conditioned_1000iter.rds
#     - iteration_map_1000.csv
#
#   boot/data/
#     - stk_ane9aS.rds
#
#   data/mse/OM_validation_C0/
#     - OM_validation_C0_block_01.rds ... block_10.rds
#     - OM_validation_C0_block_01_map.csv ... block_10_map.csv
#
#   data/mse/OM_validation_C7000/
#     - OM_validation_C7000_block_01.rds ... block_10.rds
#     - OM_validation_C7000_block_01_map.csv ... block_10_map.csv
#
# OUTPUTS:
#   outputs/OM/
#     - Diagnostic figures (.png and .pdf)
#
# Period:
#   Historical: 1989-2024
#   Projection: 2025-2054
#
# Notes:
#   - C0 evaluates biological OM behaviour without fishing.
#   - C7000 evaluates the coupled stock-fishery dynamics.
#   - Projection results comprise 100 conditioned OMs with
#     10 stochastic replicates each (1000 trajectories).
# ============================================================

rm(list = ls())

library(FLCore)
library(FLBEIA)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(here)

rds_dir <- here("data", "Rdata")
fig_dir <- here("outputs", "OM")
tab_dir <- here("outputs", "OM")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

hist_years <- 1989:2024
proj_years <- 2025:2054
all_years <- 1989:2054

save_report_plot <- function(plot, filename, width = 9, height = 6) {
  ggsave(file.path(fig_dir, paste0(filename, ".png")),
         plot = plot, width = width, height = height, dpi = 300)
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")),
         plot = plot, width = width, height = height)
}

flq_to_df <- function(x, value_name = "value") {
  
  d <- as.data.frame(x)
  
  # Standardise first FLQuant dimension name
  if("age" %in% names(d)) {
    names(d)[names(d) == "age"] <- "quant"
  }
  
  # Rename FLQuant value column
  if(!"data" %in% names(d)) {
    stop(
      paste0(
        "Column 'data' not found. Columns available: ",
        paste(names(d), collapse = ", ")
      )
    )
  }
  
  names(d)[names(d) == "data"] <- value_name
  
  d %>%
    mutate(
      quant = as.character(quant),
      year = as.integer(as.character(year)),
      season = as.character(season),
      iter = as.integer(as.character(iter))
    )
}


# ============================================================
# LOAD OBJECTS
# ============================================================

biol_file <- file.path(rds_dir, "biols_conditioned_1000iter.rds")
fleet_file <- file.path(rds_dir, "fleets_conditioned_1000iter.rds")
map_file <- file.path(rds_dir, "iteration_map_1000.csv")

biols <- readRDS(biol_file)
fleets <- readRDS(fleet_file)
iteration_map <- read.csv(map_file)

om_map <- iteration_map %>%
  filter(replicate == 1) %>%
  arrange(om)

om_iters <- om_map$iter1000


biol_name <- names(biols)[1]
biol <- biols[[biol_name]]

message("Using FLBiol: ", biol_name)

# ============================================================
# STRUCTURAL CONDITIONING DIAGNOSTICS
# ============================================================
# Stock weight-at-age is plotted once below together with fleet weight-at-age.

# ============================================================
# NATURAL MORTALITY
# ============================================================

m_flq <- slot(biol, "m")

m_df <- flq_to_df(m_flq, "M") %>%
  filter(iter %in% om_iters, year %in% c(hist_years, proj_years)) %>%
  mutate(
    om = match(iter, om_iters),
    age = as.numeric(sub("\\+$", "", quant)),
    season = factor(season, levels = c("1","2","3","4"),
                    labels = c("Q1","Q2","Q3","Q4")),
    period = if_else(year <= 2024, "Historical", "Projection"))

m_median <- m_df %>%
  group_by(year, age, season, period) %>%
  summarise(median = median(M, na.rm = TRUE), .groups = "drop")

fig13b_M <- ggplot(m_df,
  aes(x = year, y = M, group = om, colour = factor(om))) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", linewidth = 0.7) +
  geom_line(alpha = 0.20, linewidth = 0.7) +
  geom_line(data = m_median, aes(x = year, y = median, group = 1),
    inherit.aes = FALSE, colour = "black", linewidth = 0.5) +
  facet_grid(age ~ season, scales = "free_y") +
  guides(colour = "none") +
  labs(
    x = "Year", y = "Natural mortality (M)",
    title = "",
    subtitle = "") +
  theme_bw()

save_report_plot(fig13b_M, "Fig13b_M_historical_projection", 10, 8)

# ============================================================
# CATCHABILITY
# ============================================================

q_flq <- catch.q(fleets$SEINE@metiers$ALL@catches$ANE)

q_df <- flq_to_df(q_flq, "q") %>%
  filter(iter %in% om_iters, year %in% c(hist_years, proj_years)) %>%
  mutate(
    om = match(iter, om_iters),
    age = as.numeric(sub("\\+$", "", quant)),
    season = factor(season, levels = c("1","2","3","4"),
                    labels = c("Q1","Q2","Q3","Q4")),
    period = if_else(year <= 2024, "Historical", "Projection"))

q_median <- q_df %>%
  group_by(year, age, season, period) %>%
  summarise(median = median(q, na.rm = TRUE), .groups = "drop")

fig13c_q <- ggplot(q_df,
  aes(x = year, y = q, group = om, colour = factor(om))) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", linewidth = 0.7) +
  geom_line(alpha = 0.20, linewidth = 0.7) +
  geom_line(data = q_median,aes(x = year, y = median, group = 1),
    inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
  facet_grid(age ~ season, scales = "free_y") +
  guides(colour = "none") +
  labs(
    x = "Year", y = "Catchability (q)") +
  theme_bw()

save_report_plot(fig13c_q, "Fig13c_q_historical_projection", 10, 8)

# ============================================================
# WEIGHT-AT-AGE: STOCK
# ============================================================

wt_stock_df <- flq_to_df(biol@wt,"weight") |>
  filter(iter %in% om_iters,year %in% all_years) |>
  mutate(
    om=match(iter,om_iters),
    age=as.numeric(sub("\\+$","",quant)),
    season=factor(season,levels=c("1","2","3","4"),
                  labels=c("Q1","Q2","Q3","Q4")),
    period=if_else(year<=2024,"Historical","Projection")
  )

wt_stock_median <- wt_stock_df |>
  group_by(year,age,season,period) |>
  summarise(median=median(weight,na.rm=TRUE),.groups="drop")

fig13a_weight_stock <- ggplot(
  wt_stock_df,
  aes(year,weight,group=om,colour=factor(om))
) +
  geom_vline(xintercept=2024.5,linetype="dashed",linewidth=0.6) +
  geom_line(alpha=0.20,linewidth=0.7) +
  geom_line(
    data=wt_stock_median,
    aes(year,median,group=1),
    inherit.aes=FALSE,
    colour="black",linewidth=0.7
  ) +
  facet_grid(age~season,scales="free_y") +
  guides(colour="none") +
  labs(
    x="Year",y="Mean weight (kg)",
    title="Stock weight-at-age",
    subtitle="One trajectory per conditioned OM (n = 100); black line = annual median; dashed line = start of projection"
  ) +
  theme_bw()

save_report_plot(
  fig13a_weight_stock,
  "Fig13a_stock_weight_historical_projection",
  10,8
)


# ============================================================
# WEIGHT-AT-AGE: FLEET
# ============================================================

fleet_name <- names(fleets)[1]
flt <- fleets[[fleet_name]]
cobj <- flt@metiers[[1]]@catches[[biol_name]]

wt_fleet_df <- flq_to_df(cobj@landings.wt,"weight") |>
  filter(iter %in% om_iters,year %in% all_years) |>
  mutate(
    om=match(iter,om_iters),
    age=as.numeric(sub("\\+$","",quant)),
    season=factor(season,levels=c("1","2","3","4"),
                  labels=c("Q1","Q2","Q3","Q4")),
    period=if_else(year<=2024,"Historical","Projection"))

wt_fleet_median <- wt_fleet_df |>
  group_by(year,age,season,period) |>
  summarise(median=median(weight,na.rm=TRUE),.groups="drop")

fig13b_weight_fleet <- ggplot(
  wt_fleet_df,
  aes(year,weight,group=om,colour=factor(om))) +
  geom_vline(xintercept=2024.5,linetype="dashed",linewidth=0.6) +
  geom_line(alpha=0.20,linewidth=0.7) +
  geom_line(
    data=wt_fleet_median,
    aes(year,median,group=1),
    inherit.aes=FALSE,
    colour="black",linewidth=0.7) +
  facet_grid(age~season,scales="free_y") +
  guides(colour="none") +
  labs(
    x="Year",y="Mean weight (kg)",
    title="Fleet weight-at-age",
    subtitle="Landings weight; one trajectory per conditioned OM (n = 100); black line = annual median"
  ) +
  theme_bw()

save_report_plot(
  fig13b_weight_fleet,
  "Fig13b_fleet_weight_historical_projection",
  10,8
)



# ============================================================
# DYNAMIC OPEN-LOOP OM OUTPUTS
#
# C0    : biological OM validation
# C7000 : biological + fishery OM validation
#
# Historical: 100 conditioned OMs (one replicate per OM)
# Projection: 100 OMs x 10 replicates = 1000 trajectories
# Recruitment = age 0, Q3
# SSB = Q2
# F = FLBEIA bioSum indicator "f"
# Realised catch = FLCatch@landings
# ============================================================

scenario_dirs <- c(
  C0=here("data","mse","OM_validation_C0"),
  C7000=here("data","mse","OM_validation_C7000"))

block_ids <- 1:10
Blim <- 4721
Bpa <- 6561

# ============================================================
# HELPERS
# ============================================================

add_global_map <- function(dat,map_b) {
  out <- dat %>% mutate(iter_local=as.integer(as.character(iter))) %>% left_join(map_b %>% dplyr::select(iter_local,iter1000,om,replicate),by="iter_local")
  stopifnot(!anyNA(out$iter1000),!anyNA(out$om),!anyNA(out$replicate))
  out
}

standardise_dynamic <- function(x,map_b) {
  as.data.frame(x) %>% rename(value=data) %>% mutate(year=as.integer(as.character(year))) %>% add_global_map(map_b) %>% transmute(year,iter1000,om,replicate,value)
}

standardise_biosum <- function(x,map_b) {
  x %>% ungroup() %>% mutate(year=as.integer(as.character(year)),iter=as.integer(as.character(iter))) %>% add_global_map(map_b) %>% ungroup() %>% transmute(year,iter1000,om,replicate,value)
}
# ============================================================
# READ AND EXTRACT 10 BLOCKS
# ============================================================

read_validation_scenario <- function(scenario_id) {
  
  stopifnot(scenario_id %in% names(scenario_dirs))
  
  om_out_dir <- scenario_dirs[[scenario_id]]
  om_files <- file.path(om_out_dir,sprintf("OM_validation_%s_block_%02d.rds",scenario_id,block_ids))
  map_files <- file.path(om_out_dir,sprintf("OM_validation_%s_block_%02d_map.csv",scenario_id,block_ids))
  
  stopifnot(length(om_files)==10,length(map_files)==10,all(file.exists(om_files)),all(file.exists(map_files)))
  
  dynamic_list <- map(block_ids,function(b) {
    
    message("Reading ",scenario_id," block ",b," / ",length(block_ids))
    
    OM_b <- readRDS(om_files[b])
    map_b <- read.csv(map_files[b])
    
    biol <- OM_b$biols$ANE
    cobj <- OM_b$fleets$SEINE@metiers$ALL@catches$ANE
    
    N_b <- as.data.frame(biol@n) %>% rename(value=data) %>% 
      mutate(year=as.integer(as.character(year)),age=as.integer(as.character(age)),season=as.integer(as.character(season))) %>% 
      add_global_map(map_b) %>% transmute(age,year,season,iter1000,om,replicate,value)
    
    M_b <- as.data.frame(biol@m) %>% rename(value=data) %>%
      mutate(year=as.integer(as.character(year)),age=as.integer(as.character(age)),season=as.integer(as.character(season))) %>% 
      add_global_map(map_b) %>% transmute(age,year,season,iter1000,om,replicate,value)
    
    W_b <- as.data.frame(biol@wt) %>% rename(value=data) %>% 
      mutate(year=as.integer(as.character(year)),age=as.integer(as.character(age)),season=as.integer(as.character(season))) %>% 
      add_global_map(map_b) %>% transmute(age,year,season,iter1000,om,replicate,value)
    
    fleet <- OM_b$fleets$SEINE
    
    Q_b <- as.data.frame(OM_b$fleets$SEINE@metiers$ALL@catches$ANE@catch.q) %>%
      rename(value=data) %>%
      mutate(year=as.integer(as.character(year)),age=as.integer(as.character(age)),season=as.integer(as.character(season))) %>%
      add_global_map(map_b) %>%
      transmute(age,year,season,iter1000,om,replicate,value)
    
    # Recruitment: age 0, Q3
    R_b <- standardise_dynamic(biol@n["0",,,"3",,],map_b)
    
    # SSB: Q2
    Mat <- predict(biol@mat)
    SSB_q2 <- quantSums(biol@n[,,,"2",,]*biol@wt[,,,"2",,]*Mat[,,,"2",,])
    B_b <- standardise_dynamic(SSB_q2,map_b)
    
    # Fishing mortality: FLBEIA standard biological summary
    bio_b <- bioSum(OM_b,scenario=scenario_id,long=TRUE,byyear=TRUE,ssb_season=2)
    F_b <- bio_b %>% filter(indicator=="f") %>% standardise_biosum(map_b)
    
    # Realised annual catch from the fleet object
    C_ann <- seasonSums(quantSums(cobj@landings))
    C_b <- standardise_dynamic(C_ann,map_b)
    
    list(N=N_b,M=M_b,W=W_b,Q=Q_b,R=R_b,B=B_b,F=F_b,C=C_b)
  })
  
  out <- list(
    N=map_dfr(dynamic_list,"N") %>% mutate(scenario=scenario_id),
    M=map_dfr(dynamic_list,"M") %>% mutate(scenario=scenario_id),
    W=map_dfr(dynamic_list,"W") %>% mutate(scenario=scenario_id),
    Q=map_dfr(dynamic_list,"Q") %>% mutate(scenario=scenario_id),
    R=map_dfr(dynamic_list,"R") %>% mutate(scenario=scenario_id),
    B=map_dfr(dynamic_list,"B") %>% mutate(scenario=scenario_id),
    F=map_dfr(dynamic_list,"F") %>% mutate(scenario=scenario_id),
    C=map_dfr(dynamic_list,"C") %>% mutate(scenario=scenario_id)
  )
  
  stopifnot(
    n_distinct(out$R$iter1000)==1000,n_distinct(out$B$iter1000)==1000,n_distinct(out$F$iter1000)==1000,n_distinct(out$C$iter1000)==1000,
    n_distinct(out$R$om)==100,n_distinct(out$B$om)==100,n_distinct(out$F$om)==100,n_distinct(out$C$om)==100,
    all(is.finite(out$R$value)),all(is.finite(out$B$value)),all(is.finite(out$F$value)),all(is.finite(out$C$value))
  )
  
  message(scenario_id,": ",n_distinct(out$B$om)," OMs, ",n_distinct(out$B$iter1000)," trajectories")
  out
}

# ============================================================
# READ BOTH VALIDATION SCENARIOS
# ============================================================

val_C0 <- read_validation_scenario("C0")
val_C7000 <- read_validation_scenario("C7000")

# ============================================================
# 2.4.1 STRUCTURAL AND NUMERICAL VALIDATION
# ============================================================

validate_outputs <- function(x,scenario_id) {
  bind_rows(
    x$R %>% mutate(variable="Recruitment"),
    x$B %>% mutate(variable="SSB"),
    x$F %>% mutate(variable="F"),
    x$C %>% mutate(variable="Catch")
  ) %>%
    group_by(variable) %>%
    summarise(
      scenario=scenario_id,
      n_OM=n_distinct(om),
      n_trajectories=n_distinct(iter1000),
      first_year=min(year),
      last_year=max(year),
      n_NA=sum(is.na(value)),
      n_nonfinite=sum(!is.finite(value)),
      n_negative=sum(value<0,na.rm=TRUE),
      min=min(value,na.rm=TRUE),
      max=max(value,na.rm=TRUE),
      .groups="drop")
}

validation_structure <- bind_rows(
  validate_outputs(val_C0,"C0"),
  validate_outputs(val_C7000,"C7000")
)

validation_structure

# ============================================================
# CHECK PROJECTION SCENARIO IMPLEMENTATION
# ============================================================

validation_projection <- bind_rows(
  val_C0$F %>% filter(year>=2025) %>% summarise(scenario="C0",variable="F",min=min(value),max=max(value)),
  val_C0$C %>% filter(year>=2025) %>% summarise(scenario="C0",variable="Catch",min=min(value),max=max(value)),
  val_C7000$F %>% filter(year>=2025) %>% summarise(scenario="C7000",variable="F",min=min(value),max=max(value)),
  val_C7000$C %>% filter(year>=2025) %>% summarise(scenario="C7000",variable="Catch",min=min(value),max=max(value))
)

validation_projection


# ============================================================
# TABLE 2.4.1 STRUCTURAL AND NUMERICAL VALIDATION
# ============================================================

tab_structure <- validation_structure %>%
  mutate(
    range=paste0(signif(min,4)," – ",signif(max,4)),
    complete=if_else(n_NA==0 & n_nonfinite==0 & n_negative==0,"Yes","No")) %>%
  dplyr::select(scenario,variable,n_OM,n_trajectories,first_year,last_year,n_NA,n_nonfinite,n_negative,range,complete)

tab_structure

tab_scenario_check <- validation_projection %>%
  mutate(
    expected=case_when(
      scenario=="C0" & variable=="F" ~ "0",
      scenario=="C0" & variable=="Catch" ~ "0 t",
      scenario=="C7000" & variable=="F" ~ "> 0",
      scenario=="C7000" & variable=="Catch" ~ "≤ 7,000 t"),
    range=paste0(signif(min,4)," – ",signif(max,4))) %>%
  dplyr::select(scenario,variable,expected,range)

tab_scenario_check

# ============================================================
# JOINT OPEN-LOOP VALIDATION FIGURE: C0 VS C7000
# ============================================================
# Historical uncertainty: 100 conditioned OMs (one replicate per OM).
# Projection uncertainty: 100 OMs x 10 stochastic replicates.
# Lines = annual medians; ribbons = 5th-95th percentile intervals.
# ============================================================

summarise_validation_series <- function(dat) {
  historical <- dat %>%
    filter(year <= 2024, replicate == 1) %>%
    group_by(scenario, year) %>%
    summarise(median = median(value, na.rm = TRUE), q05 = quantile(value, 0.05, na.rm = TRUE), q95 = quantile(value, 0.95, na.rm = TRUE), .groups = "drop") %>%
    mutate(period = "Historical")

  projection <- dat %>%
    filter(year >= 2025) %>%
    group_by(scenario, year) %>%
    summarise(median = median(value, na.rm = TRUE), q05 = quantile(value, 0.05, na.rm = TRUE), q95 = quantile(value, 0.95, na.rm = TRUE), .groups = "drop") %>%
    mutate(period = "Projection")

  bind_rows(historical, projection)
}

validation_C <- bind_rows(val_C0$C, val_C7000$C) %>% summarise_validation_series()
validation_F <- bind_rows(val_C0$F, val_C7000$F) %>% summarise_validation_series()
validation_R <- bind_rows(val_C0$R, val_C7000$R) %>% summarise_validation_series()
validation_SSB <- bind_rows(val_C0$B, val_C7000$B) %>% summarise_validation_series()

stopifnot(all(c("scenario", "year", "median", "q05", "q95", "period") %in% names(validation_C)))
stopifnot(all(c("scenario", "year", "median", "q05", "q95", "period") %in% names(validation_F)))
stopifnot(all(c("scenario", "year", "median", "q05", "q95", "period") %in% names(validation_R)))
stopifnot(all(c("scenario", "year", "median", "q05", "q95", "period") %in% names(validation_SSB)))

p_validation_C <- ggplot(validation_C, aes(x = year, y = median, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 7000, linetype = "dotted", linewidth = 0.6) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", linewidth = 0.6) +
  labs(x = NULL, y = "Catch (t)", colour = "Scenario", fill = "Scenario", tag = "A") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_validation_F <- ggplot(validation_F, aes(x = year, y = median, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", linewidth = 0.6) +
  labs(x = NULL, y = "Fishing mortality (F)", colour = "Scenario", fill = "Scenario", tag = "B") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_validation_R <- ggplot(validation_R, aes(x = year, y = median / 1e6, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q05 / 1e6, ymax = q95 / 1e6), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", linewidth = 0.6) +
  labs(x = NULL, y = "Recruitment (millions)", colour = "Scenario", fill = "Scenario", tag = "C") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p_validation_SSB <- ggplot(validation_SSB, aes(x = year, y = median, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = Blim, linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = Bpa, linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = 2024.5, linetype = "dashed", linewidth = 0.6) +
  labs(x = "Year", y = "SSB (t)", colour = "Scenario", fill = "Scenario", tag = "D") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

fig_validation_reference_OM <- (p_validation_C / p_validation_F / p_validation_R / p_validation_SSB) +
  plot_layout(guides = "collect", heights = c(1, 1, 1, 1)) &
  theme(legend.position = "bottom", plot.tag = element_text(size = 11, face = "bold"), plot.margin = margin(5, 8, 5, 8))

save_report_plot(fig_validation_reference_OM, "Fig_reference_OM_validation_C0_C7000", width = 10, height = 12)

# ============================================================
# 2.4.2 HISTORICAL-TO-PROJECTION CONTINUITY
# ============================================================

# ------------------------------------------------------------
# 2.4.2.1 SSB AND RECRUITMENT AROUND THE PROJECTION BOUNDARY
# ------------------------------------------------------------
continuity_bio <- bind_rows(
  val_C7000$B %>% mutate(variable="SSB"),
  val_C7000$R %>% mutate(variable="Recruitment")
) %>%
  filter(year %in% c(2024,2025))


continuity_bio_clean <- continuity_bio %>%
  filter((year==2025 & variable=="Recruitment") | replicate==1)


continuity_bio_summary <- continuity_bio_clean %>%
  group_by(variable,year) %>%
  summarise(
    n=dplyr::n(),
    median=median(value),
    q05=quantile(value,0.05),
    q25=quantile(value,0.25),
    q75=quantile(value,0.75),
    q95=quantile(value,0.95),
    min=min(value),
    max=max(value),
    .groups="drop")

continuity_bio_summary


# ------------------------------------------------------------
# 2.4.2.2 HISTORICAL 2024Q4 -> 2025Q1 COHORT TRANSITION
# ------------------------------------------------------------
# ------------------------------------------------------------
# Cohort propagation: 2024 Q4 -> 2025 Q1
# ------------------------------------------------------------

N_2024_Q4 <- val_C0$N %>%
  filter(year==2024,season==4,age %in% 0:1,replicate==1) %>%
  transmute(
    om,
    cohort=if_else(age==0,"age0_2024Q4_to_age1_2025Q1","age1_2024Q4_to_age2_2025Q1"),
    N_2024Q4=value)

N_2025_Q1 <- val_C0$N %>%
  filter(year==2025,season==1,age %in% 1:2,replicate==1) %>%
  transmute(
    om,
    cohort=if_else(age==1,"age0_2024Q4_to_age1_2025Q1","age1_2024Q4_to_age2_2025Q1"),
    N_2025Q1=value)

N_cohort <- N_2024_Q4 %>%
  left_join(N_2025_Q1,by=c("om","cohort")) %>%
  mutate(survival=N_2025Q1/N_2024Q4)


M_Q4_2024 <- val_C0$M %>%
  filter(year==2024,season==4,age %in% 0:1,replicate==1) %>%
  transmute(
    om,
    cohort=if_else(age==0,"age0_to_age1","age1_to_age2"),
    M=value)

cohort_transition_2024_2025 <- N_cohort %>%
  mutate(cohort=recode(cohort,
                       "age0_2024Q4_to_age1_2025Q1"="age0_to_age1",
                       "age1_2024Q4_to_age2_2025Q1"="age1_to_age2")) %>%
  left_join(M_Q4_2024,by=c("om","cohort")) %>%
  mutate(
    Z_implied=-log(survival),
    F_implied=Z_implied-M)

transition_2024_2025_summary <- cohort_transition_2024_2025 %>%
  group_by(cohort) %>%
  summarise(
    n=dplyr::n(),
    median_survival=median(survival),
    median_Z=median(Z_implied),
    median_M=median(M),
    median_F=median(F_implied),
    min_F=min(F_implied),
    max_F=max(F_implied),
    .groups="drop")

transition_2024_2025_summary


# ------------------------------------------------------------
# 2.4.2.3 COHORT PROPAGATION UNDER C0: 2025Q4 -> 2026Q1
# ------------------------------------------------------------

N_transition_C0 <- val_C0$N %>%
  filter((year==2025 & season==4 & age %in% 0:1) | (year==2026 & season==1 & age %in% 1:2)) %>%
  transmute(
    om,replicate,
    cohort=case_when(
      year==2025 & age==0 ~ "age0_to_age1",
      year==2025 & age==1 ~ "age1_to_age2",
      year==2026 & age==1 ~ "age0_to_age1",
      year==2026 & age==2 ~ "age1_to_age2"),
    period=if_else(year==2025,"N_start","N_end"),
    value) %>%
  pivot_wider(names_from=period,values_from=value) %>%
  mutate(observed_survival=N_end/N_start)

M_Q4_2025 <- val_C0$M %>%
  filter(year==2025,season==4,age %in% 0:1) %>%
  transmute(
    om,replicate,
    cohort=if_else(age==0,"age0_to_age1","age1_to_age2"),
    M=value)

survival_C0 <- N_transition_C0 %>%
  left_join(M_Q4_2025,by=c("om","replicate","cohort")) %>%
  mutate(expected_survival=exp(-M),difference=observed_survival-expected_survival)

survival_C0_summary <- survival_C0 %>%
  group_by(cohort) %>%
  summarise(
    n=dplyr::n(),
    observed=median(observed_survival),
    expected=median(expected_survival),
    min_diff=min(difference),
    median_diff=median(difference),
    max_diff=max(difference),
    .groups="drop")

survival_C0_summary
# ------------------------------------------------------------
# 2.4.2.4 WEIGHT-AT-AGE TRANSITION
# ------------------------------------------------------------

weight_transition <- val_C0$W %>%
  filter(year %in% 2022:2025,season==2,replicate==1) %>%
  group_by(age,year) %>%
  summarise(
    n=dplyr::n(),
    median=median(value),
    q05=quantile(value,0.05),
    q95=quantile(value,0.95),
    min=min(value),
    max=max(value),
    .groups="drop")

weight_transition

weight_expected_2025 <- val_C0$W %>%
  filter(year %in% 2022:2024,season==2,replicate==1) %>%
  group_by(om,age) %>%
  summarise(expected_2025=mean(value),.groups="drop")

weight_observed_2025 <- val_C0$W %>%
  filter(year==2025,season==2,replicate==1) %>%
  dplyr::select(om,age,observed_2025=value)

weight_check <- weight_observed_2025 %>%
  left_join(weight_expected_2025,by=c("om","age")) %>%
  mutate(difference=observed_2025-expected_2025)

weight_check_summary <- weight_check %>%
  group_by(age) %>%
  summarise(
    n=dplyr::n(),
    observed=median(observed_2025),
    expected=median(expected_2025),
    min_diff=min(difference),
    median_diff=median(difference),
    max_diff=max(difference),
    .groups="drop")

weight_check_summary

# ------------------------------------------------------------
# 2.4.2.5 NATURAL MORTALITY TRANSITION
# ------------------------------------------------------------

M_transition <- val_C0$M %>%
  filter(year %in% 2022:2025,replicate==1) %>%
  group_by(age,season,year) %>%
  summarise(
    n=dplyr::n(),
    median=median(value),
    q05=quantile(value,0.05),
    q95=quantile(value,0.95),
    min=min(value),
    max=max(value),
    .groups="drop")

M_transition

# ------------------------------------------------------------
# 2.4.2.5 NATURAL MORTALITY TRANSITION
# Lorenzen weight-scaling validation
# ------------------------------------------------------------

lorenzen_b <- 0.288
M_base_season <- c("0"=2.97,"1"=1.60,"2"=2.48,"3"=2.48)/4

stk_base <- readRDS(file.path(getwd(),"boot","data","stk_ane9aS.rds"))
wt_base_proj <- yearMeans(stock.wt(stk_base[,ac(2022:2024)]))

W_base <- as.data.frame(wt_base_proj) %>%
  transmute(age=as.integer(as.character(age)),season=as.integer(as.character(season)),W_base=data)

M_expected_2025 <- val_C0$W %>%
  filter(year==2025,replicate==1) %>%
  dplyr::select(om,age,season,W_boot=value) %>%
  left_join(W_base,by=c("age","season")) %>%
  mutate(M_base=M_base_season[as.character(age)],
         expected_M=if_else(is.finite(W_boot) & is.finite(W_base) & W_boot>0 & W_base>0,M_base*(W_base/W_boot)^lorenzen_b,M_base))

M_observed_2025 <- val_C0$M %>%
  filter(year==2025,replicate==1) %>%
  dplyr::select(om,age,season,observed_M=value)

M_check <- M_observed_2025 %>%
  left_join(M_expected_2025 %>% dplyr::select(om,age,season,expected_M),by=c("om","age","season")) %>%
  mutate(difference=observed_M-expected_M)

M_check_summary <- M_check %>%
  group_by(age,season) %>%
  summarise(n=dplyr::n(),observed=median(observed_M),expected=median(expected_M),min_diff=min(difference),median_diff=median(difference),max_diff=max(difference),.groups="drop")

M_check_summary

# ------------------------------------------------------------
# 2.4.2.6 TEMPORAL CONSTANCY OF PROJECTED M
# ------------------------------------------------------------

M_projection_constancy <- val_C0$M %>%
  filter(year %in% 2025:2054,replicate==1) %>%
  group_by(om,age,season) %>%
  summarise(
    min_M=min(value),
    max_M=max(value),
    difference=max_M-min_M,
    .groups="drop")

M_projection_constancy_summary <- M_projection_constancy %>%
  summarise(
    n=dplyr::n(),
    max_abs_diff=max(abs(difference)),
    n_failed=sum(abs(difference)>1e-10))

M_projection_constancy_summary

# ------------------------------------------------------------
# 2.4.2.7 TEMPORAL CONSTANCY OF PROJECTED CATCHABILITY
# ------------------------------------------------------------

Q_projection_constancy <- val_C0$Q %>%
  filter(year %in% 2025:2054,replicate==1) %>%
  group_by(om,age,season) %>%
  summarise(
    min_Q=min(value),
    max_Q=max(value),
    difference=max_Q-min_Q,
    .groups="drop")

Q_projection_constancy_summary <- Q_projection_constancy %>%
  summarise(
    n=dplyr::n(),
    max_abs_diff=max(abs(difference)),
    n_failed=sum(abs(difference)>1e-10))

Q_projection_constancy_summary

# ------------------------------------------------------------
# OM-SPECIFIC VARIABILITY OF PROJECTED CATCHABILITY
# ------------------------------------------------------------

Q_projection_variability <- val_C0$Q %>%
  filter(year==2025,replicate==1) %>%
  group_by(age,season) %>%
  summarise(
    n_om=n_distinct(om),
    n_unique=n_distinct(value),
    median=median(value),
    q05=quantile(value,0.05),
    q95=quantile(value,0.95),
    min=min(value),
    max=max(value),
    .groups="drop")

Q_projection_variability

# ------------------------------------------------------------
# 2.4.2.8 SEASONAL CATCH SHARE
# ------------------------------------------------------------

OM_test <- readRDS(file.path(here(),"data","mse","OM_validation_C0","OM_validation_C0_block_01.rds"))

iteration_map_1000 <- read.csv(file.path(here(),"data","Rdata","iteration_map_1000.csv"))

seasonal_share <- as.data.frame(OM_test$fleets.ctrl$seasonal.share$ANE) %>%
  mutate(
    year=as.integer(as.character(year)),
    season=as.integer(as.character(season)),
    iter1000=as.integer(as.character(iter))
  ) %>%
  left_join(iteration_map_1000 %>% dplyr::select(iter1000,om,replicate),by="iter1000")


seasonal_share_expected <- seasonal_share %>%
  filter(year %in% 2015:2024) %>%
  group_by(iter1000,om,replicate,season) %>%
  summarise(expected_2025=mean(data),.groups="drop")

seasonal_share_observed <- seasonal_share %>%
  filter(year==2025) %>%
  dplyr::select(iter1000,om,replicate,season,observed_2025=data)

seasonal_share_check <- seasonal_share_observed %>%
  left_join(seasonal_share_expected,by=c("iter1000","om","replicate","season")) %>%
  mutate(difference=observed_2025-expected_2025)

seasonal_share_check_summary <- seasonal_share_check %>%
  group_by(season) %>%
  summarise(
    n=dplyr::n(),
    observed=median(observed_2025),
    expected=median(expected_2025),
    min_diff=min(difference),
    median_diff=median(difference),
    max_diff=max(difference),
    .groups="drop")

seasonal_share_check_summary

seasonal_share_check %>%
  summarise(
    max_abs_diff=max(abs(difference)),
    n_failed=sum(abs(difference)>1e-10))

seasonal_share_sum <- seasonal_share %>%
  filter(year %in% 2025:2054) %>%
  group_by(iter1000,year) %>%
  summarise(total_share=sum(data),.groups="drop")

seasonal_share_sum %>%
  summarise(
    min_sum=min(total_share),
    median_sum=median(total_share),
    max_sum=max(total_share),
    max_abs_diff=max(abs(total_share-1)),
    n_failed=sum(abs(total_share-1)>1e-10))

range(as.numeric(OM_test$fleets$SEINE@metiers$ALL@catches$ANE@landings.sel[,as.character(2025:2054),,,,]),na.rm=TRUE)

range(as.numeric(OM_test$fleets$SEINE@metiers$ALL@catches$ANE@discards.sel[,as.character(2025:2054),,,,]),na.rm=TRUE)


# ============================================================
# FIGURE X. HISTORICAL-TO-PROJECTION CONTINUITY
# ============================================================

# ------------------------------------------------------------
# SSB
# ------------------------------------------------------------

fig_SSB_dat <- val_C0$B %>%
  filter(year %in% 2020:2027) %>%
  group_by(year) %>%
  summarise(median=median(value),q05=quantile(value,0.05),q95=quantile(value,0.95),.groups="drop")

p_SSB <- ggplot(fig_SSB_dat,aes(year,median)) +
  geom_ribbon(aes(ymin=q05,ymax=q95),alpha=0.20) +
  geom_line(linewidth=0.8) +
  geom_point(size=1.8) +
  geom_vline(xintercept=2024.5,linetype=2) +
  labs(x=NULL,y="SSB (t)",tag="A") +
  theme_bw() +
  theme(panel.grid.minor=element_blank())

# ------------------------------------------------------------
# RECRUITMENT
# ------------------------------------------------------------

fig_R_dat <- val_C0$R %>%
  filter(year %in% 2020:2027) %>%
  group_by(year) %>%
  summarise(median=median(value),q05=quantile(value,0.05),q95=quantile(value,0.95),.groups="drop")

p_R <- ggplot(fig_R_dat,aes(year,median/1e6)) +
  geom_ribbon(aes(ymin=q05/1e6,ymax=q95/1e6),alpha=0.20) +
  geom_line(linewidth=0.8) +
  geom_point(size=1.8) +
  geom_vline(xintercept=2024.5,linetype=2) +
  labs(x=NULL,y="Recruitment (millions)",tag="B") +
  theme_bw() +
  theme(panel.grid.minor=element_blank())


# ------------------------------------------------------------
# WEIGHT-AT-AGE BY AGE AND SEASON
# ------------------------------------------------------------
season_cols <- c("Q1"="#1b9e77","Q2"="#d95f02","Q3"="#7570b3","Q4"="#e7298a")

fig_W_dat <- val_C0$W %>%
  filter(year %in% 2020:2027,replicate==1) %>%
  group_by(age,season,year) %>%
  summarise(median=median(value),q05=quantile(value,0.05),q95=quantile(value,0.95),.groups="drop") %>%
  mutate(age=factor(age),season=factor(season,levels=1:4,labels=paste0("Q",1:4)))

p_W <- ggplot(fig_W_dat,aes(year,median,group=season,color=season)) +
  geom_line(linewidth=0.75) +
  geom_point(size=1.5) +
  geom_vline(xintercept=2024.5,linetype=2) +
  facet_wrap(~age,nrow=1,scales="free_y",labeller=labeller(age=function(x) paste("Age",x))) +
  scale_color_manual(values=season_cols) +
  labs(x="Year",y="Weight-at-age (kg)",color="Season",tag="C") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),strip.background=element_blank())
# ------------------------------------------------------------
# CATCHABILITY BY AGE AND SEASON
# ------------------------------------------------------------

fig_Q_dat <- val_C0$Q %>%
  filter(year %in% 2020:2027,replicate==1) %>%
  group_by(age,season,year) %>%
  summarise(median=median(value),q05=quantile(value,0.05),q95=quantile(value,0.95),.groups="drop") %>%
  mutate(age=factor(age),season=factor(season,levels=1:4,labels=paste0("Q",1:4)))

p_Q <- ggplot(fig_Q_dat,aes(year,median,group=season,color=season)) +
  geom_line(linewidth=0.75) +
  geom_point(size=1.5) +
  geom_vline(xintercept=2024.5,linetype=2) +
  facet_wrap(~age,nrow=1,scales="free_y",labeller=labeller(age=function(x) paste("Age",x))) +
  scale_color_manual(values=season_cols) +
  labs(x="Year",y="Catchability (q)",color="Season",tag="D") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),strip.background=element_blank())
# ------------------------------------------------------------
# COMBINE
# ------------------------------------------------------------

fig_continuity <- (p_SSB ) / (p_R) +
  plot_layout(heights=c(1,1))

fig_continuity

fig_continuity2 <- (p_W / p_Q) +
  plot_layout(heights=c(1,1),guides="collect") &
  theme(legend.position="bottom",legend.margin=margin(t=2),plot.margin=margin(5,5,5,5))
fig_continuity2

# ============================================================
# TABLE X. HISTORICAL-TO-PROJECTION CONTINUITY CHECKS
# ============================================================

R_2024_median <- continuity_bio_summary %>%
  filter(variable == "Recruitment", year == 2024) %>%
  pull(median)

R_2025_median <- continuity_bio_summary %>%
  filter(variable == "Recruitment", year == 2025) %>%
  pull(median)

SSB_2024_median <- continuity_bio_summary %>%
  filter(variable == "SSB", year == 2024) %>%
  pull(median)

SSB_2025_median <- continuity_bio_summary %>%
  filter(variable == "SSB", year == 2025) %>%
  pull(median)

survival_max_abs_diff <- max(abs(survival_C0$difference))

weight_max_abs_diff <- max(abs(weight_check$difference))

M_scaling_max_abs_diff <- max(abs(M_check$difference))

M_constancy_max_abs_diff <- max(abs(M_projection_constancy$difference))

Q_constancy_max_abs_diff <- max(abs(Q_projection_constancy$difference))

landings_sel_range <- range(
  as.numeric(
    OM_test$fleets$SEINE@metiers$ALL@catches$ANE@landings.sel[
      , as.character(proj_years), , , ,
    ]
  ),
  na.rm = TRUE
)

discards_sel_range <- range(
  as.numeric(
    OM_test$fleets$SEINE@metiers$ALL@catches$ANE@discards.sel[
      , as.character(proj_years), , , ,
    ]
  ),
  na.rm = TRUE
)

seasonal_share_max_abs_diff <- max(abs(seasonal_share_check$difference))

seasonal_sum_max_abs_diff <- max(abs(seasonal_share_sum$total_share - 1))

table_continuity <- tibble(
  Component = c(
    "Recruitment",
    "SSB",
    "Numbers-at-age",
    "Weight-at-age",
    "Natural mortality",
    "Catchability",
    "Landings selectivity",
    "Discards selectivity",
    "Seasonal catch allocation"
  ),
  `Projection assumption` = c(
    "BH stock-recruit relationship with stochastic process error from 2025",
    "Population propagated dynamically; SSB evaluated in Q2",
    "Ageing and survival determined by quarterly mortality",
    "OM-specific mean of 2022-2024",
    "Lorenzen scaling based on projected weight-at-age; constant after 2025",
    "OM-specific estimate based on 2020-2024; constant after 2025",
    "Fixed at 1",
    "Fixed at 0",
    "OM-specific mean seasonal proportions from 2015-2024"
  ),
  `Validation result` = c(
    sprintf(
      "Median recruitment: %.2f million in 2024 and %.2f million in 2025; stochastic process error applied from 2025",
      R_2024_median / 1e6,
      R_2025_median / 1e6
    ),
    sprintf(
      "Median SSB: %s t in 2024 and %s t in 2025; transition generated by projected population dynamics",
      format(round(SSB_2024_median), big.mark = ",", scientific = FALSE),
      format(round(SSB_2025_median), big.mark = ",", scientific = FALSE)
    ),
    sprintf(
      "Under C0, 2025 Q4-2026 Q1 cohort survival matched exp(-M); max. absolute difference = %.3g",
      survival_max_abs_diff
    ),
    sprintf(
      "Projected 2025 values matched OM-specific 2022-2024 means; max. absolute difference = %.3g",
      weight_max_abs_diff
    ),
    sprintf(
      "Implemented M matched Lorenzen-scaled values (max. absolute difference = %.3g) and remained constant during 2025-2054 (max. temporal difference = %.3g)",
      M_scaling_max_abs_diff,
      M_constancy_max_abs_diff
    ),
    sprintf(
      "Catchability remained constant within each OM during 2025-2054; max. temporal difference = %.3g; OM-specific variability retained",
      Q_constancy_max_abs_diff
    ),
    sprintf(
      "Projected range: %.3g-%.3g",
      landings_sel_range[1],
      landings_sel_range[2]
    ),
    sprintf(
      "Projected range: %.3g-%.3g",
      discards_sel_range[1],
      discards_sel_range[2]
    ),
    sprintf(
      "Projected values matched OM-specific 2015-2024 means (max. absolute difference = %.3g) and seasonal proportions summed to 1 (max. absolute deviation = %.3g)",
      seasonal_share_max_abs_diff,
      seasonal_sum_max_abs_diff
    )
  )
)

knitr::kable(
  table_continuity,
  format = "pipe",
  align = c("l", "l", "l"),
  caption = paste(
    "Summary of historical-to-projection continuity checks for the Operating Model.",
    "Validation assessed consistency with predefined OM-specific projection assumptions",
    "rather than requiring numerical equality between 2024 and 2025."
  )
)


# ============================================================
# 2.4.3 UNFISHED BIOLOGICAL DIAGNOSTIC
# ============================================================
# C0 is used as a biological control scenario to evaluate
# projected population dynamics in the absence of fishing.
# ============================================================


# ------------------------------------------------------------
# 2.4.3.1 PROJECTED SSB UNDER C0
# ------------------------------------------------------------

C0_SSB_year <- val_C0$B %>%
  filter(year %in% proj_years) %>%
  group_by(year) %>%
  summarise(
    n = dplyr::n(),
    median = median(value),
    q05 = quantile(value, 0.05),
    q95 = quantile(value, 0.95),
    min = min(value),
    max = max(value),
    p_blim = mean(value < Blim),
    p_bpa = mean(value < Bpa),
    .groups = "drop"
  )

C0_SSB_summary <- C0_SSB_year %>%
  summarise(
    min_annual_median = min(median),
    max_annual_median = max(median),
    min_SSB = min(min),
    max_SSB = max(max),
    max_annual_p_blim = max(p_blim),
    max_annual_p_bpa = max(p_bpa)
  )

C0_SSB_summary


# ------------------------------------------------------------
# 2.4.3.2 PROJECTED RECRUITMENT UNDER C0
# ------------------------------------------------------------

C0_R_year <- val_C0$R %>%
  filter(year %in% proj_years) %>%
  group_by(year) %>%
  summarise(
    n = dplyr::n(),
    median = median(value),
    q05 = quantile(value, 0.05),
    q95 = quantile(value, 0.95),
    min = min(value),
    max = max(value),
    .groups = "drop"
  )

C0_R_summary <- C0_R_year %>%
  summarise(
    min_annual_median = min(median),
    max_annual_median = max(median),
    min_R = min(min),
    max_R = max(max)
  )

C0_R_summary


# ------------------------------------------------------------
# 2.4.3.3 BEVERTON-HOLT RECRUITMENT PROCESS CHECK
# ------------------------------------------------------------
# Recruitment generated during the projection is compared with
# the deterministic OM-specific Beverton-Holt expectation.
#
# log(R / R_BH) should reproduce the lognormal process-error
# distribution parameterised by sigmaR for each conditioned OM.
# ------------------------------------------------------------

BH_check <- val_C0$B %>%
  filter(year %in% proj_years) %>%
  dplyr::select(year, iter1000, om, replicate, SSB = value) %>%
  left_join(
    val_C0$R %>%
      filter(year %in% proj_years) %>%
      dplyr::select(year, iter1000, R = value),
    by = c("year", "iter1000")
  ) %>%
  left_join(
    iteration_map_1000 %>%
      dplyr::select(iter1000, a, b, sigmaR),
    by = "iter1000"
  ) %>%
  mutate(
    R_BH = a * SSB / (b + SSB),
    log_resid = log(R / R_BH)
  )

BH_check_OM <- BH_check %>%
  group_by(om) %>%
  summarise(
    sigmaR = first(sigmaR),
    observed_mean = mean(log_resid),
    observed_sd = sd(log_resid),
    expected_mean = -0.5 * sigmaR^2,
    expected_sd = sigmaR,
    mean_difference = observed_mean - expected_mean,
    sd_difference = observed_sd - expected_sd,
    .groups = "drop"
  )

BH_process_summary <- BH_check_OM %>%
  summarise(
    n_OM = dplyr::n(),
    median_mean_difference = median(mean_difference),
    q05_mean_difference = quantile(mean_difference, 0.05),
    q95_mean_difference = quantile(mean_difference, 0.95),
    median_sd_difference = median(sd_difference),
    q05_sd_difference = quantile(sd_difference, 0.05),
    q95_sd_difference = quantile(sd_difference, 0.95),
    max_abs_mean_difference = max(abs(mean_difference)),
    max_abs_sd_difference = max(abs(sd_difference))
  )

BH_check %>%
  summarise(
    n = dplyr::n(),
    n_NA = sum(is.na(log_resid)),
    n_nonfinite = sum(!is.finite(log_resid))
  )

BH_process_summary %>%
  dplyr::select(
    max_abs_mean_difference,
    max_abs_sd_difference
  )

print(C0_SSB_summary)

print(C0_R_summary)

print(BH_process_summary)

# ============================================================
# 2.4.4 CONTROLLED-FISHING DIAGNOSTIC
# ============================================================
# C7000 is used as a controlled fishing scenario to evaluate
# the coupled response of the conditioned stock and fishery.
#
# A constant annual catch advice of 7,000 t is applied during
# 2025-2054. This scenario is used as a technical perturbation
# of the OM and is not interpreted as a management strategy.
# ============================================================


# ------------------------------------------------------------
# 2.4.4.1 CATCH AND FISHING-MORTALITY RESPONSE
# ------------------------------------------------------------
# Realised catch is compared with the prescribed annual catch
# advice, while fishing mortality is allowed to emerge
# dynamically from the stock-fishery interaction.
# ------------------------------------------------------------

C7000_C_year <- val_C7000$C %>%
  filter(year %in% proj_years) %>%
  group_by(year) %>%
  summarise(
    n = dplyr::n(),
    median_Catch = median(value),
    q05_Catch = quantile(value, 0.05),
    q95_Catch = quantile(value, 0.95),
    min_Catch = min(value),
    max_Catch = max(value),
    p_full = mean(value >= 6999),
    .groups = "drop"
  )

C7000_F_year <- val_C7000$F %>%
  filter(year %in% proj_years) %>%
  group_by(year) %>%
  summarise(
    n = dplyr::n(),
    median_F = median(value),
    q05_F = quantile(value, 0.05),
    q95_F = quantile(value, 0.95),
    min_F = min(value),
    max_F = max(value),
    .groups = "drop"
  )

C7000_C_summary <- C7000_C_year %>%
  summarise(
    min_annual_median_Catch = min(median_Catch),
    max_annual_median_Catch = max(median_Catch),
    min_Catch = min(min_Catch),
    max_Catch = max(max_Catch),
    min_annual_p_full = min(p_full)
  )

C7000_F_summary <- C7000_F_year %>%
  summarise(
    min_annual_median_F = min(median_F),
    max_annual_median_F = max(median_F),
    min_F = min(min_F),
    max_F = max(max_F)
  )

C7000_C_summary

C7000_F_summary


# ------------------------------------------------------------
# 2.4.4.2 CATCH ATTAINMENT AND STOCK AVAILABILITY
# ------------------------------------------------------------
# Catch attainment is evaluated jointly with realised fishing
# mortality and SSB to determine whether failure to realise the
# prescribed catch occurs under reduced stock availability.
# ------------------------------------------------------------

C7000_CFB <- val_C7000$C %>%
  filter(year %in% proj_years) %>%
  dplyr::select(
    year,
    iter1000,
    om,
    replicate,
    Catch = value
  ) %>%
  left_join(
    val_C7000$F %>%
      filter(year %in% proj_years) %>%
      dplyr::select(
        year,
        iter1000,
        F = value
      ),
    by = c("year", "iter1000")
  ) %>%
  left_join(
    val_C7000$B %>%
      filter(year %in% proj_years) %>%
      dplyr::select(
        year,
        iter1000,
        SSB = value
      ),
    by = c("year", "iter1000")
  ) %>%
  mutate(
    attainment = Catch / 7000,
    catch_class = case_when(
      attainment >= 0.99 ~ ">=99%",
      attainment >= 0.90 ~ "90-99%",
      attainment >= 0.50 ~ "50-90%",
      TRUE ~ "<50%"
    ),
    catch_class = factor(
      catch_class,
      levels = c(">=99%", "90-99%", "50-90%", "<50%")
    )
  )

C7000_attainment_summary <- C7000_CFB %>%
  summarise(
    n = dplyr::n(),
    n_below_99 = sum(attainment < 0.99),
    p_below_99 = mean(attainment < 0.99),
    n_below_90 = sum(attainment < 0.90),
    p_below_90 = mean(attainment < 0.90),
    n_below_50 = sum(attainment < 0.50),
    p_below_50 = mean(attainment < 0.50)
  )

C7000_attainment_by_class <- C7000_CFB %>%
  group_by(catch_class) %>%
  summarise(
    n = dplyr::n(),
    median_attainment = median(attainment),
    median_F = median(F),
    q05_F = quantile(F, 0.05),
    q95_F = quantile(F, 0.95),
    median_SSB = median(SSB),
    q05_SSB = quantile(SSB, 0.05),
    q95_SSB = quantile(SSB, 0.95),
    p_blim = mean(SSB < Blim),
    p_bpa = mean(SSB < Bpa),
    .groups = "drop"
  )

C7000_attainment_summary

C7000_attainment_by_class


# ------------------------------------------------------------
# 2.4.4.3 POPULATION RESPONSE UNDER CONTROLLED FISHING
# ------------------------------------------------------------
# Annual SSB distributions are used to characterise the
# biological response of the OM to the controlled fishing
# perturbation.
#
# Probabilities below Blim and Bpa are diagnostics of this
# validation scenario and are not interpreted as MP
# performance metrics.
# ------------------------------------------------------------

C7000_SSB_year <- val_C7000$B %>%
  filter(year %in% proj_years) %>%
  group_by(year) %>%
  summarise(
    n = dplyr::n(),
    median_SSB = median(value),
    q05_SSB = quantile(value, 0.05),
    q95_SSB = quantile(value, 0.95),
    min_SSB = min(value),
    max_SSB = max(value),
    p_blim = mean(value < Blim),
    p_bpa = mean(value < Bpa),
    .groups = "drop"
  )

C7000_SSB_summary <- C7000_SSB_year %>%
  summarise(
    min_annual_median_SSB = min(median_SSB),
    max_annual_median_SSB = max(median_SSB),
    min_SSB = min(min_SSB),
    max_SSB = max(max_SSB),
    max_annual_p_blim = max(p_blim),
    year_max_p_blim = year[which.max(p_blim)],
    max_annual_p_bpa = max(p_bpa),
    year_max_p_bpa = year[which.max(p_bpa)]
  )

C7000_SSB_summary


# ------------------------------------------------------------
# 2.4.4.4 DEPLETION AND CROSS-SCENARIO DIAGNOSTIC
# ------------------------------------------------------------
# Trajectory-level diagnostics identify whether low-biomass
# events are isolated or persistent.
#
# Extremely depleted trajectories under C7000 are then paired
# with the same OM and stochastic replicate under C0. This
# provides a direct check of whether the extreme response is
# also generated by the biological OM in the absence of
# fishing.
# ------------------------------------------------------------

C7000_SSB_traj <- val_C7000$B %>%
  filter(year %in% proj_years) %>%
  group_by(iter1000, om, replicate) %>%
  summarise(
    min_SSB = min(value),
    n_blim = sum(value < Blim),
    n_bpa = sum(value < Bpa),
    ever_blim = any(value < Blim),
    ever_bpa = any(value < Bpa),
    ever_1000 = any(value < 1000),
    ever_100 = any(value < 100),
    .groups = "drop"
  )

C7000_depletion_summary <- C7000_SSB_traj %>%
  summarise(
    n_trajectories = dplyr::n(),
    n_ever_blim = sum(ever_blim),
    p_ever_blim = mean(ever_blim),
    n_ever_bpa = sum(ever_bpa),
    p_ever_bpa = mean(ever_bpa),
    n_ever_1000 = sum(ever_1000),
    p_ever_1000 = mean(ever_1000),
    n_ever_100 = sum(ever_100),
    p_ever_100 = mean(ever_100),
    max_years_blim = max(n_blim),
    max_years_bpa = max(n_bpa),
    n_OM_blim = n_distinct(om[ever_blim]),
    n_OM_1000 = n_distinct(om[ever_1000])
  )

C7000_depletion_summary


# ------------------------------------------------------------
# SAME STOCHASTIC TRAJECTORIES UNDER C0
# ------------------------------------------------------------

extreme_traj <- C7000_SSB_traj %>%
  filter(ever_1000) %>%
  dplyr::select(
    iter1000,
    om,
    replicate,
    min_SSB_C7000 = min_SSB
  )

C0_extreme_pair <- val_C0$B %>%
  filter(
    year %in% proj_years,
    iter1000 %in% extreme_traj$iter1000
  ) %>%
  group_by(iter1000, om, replicate) %>%
  summarise(
    min_SSB_C0 = min(value),
    median_SSB_C0 = median(value),
    ever_blim_C0 = any(value < Blim),
    ever_bpa_C0 = any(value < Bpa),
    .groups = "drop"
  ) %>%
  left_join(
    extreme_traj,
    by = c("iter1000", "om", "replicate")
  ) %>%
  arrange(min_SSB_C7000)

C0_extreme_summary <- C0_extreme_pair %>%
  summarise(
    n_extreme = dplyr::n(),
    min_of_min_SSB_C0 = min(min_SSB_C0),
    max_of_min_SSB_C0 = max(min_SSB_C0),
    n_blim_C0 = sum(ever_blim_C0),
    n_bpa_C0 = sum(ever_bpa_C0)
  )

# ============================================================
# FINAL VALIDATION SUMMARIES
# ============================================================

C0_SSB_summary

C0_R_summary

BH_process_summary

C7000_C_summary

C7000_F_summary

C7000_attainment_summary

C7000_attainment_by_class

C7000_SSB_summary

C7000_depletion_summary

C0_extreme_summary

message(
  "Reference OM open-loop validation completed successfully: ",
  "C0 and C7000, 100 OMs and 1000 trajectories per scenario."
)

# ============================================================
# SAVE WD VALIDATION TABLES
# ============================================================

wd_tab_dir <- here("outputs", "OM", "tables")

dir.create(wd_tab_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(table_continuity, file.path(wd_tab_dir, "reference_OM_continuity.csv"), row.names = FALSE)

write.csv(C0_SSB_year, file.path(wd_tab_dir, "C0_SSB_by_year.csv"), row.names = FALSE)

write.csv(C0_SSB_summary, file.path(wd_tab_dir, "C0_SSB_summary.csv"), row.names = FALSE)

write.csv(C0_R_year, file.path(wd_tab_dir, "C0_recruitment_by_year.csv"), row.names = FALSE)

write.csv(C0_R_summary, file.path(wd_tab_dir, "C0_recruitment_summary.csv"), row.names = FALSE)

write.csv(BH_process_summary, file.path(wd_tab_dir, "BH_process_error_summary.csv"), row.names = FALSE)

write.csv(C7000_C_year, file.path(wd_tab_dir, "C7000_catch_by_year.csv"), row.names = FALSE)

write.csv(C7000_F_year, file.path(wd_tab_dir, "C7000_F_by_year.csv"), row.names = FALSE)

write.csv(C7000_attainment_summary, file.path(wd_tab_dir, "C7000_catch_attainment_summary.csv"), row.names = FALSE)

write.csv(C7000_attainment_by_class, file.path(wd_tab_dir, "C7000_catch_attainment_by_class.csv"), row.names = FALSE)

write.csv(C7000_SSB_year, file.path(wd_tab_dir, "C7000_SSB_by_year.csv"), row.names = FALSE)

write.csv(C7000_SSB_summary, file.path(wd_tab_dir, "C7000_SSB_summary.csv"), row.names = FALSE)

write.csv(C7000_depletion_summary, file.path(wd_tab_dir, "C7000_depletion_summary.csv"), row.names = FALSE)

write.csv(C0_extreme_pair, file.path(wd_tab_dir, "C7000_extreme_trajectories_paired_C0.csv"), row.names = FALSE)
write.csv(C0_extreme_summary, file.path(wd_tab_dir, "C7000_extreme_trajectories_C0_summary.csv"), row.names = FALSE)
