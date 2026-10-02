# ============================================================
# 01_bootstrap_weight_at_age.R
# Parametric bootstrap of weight-at-age
#
# Purpose:
#   Generate alternative historical weight-at-age trajectories
#   for the SS3 bootstrap conditioning of the operating model.
#
#   Mixed-effects models are fitted separately by quarter to
#   historical weight-at-age observations (1989-2024).
#   Parametric bootstrap samples are then generated and
#   converted to SS3 wtatage.ss format.
#
# INPUTS:
#   boot/data/Watege_fun/
#     - wage_seine_2.csv
#     - wage_pela.csv
#     - wage_eco.csv
#     - wage_ecoR.csv
#
# OUTPUTS:
#   data/bootstrap/wtatage_boot/
#     - wtatage_boot_001.ss ... wtatage_boot_150.ss
#
#   boot/data/Watege_fun/
#     - weight_at_age_bootstrap.rds
#
# Notes:
#   - Historical period: 1989-2024.
#   - Mixed models: log(weight) ~ age + (1 | year).
#   - Bootstrap is parametric and reproducible using fixed seeds.
#   - Weight predictions are converted from grams to kg.
# ============================================================


# Analysis of Mean Weights by Age and Quarter Using Mixed Models in Fisheries Data

rm(list=ls())
# Load packages -----------------------------------------------------------

library(icesTAF)
library(r4ss)
library(tidyverse)
library(lubridate)
library(ggpubr)
# library(ss3diags)
library(flextable)
library(reshape2)
library(lme4)
library(here)

# Working directory
wd <- getwd()


data_newyear<-here("boot","data","Watege_fun")

# IMPORTANTE:
# modelo SS3 histórico termina en 2024
firstyear <- 1989
lastyear  <- 2024

wage_seine <- read.csv(paste0(data_newyear,"/wage_seine_2.csv")) 
wage_pela <- read.csv(paste0(data_newyear,"/wage_pela.csv"))
wage_eco <- read.csv(paste0(data_newyear,"/wage_eco.csv"))
wage_ecoR <- read.csv(paste0(data_newyear,"/wage_ecoR.csv"))

#'======================================================================================================
# COMBINE DATA ----
datasurveys <- rbind(wage_pela, wage_eco ,wage_ecoR ) %>% melt(id.vars=c("year","step","type"))
datasurveys$variable <- gsub("^X", "", datasurveys$variable)
names(datasurveys)<-c("year","step","type","age","weight")

dataall1 <- rbind(wage_seine, wage_pela, wage_eco ,wage_ecoR)%>% melt(id.vars=c("year","step","type"))
dataall1$variable <- gsub("^X", "", dataall1$variable)
names(dataall1)<-c("year","step","type","age","weight")

#### Bootstrap

# ============================================================
# PREPARAR DATOS PARA BOOTSTRAP
# ============================================================

data_filtered <- dataall1 %>% filter(year >= firstyear,year <= lastyear,!is.na(weight), weight > 0) %>%
  mutate(log_weight = log(weight), year = factor(year), age  = factor(age))

data_quarter1 <- data_filtered %>% filter(step == 1, age %in% c("1", "2", "3")) %>%
  select(year, age, log_weight)
data_quarter2 <- data_filtered %>% filter(step == 2, age %in% c("1", "2", "3")) %>%
  select(year, age, log_weight)
data_quarter3 <- data_filtered %>% filter(step == 3, age %in% c("0", "1", "2", "3")) %>%
  select(year, age, log_weight)
data_quarter4 <- data_filtered %>% filter(step == 4, age %in% c("0", "1", "2", "3")) %>%
  filter(!(age == "3" & log_weight < 3)) %>% select(year, age, log_weight)


# Se ajustan los modelos 


modelo_mixto1 <- lmer(log_weight ~ age + (1 | year), data = data_quarter1, REML = TRUE)
modelo_mixto2 <- lmer(log_weight ~ age + (1 | year), data = data_quarter2, REML = TRUE)
modelo_mixto3 <- lmer(log_weight ~ age + (1 | year), data = data_quarter3, REML = TRUE)
modelo_mixto4 <- lmer(log_weight ~ age + (1 | year), data = data_quarter4, REML = TRUE)


# Matrices para predecir 
new_q1 <- expand.grid(year = factor(firstyear:lastyear), age  = factor(1:3))
new_q2 <- expand.grid(year = factor(firstyear:lastyear), age  = factor(1:3))
new_q3 <- expand.grid(year = factor(firstyear:lastyear), age  = factor(0:3))
new_q4 <- expand.grid(year = factor(firstyear:lastyear), age  = factor(0:3))

# Función de bootstrap
# ============================================================
# FUNCIÓN PARA OBTENER PESOS PREDICHOS
# ============================================================

make_pred_fun <- function(newdata) {
  function(model) {
    pred_log <- predict(model, newdata = newdata, re.form = NULL, allow.new.levels = TRUE)
    # Misma transformación usada en el modelo original
    exp(pred_log) / 1000
  }
}

# ============================================================
# PARAMETRIC BOOTSTRAP
# ============================================================

nboot <- 150

set.seed(20260808)

boot_q1 <- bootMer(modelo_mixto1, FUN = make_pred_fun(new_q1), nsim = nboot, type = "parametric", use.u = FALSE, seed = 101)
boot_q2 <- bootMer(modelo_mixto2, FUN = make_pred_fun(new_q2), nsim = nboot, type = "parametric", use.u = FALSE, seed = 102)
boot_q3 <- bootMer(modelo_mixto3, FUN = make_pred_fun(new_q3), nsim = nboot, type = "parametric", use.u = FALSE, seed = 103)
boot_q4 <- bootMer(modelo_mixto4, FUN = make_pred_fun(new_q4), nsim = nboot, type = "parametric", use.u = FALSE, seed = 104)

#Ahora lo pasamos a formato largo 
boot_to_long <- function(boot_object, newdata, quarter) {
  map_dfr(seq_len(nrow(boot_object$t)), function(b) {
    newdata %>%
      mutate(bootstrap = b, step = quarter, weight = as.numeric(boot_object$t[b, ])) %>%
      transmute(bootstrap, year = as.integer(as.character(year)), step, age = as.integer(as.character(age)), weight)
  })}

wt_boot_long <- bind_rows(
  boot_to_long(boot_q1, new_q1, 1),
  boot_to_long(boot_q2, new_q2, 2),
  boot_to_long(boot_q3, new_q3, 3),
  boot_to_long(boot_q4, new_q4, 4))

# ============================================================
# 10. FORMATO WIDE
# ============================================================

wt_boot_wide <- wt_boot_long %>% 
                mutate(age = paste0("age_",age)) %>%
                pivot_wider(names_from = age,values_from = weight) %>%
                mutate(age_0 = replace_na(age_0,0)) %>%
                select(bootstrap,year,step,age_0,age_1,age_2,age_3) %>%
                arrange(bootstrap,year,step)

wt_boot_wide %>% filter(bootstrap == 1) %>% head()

# ============================================================
# CONSTRUIR wtatage.ss PARA CADA BOOTSTRAP
# ============================================================

make_wtatage_ss <- function(dat) {
  
  fleets <- c(-2, -1, 0:8)
  
  map_dfr(fleets,
    function(f) {
      
      x0 <- if (f == -2) {rep(0, nrow(dat))} else { dat$age_0}
      
      tibble(
        year        = dat$year,
        seas        = dat$step,
        sex         = 1,
        bio_pattern = 1,
        birthseas   = 1,
        fleet       = f,
        X0          = x0,
        X1          = dat$age_1,
        X2          = dat$age_2,
        X3          = dat$age_3)
    }) %>%
    arrange(year, fleet, seas)
}

# ============================================================
# CONSTRUIR LOS 152 wtatage
# ============================================================

watageSS_boot <- wt_boot_wide %>% split(.$bootstrap) %>%
  map(function(x) {
    x %>% select(year,step,age_0, age_1, age_2,age_3) %>%
      make_wtatage_ss()
  })

names(watageSS_boot) <- sprintf("boot_%03d",seq_along(watageSS_boot))


# ============================================================
# FUNCIÓN PARA ESCRIBIR wtatage.ss
# ============================================================

write_wtatage_ss <- function(dat, file) {
  
  header <- c(
    "3 # maxage",
    "# if Yr is negative, then fill remaining years for that Seas, growpattern, Bio_Pattern, Fleet",
    "# if season is negative, then fill remaining fleets for that Seas, Bio_Pattern, Sex, Fleet",
    "# will fill through forecast years, so be careful",
    "# fleet 0 contains begin season pop WT",
    "# fleet -1 contains mid season pop WT",
    "# fleet -2 contains maturity*fecundity",
    "#year seas sex bio_pattern birthseas fleet X0 X1 X2 X3")
  
  writeLines(header,con = file)
  
  # Datos 1989-2024
  write.table(dat,file = file,append = TRUE,row.names = FALSE,col.names = FALSE,quote = FALSE,sep = " ")
  
  # ----------------------------------------------------------
  # FILA TERMINAL SS3
  # ----------------------------------------------------------
  
  terminal_row <- tibble(
    year        = -9999,
    seas        = 1,
    sex         = 1,
    bio_pattern = 1,
    birthseas   = 1,
    fleet       = -2,
    X0          = 0,
    X1          = dat$X1[1],
    X2          = dat$X2[1],
    X3          = dat$X3[1])
  
  write.table(terminal_row,file = file,append = TRUE,row.names = FALSE,col.names = FALSE,quote = FALSE,sep = " ")
}


# ============================================================
# DIRECTORIO DE SALIDA
# ============================================================

wt_boot_dir <- file.path("data","bootstrap","wtatage_boot")

dir.create(wt_boot_dir,recursive = TRUE,showWarnings = FALSE)


# ============================================================
# ESCRIBIR TODOS LOS ARCHIVOS
# ============================================================

iwalk(watageSS_boot,
  function(dat, boot_name) {
    write_wtatage_ss(dat,file.path(wt_boot_dir,paste0("wtatage_",boot_name,".ss")))
  }
)

# validar salida
wt_files <- list.files(wt_boot_dir,pattern = "^wtatage_boot_[0-9]{3}\\.ss$",full.names = TRUE)

length(wt_files)

# ============================================================
# COMPARAR PESOS BASE vs BOOTSTRAP
# ============================================================

# Load the reference weight-at-age object for bootstrap diagnostics.
weight_env <- new.env()
load(file.path(data_newyear, "watage_mixtos.RData"), envir = weight_env)
stopifnot(exists("watage_mixtos", envir = weight_env, inherits = FALSE))
watage_mixtos <- weight_env$watage_mixtos

# ---- 1. Base en formato largo ----

wt_base_long <- watage_mixtos %>%
  pivot_longer(cols = starts_with("age_"), names_to = "age", values_to = "weight_base") %>%
  mutate(year = as.integer(as.character(year)), age = as.integer(gsub("age_", "", age))) %>%
  filter(weight_base > 0)


# ---- 2. Resumen de los bootstrap ----

wt_boot_summary <- wt_boot_long %>% group_by(year, step, age) %>%
  summarise(median = median(weight, na.rm = TRUE),
            p05 = quantile(weight, 0.05, na.rm = TRUE),
            p25 = quantile(weight, 0.25, na.rm = TRUE),
            p75 = quantile(weight, 0.75, na.rm = TRUE),
            p95 = quantile(weight, 0.95, na.rm = TRUE), .groups = "drop")

# ---- 3. Unir con modelo base ----

wt_compare <- wt_boot_summary %>% left_join(wt_base_long, by = c("year", "step", "age"))

#plot 

wt_observed <- data_filtered %>%
  transmute(
    year = as.integer(as.character(year)),
    step = as.integer(step),
    age = as.integer(as.character(age)),
    type = type,
    weight_obs = weight / 1000)


set.seed(123)

boots_show <- sample(unique(wt_boot_long$bootstrap),20)

ggplot() +
  geom_line(data = wt_boot_long %>% filter(bootstrap %in% boots_show),
    aes(x = year, y = weight, group = bootstrap), color = "steelblue", alpha = 0.1, linewidth = 0.6) +  
  geom_point(data = wt_boot_long %>% filter(bootstrap %in% boots_show),
            aes(x = year, y = weight, group = bootstrap), color = "steelblue", alpha = 0.25, linewidth = 0.3) +
  geom_line(data = wt_base_long, aes(x = year, y = weight_base), linetype = 1, linewidth = 0.9) +
  facet_grid(age ~ step, labeller = labeller(
      step = c(`1` = "Q1", `2` = "Q2", `3` = "Q3", `4` = "Q4"), age = label_both), scales = "free_y") +
  labs(x = "Year", y = "Mean weight (kg)") +
  theme_bw()


fig_wt_boot_obs <- ggplot() +
  # Banda bootstrap 5-95%
  geom_ribbon(data = wt_compare,aes(x = year,ymin = p05,ymax = p95),alpha = 0.20) +
  # Mediana bootstrap
  geom_line(data = wt_compare,aes(x = year,y = median),linewidth = 0.7) +
  # Serie base
  geom_line(data = wt_compare,aes(x = year,y = weight_base),linewidth = 0.8,linetype = 2) +
  # Puntos observados
  geom_point(data = wt_observed,aes(x = year,y = weight_obs),shape = 1,size = 1.5,alpha = 0.7) + 
  geom_point(data = wt_observed,aes(x = year,y = weight_obs,shape = type),size = 1.7,alpha = 0.7)+
  facet_grid(age ~ step,labeller = labeller(step = c(`1`="Q1",`2`="Q2",`3`="Q3",`4`="Q4"),age = label_both),scales = "free_y") +
  labs(x = "Year",y = "Mean weight (kg)") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey95"),
    axis.text.x = element_text(angle = 45,hjust = 1))

fig_wt_boot_obs

#
cv_base <- wt_base_long %>% group_by(step, age) %>%
  summarise(mean = mean(weight_base), sd = sd(weight_base), CV = sd / mean, .groups = "drop")

cv_boot <- wt_boot_long %>% group_by(bootstrap, step, age) %>%
  summarise(mean = mean(weight), sd = sd(weight), CV = sd / mean, .groups = "drop")

cv_compare <- cv_boot %>% group_by(step, age) %>%
  summarise(CV_boot_median = median(CV),
            CV_boot_p05 = quantile(CV, 0.05),
            CV_boot_p95 = quantile(CV, 0.95), .groups = "drop") %>%
  left_join(cv_base %>% select(step, age, CV_base = CV), by = c("step", "age"))

print(cv_compare, n = Inf)


# Guardamos la incertidumbre que hemos introducido
weight_boot_diagnostics <- tibble(step = 1:4,
  sd_year = c(0.089606, 0.10557, 0.082151, 0.11345),
  sd_residual = c(0.257516, 0.18510, 0.307423, 0.26282)) %>%
  mutate(ratio = sd_year / sd_residual) %>%
  left_join(cv_compare, by = "step")


# Para watage
wt_boot_wide <- wt_boot_long %>%
  mutate(age = paste0("age_", age)) %>%
  pivot_wider(names_from = age,values_from = weight) %>%
  mutate(age_0 = replace_na(age_0, 0)) %>%
  arrange(bootstrap,year,step)


saveRDS(list(weights = wt_boot_long,
             watageSS = watageSS_boot,
             diagnostics = weight_boot_diagnostics),
        file = file.path(data_newyear, "weight_at_age_bootstrap.rds"))


