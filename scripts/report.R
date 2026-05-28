##==============================================================================
## FIGURAS MÍNIMAS PARA ESCENARIOS FCAP
##==============================================================================

rm(list = ls())
graphics.off()

library(dplyr)
library(tidyr)
library(stringr)
library(scales)
library(ggplot2)
library(here)

wd <- here()
setwd(wd)

dir.create("report/figs", showWarnings = FALSE, recursive = TRUE)

S <- readRDS("output/summary_flbeia_Fcap_it50.rds")
# Cambia 1000 si tu archivo tiene otro número de iteraciones

bioQ    <- S$bioQ %>% mutate(year = as.numeric(year))
risk    <- S$risk %>% mutate(year = as.numeric(year))
proj.yr <- 2025
Blim    <- S$meta$refs$Blim
Bpa     <- S$meta$refs$Bpa

parse_fcap <- function(df) {
  df %>%
    mutate(
      Fcap_num = as.numeric(
        str_replace(
          str_remove(str_extract(scenario, "Fcap_[0-9p]+"), "Fcap_"),"p", ".")),
      HCR = paste0("F = ", Fcap_num))
}

##==============================================================================
## FIGURA 1: SSB, recruitment, catch y F
##==============================================================================

inds <- c("ssb", "rec", "catch", "f")

bioQ_fcap <- bioQ %>%
  parse_fcap() %>%
  filter(indicator %in% inds)

ref_lines_ssb <- data.frame(
  indicator = "ssb",
  yint = c(Blim, Bpa),
  Ref = c("Blim", "Bpa"))

fig_core <- ggplot(bioQ_fcap,
  aes(x = year, y = q50, color = HCR, fill = HCR)) +
  geom_ribbon(aes(ymin = q05, ymax = q95),alpha = 0.15,colour = NA) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ indicator, scales = "free", ncol = 2) +
  geom_vline(xintercept = proj.yr, linetype = "longdash") +
  geom_hline(
    data = ref_lines_ssb,
    aes(yintercept = yint, linetype = Ref),
    inherit.aes = FALSE,
    linewidth = 0.5,
    color = "black") +
  theme_bw(base_size = 9) +
  theme(legend.position = "right",strip.text = element_text(face = "bold")) +
  labs(x = NULL,y = NULL,color = "Fcap",fill = "Fcap",
    title = "")

ggsave("report/figs/fig1_Fcap_core.png", fig_core,width = 9,height = 5,dpi = 300)
print(fig_core)

##==============================================================================
## FIGURA 2: Riesgo bajo Blim y Bpa
##==============================================================================
risk_fcap <- risk %>%
  parse_fcap() %>%
  mutate(
    year = as.numeric(year),
    indicator = recode(
      indicator,
      pBlim = "P(SSB < Blim)",
      pBpa  = "P(SSB < Bpa)"))

fig_risk <- ggplot(risk_fcap,
  aes(x = year, y = value, color = HCR)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_wrap(~ indicator, ncol = 1) +
  geom_vline(xintercept = proj.yr, linetype = "longdash") +
  geom_hline(yintercept = 0.05, linetype = "longdash") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw(base_size = 9) +
  theme(legend.position = "right",strip.text = element_text(face = "bold")) +
  labs(x = NULL,y = "Probability",color = "Fcap",title = "Bpa = ; Blim =")

ggsave("report/figs/fig2_Fcap_risk.png",fig_risk,width = 7,height = 5,dpi = 300)
print(fig_risk)



