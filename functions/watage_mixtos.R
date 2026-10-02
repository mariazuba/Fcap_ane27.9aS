
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

# Working directory
wd <- getwd()


data_newyear<-paste0(getwd(),"/boot/data/Watege_fun")
lastyear<-2025
 # load("boot/initial/data/update_NewYear/inputbase.RData")
#' 
#' 
#' #'======================================================================================================
#' # Average weight data by age
#' #'======================================================================================================
#' ## Commercial fleet SEINE----
#' wage_seine$weight[wage_seine$weight==0] <- NA
#' wage_seine0 <- wage_seine %>%
#'   filter(!(age %in% c(2, 3) & weight <= 15)) %>%
#'   filter(!(step %in% c(1, 2) & age == 0)) %>%
#'   filter(!(age == 0 & weight <= 2)) %>%
#'   filter(year >= 1989 & weight <= 40)
#' 
#' wage_seine0 <- wage_seine[,c(1,2,4,5)] %>%
#'   pivot_wider(
#'     names_from = age, 
#'     values_from = weight
#'   )%>% mutate(type = "SEINE")
#' #wage_seine0
#' #================================================================
#' ## Pelago spring survey ----
# wage_pela$weight[wage_pela$weight==0] <- NA
# wage_Pelago0 <- wage_pela[,c(1,2,4,5)] %>%
#   mutate(step = case_when(
#     step %in% 1:3   ~ 1,
#     step %in% 4:6   ~ 2,
#     step %in% 7:9   ~ 3,
#     step %in% 10:12 ~ 4
#   ))
# #wage_Pelago0
# wage_Pelago0.1 <- wage_Pelago0 %>%
#   pivot_wider(
#     names_from = age,
#     values_from = weight
#   )%>% mutate(type = "Pelago")
#' #================================================================
#' ## Ecocadiz summer survey ----
#' wage_eco$weight[wage_eco$weight == 0] <- NA
#' wage_Ecocadiz0 <- wage_eco[,c(1,2,4,5)] %>%
#'   mutate(step = case_when(
#'     step %in% 1:3   ~ 1,
#'     step %in% 4:6   ~ 2,
#'     step %in% 7:9   ~ 3,
#'     step %in% 10:12 ~ 4
#'   ))
#' #wage_Ecocadiz0
#' wage_Ecocadiz0.1 <- wage_Ecocadiz0  %>%
#'   pivot_wider(
#'     names_from = age, 
#'     values_from = weight
#'   )%>% mutate(type = "Ecocadiz")
#' 
#' #================================================================
#' ## EcocadizReclutas fall survey----
#' 
#' wage_ecoR$weight[wage_ecoR$weight == 0] <- NA
#' wage_EcocadizRec0 <- wage_ecoR[,c(1,2,4,5)] %>%
#'   mutate(step = case_when(
#'     step %in% 1:3   ~ 1,
#'     step %in% 4:6   ~ 2,
#'     step %in% 7:9   ~ 3,
#'     step %in% 10:12 ~ 4
#'   ))
#' #wage_EcocadizRec0
#' wage_EcocadizRec0.1 <- wage_EcocadizRec0 %>%
#'   pivot_wider(
#'     names_from = age, 
#'     values_from = weight
#'   )%>% mutate(type = "EcocadizRec")
#' 
#' 
#' 
# write.taf(list(wage_eco=wage_Ecocadiz0.1,
#                wage_ecoR=wage_EcocadizRec0.1,
#                wage_pela=wage_Pelago0.1,
#                wage_seine=wage_seine0),dir=data_newyear)



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
#'======================================================================================================
# CLEAN DATA ----
# Filter data with non-NA weight
data_filtered <- dataall1[!is.na(dataall1$weight), ]
# Transform weight column to natural log
data_filtered$log_weight <- log(data_filtered$weight)

data_quarter1 <- data_filtered %>% filter(step == 1) %>% select(year, age, log_weight)%>%
  filter(!is.infinite(log_weight) & !is.na(log_weight))
data_quarter2 <- data_filtered %>% filter(step == 2) %>% select(year, age, log_weight)%>%
  filter(!is.infinite(log_weight) & !is.na(log_weight))
data_quarter3 <- data_filtered %>% filter(step == 3) %>% select(year, age, log_weight)%>%
  filter(!is.infinite(log_weight) & !is.na(log_weight))
data_quarter4 <- data_filtered %>% filter(step == 4) %>% select(year, age, log_weight)%>%
  filter(!is.infinite(log_weight) & !is.na(log_weight))%>%
  filter(!is.infinite(log_weight) & !is.na(log_weight))
data_quarter4 <- data_quarter4[-which(data_quarter4$age == 3 & data_quarter4$log_weight < 3), ]
#'======================================================================================================
# APPLY MIXED MODELS ----

# Fit the mixed linear model
modelo_mixto1 <- lmer(log_weight ~ age + (1 | year), data = data_quarter1)
modelo_mixto2 <- lmer(log_weight ~ age + (1 | year), data = data_quarter2)
modelo_mixto3 <- lmer(log_weight ~ age + (1 | year), data = data_quarter3)
modelo_mixto4 <- lmer(log_weight ~ age + (1 | year), data = data_quarter4)

# View model summary
summary(modelo_mixto1)

#'======================================================================================================
# MIXED MODEL PREDICTIONS ----
# Create a new data frame for predictions

nuevos_datos1 <- expand.grid(year = 1989:lastyear, age = 1:3)
nuevos_datos2 <- expand.grid(year = 1989:lastyear, age = 1:3)
nuevos_datos3 <- expand.grid(year = 1989:lastyear, age = 0:3)
nuevos_datos4 <- expand.grid(year = 1989:lastyear, age = 0:3)

# Convert 'year' and 'age' columns to factors
nuevos_datos1$year <- as.factor(nuevos_datos1$year)
nuevos_datos1$age <- as.factor(nuevos_datos1$age)
nuevos_datos2$year <- as.factor(nuevos_datos2$year)
nuevos_datos2$age <- as.factor(nuevos_datos2$age)
nuevos_datos3$year <- as.factor(nuevos_datos3$year)
nuevos_datos3$age <- as.factor(nuevos_datos3$age)
nuevos_datos4$year <- as.factor(nuevos_datos4$year)
nuevos_datos4$age <- as.factor(nuevos_datos4$age)

# Make predictions using the mixed model
nuevos_datos1$predictions <- predict(modelo_mixto1, newdata = nuevos_datos1, allow.new.levels = TRUE)
nuevos_datos2$predictions <- predict(modelo_mixto2, newdata = nuevos_datos2, allow.new.levels = TRUE)
nuevos_datos3$predictions <- predict(modelo_mixto3, newdata = nuevos_datos3, allow.new.levels = TRUE)
nuevos_datos4$predictions <- predict(modelo_mixto4, newdata = nuevos_datos4, allow.new.levels = TRUE)

# Merge data by 'year' and 'age' columns
datos_unidos1 <- merge(data_quarter1, nuevos_datos1, by = c("year", "age"), all = TRUE)
datos_unidos2 <- merge(data_quarter2, nuevos_datos2, by = c("year", "age"), all = TRUE)
datos_unidos3 <- merge(data_quarter3, nuevos_datos3, by = c("year", "age"), all = TRUE)
datos_unidos4 <- merge(data_quarter4, nuevos_datos4, by = c("year", "age"), all = TRUE)

# Add a column identifying the quarter
datos_unidos1$quarter <- "Q1"
datos_unidos2$quarter <- "Q2"
datos_unidos3$quarter <- "Q3"
datos_unidos4$quarter <- "Q4"

# Combine four datasets into one
datos_combinados <- rbind(datos_unidos1, datos_unidos2, datos_unidos3, datos_unidos4)

# Plot observed and predicted values using facet_wrap to separate by quarter
datos_combinados$year <- as.numeric(datos_combinados$year)

# Format predictions dataset
# Manually assign the step values to each dataset
nuevos_datos1$step <- 1
nuevos_datos2$step <- 2
nuevos_datos3$step <- 3
nuevos_datos4$step <- 4

dataall_end <- rbind(nuevos_datos1, nuevos_datos2, nuevos_datos3, nuevos_datos4)
dataall_end$predictions <- exp(dataall_end$predictions) / 1000

watage_mixtos <- dataall_end %>%
  pivot_wider(names_from = age, values_from = predictions, names_prefix = "age_")
watage_mixtos[is.na(watage_mixtos)] <- 0
watage_mixtos <- watage_mixtos %>% select(year, step, age_0, age_1, age_2, age_3)

# data_esc <- "boot/data"
# write.taf(list(wtage_mixmodel = watage_mixtos), dir = data_esc)

#'======================================================================================================
#Plots ----
#'======================================================================================================
#'*FLEETS*
fig_fleet <- ggplot(dataall1 %>% filter(type=="SEINE"), aes(x = year, y = weight / 1000, color = factor(age)), shape = 1, size = 2) +
  geom_point() +
  geom_line() +
  labs(title = "Commercial fleet - SEINE",
       x = "Year",
       y = "Weight mean (Kg)",
       color = "Age") +
  facet_wrap(~ step, ncol = 2, as.table = TRUE, strip.position = "top",
             labeller = labeller(step = c("1" = "Q1", 
                                          "2" = "Q2",
                                          "3" = "Q3", 
                                          "4" = "Q4"))) +
  theme(panel.grid = element_line(color = NA)) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 6),
        panel.background = element_rect(colour = "gray", fill = "gray99"),
        strip.background = element_rect(colour = "gray", fill = "gray99"),
        legend.title = element_text(size = 6, face = "bold"), 
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate X-axis labels
        axis.text.y = element_text(size = 10),  # Y-axis label size
        axis.title.x = element_text(margin = margin(t = 10)),  # Margin for X-axis title
        axis.title.y = element_text(margin = margin(r = 10))) 

#ggsave(file.path(paste0(path_rep, "/fig_weight_by_quarters_SEINE.png")), fig0,  width = 5, height = 5)

#'======================================================================================================
#'*SURVEYS*
fig_surveys <- ggplot(datasurveys, aes(x = year, y = weight / 1000, color = factor(age)), shape = 1, size = 2) +
  geom_point() +
  geom_line() +
  labs(title = "SURVEYS",
       x = "Year",
       y = "Weight mean (Kg)",
       color = "Age") +
  facet_wrap(~ type, ncol = 2, as.table = TRUE, strip.position = "top",
             labeller = labeller(type = c("Pelago" = "1.PELAGO", 
                                          "Ecocadiz" = "2.ECOCADIZ",
                                          "EcocadizRec" = "3.ECOCADIZ-RECRUITS"))) +
  theme(panel.grid = element_line(color = NA)) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 6),
        panel.background = element_rect(colour = "gray", fill = "gray99"),
        strip.background = element_rect(colour = "gray", fill = "gray99"),
        legend.title = element_text(size = 6, face = "bold"), 
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate X-axis labels
        axis.text.y = element_text(size = 10),  # Y-axis label size
        axis.title.x = element_text(margin = margin(t = 10)),  # Margin for X-axis title
        axis.title.y = element_text(margin = margin(r = 10))) 
#ggsave(file.path(paste0(path_rep, "/fig_weight_by_quarters_SURVEY.png")), fig0.1,  width = 5, height = 5)
#'======================================================================================================

# Plot the improved graph
fig_watage <- ggplot(datos_combinados, aes(x = year)) +
  geom_point(aes(y = exp(log_weight) / 1000, color = factor(age)), shape = 1, size = 2) +
  geom_line(aes(y = exp(predictions) / 1000, color = factor(age), group = age), size = 1) +
  labs(title = "",
       x = "Year",
       y = "Weight mean (Kg)",
       color = "Age") +
  facet_wrap(~ quarter) +
  theme(panel.grid = element_line(color = NA)) +
  theme(plot.title = element_text(size = 5),
        axis.title = element_text(size = 6),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 6),
        panel.background = element_rect(colour = "gray", fill = "gray99"),
        strip.background = element_rect(colour = "gray", fill = "gray99"),
        legend.title = element_text(size = 6, face = "bold"), 
        legend.text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate X-axis labels
        axis.text.y = element_text(size = 10),  # Y-axis label size
        axis.title.x = element_text(margin = margin(t = 10)),  # Margin for X-axis title
        axis.title.y = element_text(margin = margin(r = 10))) + # Margin for Y-axis title
  scale_x_continuous(breaks = seq(min(datos_combinados$year), max(datos_combinados$year), by = 5)) +  # Clearer year interval
  scale_y_continuous(limits = c(0, 50 / 1000))  # Adjust Y-axis limits

#ggsave(file.path(paste0(path_rep, "/fig_weight_by_quarters_obs_est.png")), fig1,  width = 5, height = 5)

# Write .RData ----

# Create data frames for each quarter in a list
watage_fleets <- lapply(c(-2, -1, 0:8), function(f) {
  data.frame(
    year = watage_mixtos$year,
    seas = watage_mixtos$step,
    sex = 1,
    bio_pattern = 1,
    birthseas = 1,
    fleet = f,
    "0" = if (f == -2) watage_mixtos$age_0 * 0 else watage_mixtos$age_0,
    "1" = watage_mixtos$age_1,
    "2" = watage_mixtos$age_2,
    "3" = watage_mixtos$age_3
  )
})

# Combine all data frames in the list into a single data frame and sort by year
watageSS <- bind_rows(watage_fleets) %>% arrange(year)

# Rename columns to match the required format
colnames(watageSS) <- c("year", "seas", "sex", "bio_pattern", "birthseas", "fleet", "0", "1", "2", "3")

# Convert the year column to numeric to avoid errors
watageSS$year <- as.numeric(as.character(watageSS$year))



write.taf(list(watageSS=watageSS),dir=data_newyear)

save(watage_mixtos,fig_fleet,fig_surveys,fig_watage,
     file=paste0(data_newyear,"/watage_mixtos.RData"))





