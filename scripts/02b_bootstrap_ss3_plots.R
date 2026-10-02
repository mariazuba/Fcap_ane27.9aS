# ============================================================
# 02B_bootstrap_ss3_plots.R
# Diagnostic plots for SS3 bootstrap conditioning
#
# Purpose:
#   Generate diagnostic figures comparing the base SS3 model
#   with the 100 validated bootstrap historical conditions
#   retained for operating-model conditioning.
#
# Diagnostics include:
#   - Survey indices
#   - Catch by quarter
#   - Age compositions
#   - Numbers-at-age
#   - Spawning stock biomass (SSB)
#   - Recruitment
#   - Fishing mortality
#   - Selectivity
#   - Survey catchability
#   - Weight-at-age
#
# INPUTS:
#   data/Rdata/
#     - bootstrap_plot_data.rds
#     - iteration_map.csv
#
#   boot/data/ss3_model/
#     - wtatage.ss
#
#   data/bootstrap/wtatage_boot/
#     - wtatage_boot_*.ss
#
# OUTPUTS:
#   outputs/bootstrap/
#     - Bootstrap diagnostic figures (.png)
#
# Notes:
#   - Bootstrap trajectories correspond to the 100 validated
#     SS3 runs selected during bootstrap conditioning.
#   - Black lines/points represent the base SS3 model.
#   - Bootstrap trajectories represent uncertainty in the
#     historical operating-model conditioning.
# ============================================================

rm(list=ls())

library(tidyverse)
library(here)

rds_boot_dir <- here("data","Rdata")
output_boot_dir <- here("outputs","bootstrap")
dir.create(output_boot_dir,recursive=TRUE,showWarnings=FALSE)

dat <- readRDS(file.path(rds_boot_dir,"bootstrap_plot_data.rds"))
list2env(dat,envir=.GlobalEnv)
rm(dat)

# ============================================================
# INDICES
# ============================================================

fig_indices <- ggplot() +
  geom_line(data=boot_indices,aes(year,obs,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
  geom_line(data=base_indices,aes(year,exp),color="black",linewidth=0.5) +
  geom_point(data=base_indices,aes(year,exp),color="black",size=0.7) +
  facet_wrap(~index_name,scales="free_y",ncol=2) +
  labs(x="Year",y="Biomass Survey index",title="") +
  theme_bw() +
  theme(strip.text=element_text(face="bold"),plot.title=element_text(face="bold"))+
  theme(legend.position="none")

fig_indices
ggsave(file.path(output_boot_dir,"bootstrap_indices.png"),fig_indices,width=10,height=7,dpi=300)

# ============================================================
# CATCH
# ============================================================
boot_catch_plot <- boot_catch %>% filter(fleet==seas)
base_catch_plot <- base_catch %>% filter(fleet==seas)

fig_catch <- ggplot() +
  geom_line(data=boot_catch_plot,aes(year,catch,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
   geom_line(data=base_catch_plot,aes(year,obs),color="black",linewidth=0.5) +
  geom_point(data=base_catch_plot,aes(year,obs),color="black",size=0.7) +
  facet_wrap(~seas,scales="free_y",ncol=2,labeller=labeller(seas=c(`1`="Q1",`2`="Q2",`3`="Q3",`4`="Q4"))) +
  labs(x="Year",y="Catch (t)",title="") +
  theme_bw() +
  theme(strip.text=element_text(face="bold"),plot.title=element_text(face="bold"))+
  theme(legend.position="none")


fig_catch
ggsave(file.path(output_boot_dir,"bootstrap_catch.png"),fig_catch,width=10,height=7,dpi=300)

# ============================================================
# AGE COMPOSITIONS
# ============================================================

fleet_names <- c("SEINE_Q1","SEINE_Q2","SEINE_Q3","SEINE_Q4","PELAGO","ECOCADIZ")

agecomp_plots <- setNames(lapply(fleet_names,\(flt) ggplot() +
                                   geom_line(data=filter(boot_agecomp,fleet_name==flt),aes(age,prop,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
                                   geom_line(data=filter(base_agecomp,fleet_name==flt),aes(age,expected,group=year),color="black",linewidth=0.5) +
                                   geom_point(data=filter(base_agecomp,fleet_name==flt),aes(age,expected),color="black",size=0.7) +
                                   facet_wrap(~year,ncol=6) +
                                   scale_y_continuous(limits=c(0,1)) +
                                   labs(x="Age",y="Proportion",title=flt) +
                                   theme_bw()+
                                   theme(strip.text=element_text(face="bold"),plot.title=element_text(face="bold"))+
                                   theme(legend.position="none")),fleet_names)

agecomp_plots$SEINE_Q1
agecomp_plots$SEINE_Q
agecomp_plots$SEINE_Q3
agecomp_plots$SEINE_Q4
agecomp_plots$PELAGO
agecomp_plots$ECOCADIZ

walk(fleet_names,\(flt) ggsave(file.path(output_boot_dir,paste0("bootstrap_agecomp_",flt,".png")),agecomp_plots[[flt]],width=9,height=7,dpi=300))

# ============================================================
# NUMBERS AT AGE
# ============================================================

natage_plots <- setNames(lapply(1:4,\(s) ggplot() +
                                  geom_line(data=filter(boot_natage,Seas==s),aes(Yr,N,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
                                  geom_line(data=filter(base_natage,Seas==s),aes(Yr,N),color="black",linewidth=0.5) +
                                  geom_point(data=filter(base_natage,Seas==s),aes(Yr,N),color="black",size=0.7) +
                                  facet_wrap(~age,scales="free_y",ncol=4,labeller=labeller(age=c(`0`="age-0",`1`="age-1",`2`="age-2",`3`="age-3+"))) +
                                  labs(x="Year",y="Numbers",title=paste("N-at-age  Q",s)) +
                                  theme_bw()+
                                  theme(strip.text=element_text(face="bold"),plot.title=element_text(face="bold"))+
                                  theme(legend.position="none")),paste0("Season_",1:4))

natage_plots$Season_1
natage_plots$Season_2
natage_plots$Season_3
natage_plots$Season_4

walk(1:4,\(s) ggsave(file.path(output_boot_dir,paste0("bootstrap_natage_season",s,".png")),natage_plots[[s]],width=10,height=4,dpi=300))

# ============================================================
# SPAWNING STOCK BIOMASS
# ============================================================

fig_ssb <- ggplot() +
  geom_line(data=boot_ssb,aes(Yr,SSB,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
  geom_line(data=base_ssb,aes(Yr,SSB),color="black",linewidth=0.5) +
  geom_point(data=base_ssb,aes(Yr,SSB),color="black",size=0.7) +
  labs(x="Year",y="SSB (t)",title="") +
  theme_bw() +
  theme(plot.title=element_text(face="bold"),legend.position="none")

fig_ssb
ggsave(file.path(output_boot_dir,"bootstrap_ssb.png"),fig_ssb,width=7,height=4,dpi=300)

# ============================================================
# RECRUITMENT
# ============================================================

fig_rec <- ggplot() +
  geom_line(data=boot_rec,aes(Yr,R,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
  geom_line(data=base_rec,aes(Yr,R),color="black",linewidth=0.5) +
  geom_point(data=base_rec,aes(Yr,R),color="black",size=0.7) +
  labs(x="Year",y="Recruitment") +
  theme_bw() +
  theme(legend.position="none")

fig_rec
ggsave(file.path(output_boot_dir,"bootstrap_recruitment.png"),fig_rec,width=7,height=4,dpi=300)

# ============================================================
# FISHING MORTALITY
# ============================================================

fig_F <- ggplot() +
  geom_line(data=boot_F,aes(Yr,F,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
  geom_line(data=base_F,aes(Yr,F),color="black",linewidth=0.5) +
  geom_point(data=base_F,aes(Yr,F),color="black",size=0.7) +
  facet_wrap(~fleet,ncol=2,scales="free_y") +
  labs(x="Year",y="Fishing mortality",title="") +
  theme_bw()+
  theme(legend.position="none")


fig_F
ggsave(file.path(output_boot_dir,"bootstrap_F.png"),fig_F,width=10,height=7,dpi=300)

# ============================================================
# SELECTIVITY
# ============================================================

survey_names <- c("PELAGO","ECOCADIZ","BOCADEVA","ECORECLUTAS")
fleet_names <- c("SEINE_Q1","SEINE_Q2","SEINE_Q3","SEINE_Q4")

# Years with survey observations
survey_years <- base_indices %>% distinct(index_name,year)

# Restrict selectivity to relevant years
base_sel_survey <- base_sel %>% filter(fleet_name %in% survey_names) %>% inner_join(survey_years,by=c("fleet_name"="index_name","Yr"="year"))
boot_sel_survey <- boot_sel %>% filter(fleet_name %in% survey_names) %>% inner_join(survey_years,by=c("fleet_name"="index_name","Yr"="year"))

base_sel_fleet <- base_sel %>% filter(fleet_name %in% fleet_names,Yr>=1989,Yr<=2024)
boot_sel_fleet <- boot_sel %>% filter(fleet_name %in% fleet_names,Yr>=1989,Yr<=2024)

base_sel_plot <- bind_rows(base_sel_fleet,base_sel_survey)
boot_sel_plot <- bind_rows(boot_sel_fleet,boot_sel_survey)

# Identify unique selectivity profiles within relevant years
base_sel_profiles <- base_sel_plot %>% group_by(fleet_name,Yr) %>% summarise(profile=paste(round(Sel,8),collapse="_"),.groups="drop")

sel_blocks <- base_sel_profiles %>% group_by(fleet_name,profile) %>% summarise(Yr_rep=first(Yr),Yr_ini=min(Yr),Yr_fin=max(Yr),.groups="drop") %>% group_by(fleet_name) %>% arrange(Yr_ini,.by_group=TRUE) %>% mutate(block=row_number(),block_lab=ifelse(Yr_ini==Yr_fin,as.character(Yr_ini),paste0(Yr_ini,"-",Yr_fin))) %>% ungroup()

# Keep one representative year per selectivity profile
base_sel_plot <- base_sel_plot %>% inner_join(sel_blocks %>% select(fleet_name,Yr=Yr_rep,block,block_lab),by=c("fleet_name","Yr"))
boot_sel_plot <- boot_sel_plot %>% inner_join(sel_blocks %>% select(fleet_name,Yr=Yr_rep,block,block_lab),by=c("fleet_name","Yr"))

fleet_names_sel <- unique(base_sel_plot$fleet_name)

# Check survey years and selectivity blocks
base_indices %>% group_by(index_name) %>% summarise(years=paste(sort(unique(year)),collapse=", "),.groups="drop")
sel_blocks

# Plots
sel_seine <- ggplot() +
  geom_line(data=filter(boot_sel_plot,fleet_name %in% fleet_names),aes(age,Sel,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
  geom_line(data=filter(base_sel_plot,fleet_name %in% fleet_names),aes(age,Sel,group=fleet_name),color="black",linewidth=0.6) +
  geom_point(data=filter(base_sel_plot,fleet_name %in% fleet_names),aes(age,Sel),color="black",size=0.7) +
  facet_wrap(~fleet_name,nrow=1) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="Age",y="Selectivity") +
  theme_bw() +
  theme(strip.text=element_text(face="bold"),legend.position="none")

sel_seine
ggsave(file.path(output_boot_dir,"bootstrap_selectivity_SEINE.png"),sel_seine,width=10,height=3,dpi=300)

sel_survey_comp <- bind_rows(
  boot_sel %>% filter(fleet_name=="PELAGO",Yr>=1999,Yr<=2024) %>% mutate(panel="PELAGO 1999-2024"),
  boot_sel %>% filter(fleet_name=="ECOCADIZ",Yr>=2004,Yr<=2014) %>% mutate(panel="ECOCADIZ 2004-2014"),
  boot_sel %>% filter(fleet_name=="ECOCADIZ",Yr>=2015,Yr<=2023) %>% mutate(panel="ECOCADIZ 2015-2023")
)

base_survey_comp <- bind_rows(
  base_sel %>% filter(fleet_name=="PELAGO",Yr>=1999,Yr<=2024) %>% mutate(panel="PELAGO 1999-2024"),
  base_sel %>% filter(fleet_name=="ECOCADIZ",Yr>=2004,Yr<=2014) %>% mutate(panel="ECOCADIZ 2004-2014"),
  base_sel %>% filter(fleet_name=="ECOCADIZ",Yr>=2015,Yr<=2023) %>% mutate(panel="ECOCADIZ 2015-2023")
)

sel_survey_comp <- sel_survey_comp %>% group_by(bootstrap,panel,age) %>% slice(1) %>% ungroup()
base_survey_comp <- base_survey_comp %>% group_by(panel,age) %>% slice(1) %>% ungroup()

fig_sel_survey <- ggplot() +
  geom_line(data=sel_survey_comp,aes(age,Sel,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
  geom_line(data=base_survey_comp,aes(age,Sel,group=panel),color="black",linewidth=0.6) +
  geom_point(data=base_survey_comp,aes(age,Sel),color="black",size=0.7) +
  facet_wrap(~panel,nrow=1) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="Age",y="Selectivity") +
  theme_bw() +
  theme(strip.text=element_text(face="bold"),legend.position="none")

fig_sel_survey
ggsave(file.path(output_boot_dir,"bootstrap_selectivity_surveys.png"),fig_sel_survey,width=9,height=3,dpi=300)
# ============================================================
# CATCHABILITY
# ============================================================

fig_q <- ggplot(boot_q,aes(index_name,q,colour=bootstrap,)) +
  geom_jitter(width=0.12,alpha=0.52,size=1) +
  geom_point(data=base_q,aes(index_name,q),color="black",size=3) +
  labs(x=NULL,y="Catchability (q)",title="Catchability estimates") +
  theme_bw()+
  theme(legend.position="none")


fig_q
ggsave(file.path(output_boot_dir,"bootstrap_catchability.png"),fig_q,width=6,height=4,dpi=300)


# Pesos medios

iteration_map <- read.csv(file.path(rds_boot_dir,"iteration_map.csv"))
# Pesos medios
wt_boot_dir <- here("data","bootstrap","wtatage_boot")
wt_files_100 <- file.path(wt_boot_dir,paste0("wtatage_",iteration_map$bootstrap,".ss"))

wt_base_file <- here("boot","data","ss3_model","wtatage.ss")

wt_base <- read.table(wt_base_file,skip=8,col.names=c("year","seas","sex","bio_pattern","birthseas","fleet","X0","X1","X2","X3")) %>%
  filter(year>=1989,fleet==0) %>%
  pivot_longer(X0:X3,names_to="age",values_to="weight") %>%
  mutate(age=as.integer(sub("X","",age))) %>%
  filter(weight>0)

wt_boot <- map2_dfr(wt_files_100,iteration_map$bootstrap,\(f,id) read.table(f,skip=8,col.names=c("year","seas","sex","bio_pattern","birthseas","fleet","X0","X1","X2","X3")) %>%
                      filter(year>=1989,fleet==0) %>%
                      pivot_longer(X0:X3,names_to="age",values_to="weight") %>%
                      mutate(bootstrap=id,age=as.integer(sub("X","",age))) %>%
                      filter(weight>0))

fig_wt <- ggplot() +
  geom_line(data=wt_boot,aes(year,weight,group=bootstrap,colour=bootstrap),alpha=0.52,linewidth=0.5) +
  geom_line(data=wt_base,aes(year,weight),color="black",linewidth=0.7) +
  geom_point(data=wt_base,aes(year,weight),color="black",size=0.7) +
  facet_grid(age~seas,labeller=labeller(seas=c(`1`="Q1",`2`="Q2",`3`="Q3",`4`="Q4"),age=label_both),scales="free_y") +
  labs(x="Year",y="Mean weight (kg)") +
  theme_bw() +
  theme(strip.text=element_text(face="bold"),legend.position="none")

fig_wt
ggsave(file.path(output_boot_dir,"bootstrap_weight_at_age.png"),fig_wt,width=10,height=8,dpi=300)
