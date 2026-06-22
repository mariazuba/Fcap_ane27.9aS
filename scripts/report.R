################################################################################
# report.R
# Generate capelin-style diagnostic plots for Fcap x Besc grid
# Anchovy 9aS MSE
################################################################################

rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(here)
  library(stringr)
  library(purrr)
  library(scales)
})

setwd(here())

#===============================================================================
# 1. User options
#===============================================================================

experiment_name <- "Fcap_Besc_grid"

# Choose one SR scenario for SR-specific plots.
# Options expected: "hist", "SS3", "s05", "s07", "extreme"
SR_to_plot <- "s07" #"extreme" #"hist"
CP_cluster_to_plot <- "CP_cluster_1" #"base", "CP_cluster_1","CP_cluster_2","CP_cluster_3"

risk_threshold <- 0.05

# Select HCRs to show in trajectory plots. NULL uses all for the chosen SR.
selected_Fcap <- NULL
selected_Besc <- NULL
# Example:
# selected_Fcap <- c(0.75, 1.00, 1.25, 1.50, 2.00)
# selected_Besc <- c(5000, 6561, 8000, 10000)

#===============================================================================
# 2. Paths and read outputs
#===============================================================================

out_dir    <- here("outputs", "Fcap_results", experiment_name)
report_dir <- file.path(out_dir, "report")
plot_dir   <- file.path(report_dir, "plots")
table_dir  <- file.path(report_dir, "tables")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

obj_file <- file.path(out_dir, "report_outputs.rds")


if (!file.exists(obj_file)) {
  obj_file <- file.path(out_dir, paste0(experiment_name, "_report_outputs.rds"))}
if (!file.exists(obj_file)) {stop("Cannot find report output object in: ", out_dir)}

obj <- readRDS(obj_file)

annual_all        <- obj$annual_all
annual_summary    <- obj$annual_summary
ssb_hist_all      <- obj$ssb_hist_all
risk_by_year      <- obj$risk_by_year
performance_table <- obj$performance_table
candidate_hcrs    <- obj$candidate_hcrs
risk3_by_hcr      <- obj$risk3_by_hcr
settings          <- obj$settings

if (!is.null(settings$risk_threshold)) risk_threshold <- settings$risk_threshold

available_SR <- sort(unique(performance_table$recruitment_scenario))

if (!SR_to_plot %in% available_SR) {
  stop("SR_to_plot not found. Available options: ", paste(available_SR, collapse = ", "))}

available_clusters <- sort(unique(performance_table$CP_cluster))


if (!CP_cluster_to_plot %in% available_clusters) {
  stop(
    "CP_cluster_to_plot not found. Available options: ",
    paste(available_clusters, collapse = ", ")
  )
}

#===============================================================================
# 3. Plot helper functions
#===============================================================================

escape_palette <- function(x) {
  vals <- sort(unique(x))
  pal <- c("#4B0055", "#2C6E91", "#2DBE7E", "#E0AC27", "#777777", "#C65A9A")
  setNames(pal[seq_along(vals)], vals)
}

hcr_label <- function(Fcap, Bescapement) {paste0("Fcap ", Fcap, ", Besc ", Bescapement)}

clean_theme <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      axis.title = element_text(face = "plain"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.9))
    )
}

save_plot <- function(plot, filename, width = 10, height = 7) {
  ggsave(file.path(plot_dir, filename), plot, width = width, height = height, dpi = 300)
}

#===============================================================================
# 4. Filter data for selected SR
#===============================================================================

filter_hcr_selection <- function(df) {
  out <- df |> filter(recruitment_scenario == SR_to_plot,CP_cluster == CP_cluster_to_plot)
  
  if (!is.null(selected_Fcap)) out <- out |> filter(Fcap %in% selected_Fcap)
  if (!is.null(selected_Besc)) out <- out |> filter(Bescapement %in% selected_Besc)
  
  out
}

annual_all_sr     <- filter_hcr_selection(annual_all)
annual_summary_sr <- filter_hcr_selection(annual_summary)
risk_by_year_sr   <- filter_hcr_selection(risk_by_year)
performance_sr    <- filter_hcr_selection(performance_table)


annual_summary_sr <- annual_summary_sr |> mutate(
                      Besc_lab = paste0(Bescapement, " escapement"),
                      Fcap_lab = paste0("Fcap = ", Fcap),
                      hcr_lab = hcr_label(Fcap, Bescapement))

risk_by_year_sr <- risk_by_year_sr |> mutate(
                    Besc_lab = paste0(Bescapement, " escapement"),
                    Fcap_lab = paste0("Fcap = ", Fcap),
                    hcr_lab = hcr_label(Fcap, Bescapement))

performance_sr <- performance_sr |> mutate(
                    Besc_lab = paste0(Bescapement, " escapement"),
                    Fcap_lab = paste0("Fcap = ", Fcap),
                    hcr_lab = hcr_label(Fcap, Bescapement))

cols_besc <- escape_palette(performance_table$Bescapement)

write_csv(performance_sr,file.path(table_dir, paste0("performance_", SR_to_plot, "_", CP_cluster_to_plot, ".csv")))
write_csv(candidate_hcrs, file.path(table_dir, paste0("candidate_HCRs_", SR_to_plot, "_", CP_cluster_to_plot, ".csv")))

#===============================================================================
# 5. Figure 1: SSB trajectories by escapement, chosen SR
#===============================================================================

ssb_hist_sr <- filter_hcr_selection(ssb_hist_all) |>
  mutate(
    Besc_lab = paste0(Bescapement, " escapement"),
    Fcap_lab = paste0("Fcap = ", Fcap),
    hcr_lab  = paste0(CP_cluster, " | Fcap = ", Fcap, ", Besc = ", Bescapement),
    period = "Historical")

ssb_proj_sr <- annual_summary_sr |>
  mutate(
    Besc_lab = paste0(Bescapement, " escapement"),
    Fcap_lab = paste0("Fcap = ", Fcap),
    hcr_lab  = paste0(CP_cluster, " | Fcap = ", Fcap, ", Besc = ", Bescapement),
    period = "Projection")

ssb_plot_df <- bind_rows(ssb_hist_sr, ssb_proj_sr)

# Panel "All compared"
ssb_all <- ssb_plot_df |> mutate(panel = "All compared")
# Panel by Bescapement
ssb_by_besc <- ssb_plot_df |> mutate(panel = Besc_lab)

ssb_plot_long <- bind_rows(ssb_all, ssb_by_besc) |>
                 mutate(panel = factor(panel, levels = c("All compared", sort(unique(Besc_lab)))))

besc_levels <- c("All compared","6561 escapement","5000 escapement", "8000 escapement", "10000 escapement")

ssb_plot_long <- ssb_plot_long |> mutate(panel = factor(panel,levels = besc_levels))

p_ssb <- ggplot(ssb_plot_long, aes(x = year)) +
                 geom_ribbon(aes(ymin = p05_SSB, ymax = p95_SSB, group = hcr_lab,
                                 fill =factor(Fcap)), alpha = 0.15, colour = NA) +
                 geom_line(aes(y = p05_SSB, group = hcr_lab), colour = "black", linetype = "dotted", linewidth = 0.35) +
                 geom_line(aes(y = p95_SSB, group = hcr_lab), colour = "black", linetype = "dotted", linewidth = 0.35) +
                 geom_line(aes(y = median_SSB, colour = factor(Fcap), group = hcr_lab), linewidth = 0.8) +
                 geom_hline(aes(yintercept = Blim), colour = "red", linewidth = 0.5) +
                 facet_wrap(~ panel, scales = "free_y", nrow = 2) +
                 labs(x = NULL, y = "SSB", colour = "Fcap",fill = "Fcap", title = paste("Projected SSB trajectories - SR scenario:", SR_to_plot," - Cluster",CP_cluster_to_plot),
                 subtitle = "Solid line: median; grey band and dotted lines: 5th–95th percentiles; red line: Blim") +
                 clean_theme(base_size = 11) +
                 theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 11))

save_plot(p_ssb, paste0("Fig1_SSB_trajectories_capelin_style_",SR_to_plot,"_",CP_cluster_to_plot,".png"), width = 15, height = 5)
#===============================================================================
# 6. Figure 2: Annual probabilities of opening and risk
#===============================================================================

risk_open_long <- risk_by_year_sr |> select(recruitment_scenario, Fcap, Bescapement, Besc_lab, Fcap_lab, hcr_lab,
                                             year, prob_fishery_open, prob_fishery_closed, prob_SSB_below_Blim, prob_SSB_below_Bescapement) |>
                  pivot_longer(cols = c(prob_fishery_closed, prob_fishery_open), names_to = "fishery_status", values_to = "probability") |>
                  mutate(fishery_status = recode(fishery_status, 
                                                 prob_fishery_closed = "Fishery closed", 
                                                 prob_fishery_open   = "Fishery opened"),
                         fishery_status = factor(fishery_status, levels = c("Fishery closed", "Fishery opened")))

p_open_risk <- ggplot(risk_open_long, aes(x = year)) +
                geom_col(aes(y = probability, fill = fishery_status), position = "stack", width = 0.9, alpha = 0.85) +
                geom_line(aes(y = prob_SSB_below_Blim, group = hcr_lab, colour = "< Blim"), linewidth = 0.65) +
                geom_point(aes(y = prob_SSB_below_Blim, colour = "< Blim"), size = 1.2) +
                geom_hline(aes(yintercept = risk_threshold, linetype = "5% threshold"), colour = "black") +
                facet_grid(Fcap_lab ~ Besc_lab) +
                coord_cartesian(ylim = c(0, 1)) +
                scale_fill_manual(name = NULL, values = c("Fishery closed" = "#CC79A7", "Fishery opened" = "#009E73")) +
                scale_colour_manual(name = NULL, values = c("< Blim" = "black")) +
                scale_linetype_manual(name = NULL, values = c("5% threshold" = "dashed")) +
                labs(x = NULL, y = "Proportion of iterations", 
                     title = paste("Annual fishery opening and risk - SR scenario:", SR_to_plot," - Cluster",CP_cluster_to_plot),
                subtitle = "Bars: fishery opened/closed; line: P(SSB < Blim)") +
                clean_theme(base_size = 11) +
                theme(legend.position = "bottom", legend.box = "horizontal")

save_plot(p_open_risk, paste0("Fig2_annual_opening_risk_", SR_to_plot,"_",CP_cluster_to_plot,".png"), width = 14, height = 10)
#===============================================================================
# 7. Figure 3: Catch-SSB trade-off, chosen SR
#===============================================================================
#==============================================================================
# Fig3. Catch-SSB trade-off, capelin style
#==============================================================================

tradeoff_all  <- performance_sr |> mutate(panel = "All compared")
tradeoff_besc <- performance_sr |> mutate(panel = case_when(Bescapement == 6561 ~ "6561 escapement (Bpa)", TRUE ~ paste0(Bescapement, " escapement")))

tradeoff_plot_df <- bind_rows(tradeoff_all, tradeoff_besc) |>
                    mutate(panel = factor(panel,
                    levels = c("All compared", "6561 escapement (Bpa)", "5000 escapement", "8000 escapement", "10000 escapement")),
                    Fcap_factor = factor(Fcap, levels = c(0.75, 1, 1.25, 1.5, 2)))

cols_fcap <- c("0.75" = "#440154", "1"= "#31688E", "1.25" = "#35B779", "1.5" = "#FDE725", "2" = "#E69F00")

p_tradeoff <- ggplot(tradeoff_plot_df,
              aes(x = median_SSB_open_years, y = median_Catch_open_years, colour = Fcap_factor)) +
              geom_segment(aes(x = p05_SSB_open_years, xend = p95_SSB_open_years, y = median_Catch_open_years, yend = median_Catch_open_years),
              linewidth = 0.7, alpha = 0.7) +
              geom_segment(aes(x = median_SSB_open_years, xend = median_SSB_open_years, y = p05_Catch_open_years, yend = p95_Catch_open_years),
              linewidth = 0.7, alpha = 0.7) +
              geom_point(size = 3.6) +
              facet_wrap(~ panel, nrow = 1) +
              scale_colour_manual(values = cols_fcap, name = "Fcap") +
              labs(x = "Median SSB, open years only", y = "Median catch, open years only",
              title = paste("Catch-SSB trade-off - SR scenario:", SR_to_plot),
              subtitle = "Points: medians; error bars: 5th-95th percentiles across iterations") +
              clean_theme(base_size = 13) +
              theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

save_plot(p_tradeoff, paste0("Fig3_tradeoff_Catch_SSB_", SR_to_plot, "_",CP_cluster_to_plot,".png"),width = 15, height = 5)
#===============================================================================
# 8. Figure 4: Probability of collapse through time
#===============================================================================
risk_by_year_sr <- risk_by_year_sr |>
  mutate(
    Besc_lab = case_when(
      Bescapement == 6561 ~ "6561 escapement (Bpa)",
      TRUE ~ paste0(Bescapement, " escapement")
    ),
    Besc_lab = factor(
      Besc_lab,
      levels = c(
        "5000 escapement",
        "6561 escapement (Bpa)",
        "8000 escapement",
        "10000 escapement"
      )
    )
  )

cols_besc <- c(
  "5000 escapement" = "#440154",
  "6561 escapement (Bpa)" = "#31688E",
  "8000 escapement" = "#35B779",
  "10000 escapement" = "#E69F00"
)

p_collapse <- ggplot(
  risk_by_year_sr,
  aes(
    x = year,
    y = prob_SSB_below_Blim,
    colour = Besc_lab,
    group = hcr_lab
  )
) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_hline(yintercept = risk_threshold, linetype = "dashed") +
  facet_wrap(~ Fcap_lab) +
  scale_colour_manual(values = cols_besc, name = "Bescapement", drop = FALSE) +
  labs(
    x = NULL,
    y = "P(SSB < Blim)",
    title = paste("Annual probability of falling below Blim - SR scenario:", SR_to_plot),
    subtitle = "Dashed line: 5% risk threshold"
  ) +
  clean_theme(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

save_plot(p_collapse, paste0("Fig4_annual_collapse_probability_", SR_to_plot, "_", CP_cluster_to_plot, ".png"), width = 12, height = 7)

#===============================================================================
# 9. Figure 5: Summary performance panel, all SR scenarios
#===============================================================================

panel_df <- performance_table |> mutate(SR_lab      = recruitment_scenario,
                                        Besc_factor = factor(Bescapement),
                                        Besc_lab    = paste0(Bescapement, " escape"),
                                        Fcap_factor = factor(Fcap),
                                        hcr_x       = as.numeric(factor(Bescapement, levels = sort(unique(Bescapement)))))

metric_points <- panel_df |> transmute(SR_lab, Fcap_factor, Besc_factor, Besc_lab, hcr_x,
                    `Max risk across years, % iters <= Blim` = 100 * max_risk3_Blim,
                    `# years > 5%, % iters <= Blim` = n_years_risk3_Blim_gt_5,
                    `Max risk across years, % iters <= Bescapement` = 100 * max_risk3_Bescapement,
                    `# years > 5%, % iters <= Bescapement` = n_years_risk3_Bescapement_gt_5,
                    `# years fishery closed` = n_years_closed_mean,
                    `Median catch in tonnes, open years only` = median_Catch_open_years,
                    `Median SSB in tonnes, all years` = median_SSB_all_years,
                    `Median SSB in tonnes, open years only` = median_SSB_open_years) |>
                pivot_longer(cols = -c(SR_lab, Fcap_factor, Besc_factor, Besc_lab, hcr_x), names_to = "metric", values_to = "value")

metric_intervals <- panel_df |> transmute(SR_lab, Fcap_factor, Besc_factor, Besc_lab, hcr_x,
                        `# years fishery closed_low` = n_years_closed_p05,
                        `# years fishery closed_high` = n_years_closed_p95,
                        `Median catch in tonnes, open years only_low` = p05_Catch_open_years,
                        `Median catch in tonnes, open years only_high` = p95_Catch_open_years,
                        `Median SSB in tonnes, all years_low` = p05_SSB_all_years,
                        `Median SSB in tonnes, all years_high` = p95_SSB_all_years,
                        `Median SSB in tonnes, open years only_low` = p05_SSB_open_years,
                        `Median SSB in tonnes, open years only_high` = p95_SSB_open_years) |>
                    pivot_longer(cols = -c(SR_lab, Fcap_factor, Besc_factor, Besc_lab, hcr_x), names_to = "metric_stat", values_to = "interval_value") |>
                    mutate(bound = if_else(str_detect(metric_stat, "_low$"), "low", "high"), metric = str_remove(metric_stat, "_(low|high)$")) |>
                    select(-metric_stat) |> pivot_wider(names_from = bound, values_from = interval_value)

summary_panel_df <- metric_points |> left_join(metric_intervals, by = c("SR_lab", "Fcap_factor", "Besc_factor", "Besc_lab", "hcr_x", "metric")) |>
                    mutate(metric = factor(metric, levels = c("Max risk across years, % iters <= Blim",
                                                              "# years > 5%, % iters <= Blim",
                                                              "Max risk across years, % iters <= Bescapement",
                                                              "# years > 5%, % iters <= Bescapement",
                                                              "# years fishery closed",
                                                              "Median catch in tonnes, open years only",
                                                              "Median SSB in tonnes, all years",
                                                              "Median SSB in tonnes, open years only")))

p_summary_panel <- ggplot(summary_panel_df, aes(x = Fcap_factor, y = value, colour = Besc_factor)) +
                    geom_hline(data = tibble(metric = factor("Max risk across years, % iters <= Blim", 
                                             levels = levels(summary_panel_df$metric)), yint = 5),
                               aes(yintercept = yint), inherit.aes = FALSE, linetype = "dashed", colour = "red") +
                    geom_errorbar(data = summary_panel_df |> filter(!is.na(low), !is.na(high)),
                               aes(ymin = low, ymax = high), width = 0.08, linewidth = 0.45, alpha = 0.8, position = position_dodge(width = 0.45)) +
                    geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
                    facet_grid(metric ~ SR_lab, scales = "free_y") +
                    scale_colour_manual(values = cols_besc, name = "Bescapement") +
                    labs(x = "Fcap", y = NULL, 
                         title = "Performance statistics by SR scenario",
                         subtitle = "Points: performance value; error bars: 5th-95th percentiles where applicable; red dashed line: 5% risk threshold") +
                    clean_theme(base_size = 10) +
                    theme(strip.text.y = element_text(angle = 0, hjust = 0, face = "bold", colour = "darkred"),
                          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

save_plot(p_summary_panel, "Fig5_summary_performance_panel_all_SR.png", width = 15, height = 14)

#===============================================================================
# 10. Save report objects
#===============================================================================

report_selected_outputs <- list(
                            SR_to_plot        = SR_to_plot,
                            performance_sr    = performance_sr,
                            annual_summary_sr = annual_summary_sr,
                            risk_by_year_sr   = risk_by_year_sr,
                            summary_panel_df  = summary_panel_df,
                            plots             = list(ssb           = p_ssb,
                                                     opening_risk  = p_open_risk,
                                                     tradeoff      = p_tradeoff,
                                                     collapse      = p_collapse,
                                                     summary_panel = p_summary_panel))

saveRDS(report_selected_outputs, file.path(report_dir, paste0("report_outputs_", SR_to_plot, "_",CP_cluster_to_plot,".rds")))

message("Report saved in: ", report_dir)

