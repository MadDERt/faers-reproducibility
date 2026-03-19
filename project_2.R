# ==============================================================================
# Script: Association of Antibiotic Use with Secondary Primary Malignancies (SPM)
#         in CAR-T Recipients – A FAERS-Based Disproportionality and TTO Analysis
# Description: This script analyzes the impact of concomitant antibiotic use on
#              the occurrence of Secondary Primary Malignancies (SPM) in patients
#              treated with CAR-T therapies (CD19 and BCMA targeted). It includes
#              signal detection, forest plot of RORs, time-to-onset (TTO) analysis,
#              and demographic summaries.
# Author: [Zhangyu Wang]
# Date: [2026/3/18]
# ==============================================================================

#1. Load Required Libraries ----------------------------------------------
library(faers)
library(dplyr)
library(ggplot2)
library(data.table)
library(survminer)
library(survival)
library(forestploter)
library(grid)
library(cowplot)
library(patchwork)

#2. Data Acquisition and Preprocessing -----------------------------------
# 2.1 Download FAERS data (2017 Q2 – 2024 Q1)
datafaers <- faers(
  years = c(
    2017, 2017, 2017,        # 2017: q2, q3, q4
    2018, 2018, 2018, 2018,  # 2018: q1-q4
    2019, 2019, 2019, 2019,  # 2019: q1-q4
    2020, 2020, 2020, 2020,  # 2020: q1-q4
    2021, 2021, 2021, 2021,  # 2021: q1-q4
    2022, 2022, 2022, 2022,  # 2022: q1-q4
    2023, 2023, 2023, 2023,  # 2023: q1-q4
    2024                       # 2024: q1
  ),
  quarters = c(
    "q2", "q3", "q4",        # 2017
    "q1", "q2", "q3", "q4",  # 2018
    "q1", "q2", "q3", "q4",  # 2019
    "q1", "q2", "q3", "q4",  # 2020
    "q1", "q2", "q3", "q4",  # 2021
    "q1", "q2", "q3", "q4",  # 2022
    "q1", "q2", "q3", "q4",  # 2023
    "q1"                       # 2024
  ),
  dir = "~/projects/test/data",
  compress_dir = "~/projects/TCGA/compress"
)

# Standardize drug names
datafaers_stand <- faers_standardize(datafaers, "~/projects/test/meddra_26.0")
datafaers_dedup <- faers_dedup(datafaers_stand)

#3. Define CAR-T Therapies and Create Sub‑Cohorts ----------------------
# 3.1 Anti‑CD19 CAR‑T
cd19_generic <- "axicabtagene|brexucabtagene|lisocabtagene|tisagenlecleucel"
cd19_brand   <- "yescarta|tecartus|breyanzi|kymriah"
CD19_pattern <- paste0(cd19_generic, "|", cd19_brand)

CD19_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(CD19_pattern, x$drugname, ignore.case = TRUE) |
              grepl(CD19_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

# 3.2 Anti‑BCMA CAR‑T
bcma_generic <- "idecabtagene|ciltacabtagene"
bcma_brand   <- "abecma|carvykti"
BCMA_pattern <- paste0(bcma_generic, "|", bcma_brand)

BCMA_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(BCMA_pattern, x$drugname, ignore.case = TRUE) |
              grepl(BCMA_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

# 3.3 Combined CAR‑T (any target)
CART_pattern <- paste0(CD19_pattern, "|", BCMA_pattern)
CART_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(CART_pattern, x$drugname, ignore.case = TRUE) |
              grepl(CART_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

# 3.4 Individual drugs (for completeness, not used later)
axi_pattern <- "axicabtagene"
axi_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(axi_pattern, x$drugname, ignore.case = TRUE) |
              grepl(axi_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

tisa_pattern <- "tisagenlecleucel"
tisa_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(tisa_pattern, x$drugname, ignore.case = TRUE) |
              grepl(tisa_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

cilta_pattern <- "ciltacabtagene"
cilta_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(cilta_pattern, x$drugname, ignore.case = TRUE) |
              grepl(cilta_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

#4. Signal Detection for Secondary Primary Malignancies (SPM) # --------
# CD19
result_SPMs_CD19 <- faers_phv_signal(
  CD19_data,
  .full = datafaers_dedup,
  .events = "hlgt_name",
  .fn = ~ .x[soc_code == "10029104"]
)

# BCMA
result_SPMs_BCMA <- faers_phv_signal(
  BCMA_data,
  .full = datafaers_dedup,
  .events = "hlgt_name",
  .fn = ~ .x[soc_code == "10029104"]
)

# Combined CAR-T (HLGT level)
result_SPMs_CART <- faers_phv_signal(
  CART_data,
  .full = datafaers_dedup,
  .events = "hlgt_name",
  .fn = ~ .x[soc_code == "10029104"]
)

# Combined CAR-T (LLT level, for finer detail)
result_SPMs_CART_LLT <- faers_phv_signal(
  CART_data,
  .full = datafaers_dedup,
  .events = "llt_name",
  .fn = ~ .x[soc_code == "10029104"]
)

#5. Define Antibiotic Drugs (generic and brand names) --------------------
abx_vector <- c(
  "penicillin", "amoxicillin", "ampicillin", "piperacillin", "oxacillin",
  "dicloxacillin", "nafcillin", "ticarcillin",
  "cef", "cephalexin", "cefazolin", "ceftriaxone", "cefotaxime", "ceftazidime",
  "cefepime", "cefuroxime", "cefdinir", "cefixime", "cefpodoxime",
  "meropenem", "imipenem", "ertapenem", "doripenem",
  "aztreonam",
  "clavulanate", "sulbactam", "tazobactam", "avibactam", "vaborbactam",
  "azithromycin", "erythromycin", "clarithromycin", "telithromycin",
  "tetracycline", "doxycycline", "minocycline", "tigecycline",
  "gentamicin", "tobramycin", "amikacin", "streptomycin", "neomycin",
  "ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin", "norfloxacin",
  "gemifloxacin", "delafloxacin",
  "vancomycin", "teicoplanin", "telavancin", "dalbavancin", "oritavancin",
  "clindamycin", "lincomycin",
  "metronidazole", "tinidazole", "ornidazole",
  "sulfamethoxazole", "trimethoprim", "sulfadiazine", "sulfasalazine",
  "nitrofurantoin", "furazolidone",
  "colistin", "polymyxin",
  "linezolid", "tedizolid",
  "daptomycin", "fosfomycin",
  "isoniazid", "rifampin", "rifampicin", "rifabutin", "rifapentine",
  "ethambutol", "pyrazinamide", "bedaquiline", "delamanid",
  "atovaquone", "pentamidine", "ivermectin", "albendazole", "mebendazole",
  "augmentin", "zithromax", "bactrim", "septra", "cipro", "levaquin",
  "avelox", "vancocin", "flagyl", "zyvox", "cubicin"
)

antibiotic_pattern <- paste0("\\b(", paste(abx_vector, collapse = "|"), ")\\b")

#6. Function to Analyze Association Between Antibiotics and SPM  --------
analyze_cart_abx_spm <- function(cart_data,
                                 abx_pattern = antibiotic_pattern,
                                 target_name = "BCMA") {

  message(paste0("Analyzing ", target_name, " target data..."))

  # Extract drug table
  drug_data <- as.data.table(faers_get(cart_data, "drug"))

  # Identify patients exposed to antibiotics (role_cod == "C" = concomitant)
  with_abx_ids <- drug_data[
    grepl(abx_pattern, drugname, ignore.case = TRUE) &
      role_cod == "C",
    unique(primaryid)
  ]

  # Control group (no antibiotic exposure)
  all_ids <- drug_data[, unique(primaryid)]
  without_abx_ids <- setdiff(all_ids, with_abx_ids)

  # Create FAERS subsets
  cart_with_abx <- faers_keep(cart_data, primaryid = with_abx_ids)
  cart_without_abx <- faers_keep(cart_data, primaryid = without_abx_ids)

  # Identify SPM patients (SOC 10029104)
  reac_with_abx <- as.data.table(faers_get(cart_with_abx, "reac"))
  reac_without_abx <- as.data.table(faers_get(cart_without_abx, "reac"))

  spm_with_ids <- reac_with_abx[soc_code == "10029104", unique(primaryid)]
  spm_without_ids <- reac_without_abx[soc_code == "10029104", unique(primaryid)]

  # Count totals
  n_total_with   <- length(faers_primaryid(cart_with_abx))
  n_total_without <- length(faers_primaryid(cart_without_abx))
  n_spm_with     <- length(spm_with_ids)
  n_spm_without  <- length(spm_without_ids)

  # Build summary table
  res_table <- data.table(
    Target = target_name,
    Group = c("With Antibiotics", "Without Antibiotics"),
    Total_Patients = c(n_total_with, n_total_without),
    SPM_Patients = c(n_spm_with, n_spm_without),
    Non_SPM_Patients = c(n_total_with - n_spm_with, n_total_without - n_spm_without),
    Incidence = c(
      paste0(round(n_spm_with / n_total_with * 100, 2), "%"),
      paste0(round(n_spm_without / n_total_without * 100, 2), "%")
    )
  )

  # Fisher's exact test
  mat <- matrix(c(n_spm_with, n_spm_without,
                  res_table[1, Non_SPM_Patients], res_table[2, Non_SPM_Patients]),
                nrow = 2)
  fisher_res <- fisher.test(mat)

  # ROR (Reporting Odds Ratio)
  ror <- (n_spm_with / (n_total_with - n_spm_with)) /
    (n_spm_without / (n_total_without - n_spm_without))

  message(paste0("--- ", target_name, " analysis complete ---"))
  message(paste0("Antibiotic group incidence: ", res_table[1, Incidence]))
  message(paste0("Non-antibiotic group incidence: ", res_table[2, Incidence]))
  message(paste0("ROR estimate: ", round(ror, 2)))

  return(list(
    table = res_table,
    fisher = fisher_res,
    ror = ror,
    data_with_abx = cart_with_abx,
    data_without_abx = cart_without_abx
  ))
}

# 6.1 Apply the function to each CAR-T target
CART_list <- list(
  BCMA = BCMA_data,
  CD19 = CD19_data,
  CART = CART_data
)
CART_anti_results <- lapply(names(CART_list), function(name) {
  analyze_cart_abx_spm(cart_data = CART_list[[name]], target_name = name)
})
names(CART_anti_results) <- names(CART_list)

# Combine results into a single data table
dt_with <- rbindlist(lapply(CART_anti_results, function(x) x$table[Group == "With Antibiotics",
                                                                   .(Target, SPM_With = SPM_Patients, Non_SPM_With = Non_SPM_Patients, Total_With = Total_Patients)]))
dt_without <- rbindlist(lapply(CART_anti_results, function(x) x$table[Group == "Without Antibiotics",
                                                                      .(Target, SPM_Without = SPM_Patients, Non_SPM_Without = Non_SPM_Patients, Total_Without = Total_Patients)]))

dt_merge <- merge(dt_with, dt_without, by = "Target")

dt_merge[, `:=`(
  ror = (SPM_With / Non_SPM_With) / (SPM_Without / Non_SPM_Without),
  se = sqrt(1/SPM_With + 1/SPM_Without + 1/Non_SPM_With + 1/Non_SPM_Without),
  p_val = sapply(1:.N, function(i) {
    m <- matrix(c(SPM_With[i], Non_SPM_With[i], SPM_Without[i], Non_SPM_Without[i]), nrow = 2)
    fisher.test(m)$p.value
  })
)]

dt_merge[, `:=`(
  lower = exp(log(ror) - 1.96 * se),
  upper = exp(log(ror) + 1.96 * se)
)]

dt_merge[, `:=`(
  Subgroup = ifelse(Target == "CART", "Overall CAR-T", paste0(Target, " Group")),
  `Incidence (With vs Without)` = paste0(SPM_With, "/", Total_With, " vs ",
                                         SPM_Without, "/", Total_Without),
  `ROR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", ror, lower, upper),
  `P Value` = ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)),
  ` ` = paste(rep(" ", 20), collapse = " ")
)]

dt_merge <- dt_merge[order(Target == "CART", decreasing = TRUE)]

#7. Forest Plot for Antibiotic‑SPM Association (Figure Panel A) ---------
GLOBAL_FONT <- "serif"
BASE_SIZE <- 12

theme_forest_unified <- forest_theme(
  base_size = BASE_SIZE,
  base_family = GLOBAL_FONT,
  ci_pch = 15,
  ci_col = "#377EB8",
  ci_fill = "#377EB8",
  ci_lwd = 1.8,
  ci_Theight = 0.2,
  refline_lty = "dashed",
  refline_col = "grey20",
  core = list(padding = unit(c(4, 5), "mm")),
  header = list(fg_params = list(fontface = "bold", size = BASE_SIZE + 1))
)
theme_forest_unified$xlab_hjust <- 0.2

dt_plot <- copy(dt_merge)
dt_plot[, Subgroup := paste0("  ", Subgroup)]

p_forest <- forest(
  data = dt_plot[, .(Subgroup, `Incidence (With vs Without)`, `ROR (95% CI)`, `P Value`, ` `)],
  est = dt_plot$ror,
  lower = dt_plot$lower,
  upper = dt_plot$upper,
  ci_column = 5,
  x_trans = "log10",
  xlim = c(0.1, 15),
  ticks_at = c(0.1, 0.5, 1, 2, 5, 10, 15),
  xlab = "Reporting Odds Ratio (Log Scale)",
  theme = theme_forest_unified
)

p_forest <- add_border(p_forest, part = "header", where = "top", gp = gpar(lwd = 2))
p_forest <- add_border(p_forest, part = "header", where = "bottom", gp = gpar(lwd = 1.5))
p_forest <- add_border(p_forest, part = "body", row = nrow(dt_plot), where = "bottom", gp = gpar(lwd = 2))

#8. Time‑to‑Onset (TTO) Analysis for SPM (with/without antibioti --------
# 8.1 Helper: robust date conversion (handles YYYYMMDD, YYYYMM, YYYY)
to_date_robust <- function(x) {
  x <- as.character(x)
  d <- as.Date(rep(NA, length(x)))

  # 8-digit complete date (YYYYMMDD)
  idx8 <- which(nchar(x) == 8)
  if(length(idx8) > 0) d[idx8] <- as.Date(x[idx8], format = "%Y%m%d")

  # 6-digit (YYYYMM) -> assume first day of month
  idx6 <- which(is.na(d) & nchar(x) == 6)
  if(length(idx6) > 0) d[idx6] <- as.Date(paste0(x[idx6], "01"), format = "%Y%m%d")

  # 4-digit (YYYY) -> assume January 1st
  idx4 <- which(is.na(d) & nchar(x) == 4)
  if(length(idx4) > 0) d[idx4] <- as.Date(paste0(x[idx4], "0101"), format = "%Y%m%d")

  return(d)
}

# 8.2 Main TTO function
run_tto_analysis_final <- function(cart_data, target_name = "BCMA", x_limit = 2000) {

  message(paste0(">>> Running TTO analysis for: ", target_name))

  # Extract data tables
  drug_dt <- as.data.table(faers_get(cart_data, "drug"))
  reac_dt <- as.data.table(faers_get(cart_data, "reac"))
  ther_dt <- as.data.table(faers_get(cart_data, "ther"))
  demo_dt <- as.data.table(faers_get(cart_data, "demo"))

  # Identify antibiotic-exposed patients (same logic as in ROR analysis)
  abx_ids <- drug_dt[
    grepl(antibiotic_pattern, drugname, ignore.case = TRUE) & role_cod == "C",
    unique(primaryid)
  ]

  # ----------------------------------------------------------------------------
  # Retrieve CAR-T infusion date (robust recovery logic)
  # ----------------------------------------------------------------------------
  cart_regex <- "axicabtagene|brexucabtagene|lisocabtagene|tisagenlecleucel|yescarta|tecartus|breyanzi|kymriah|idecabtagene|ciltacabtagene|abecma|carvykti"
  cart_drugs <- drug_dt[grepl(cart_regex, drugname, ignore.case = TRUE)]

  # Dates from drug table (if start_dt exists)
  d_dates <- data.table(primaryid = character(), t_drug = as.Date(character()))
  if ("start_dt" %in% names(cart_drugs)) {
    d_dates <- cart_drugs[!is.na(start_dt) & start_dt != "",
                          .(t_drug = min(to_date_robust(start_dt), na.rm = TRUE)),
                          by = primaryid]
  }

  # Dates from therapy table (ther)
  if ("dsg_drug_seq" %in% names(ther_dt)) {
    setnames(ther_dt, "dsg_drug_seq", "drug_seq", skip_absent = TRUE)
  }
  t_dates <- merge(ther_dt, cart_drugs[, .(primaryid, drug_seq)], by = c("primaryid", "drug_seq"))
  t_dates <- t_dates[!is.na(start_dt) & start_dt != "",
                     .(t_ther = min(to_date_robust(start_dt), na.rm = TRUE)),
                     by = primaryid]

  # Combine and take the earliest available date
  start_dt_final <- merge(d_dates, t_dates, by = "primaryid", all = TRUE)
  start_dt_final[, start_date := pmin(t_drug, t_ther, na.rm = TRUE)]
  start_dt_final <- start_dt_final[, .(primaryid, start_date)]

  # ----------------------------------------------------------------------------
  # Retrieve SPM event date (event_dt, init_fda_dt, mfr_dt – priority order)
  # ----------------------------------------------------------------------------
  spm_ids <- reac_dt[soc_code == "10029104", unique(primaryid)]

  event_dt_df <- demo_dt[primaryid %in% spm_ids,
                         .(primaryid,
                           d_ev  = to_date_robust(event_dt),
                           d_fda = to_date_robust(init_fda_dt),
                           d_mfr = to_date_robust(mfr_dt))]

  # Take the earliest available date among the three
  event_dt_df[, event_date := pmin(d_ev, d_fda, d_mfr, na.rm = TRUE)]
  event_dt_df <- event_dt_df[!is.na(event_date), .(event_date = min(event_date)), by = primaryid]

  # ----------------------------------------------------------------------------
  # Compute TTO and filter valid cases (≥0 and ≤ x_limit)
  # ----------------------------------------------------------------------------
  tto_df <- merge(start_dt_final, event_dt_df, by = "primaryid")
  tto_df[, tto := as.numeric(event_date - start_date)]
  tto_df <- tto_df[tto >= 0 & tto <= x_limit]

  tto_df[, group := ifelse(primaryid %in% abx_ids, "With Antibiotics", "Without Antibiotics")]
  tto_df[, status := 1]  # event indicator

  # ----------------------------------------------------------------------------
  # If too few cases, skip
  # ----------------------------------------------------------------------------
  if (nrow(tto_df) < 2) {
    message(paste0("!!! ", target_name, " insufficient sample, skipping."))
    return(NULL)
  }

  # Fit cumulative hazard (Nelson-Aalen) model
  fit <- survfit(Surv(tto, status) ~ group, data = tto_df)

  # ----------------------------------------------------------------------------
  # Create plot with risk table
  # ----------------------------------------------------------------------------
  theme_clinical_bold <- function() {
    theme_classic() +
      theme(
        text = element_text(family = "serif"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
        axis.text = element_text(face = "bold", color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_rect(colour = "black", fill = "grey90"),
        strip.text = element_text(face = "bold")
      )
  }

  p <- ggsurvplot(
    fit,
    data = tto_df,
    fun = "cumhaz",
    pval = TRUE,
    pval.coord = c(40, 0.15),
    palette = c("#E41A1C", "#377EB8"),
    xlab = "Days from CAR-T Infusion",
    ylab = "Cumulative Hazard",
    font.x = c(10, "bold", "black"),
    font.y = c(10, "bold", "black"),
    legend.labs = c("With Antibiotics", "Without Antibiotics"),
    legend.title = "Exposure Group",
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.18,
    risk.table.y.text = FALSE,
    risk.table.fontsize = 2.8,
    ggtheme = theme_clinical_bold()
  )

  # Fine-tune plot aesthetics
  p$plot <- p$plot +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.margin = margin(t = 20, r = 15, b = 15, l = 10),
      legend.text = element_text(family = "serif", size = 6.5),
      legend.margin = margin(t = -12, b = 0, r = 0, l = 0),
      legend.spacing.x = unit(0.1, "cm")
    )

  p$table <- p$table +
    theme_clinical_bold() +
    theme(
      plot.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line.x = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.4),
      panel.border = element_blank(),
      plot.margin = margin(t = -15, r = 10, b = 5, l = 5)
    )

  # Summary statistics
  stats <- tto_df[, .(Count = .N, Median = median(tto), Q1 = quantile(tto, 0.25), Q3 = quantile(tto, 0.75)), by = group]

  return(list(data = tto_df, plot = p, stats = stats))
}

# 8.3 Run TTO analysis for all three cohorts
tto_results <- lapply(names(CART_list), function(name) {
  run_tto_analysis_final(cart_data = CART_list[[name]], target_name = name)
})
names(tto_results) <- names(CART_list)

#9. Combine Forest Plot and TTO Plots into a Multi‑Panel Figure ( --------
# 9.1 Convert forest plot to grob
forest_grob <- grid::grid.grabExpr(print(p_forest))
# 9.2 Convert each TTO plot to grob (only if not NULL)
grab_styled_surv <- function(fit_obj) {
  if (is.null(fit_obj)) return(NULL)
  grid::grid.grabExpr(print(fit_obj))
}
tto_grobs <- list(
  CART = grab_styled_surv(tto_results[["CART"]][["plot"]]),
  CD19  = grab_styled_surv(tto_results[["CD19"]][["plot"]]),
  BCMA  = grab_styled_surv(tto_results[["BCMA"]][["plot"]])
)

# 9.3 Arrange bottom row (three TTO plots side by side)
bottom_row <- plot_grid(
  tto_grobs[["CART"]], tto_grobs[["CD19"]], tto_grobs[["BCMA"]],
  labels = c("B", "C", "D"),
  label_size = 22,
  label_fontfamily = GLOBAL_FONT,
  ncol = 3,
  align = 'h',
  axis = 'bt'
)

#9.4 Combine forest plot (top) and bottom row
final_manuscript_grid <- plot_grid(
  forest_grob,
  bottom_row,
  labels = c("A", ""),
  label_size = 22,
  label_fontfamily = GLOBAL_FONT,
  ncol = 1,
  rel_heights = c(1, 1.3)
)

# Display and save
grid.newpage()
grid.draw(final_manuscript_grid)

ggsave("Figure4.pdf",
       plot = final_manuscript_grid,
       width = 8, height = 6,
       device = cairo_pdf)

#10. Demographics of the Overall CAR-T Cohort ----------------------------
demo_dt <- as.data.table(CART_data@data$demo)

# Basic age stats
mean_age   <- mean(demo_dt$age, na.rm = TRUE)
median_age <- median(demo_dt$age, na.rm = TRUE)
sd_age     <- sd(demo_dt$age, na.rm = TRUE)
range_age  <- range(demo_dt$age, na.rm = TRUE)

#10.1 Counts by drug target (hard-coded from prior knowledge; adjust if available)
drug_summary <- data.table(
  Target = c("anti_BCMA", "anti_CD19"),
  count = c(1681, 9560)
)

#10.2 Gender summary
sex_summary <- demo_dt[, .(count = .N), by = sex]
sex_summary[sex == "Male",   label := "Male"]
sex_summary[sex == "Female", label := "Female"]
sex_summary[is.na(label),    label := "Unknown/Other"]
sex_summary <- sex_summary[, .(count = sum(count)), by = label]
#10.3 Reporting country (top 8, rest as Others)
country_summary <- demo_dt[, .(count = .N), by = .(Country = reporter_country)][order(-count)]
top_countries <- country_summary[1:8]
others_count <- country_summary[9:.N, sum(count)]
country_summary_plot <- rbind(top_countries, data.table(Country = "Others", count = others_count))

#10.4 Yearly reporting trend
demo_dt[, year := floor(event_dt / 10000)]  # extract year from event_dt (YYYYMMDD)
year_summary <- demo_dt[, .(count = .N), by = year][!is.na(year)][order(year)]

#10.5 Plotting function (unified theme)
theme_baseline <- function() {
  theme_classic() +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text = element_text(color = "black", face = "bold", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)
    )
}

# A. Cases by CAR-T target
p1 <- ggplot(drug_summary, aes(x = reorder(Target, count), y = count, fill = Target)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("#377EB8", "#4DAF4A")) +
  labs(title = "Cases by CAR-T Target", x = "", y = "Case Count") +
  theme_baseline() + theme(legend.position = "none")

# B. Yearly reporting trend
p2 <- ggplot(year_summary, aes(x = year, y = count, group = 1)) +
  geom_line(color = "#E41A1C", linewidth = 1.2) +
  geom_point(size = 3, color = "#E41A1C") +
  labs(title = "Yearly Reporting Trend", x = "Year", y = "N") +
  theme_baseline()

# C. Gender distribution
p3 <- ggplot(sex_summary, aes(x = label, y = count, fill = label)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "grey70")) +
  labs(title = "Gender Distribution", x = "", y = "N") +
  theme_baseline() + theme(legend.position = "none")

# D. Top reporting countries
p4 <- ggplot(country_summary_plot, aes(x = reorder(Country, count), y = count)) +
  geom_bar(stat = "identity", fill = "grey40", color = "black", width = 0.7) +
  coord_flip() +
  labs(title = "Top Reporting Countries", x = "", y = "N") +
  theme_baseline()

#10.6 Combine demographic plots
figure_1 <- (p1 | p2) / (p3 | p4) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold", family = "serif", size = 16))

print(figure_1)
ggsave("CART_Demographics.pdf", figure_1, width = 12, height = 10, device = cairo_pdf)
#11. Save Workspace ------------------------------------------------------
save(list = ls(all.names = TRUE), file = "CAT.RData")
# load("~/projects/TCGA/Rdata/CAT.RData")  # Uncomment to reload

# ==============================================================================
# End of Script
# ==============================================================================
