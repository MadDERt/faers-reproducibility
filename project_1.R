# ==============================================================================
# Script: Comprehensive Analysis of ICI-Associated Severe Cardiac Adverse Events
# Description: FAERS data mining to identify and characterize severe cardiac AEs
#              related to immune checkpoint inhibitors (ICIs).
# Author: [Zhangyu Wang]
# Date: [2026/3/18]
# ==============================================================================


# 1. Load Required Libraries ----------------------------------------------
library(future.apply)
library(purrr)
library(UCSCXenaTools)
library(dplyr)
library(stringr)
library(tidyverse)
library(biomaRt)
library(data.table)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggsci)
library(scales)
library(faers)

# 2. Data Acquisition and Preprocessing -----------------------------------
# Define time range: 2014-2022, all quarters
years    <- rep(2014:2022, each = 4)
quarters <- rep(c("q1", "q2", "q3", "q4"), times = 9)

# Download FAERS data
datafaers <- faers(
  years        = years,
  quarters     = quarters,
  format       = "ascii",
  dir          = "~/projects/test/data",
  compress_dir = "~/projects/TCGA/compress"
)

# Standardize drug names using MedDRA v26.0
datafaers_stand <- faers_standardize(datafaers, "~/projects/test/meddra_26.0")

# Remove duplicate reports
datafaers_dedup <- faers_dedup(datafaers_stand)

# 3. Define Immune Checkpoint Inhibitors (ICI) and Build Cohort -----------
# 3.1 ICI drug names (primary suspect drugs)
Immunity_names <- "Nivolumab|Pembrolizumab|Cemiplimab|Atezolizumab|Avelumab|Durvalumab|Dostarlimab"
Immunity_pattern <- paste(Immunity_names, collapse = "|")

# Expand with FDA drug names (brand/generic)
fda_Immunity <- fda_drugs()[grepl(Immunity_pattern, ActiveIngredient, ignore.case = TRUE)]
Immunity_pattern <- paste0(
  unique(tolower(c(Immunity_names, fda_Immunity$DrugName))),
  collapse = "|"
)

# Filter reports where ICI is the primary suspect (role_cod == "PS")
Immunity_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(Immunity_pattern, x$drugname, ignore.case = TRUE) |
              grepl(Immunity_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

# 3.2 Apply inclusion/exclusion criteria

# Include only patients aged ≥18 years
Immunity_data <- faers_filter(
  Immunity_data,
  .fn = function(x) { x[age >= 18, primaryid] },
  .field = "demo"
)

# Exclude records with missing gender or age
Immunity_data <- faers_filter(
  Immunity_data,
  .fn = function(x) {
    idx <- !is.na(x$gender) & !is.na(x$age_in_years)
    x[idx, unique(primaryid)]
  },
  .field = "demo"
)

# Exclude reports with cardiac tumor indications
excluded_indications <- c(
  "benign cardiac neoplasm",
  "cardiac neoplasm malignant",
  "cardiac myxoma",
  "pericarditis malignant",
  "pericardial mesothelioma malignant",
  "pericardial mesothelioma malignant recurrent",
  "pericardial effusion malignant",
  "metastases to heart",
  "cardiac neoplasm unspecified",
  "pericardial neoplasm",
  "malignant pericardial neoplasm",
  "benign pericardium neoplasm",
  "cardiac fibroma",
  "cardiac haemangioma benign",
  "cardiac neurofibroma",
  "cardiac teratoma",
  "cardiac valve fibroelastoma",
  "primary cardiac lymphoma",
  "leukemic cardiac infiltration",
  "pericardial lipoma",
  "cardiac lipoma"
)
excluded_pattern <- paste(excluded_indications, collapse = "|")
Immunity_data <- faers_filter(
  Immunity_data,
  .fn = ~ .x[grepl(excluded_pattern, indi_pt, ignore.case = TRUE), primaryid],
  .field = "indi",
  .invert = TRUE
)

# Exclude concomitant drugs known to cause cardiac toxicity
excluded_concomitant_drugs <- paste(
  "doxorubicin", "epirubicin", "afatinib", "alectinib", "brigatinib",
  "ceritinib", "crizotinib", "erlotinib", "gefitinib", "lorlatinib",
  "osimertinib", sep = "|"
)
Immunity_data <- faers_filter(
  Immunity_data,
  .fn = ~ .x[grepl(excluded_concomitant_drugs, drugname, ignore.case = TRUE) &
               role_cod == "C", primaryid],
  .field = "drug",
  .invert = TRUE
)

# Exclude concomitant anti-CTLA-4 agents (ipilimumab, tremelimumab)
excluded_anti_ctla4 <- "ipilimumab|tremelimumab"
Immunity_data <- faers_filter(
  Immunity_data,
  .fn = ~ .x[grepl(excluded_anti_ctla4, drugname, ignore.case = TRUE), primaryid],
  .field = "drug",
  .invert = TRUE
)

# Remove records missing reporter country or occupation
Immunity_data <- faers_filter(
  Immunity_data,
  .fn = function(x) {
    x[!is.na(reporter_country) & !is.na(occp_cod), primaryid]
  },
  .field = "demo"
)

# 4. Identify Cardiac Adverse Events -------------------------------------

# 4.1 Load MedDRA and extract cardiac PTs (Primary SOC)
meddra_obj <- meddra(path = "~/projects/test/meddra_26.0", primary_soc = TRUE)
hierarchy <- meddra_hierarchy(meddra_obj)
cardiac_pts <- hierarchy[
  soc_name == "Cardiac disorders" & primary_soc_fg == "Y",
  .(pt_code, pt_name)
]
cardiac_pt_names <- unique(cardiac_pts$pt_name)

# 4.2 Disproportionality analysis (PRR) for cardiac events

# Composite signal detection (all cardiac PTs)
results_cd <- faers_phv_signal_composite(
  Immunity_data,
  .event_set = cardiac_pt_names,
  .full      = datafaers_dedup,
  .methods   = "prr"
)

# Signal detection at individual PT level
result_cAEs <- faers_phv_signal(
  Immunity_data,
  .events   = "pt_name",
  .full     = datafaers_dedup,
  .methods  = "prr"
)

# 4.3 Define severe cardiac AEs (scAEs): cases with PRR lower CI > 1 and a ≥ 3
result_scAEs <- result_cAEs[
  pt_name %in% cardiac_pt_names &
    a >= 3 &
    prr_ci_low > 1
]
scAEs_pt <- result_scAEs$pt_name

# Create FAERS object for scAEs reports
scAEs_data <- faers_filter(
  Immunity_data,
  .fn = function(reac_data) {
    reac_data[pt %in% scAEs_pt, primaryid]
  },
  .field = "reac"
)

# 5. Descriptive Statistics of scAEs Cohort -------------------------------
# Age median
median_age <- median(scAEs_data@data$demo$age_in_years, na.rm = TRUE)

# Male proportion
male_count <- sum(scAEs_data@data$demo$sex == "Male", na.rm = TRUE)
total_count <- sum(!is.na(scAEs_data@data$demo$sex))
male_proportion <- male_count / total_count

# Lung cancer indication (subset)
lung_cancer_pts <- c(
  "Adenosquamous cell lung cancer stage III",
  "Large cell lung cancer metastatic",
  "Lung adenocarcinoma",
  "Lung adenocarcinoma metastatic",
  "Lung adenocarcinoma stage II",
  "Lung adenocarcinoma stage III",
  "Lung adenocarcinoma stage IV",
  "Lung cancer metastatic",
  "Lung carcinoma cell type unspecified stage 0",
  "Lung carcinoma cell type unspecified stage III",
  "Lung carcinoma cell type unspecified stage IV",
  "Lung neoplasm malignant",
  "Lung squamous cell carcinoma metastatic",
  "Lung squamous cell carcinoma recurrent",
  "Lung squamous cell carcinoma stage II",
  "Lung squamous cell carcinoma stage III",
  "Lung squamous cell carcinoma stage IV",
  "Metastases to lung",
  "Non-small cell lung cancer",
  "Non-small cell lung cancer metastatic",
  "Non-small cell lung cancer stage II",
  "Non-small cell lung cancer stage III",
  "Non-small cell lung cancer stage IIIA",
  "Non-small cell lung cancer stage IIIB",
  "Non-small cell lung cancer stage IV",
  "Sarcomatoid carcinoma of the lung",
  "Small cell lung cancer",
  "Small cell lung cancer extensive stage",
  "Small cell lung cancer metastatic",
  "Squamous cell carcinoma of lung"
)

lung_cancer_data <- faers_filter(
  scAEs_data,
  .fn = function(indi_data) {
    unique(indi_data[indi_pt %in% lung_cancer_pts, primaryid])
  },
  .field = "indi"
)

# 6. Visualizations -------------------------------------------------------
# 6.1 Volcano Plot for scAEs (Top 5 Strongest Signals)
plot_data <- result_scAEs %>%
  mutate(
    log2PRR = log2(prr),
    se = (log(prr_ci_high) - log(prr_ci_low)) / (2 * 1.96),
    negLog10P = abs(log(prr) / se)  # Approximate Z-score based p-value
  )
top5_signals <- plot_data %>% slice_max(order_by = prr, n = 5)

p_volcano <- ggplot(plot_data, aes(x = log2PRR, y = negLog10P)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_point(aes(size = a, fill = log2PRR), shape = 21, color = "white", alpha = 0.7) +
  geom_point(data = top5_signals, aes(size = a), shape = 21, color = "black", fill = NA, stroke = 1.2) +
  scale_fill_gradient(low = "#FEB24C", high = "#E31A1C", name = "Log2(PRR)") +
  geom_text_repel(
    data = top5_signals,
    aes(label = pt_name),
    family = "serif",
    size = 4,
    fontface = "bold",
    box.padding = 1.0,
    point.padding = 0.5,
    segment.color = "black",
    arrow = arrow(length = unit(0.01, "npc")),
    max.overlaps = Inf
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_classic() +
  labs(
    title = "Top 5 Strongest Signals of Severe Cardiac AEs",
    subtitle = "Disproportionality analysis via Proportional Reporting Ratio (PRR)",
    x = "Signal Strength [Log2 (PRR)]",
    y = "Statistical Significance [-Log10 (P-value)]",
    size = "Case Count"
  ) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    axis.text = element_text(face = "bold", color = "black", size = 11),
    axis.title = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

# 6.2 Demographics: Age distribution and sex proportion
sex_df <- scAEs_data@data$demo %>%
  filter(sex %in% c("Male", "Female")) %>%
  group_by(sex) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

p1_age <- ggplot(scAEs_data@data$demo, aes(x = age_in_years)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 5,
                 fill = "#377EB8", color = "white", alpha = 0.7) +
  geom_density(color = "#E41A1C", linewidth = 1) +
  labs(title = "Age Distribution", x = "Age (Years)", y = "Density") +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

p1_sex <- ggplot(sex_df, aes(x = sex, y = percentage, fill = sex)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            vjust = -0.5, fontface = "bold", family = "serif") +
  scale_fill_manual(values = c("Male" = "#377EB8", "Female" = "#E41A1C")) +
  labs(title = "Gender Proportion", x = NULL, y = "Percentage (%)") +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

layout_p1 <- (p1_age | p1_sex) + plot_annotation(tag_levels = 'A')
# 6.3 Time-to-Onset Distribution (Days from ICI start to cardiac event)
demo_dt <- scAEs_data@data$demo[, .(primaryid, event_dt)]
demo_dt[, event_date := as.Date(as.character(event_dt), format = "%Y%m%d")]

drug_dt <- scAEs_data@data$ther[, .(primaryid, start_dt)]
drug_dt[, start_date := as.Date(as.character(start_dt), format = "%Y%m%d")]

onset_df <- merge(
  demo_dt,
  drug_dt[, .(start_date = min(start_date, na.rm = TRUE)), by = primaryid],
  by = "primaryid"
)
onset_df[, diff_time := as.numeric(event_date - start_date)]

onset_plot_data <- onset_df[diff_time >= 0 & diff_time <= 365]
median_val <- median(onset_plot_data$diff_time, na.rm = TRUE)

p_onset <- ggplot(onset_plot_data, aes(x = diff_time)) +
  geom_density(fill = "#4DAF4A", alpha = 0.3, color = "#4DAF4A", linewidth = 1) +
  geom_rug(alpha = 0.1, color = "grey30") +
  geom_vline(xintercept = median_val, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("label", x = median_val + 40, y = 0.01,
           label = paste("Median:", round(median_val, 0), "days"),
           family = "serif", fontface = "bold", fill = "white") +
  scale_x_continuous(breaks = seq(0, 365, 30), limits = c(0, 180)) +  # Focus on first 6 months
  labs(
    title = "Time-to-Onset Distribution of Cardiac scAEs",
    subtitle = "Days from ICI Treatment Initiation to Cardiac Event",
    x = "Days to Onset",
    y = "Density"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
    axis.text = element_text(face = "bold", color = "black", size = 11),
    axis.title = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

# 6.4 Clinical Outcomes of scAEs
outc_map <- c(
  "DE" = "Death",
  "HO" = "Hospitalization",
  "LT" = "Life-threatening",
  "DS" = "Disability",
  "OT" = "Other Serious",
  "RI" = "Required Intervention",
  "CA" = "Congenital Anomaly"
)

outc_df <- as.data.frame(table(scAEs_data@data[["outc"]]$outc_cod)) %>%
  rename(outc_cod = Var1, count = Freq) %>%
  filter(outc_cod %in% names(outc_map)) %>%
  mutate(
    Outcome = outc_map[as.character(outc_cod)],
    percentage = count / sum(count) * 100
  ) %>%
  arrange(desc(count))

p_outcome <- ggplot(outc_df, aes(x = reorder(Outcome, -count), y = count, fill = Outcome)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, percentage)),
            vjust = -0.2, family = "serif", fontface = "bold", size = 3.5) +
  scale_fill_npg() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Clinical Outcomes of Severe Cardiac AEs",
    subtitle = paste("Total serious cases analyzed: n =", sum(outc_df$count)),
    x = NULL,
    y = "Number of Cases"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30"),
    axis.text = element_text(face = "bold", color = "black", size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title.y = element_text(face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

# 6.5 Combine All Panels into a Multi-Panel Figure
p_volcano <- p_volcano + theme(plot.margin = margin(10, 10, 10, 10))
layout_p1 <- layout_p1 + theme(plot.margin = margin(10, 10, 10, 10))
p_onset   <- p_onset   + theme(plot.margin = margin(10, 10, 10, 10))
p_outcome <- p_outcome + theme(plot.margin = margin(10, 10, 10, 10))

combined_figure <- (p_volcano | layout_p1) / (p_onset | p_outcome)

final_plot <- combined_figure +
  plot_annotation(
    title = "Comprehensive Analysis of ICI-associated Severe Cardiac Adverse Events",
    subtitle = "From Signal Detection to Clinical Characterization (FAERS 2014-2021)",
    tag_levels = 'A',
    caption = "Notes: Panel A uses PRR disproportionality analysis; Panel C includes cases with complete date records."
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 20, family = "serif"),
    plot.title = element_text(face = "bold", size = 18, family = "serif", hjust = 0.5),
    plot.subtitle = element_text(size = 12, family = "serif", hjust = 0.5, color = "grey30")
  )

# Save combined figure
ggsave(
  filename = "Main_Figure_Combined_Analysis.pdf",
  plot = final_plot,
  device = cairo_pdf,
  width = 15,
  height = 13
)
ggsave(
  filename = "Main_Figure_Combined_Analysis.png",
  plot = final_plot,
  type = "cairo",
  dpi = 300,
  width = 15,
  height = 13
)


# 7. Save Workspace for Future Use ----------------------------------------

save(list = ls(all.names = TRUE), file = "faers.RData")
# load("faers.RData")  # Uncomment to reload saved environment

# ==============================================================================
# End of Script
# ==============================================================================
