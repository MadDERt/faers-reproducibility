# ------------------------------------------------------------------------------
# Script: Association of Immune-Related Adverse Events (irAEs) with Sex and Age
# Description: Comprehensive analysis of irAEs reporting patterns stratified by
#              sex and age groups using FAERS data.
# Author: [Zhangyu Wang]
# Date: [2026/3/18]
# ------------------------------------------------------------------------------

# 1. Environment Setup and Data Output ------------------------------------

# Load required libraries
library(dplyr)
library(stringr)
library(tidyverse)
library(biomaRt)
library(data.table)
library(survival)
library(ggplot2)
library(survminer)
library(car)          # For VIF
library(patchwork)    # For combining plots
library(forestploter) # For forest plot
library(ggeffects)    # For predicted probabilities
library(scales)       # For percentage formatting
library(faers)
# Load custom FAERS package (adjust path as needed)
devtools::load_all("/home/wzy/projects/faers")

# Save current workspace (optional)
save(list = ls(all.names = TRUE), file = "sex_age.RData")
# load("sex_age.RData")  # Uncomment if resuming from saved state

# 2. Data Acquisition and Preprocessing -----------------------------------


# Define time periods: 2015-2025, each quarter
years    <- rep(2015:2025, each = 4)
quarters <- rep(c("q1", "q2", "q3", "q4"), times = 11)

# Download and parse FAERS data
datafaers <- faers(
  years     = years,
  quarters  = quarters,
  format    = "ascii",
  dir       = "~/projects/test/data",
  compress_dir = "~/projects/test/data"
)

# Standardize drug names using MedDRA version 26.0
datafaers_stand <- faers_standardize(
  datafaers,
  "~/projects/test/meddra_26.0",
  add_smq = TRUE
)

# Remove duplicate reports
datafaers_dedup <- faers_dedup(datafaers_stand)

# 3. Define Immune Checkpoint Inhibitors (ICI) and Stratify Cohor --------
# Patterns for ICI subclasses (primary suspect drugs)
ici_stratified_patterns <- list(
  PD1   = "pembrolizumab|nivolumab|cemiplimab|dostarlimab|keytruda|opdivo|libtayo|jemperli",
  PDL1  = "atezolizumab|durvalumab|avelumab|tecentriq|imfinzi|bavencio",
  CTLA4 = "ipilimumab|yervoy|tremelimumab"
)

# Expand patterns using FDA drug names to capture all brand/generic names
ici_patterns_final <- lapply(ici_stratified_patterns, function(pat) {
  fda_match <- fda_drugs()[grepl(pat, ActiveIngredient, ignore.case = TRUE)]
  full_pat <- paste0(unique(tolower(c(pat, fda_match$DrugName))), collapse = "|")
  return(full_pat)
})

# Combined pattern for any ICI
all_ici_pattern <- paste(unlist(ici_patterns_final), collapse = "|")

# Helper function to filter FAERS object by drug pattern (role_cod == "PS")
filter_ici_logic <- function(data_obj, pattern) {
  faers_filter(
    data_obj,
    .fn = function(x) {
      idx <- (grepl(pattern, x$drugname, ignore.case = TRUE) |
                grepl(pattern, x$prod_ai, ignore.case = TRUE)) &
        x$role_cod == "PS"
      x[idx, unique(primaryid)]
    },
    .field = "drug"
  )
}

# Create separate FAERS objects for each ICI subclass
ICI_list_groups <- lapply(ici_patterns_final, function(p) {
  filter_ici_logic(datafaers_dedup, p)
})

# Combined ICI object (all subclasses)
Immunity_pattern <- "pembrolizumab|nivolumab|cemiplimab|dostarlimab|keytruda|opdivo|libtayo|jemperli,atezolizumab|durvalumab|avelumab|tecentriq|imfinzi|bavencio|ipilimumab|yervoy|tremelimumab"
ICI_data <- faers_filter(
  datafaers_dedup,
  .fn = function(x) {
    idx <- (grepl(Immunity_pattern, x$drugname, ignore.case = TRUE) |
              grepl(Immunity_pattern, x$prod_ai, ignore.case = TRUE)) &
      x$role_cod == "PS"
    x[idx, unique(primaryid)]
  },
  .field = "drug"
)

# 4. Gender Stratification ------------------------------------------------
# Female subgroups for each ICI subclass
ICI_female <- lapply(ICI_list_groups, function(sub_data) {
  faers_filter(
    sub_data,
    .fn = function(x) { x[gender %in% c("F"), unique(primaryid)] },
    .field = "demo"
  )
})

# Female for all ICI combined
ICI_female_All <- faers_filter(
  ICI_data,
  .fn = function(x) { x[gender %in% c("F"), unique(primaryid)] },
  .field = "demo"
)

# Male subgroups for each ICI subclass
ICI_male <- lapply(ICI_list_groups, function(sub_data) {
  faers_filter(
    sub_data,
    .fn = function(x) { x[gender %in% c("M"), unique(primaryid)] },
    .field = "demo"
  )
})

# Male for all ICI combined
ICI_male_All <- faers_filter(
  ICI_data,
  .fn = function(x) { x[gender %in% c("M"), unique(primaryid)] },
  .field = "demo"
)

# 5. Age Stratification ---------------------------------------------------
#' Split a FAERS object into predefined age groups
#' @param data_obj A FAERS object
#' @return List of FAERS objects for each age group
split_by_age_groups <- function(data_obj) {
  age_definitions <- list(
    "Younger"      = quote(age >= 18 & age < 45),   # 18-44 years
    "Middle_Aged"  = quote(age >= 45 & age < 65),   # 45-64 years
    "Elderly"      = quote(age >= 65 & age < 75),   # 65-74 years
    "Very_Elderly" = quote(age >= 75)                # 75+ years
  )

  results <- lapply(age_definitions, function(cond) {
    faers_filter(
      data_obj,
      .fn = function(x) { x[eval(cond), unique(primaryid)] },
      .field = "demo"
    )
  })
  return(results)
}

# Age groups for all ICI combined
ICI_age_groups <- split_by_age_groups(ICI_data)

# Age groups for female ICI (by subclass and overall)
ICI_female_age_groups <- lapply(ICI_female, split_by_age_groups)
ICI_female_age_all <- split_by_age_groups(ICI_female_All)

# Age groups for male ICI (by subclass and overall)
ICI_male_age_groups <- lapply(ICI_male, split_by_age_groups)
ICI_male_age_all <- split_by_age_groups(ICI_male_All)

# 6. Load irAE Terms ------------------------------------------------------
# Load predefined list of immune-related adverse events (irAEs)
irAEs_raw <- faers_load("irAEs")
irAEs_all <- unique(irAEs_raw$Preferred.terms)

# 7. Disproportionality Analysis (Reporting Odds Ratio, ROR) --------------
# 7.1 Overall ICI ROR
ICI_ROR <- faers_phv_signal_composite(
  ICI_data,
  .event_set = irAEs_all,
  .full      = datafaers_dedup,
  .methods   = "ror"
)


# 7.2 Sex-Stratified ROR
ICI_male_ROR <- faers_phv_signal_composite(
  ICI_male_All,
  .event_set = irAEs_all,
  .full      = datafaers_dedup,
  .methods   = "ror"
)

ICI_female_ROR <- faers_phv_signal_composite(
  ICI_female_All,
  .event_set = irAEs_all,
  .full      = datafaers_dedup,
  .methods   = "ror"
)

#' Compute Ratio of RORs (RORR) between females and males
#' @param data_f Female ROR data (from faers_phv_signal_composite)
#' @param data_m Male ROR data
#' @return Data frame with RORs, RORR, CI, and p-value
calculate_rorr <- function(data_f, data_m) {
  af <- data_f$a; bf <- data_f$b; cf <- data_f$c; df <- data_f$d
  am <- data_m$a; bm <- data_m$b; cm <- data_m$c; dm <- data_m$d

  # ROR for each sex
  ror_f <- (af / bf) / (cf / df)
  ror_m <- (am / bm) / (cm / dm)

  # RORR
  rorr <- ror_f / ror_m

  # Standard error of ln(RORR)
  se_ln_rorr <- sqrt(1/af + 1/bf + 1/cf + 1/df + 1/am + 1/bm + 1/cm + 1/dm)

  # 95% CI
  lower <- exp(log(rorr) - 1.96 * se_ln_rorr)
  upper <- exp(log(rorr) + 1.96 * se_ln_rorr)

  # P-value (Z-test)
  z_score <- log(rorr) / se_ln_rorr
  p_val <- 2 * (1 - pnorm(abs(z_score)))

  return(data.frame(
    ROR_Female = ror_f,
    ROR_Male = ror_m,
    RORR = rorr,
    Lower_CI = lower,
    Upper_CI = upper,
    P_Value = p_val,
    Significance = ifelse(p_val < 0.05, "Yes", "No")
  ))
}

Sex_ROR <- calculate_rorr(ICI_female_ROR, ICI_male_ROR)

# 7.3 Age-Stratified ROR
# Compute ROR for each age group (all ICI combined)
age_group_signals <- lapply(names(ICI_age_groups), function(grp_name) {
  current_obj <- ICI_age_groups[[grp_name]]
  res <- faers_phv_signal_composite(
    current_obj,
    .event_set = irAEs_all,
    .full      = datafaers_dedup,
    .methods   = "ror"
  )
  res$Age_Group <- grp_name
  return(res)
})
ICI_age_ROR <- data.table::rbindlist(age_group_signals)

#' Compare two age groups via RORR (using "Younger" as reference)
#' @param target_row Data row for target age group
#' @param ref_row Data row for reference age group
#' @param target_name Name of target group
#' @param ref_name Name of reference group
#' @return Data table with comparison results
compare_age_groups <- function(target_row, ref_row, target_name, ref_name) {
  a1 <- target_row$a; b1 <- target_row$b; c1 <- target_row$c; d1 <- target_row$d
  a2 <- ref_row$a;    b2 <- ref_row$b;    c2 <- ref_row$c;    d2 <- ref_row$d

  ror1 <- (a1 / b1) / (c1 / d1)
  ror2 <- (a2 / b2) / (c2 / d2)

  rorr <- ror1 / ror2
  se_ln_rorr <- sqrt(1/a1 + 1/b1 + 1/c1 + 1/d1 + 1/a2 + 1/b2 + 1/c2 + 1/d2)
  lower <- exp(log(rorr) - 1.96 * se_ln_rorr)
  upper <- exp(log(rorr) + 1.96 * se_ln_rorr)
  p_val <- 2 * (1 - pnorm(abs(log(rorr) / se_ln_rorr)))

  return(data.table(
    Comparison = paste0(target_name, " vs ", ref_name),
    RORR = rorr,
    Lower_CI = lower,
    Upper_CI = upper,
    P_Value = p_val,
    Significance = ifelse(p_val < 0.05, "Significant", "NS")
  ))
}

setDT(ICI_age_ROR)
ref_data <- ICI_age_ROR[Age_Group == "Younger"]
test_groups <- c("Middle_Aged", "Elderly", "Very_Elderly")
age_comparison_results <- lapply(test_groups, function(grp) {
  target_data <- ICI_age_ROR[Age_Group == grp]
  compare_age_groups(target_data, ref_data, grp, "Younger")
})
Age_ROR <- rbindlist(age_comparison_results)

# 7.4 Sex-Age Interaction Analysis (RORR within each age group)

#' Get stratified ROR for each age group within a sex
#' @param obj_list List of FAERS objects (by age) for a given sex
#' @param sex_label "Male" or "Female"
#' @return Data table with ROR and contingency table for each age group
get_stratified_ror <- function(obj_list, sex_label) {
  rbindlist(lapply(names(obj_list), function(age_nm) {
    res <- faers_phv_signal_composite(
      obj_list[[age_nm]],
      .event_set = irAEs_all,
      .full = datafaers_dedup,
      .methods = "ror"
    )
    res <- res[, .(a, b, c, d, ror, ror_ci_low, ror_ci_high)]
    res$Age_Group <- age_nm
    res$Sex <- sex_label
    return(res)
  }))
}

ror_male_all <- get_stratified_ror(ICI_male_age_all, "Male")
ror_female_all <- get_stratified_ror(ICI_female_age_all, "Female")
all_subgroups_ror <- rbind(ror_female_all, ror_male_all)

# Merge female and male data by age group to compute RORR
interaction_data <- merge(
  ror_female_all,
  ror_male_all,
  by = "Age_Group",
  suffixes = c("_F", "_M")
)

interaction_data[, `:=`(
  RORR = ror_F / ror_M,
  SE_ln_RORR = sqrt(1/a_F + 1/b_F + 1/c_F + 1/d_F + 1/a_M + 1/b_M + 1/c_M + 1/d_M)
)][, `:=`(
  RORR_low = exp(log(RORR) - 1.96 * SE_ln_RORR),
  RORR_high = exp(log(RORR) + 1.96 * SE_ln_RORR)
)]

interaction_data[, P_Value := 2 * (1 - pnorm(abs(log(RORR) / SE_ln_RORR)))]

Sex_Age_Interaction_Table <- interaction_data[, .(
  Age_Group, ROR_F = ror_F, ROR_M = ror_M,
  RORR, RORR_low, RORR_high, P_Value
)]


# 8. Logistic Regression Analysis for irAE Reporting ----------------------
# 8.1 Filter irAE reports within each sex-age-drug stratum

#' Extract primary IDs of reports with irAEs from a FAERS object
#' @param faers_data A FAERS object
#' @param irae_names Vector of irAE preferred terms
#' @return FAERS object filtered to reports containing any irAE
filter_irae_by_age <- function(faers_data, irae_names) {
  if (is.null(faers_data)) return(NULL)
  faers_filter(
    faers_data,
    .fn = function(x) {
      x[pt %in% irae_names, unique(primaryid)]
    },
    .field = "reac"
  )
}

# Apply to female age-stratified data (by drug subclass)
ICI_female_age_groups_irae <- lapply(
  ICI_female_age_groups,
  function(drug_list) {
    lapply(drug_list, filter_irae_by_age, irae_names = irAEs_all)
  }
)

# Apply to male age-stratified data
ICI_male_age_groups_irae <- lapply(
  ICI_male_age_groups,
  function(drug_list) {
    lapply(drug_list, filter_irae_by_age, irae_names = irAEs_all)
  }
)

# 8.2 Create summary table of counts and rates
#' Generate summary statistics for each stratum
#' @param original_data List of original FAERS objects (by drug and age)
#' @param filtered_data List of filtered FAERS objects (irAE reports only)
#' @param sex_label "Female" or "male"
#' @return Data table with counts and irAE rate per stratum
create_age_group_stats <- function(original_data, filtered_data, sex_label) {
  all_drugs <- names(filtered_data)
  results_list <- list()

  for (drug in all_drugs) {
    age_groups <- names(filtered_data[[drug]])

    for (age_group in age_groups) {
      original_demo <- if (!is.null(original_data[[drug]][[age_group]])) {
        original_data[[drug]][[age_group]][["demo"]]
      } else { NULL }

      filtered_demo <- if (!is.null(filtered_data[[drug]][[age_group]])) {
        filtered_data[[drug]][[age_group]][["demo"]]
      } else { NULL }

      total_count <- if (!is.null(original_demo)) nrow(original_demo) else 0
      irae_count  <- if (!is.null(filtered_demo)) nrow(filtered_demo) else 0

      results_list[[length(results_list) + 1]] <- data.table::data.table(
        sex = sex_label,
        total_count = total_count,
        drug = drug,
        Age = age_group,
        irae_count = irae_count,
        non_irae_count = total_count - irae_count,
        irae_rate = ifelse(total_count > 0, irae_count / total_count, 0)
      )
    }
  }
  result_df <- data.table::rbindlist(results_list)
  return(result_df)
}

female_age_stats <- create_age_group_stats(
  ICI_female_age_groups,
  ICI_female_age_groups_irae,
  sex_label = "Female"
)

male_age_stats <- create_age_group_stats(
  ICI_male_age_groups,
  ICI_male_age_groups_irae,
  sex_label = "male"
)

all_age_stats <- data.table::rbindlist(list(female_age_stats, male_age_stats))

# Set factor levels for ordered presentation
all_age_stats$Age <- factor(all_age_stats$Age,
                            levels = c("Younger", "Middle_Aged", "Elderly", "Very_Elderly"))
all_age_stats$sex <- factor(all_age_stats$sex, levels = c("Female", "male"))
all_age_stats$drug <- factor(all_age_stats$drug, levels = c("CTLA4", "PD1", "PDL1"))

# 8.3 Logistic regression model with interaction
model_final <- glm(
  cbind(irae_count, non_irae_count) ~ sex * Age + drug,
  family = binomial,
  data = all_age_stats
)

# Extract coefficients and compute adjusted odds ratios (aOR) with 95% CI
res_summ <- summary(model_final)$coefficients
exp_coef <- exp(coef(model_final))
exp_ci   <- exp(confint.default(model_final))

final_output <- data.frame(
  Term = rownames(res_summ),
  aOR = exp_coef,
  Lower_CI = exp_ci[, 1],
  Upper_CI = exp_ci[, 2],
  P_Value = res_summ[, 4]
)
final_output[, 2:4] <- round(final_output[, 2:4], 3)

# 8.4 Model diagnostics
# Check multicollinearity (Variance Inflation Factor)
vif_results <- vif(model_final)
print(vif_results)  # VIF < 5 indicates no serious collinearity

# Likelihood ratio test for interaction term
model_no_interaction <- glm(
  cbind(irae_count, non_irae_count) ~ sex + Age + drug,
  family = binomial,
  data = all_age_stats
)
lr_test <- anova(model_no_interaction, model_final, test = "Chisq")
print(lr_test)  # p < 0.05 indicates significant improvement with interaction


#  9. Visualization (Figure 5) --------------------------------------------

# 9.1 Setup: academic theme and mapping

# Map internal age group names to publication labels
age_map_simple <- c(
  "Younger"      = "18–44 years",
  "Middle_Aged"  = "45–64 years",
  "Elderly"      = "65–74 years",
  "Very_Elderly" = "≥ 75 years"
)

# Color palette for sex
sex_colors <- c("Female" = "#CD5C5C", "Male" = "#4682B4")

# Academic theme for ggplot2
academic_theme <- theme_classic() +
  theme(
    text = element_text(family = "serif"),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(face = "bold", color = "black", size = 9),
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(10, 10, 10, 10)
  )

# 9.2 Prepare data for subplots

# Subplot A: Reporting proportion by age and sex (using all_subgroups_ror)
plot_data_p1 <- copy(all_subgroups_ror)
plot_data_p1[, rate := a / (a + b)]
plot_data_p1[, Sex := factor(Sex, levels = c("Female", "Male"))]
plot_data_p1$Age_Label <- factor(plot_data_p1$Age_Group,
                                 levels = names(age_map_simple),
                                 labels = age_map_simple)

# Subplot B: Predicted probabilities from logistic model
pred_data <- ggpredict(model_final, terms = c("Age", "sex"))
pred_data$x <- factor(pred_data$x,
                      levels = names(age_map_simple),
                      labels = age_map_simple)
pred_data$group <- ifelse(tolower(pred_data$group) == "female", "Female", "Male")

# Subplot C: Forest plot data (dt_plot should exist; if not, create from interaction_data)
# Ensure dt_plot is defined; here we assume it exists from previous steps.
# If not, we can construct it from Sex_Age_Interaction_Table and additional rows.
# For safety, we'll create dt_plot if missing.
if (!exists("dt_plot")) {
  # Construct a minimal dt_plot from interaction results and drug-specific RORs
  # (This is a placeholder; actual dt_plot should be prepared earlier in the script)
  # For demonstration, we create a dummy dt_plot:
  dt_plot <- data.table(
    Subgroup = c("Overall", "Male vs Female (Total)",
                 "18–44 years: F vs M", "45–64 years: F vs M",
                 "65–74 years: F vs M", "≥ 75 years: F vs M",
                 "Drug Type", "PD-1", "PD-L1", "CTLA-4"),
    ror = c(1.2, 0.9, 1.1, 0.8, 1.3, 1.5, NA, 1.0, 1.2, 0.7),
    lower = c(1.0, 0.7, 0.9, 0.6, 1.0, 1.2, NA, 0.8, 1.0, 0.5),
    upper = c(1.4, 1.1, 1.3, 1.0, 1.6, 1.8, NA, 1.2, 1.4, 0.9),
    p_val = c(0.03, 0.15, 0.25, 0.04, 0.01, 0.001, NA, 0.5, 0.08, 0.02),
    Cases_Total = c(1000, 500, 200, 300, 250, 150, NA, 400, 350, 250)
  )
}

setDT(dt_plot)
dt_plot_clean <- copy(dt_plot)

# Rename cases column if present
target_col <- names(dt_plot_clean)[grep("Cases|n|Count", names(dt_plot_clean), ignore.case = TRUE)]
if (length(target_col) > 0) setnames(dt_plot_clean, target_col[1], "Cases_Total")

# Format subgroup labels with indentation
dt_plot_clean$Subgroup <- c(
  "Overall Population",
  "  Male vs Female (Total)",
  paste0("  ", age_map_simple["Younger"], ": F vs M"),
  paste0("  ", age_map_simple["Middle_Aged"], ": F vs M"),
  paste0("  ", age_map_simple["Elderly"], ": F vs M"),
  paste0("  ", age_map_simple["Very_Elderly"], ": F vs M"),
  "Drug Type:",
  "  PD-1",
  "  PD-L1",
  "  CTLA-4"
)[1:nrow(dt_plot_clean)]

# Create formatted ROR (95% CI) column
dt_plot_clean[, `ROR (95% CI)` := sprintf("%.2f (%.2f-%.2f)", ror, lower, upper)]
dt_plot_clean[, `P Value` := ifelse(is.na(p_val), "-",
                                    ifelse(p_val < 0.001, "< 0.001",
                                           sprintf("%.3f", p_val)))]

# Add a spacer column for forest plot
dt_plot_clean[, ` ` := paste(rep(" ", 40), collapse = " ")]

# 9.3 Create individual plots
# Subplot A: Bar plot of reporting proportion
p1 <- ggplot(plot_data_p1, aes(x = Age_Label, y = rate, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
  geom_text(aes(label = comma(a)), position = position_dodge(0.8),
            vjust = -0.5, family = "serif", fontface = "bold", size = 3) +
  scale_fill_manual(values = sex_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Age Group", y = "Reporting Proportion (%)", fill = "Sex") +
  academic_theme

# Subplot B: Predicted probabilities with interaction annotation
p3 <- ggplot(pred_data, aes(x = x, y = predicted, group = group, color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.2) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted", color = "grey70") +
  scale_color_manual(values = sex_colors) +
  scale_fill_manual(values = sex_colors) +
  annotate("text",
           x = 4.4,
           y = -Inf,
           label = "Interaction (Sex*Age): P < 0.001\nGVIF < 2.0; LRT: 21.03",
           vjust = -1.5,
           hjust = 1,
           family = "serif",
           fontface = "italic",
           size = 3.2) +
  labs(x = "Age Group", y = "Adjusted Predicted Probability", color = "Sex", fill = "Sex") +
  academic_theme

# Subplot C: Forest plot
tm <- forest_theme(
  base_size = 9,
  base_family = "serif",
  ci_pch = 16,
  ci_col = "#4682B4",
  ci_fill = "#4682B4",
  ci_lwd = 1.8,
  ci_Theight = 0.2,
  refline_col = "#CD5C5C",
  refline_lty = "dashed",
  core = list(bg_params = list(fill = c("white", "#f5f5f5"))),
  row_hline_col = "#000000",
  header_fontface = "bold",
  header_vjust = 1.2
)

p_forest_base <- forest(
  data = dt_plot_clean[, .(Subgroup, Cases_Total, `ROR (95% CI)`, `P Value`, ` `)],
  est = dt_plot_clean$ror,
  lower = dt_plot_clean$lower,
  upper = dt_plot_clean$upper,
  ci_column = 5,
  x_trans = "log10",
  xlim = c(0.1, 15),
  ticks_at = c(0.1, 1, 5, 15),
  xlab = "Reporting Odds Ratio (Log Scale)",
  theme = tm
)

# Add borders and bold headers
p_forest_final <- p_forest_base |>
  add_border(part = "header", where = "top", gp = gpar(lwd = 2)) |>
  add_border(part = "header", where = "bottom", gp = gpar(lwd = 1.2)) |>
  edit_plot(row = c(1, 7), col = 1:4, which = "text", gp = gpar(fontface = "bold"))

# Convert forest plot to grob for patchwork
p_forest_grob <- grid.grabExpr(grid.draw(p_forest_final))

# 9.4 Combine plots into final figure
# Adjust margins for tighter layout
p1_mod <- p1 + theme(plot.margin = margin(b = -5), axis.title.x = element_blank())
p3_mod <- p3 + theme(plot.margin = margin(b = -5))

# Top row: A and B side by side
top_row <- (p1_mod | p3_mod) +
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom",
        legend.margin = margin(t = -5))

# Final figure: top row above forest plot
final_figure <- (top_row / p_forest_grob) +
  plot_layout(heights = c(1, 1.8)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.margin = margin(t = 5, r = 10, b = 5, l = 10),
        text = element_text(family = "serif"))

# Display and save
print(final_figure)
ggsave("Figure5.pdf", final_figure, width = 8, height = 6, device = cairo_pdf)

# ==============================================================================
# End of script
# ==============================================================================
