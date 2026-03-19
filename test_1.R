# ==============================================================================
# Script: Performance Benchmarking of FAERS Data Processing Pipeline
# Description: This script benchmarks individual steps and the full workflow of
#              the FAERS data processing pipeline, including data parsing,
#              standardization, deduplication, and signal detection for
#              immune checkpoint inhibitors (ICIs). Results are visualized as
#              boxplots.
# Author: [Zhangyu Wang]
# Date: [2026/3/18]
# ==============================================================================


# 1. Load Required Libraries ----------------------------------------------
library(data.table)
library(ggplot2)
library(bench)
library(faers)               # Ensure the custom FAERS package is installed
library(ggsci)               # For scientific color palettes
library(scales)               # For axis label formatting

# 2. Define Global Parameters ---------------------------------------------
# Directory paths
compress_dir <- "~/projects/TCGA/compress"
raw_data_dir <- "~/projects/test/data"
meddra_path  <- "~/projects/test/meddra_26.0"

# 3. Define ICI Drug Patterns (for filtering) ----------------------------
# Generic and brand names of immune checkpoint inhibitors (primary suspect drugs)
Immunity_names <- "Nivolumab|Pembrolizumab|Cemiplimab|Atezolizumab|Avelumab|Durvalumab|Dostarlimab"
Immunity_pattern <- paste(Immunity_names, collapse = "|")

# Expand with FDA drug names (brand names) to ensure comprehensive capture
fda_Immunity <- fda_drugs()[grepl(Immunity_pattern, ActiveIngredient, ignore.case = TRUE)]
Immunity_pattern <- paste0(
  unique(tolower(c(Immunity_names, fda_Immunity$DrugName))),
  collapse = "|"
)

# 4. Prepare Input Objects for Benchmarking #    (Pre-computed obj --------
# Raw FAERS data for Q1–Q4 2015
pre_raw <- faers(
  years        = 2015,
  quarters     = c("q1", "q2", "q3", "q4"),
  dir          = raw_data_dir,
  compress_dir = compress_dir
)

# Standardized data (MedDRA mapping)
pre_std <- faers_standardize(pre_raw, meddra_path)

# Deduplicated data
pre_dedup <- faers_dedup(pre_std)

# Filtered data for ICI cases (primary suspect)
pre_imm <- faers_filter(
  pre_dedup,
  .fn = function(x) {
    idx <- grepl(Immunity_pattern, x$drugname, ignore.case = TRUE)
    x[idx, primaryid]
  },
  .field = "drug"
)

# 5. Run Benchmarking Compare execution times of individual  --------
results_pipeline <- bench::mark(
  # --- Individual steps (using pre‑computed inputs where appropriate) ---
  "Data Parsing" = {
    faers(
      years        = 2015,
      quarters     = c("q1", "q2", "q3", "q4"),
      dir          = raw_data_dir,
      compress_dir = compress_dir
    )
  },
  "Standardization" = {
    faers_standardize(pre_raw, meddra_path)
  },
  "Deduplication" = {
    faers_dedup(pre_std)
  },
  "Signal Detection" = {
    faers_phv_signal(pre_imm, .full = pre_dedup)
  },

  # --- Full end‑to‑end workflow (no pre‑computed inputs) ---
  "Full Workflow" = {
    raw   <- faers(
      years        = 2015,
      quarters     = c("q1", "q2", "q3", "q4"),
      dir          = raw_data_dir,
      compress_dir = compress_dir
    )
    std   <- faers_standardize(raw, meddra_path)
    dedup <- faers_dedup(std)
    imm   <- faers_filter(
      dedup,
      .fn = function(x) {
        idx <- grepl(Immunity_pattern, x$drugname, ignore.case = TRUE)
        x[idx, primaryid]
      },
      .field = "drug"
    )
    res   <- faers_phv_signal(imm, .full = dedup)
  },

  iterations = 5,      # Number of repetitions per expression
  min_time   = 2.0,    # Minimum total time per expression
  memory     = TRUE,   # Track memory usage
  check      = FALSE   # Do not check equality of results (not needed for timing)
)

# 6. Save Benchmark Results -----------------------------------------------
saveRDS(results_pipeline, "results_integrated_summary.rds")

# (Optional: reload if needed)
# results_pipeline <- readRDS("results_integrated_summary.rds")

# 7. Visualize Benchmark Results (Boxplot) --------------------------------
# 7.1 Define mapping from internal names to publication‑ready labels

name_map <- c(
  "Data Parsing"      = "Data Parsing",
  "Standardization"   = "Terminology Mapping",
  "Deduplication"     = "Data Deduplication",
  "Signal Detection"  = "Signal Detection",
  "Full Workflow"     = "Integrated Pipeline"
)

# 7.2 Extract raw timing data from the bench::mark result object
plot_data <- data.table::rbindlist(lapply(1:nrow(results_pipeline), function(i) {
  data.frame(
    expression = as.character(results_pipeline$expression[i]),
    time = as.numeric(results_pipeline$time[[i]])   # time is a list column
  )
}))

# 7.3 Set factor levels to preserve order on x‑axis
plot_data$expression <- factor(
  plot_data$expression,
  levels = names(name_map),
  labels = name_map
)

# 7.4 Create the boxplot
p_integrated <- ggplot(plot_data, aes(x = expression, y = time, fill = expression)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.75,
    color = "black",
    linewidth = 0.6
  ) +
  geom_jitter(
    width = 0.1,
    alpha = 0.5,
    size = 1.2,
    color = "grey20"
  ) +
  scale_y_continuous(
    labels = label_number(accuracy = 0.1),
    name = "Execution Time (Seconds)"
  ) +
  scale_fill_npg() +   # Nature Publishing Group color palette
  labs(x = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "serif", size = 10),
    aspect.ratio = 0.6,                              # Squeeze plot for compactness
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      color = "black",
      face = "bold",
      size = 8
    ),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(face = "bold", size = 10),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 40)
  )

# Display the plot
print(p_integrated)

# 7.5 Save the plot (optional)
ggsave(
  filename = "benchmark_boxplot.pdf",
  plot = p_integrated,
  device = cairo_pdf,
  width = 8,
  height = 6
)

# ==============================================================================
# End of Script
# ==============================================================================
