# ==============================================================================
# Script: Performance Scaling Benchmark of FAERS Processing Pipeline
# Description: This script evaluates how the FAERS data processing pipeline
#              scales with increasing data volume (number of quarters). It
#              measures execution time and memory allocation for the full
#              pipeline (parsing, standardization, deduplication, signal
#              detection) across 1, 4, 8, 16, 24, and 32 quarters. Results
#              are visualized as dual‑axis (time + memory) and throughput plots.
# Author: [Zhangyu Wang]
# Date: [2026/3/18]
# ==============================================================================


# 1. Load Required Libraries ----------------------------------------------
library(data.table)
library(ggplot2)
library(bench)
library(dplyr)
library(patchwork)
library(scales)
library(faers)
# 2. Define Global Paths and Parameters -----------------------------------
compress_dir <- "~/projects/TCGA/compress"
raw_data_dir <- "~/projects/test/data"
meddra_path  <- "~/projects/test/meddra_26.0"

# Target drug pattern for signal detection (example: anti‑PD‑1 agents)
target_drug <- "Nivolumab|Pembrolizumab"

# Data volume levels (number of quarters to test)
q_levels <- c(1, 4, 8, 16, 24, 32)

# 3. Run Benchmarks at Different Data Volumes ----------------------------
results <- bench::press(
  q_count = q_levels,
  {
    # Construct dynamic time range based on q_count
    total_years <- (q_count - 1) %/% 4
    years_vec   <- 2014 + (0:total_years)
    all_q       <- rep(c("q1", "q2", "q3", "q4"), length.out = q_count)
    all_y       <- rep(2014:(2014 + total_years), each = 4)[1:q_count]

    bench::mark(
      full_pipeline = {
        # Force garbage collection before each run for accurate memory measurement
        gc()

        # Step 1: Parse raw FAERS ASCII files
        raw <- faers(
          years        = all_y,
          quarters     = all_q,
          dir          = raw_data_dir,
          compress_dir = compress_dir
        )

        # Step 2: Standardize drug names using MedDRA
        std <- faers_standardize(raw, meddra_path)

        # Step 3: Remove duplicate reports
        dedup <- faers_dedup(std)

        # Step 4: Signal detection for target drug (e.g., nivolumab/pembrolizumab)
        target_ids <- faers_filter(
          dedup,
          .fn = function(x) {
            idx <- grepl(target_drug, x$drugname, ignore.case = TRUE)
            x[idx, primaryid]
          },
          .field = "drug"
        )
        sig_res <- faers_phv_signal(target_ids, .full = dedup)
        sig_res   # Return result (but not used, only for timing)
      },
      iterations = 3,          # Number of repetitions per q_count
      check      = FALSE,      # Do not check result equality (timing only)
      memory     = TRUE        # Record memory allocation
    )
  }
)

# 4. Save Raw Benchmark Results ------------------------------------------
saveRDS(results, "bench_results.rds")
# results <- readRDS("bench_results.rds")   # Uncomment to reload

# 5. Compute Linear Fit R² (for annotation in plot) -----------------------
model <- lm(as.numeric(median) ~ q_count, data = results)
r_squared <- round(summary(model)$r.squared, 4)

# 6. Prepare Data for Visualization ---------------------------------------
perf_data <- results %>%
  mutate(
    quarters   = as.numeric(as.character(q_count)),
    time_min   = as.numeric(median) / 60,          # Convert seconds to minutes
    memory_gb  = as.numeric(mem_alloc) / (1024^3), # Convert bytes to GB
    throughput = quarters / time_min                # Quarters processed per minute
  )

# 7. Create Academic Theme for Plots --------------------------------------
academic_theme <- theme_classic() +
  theme(
    text            = element_text(family = "serif"),
    axis.title      = element_text(face = "bold", size = 12),
    axis.text       = element_text(face = "bold", color = "black", size = 11),
    legend.position = "bottom",
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.margin     = margin(10, 20, 10, 10)
  )

# 8. Plot 1: Execution Time & Memory Allocation (Dual Y‑Axis) -------------
# Scaling factor for secondary axis (memory)
coeff <- max(perf_data$memory_gb) / max(perf_data$time_min)

p1 <- ggplot(perf_data, aes(x = quarters)) +
  # Time line (left axis)
  geom_line(aes(y = time_min, color = "Processing Time"), linewidth = 1.2) +
  geom_point(aes(y = time_min, color = "Processing Time"),
             size = 3, shape = 21, fill = "white", stroke = 1.5) +
  # Memory line (right axis, scaled)
  geom_line(aes(y = memory_gb / coeff, color = "Memory Usage"),
            linewidth = 1.2, linetype = "dashed") +
  geom_point(aes(y = memory_gb / coeff, color = "Memory Usage"),
             size = 3, shape = 21, fill = "white", stroke = 1.5) +
  # R² annotation
  annotate("label",
           x = 16,
           y = max(perf_data$time_min) * 0.85,
           label = paste0("bold(R)^bold('2') ~ bold('= ", r_squared, "')"),
           parse = TRUE,
           family = "serif",
           size = 3.88,
           color = "black",
           fill = "white",
           alpha = 0.8) +
  scale_x_continuous(breaks = c(1, 4, 8, 16, 24, 32)) +
  scale_y_continuous(
    name = "Execution Time (min)",
    sec.axis = sec_axis(~ . * coeff, name = "Memory Allocation (GB)")
  ) +
  scale_color_manual(values = c("Processing Time" = "#E41A1C",
                                "Memory Usage"    = "#377EB8")) +
  labs(
    x = "Data Volume (Quarters)",
    color = "Performance Metrics"
  ) +
  academic_theme

# 9. Plot 2: Throughput (Quarters per Minute) -----------------------------
avg_tp <- mean(perf_data$throughput)

p2 <- ggplot(perf_data, aes(x = quarters, y = throughput)) +
  geom_area(fill = "#4DAF4A", alpha = 0.1) +
  geom_line(color = "#4DAF4A", linewidth = 1.2) +
  geom_point(size = 4, shape = 21, fill = "white",
             color = "#4DAF4A", stroke = 1.5) +
  geom_hline(yintercept = avg_tp, linetype = "dotted",
             color = "grey20", linewidth = 0.8) +
  annotate("label",
           x = 16,
           y = avg_tp,
           label = paste0("Mean Throughput: ", round(avg_tp, 2), " Qtr/min"),
           family = "serif",
           fontface = "bold",
           color = "black",
           fill = "white",
           alpha = 0.8) +
  scale_x_continuous(breaks = c(1, 4, 8, 16, 24, 32)) +
  scale_y_continuous(limits = c(0, max(perf_data$throughput) * 1.5)) +
  labs(
    x = "Data Volume (Quarters)",
    y = "Throughput (Quarters/min)"
  ) +
  academic_theme

# 10. Display and Save Plots ----------------------------------------------
# Combine plots side by side (optional, if desired)
# combined <- p1 + p2 + plot_layout(ncol = 2)
# print(combined)

print(p1)
print(p2)

# Save individual plots (optional)
ggsave("scaling_time_memory.pdf", p1, width = 8, height = 6, device = cairo_pdf)
ggsave("scaling_throughput.pdf", p2, width = 8, height = 6, device = cairo_pdf)

# ==============================================================================
# End of Script
# ==============================================================================

