# faers-reproducibility: Analysis Scripts for the faers Package Study

[cite_start]This repository contains the source code, R scripts, and analysis pipelines required to reproduce the findings presented in the manuscript: **"faers: An Integrated R Package for End-to-End FDA Adverse Event Analysis"**. [cite: 1]

## 📋 Overview

[cite_start]The scripts provided here facilitate the complete reproduction of the performance benchmarks and clinical case studies discussed in the paper. [cite: 21, 154, 179] [cite_start]The core functional logic for data processing is powered by the [faers](https://github.com/WangLabCSU/faers) R package (version 1.6.0). [cite: 48]

## 📂 Repository Structure

-   [cite_start]`project_1.R`: Reproduction of the **Anti-PD-1/PD-L1-Associated Cardiotoxicity** study, covering data from Q1 2014 to Q4 2022. [cite: 180, 182]
-   [cite_start]`project_2.R`: Reproduction of the **CAR-T Cell Therapy-Associated Secondary Primary Malignancies (SPMs)** study and antibiotic interactions. [cite: 193, 194]
-   [cite_start]`project_3.R`: Exploratory **Subgroup Interaction Analysis** of Sex and Age on irAE reporting risk (2015–2025). [cite: 206, 208]
-   [cite_start]`test_1.R` & `test_2.R`: Scripts for **Performance Benchmarking** and **Scalability Testing** (processing up to 32 quarters of data). [cite: 154, 164, 165]
-   `faers_paper.Rproj`: RStudio project file for standardized environment configuration.

## 🛠️ Requirements

### Software

-   **R** (version $\ge$ 4.0.0)
-   [cite_start]**faers** R package (v1.6.0) [cite: 289]
-   [cite_start]Other dependencies: `data.table`, `ggplot2`, `patchwork`. [cite: 239]

### Data Access

[cite_start]Raw FAERS data should be obtained from the [FDA Quarterly Data Files](https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html). [cite: 290] [cite_start]The scripts utilize the automated pipeline within the `faers` package to handle data acquisition and parsing. [cite: 17, 72]

## 🚀 Reproduction Steps

1.  **Clone the repository**: \`\`\`bash git clone [git\@github.com](mailto:git@github.com){.email}:MadDERt/faers-reproducibility.git cd faers-reproducibility
