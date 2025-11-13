# EMAuxiliary Workflow: Preparing, Evaluating, and Integrating Auxiliary Variables in EMA Data

This repository accompanies the paper:

> **[Paper title, authors, year — to be updated when accepted]**  
> The workflow provides a three-step pipeline for preparing ecological momentary assessment (EMA) datasets, evaluating potential auxiliary variables, and incorporating them into multivariate multilevel models using the **EMAuxiliary** R function (which generates **Blimp** input code for auxiliary-variable modeling).

---

## 📘 Overview of the Three-Step Workflow

| Step | Purpose | Output / Next Input |
|------|----------|---------------------|
| **Step 1 – Prepare EMA Data** (`Step1_EMA_DataPrep.Rmd`) | Helps create a complete long-format EMA dataset where each scheduled prompt—completed or missed—is represented. Add lagged variables, within-person centering, and person-level summaries. | `ema_prepared.csv` |
| **Step 2 – Evaluate Auxiliary Variables** (`Step2_AuxEval.Rmd`) | Evaluate and rank potential auxiliary variables based on prediction of missingness and correlations with focal variables. | ranked table of auxiliary variables |
| **Step 3 – Fit Models with EMAuxiliary** (`Step3_ExamplesEMAuxiliary.Rmd`, `functions/EMAuxiliary.R`) | Incorporate selected auxiliaries into multilevel models in **Blimp**. Handles latent centering, variable typing, and Bayesian convergence diagnostics automatically. | full Blimp output and model estimates |

Each component is fully stand-alone and can be run independently, but the workflow is designed to flow from one step to the next.

---

## 🚀 Getting Started

### 1. Installation and Requirements
- **R (≥ 4.2)**  
- **Blimp**: [Download Blimp](https://www.appliedmissingdata.com/blimp)  
- **R interface for Blimp (rblimp)**:
  ```r
  remotes::install_github("blimp-stats/rblimp")
Required R packages:
dplyr, tidyr, lme4, broom.mixed, gt, DT, lubridate, ggplot2, plotly, data.table, purrr

### 2. Using the Tutorials
   
Step 1: Run relevant pieces of the notebook Step1_EMA_DataPrep.Rmd to reconstruct or simulate your EMA dataset and create ema_prepared.csv.

Step 2: Run Step2_AuxEval.Rmd using ema_prepared.csv as input. This will output ranked candidate auxiliaries.

Step 3: Use the EMAuxiliary() function (located in functions/EMAuxiliary.R) to generate a Blimp model that includes your selected auxiliaries. Run Step3_ExamplesEMAuxiliary.Rmd for tutorial.

## 📦 Repository Contents
EMAuxiliary-Workflow/

├── README.md

├── Step1_PrepareEMA.Rmd

├── Step2_EvaluateAux.Rmd

├── Step3_ExamplesEMAuxiliary.Rmd

├── functions/

├     └── EMAuxiliary.R

├     └── EMAuxiliary Reference Manual.pdf

    
## 📎 Citation and Archival Links
GitHub repository (latest version): https://github.com/schneids111/EMAuxiliary-Workflow

Permanent OSF snapshot (DOI): [to be added after upload]

Please cite both GitHub and OSF entries when referencing these materials.

## 🧠 Contact
For questions, please contact:
Stefan Schneider, PhD, University of Southern California
Email: schneids@usc.edu
