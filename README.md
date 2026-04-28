<p align="center">
  <img src="logo.png" width="250"/>
</p>

---

# README: ACORN-HAI Analysis Guide

## Introduction

The ACORN-HAI study, a prospective cohort conducted from September 2022 to December 2024, aims to establish a large-scale, multi-center, patient-centered surveillance network focusing on antimicrobial resistance in severe healthcare-associated infections. It also lays the groundwork for future interventional clinical trials targeting multidrug-resistant infections by building microbiology laboratory capacity and developing robust data collection and sharing platforms.

This README provides a step-by-step guide for conducting the analysis of the ACORN-HAI cohort using R. It covers key aspects such as baseline characteristics, antibiotic resistance, clinical outcomes, and antibiotic prescriptions, particularly highlighting **Carbapenem-resistant *Acinetobacter* (CRA)**, **Third-generation cephalosporin-resistant *Enterobacterales* (3GCRE)**, and **Carbapenem-resistant *Enterobacterales* (CRE)**. Throughout the guide, you will find explanations, code examples, and the implications of each component.

---

## Step-by-Step Guide

### Preparation

#### Step 1: Install R and RStudio
Make sure **R** and **RStudio** are installed on your computer:
- Download R from: [https://cran.r-project.org/](https://cran.r-project.org/)
- Download RStudio from: [https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

#### Step 2: Download the raw data
Request for raw REDCap data files from **Sujie** (jiesu@nus.edu.sg).

**Important:** Do not modify the names of the raw data files.

#### Step 3: Download and extract the code
Go to the `<> Code` section and download the ZIP folder containing the scripts. Once downloaded, unzip the folder to access the necessary code files.

---

### Cleaning raw data

#### Step 1: Set up `data/raw_data/` and place the raw data
**Create `raw_data`** folder inside the `data` directory, then place the raw data files (F01, F02, F03, F04, F07a, F07b, F07c, F07d, F07e, F07m) in the folder.

**Important:** Renaming the raw data files will change their order when imported into R, as they are automatically labeled `data[[1]]`, `data[[2]]`, etc. This can disrupt the data cleaning process. **Do not rename the raw data files.**

#### Step 2: Set up `data/clean_data_excel/` and `data/clean_data_RData/`

1. **Create directories**:
   - Inside `data/`, add `clean_data_excel` for Excel files and `clean_data_RData` for RData files.

2. **Directory structure**:
   ```
   data/
   ├── clean_data_excel/
   └── clean_data_RData/
   ```
#### Step 3: Run the data cleaning script
Open the `clean_data.Rmd` file in RStudio, and click **Run All** to process and clean the raw data for analysis.

#### Step 4: Output files
After cleaning, the following nine Excel files will be available in the `data/clean_data_excel/` folder:
  - `infection_types_index`: Infection types for each patient.
  - `baseline_outcomes_index`: Baseline and outcome-related variables.
  - `ast_all`: Antimicrobial susceptibility test (AST) results for all episodes.
  - `ast_all_index`: AST results for the index episodes.
  - `anti_treat_index`: Antibiotic usage for the index episodes.
  - `all_vap_bsi`: Relevant variables for all episodes.
  - `vap_bsi_index`: Relevant variables for the index episodes.
  - `df_ast`; `each_ast`: Prepare data for AST visualizations.

For details on specific variables, refer to the [data directory](https://docs.google.com/spreadsheets/d/1qLqACtCwm7IUfF0Fh_TJnrfE94kV-5Dq_Cn5IjIzS9c/edit?gid=766714505#gid=766714505).

**Note:** The cleaned data files are ready for analysis in SPSS, STATA, SAS, R, or other statistical software.

#### Step 5: Derive resistance variables
Run `class_organism.R` to derive resistance indicators from AST data:
- Antimicrobial resistance (AMR): resistance to ≥1 tested antimicrobial agent; intrinsic resistance excluded.
- Multidrug resistance (MDR): resistance to ≥1 agent in ≥3 antibiotic classes; classified using CMI 2012.

The updated dataset will be saved to:
- `data/clean_data_RData/ast_all_index.RData`
- `data/clean_data_excel/ast_all_index.xlsx`
---

### Demographic characteristics and antibiotic resistance profiles

#### Preparing data for visualization  
To prepare the data for plotting, run the following scripts in your R environment:
  - `descriptive_analysis/data_for_plot_1.R`
  - `descriptive_analysis/data_for_plot_2.R`

Each script will generate the data for plotting, with the output saved in `data/clean_data_RData/`.

#### Baseline characteristics 
Run `descriptive_analysis/table_baseline.R` script to generate a baseline characteristics table.

#### Proportion of infection types
Run the `descriptive_analysis/proportion_infection_types.R` script to generate proportion of infection types across countries with total index episodes.

#### Stacked charts
Run the `descriptive_analysis/stacked_charts_ast.R` script to generate stacked charts showing the proportions of  AST results by antibiotic class for the index episodes.

#### Pie charts
Run the following scripts to create pie charts displaying the proportions of AST results by antibiotics for the index episodes:
  - `descriptive_analysis/pie_charts_ast_VAP.R` for VAP.
  - `descriptive_analysis/pie_charts_ast_BSI_hosp.R` for hospital-acquired BSI.
  - `descriptive_analysis/pie_charts_ast_BSI_health.R` for healthcare-associated BSI.

#### Heatmap
Run the `descriptive_analysis/heatmap_ast.R` script to generate a heatmap of resistant organism proportions for the index episodes.

#### Antibiotic resistance profiles
Run the `descriptive_analysis/amr_profiles.R` script to visualize antibiotic resistance profiles across different infection types.

#### Prescriptions 
Run the `descriptive_analysis/sankey.R` script to illustrate the transition from empirical to definitive antibiotic prescriptions.

**Note:** Tables are saved in `descriptive_analysis/output/table/`, and figures in `descriptive_analysis/output/figure/`.

---

### Incidence and clinical outcomes

#### Preparing data for analysis
1. **Create `data/` folder** in the following directories:
   - `incidence/`
   - `outcome_attributable_mortality/`
   - `outcome_attributable_fbis/`
   - `outcome_attributable_eq5d/`

2. **Copy and paste**:
   - Copy the entire `data/clean_data_RData/` folder and paste it inside each of the directories listed above.

#### Incidence
- Open the `incidence` folder.
- Run the R scripts step by step.
- Tables and figures will be saved in `incidence/output/table/` and `incidence/output/figure/` directories.

#### Attributable mortality
- Open the `outcome_attributable_mortality` folder.
- Run the R scripts step by step.
- Tables and figures will be saved in `outcome_attributable_mortality/output/table/` and `outcome_attributable_mortality/output/figure/` directories.

#### Attributable FBIS
- Open the `outcome_attributable_fbis` folder.
- Run the R scripts step by step.
- Tables and figures will be saved in `outcome_attributable_fbis/output/table/` and `outcome_attributable_fbis/output/figure/` directories.

#### Attributable EQ-5D-3L
- Open the `outcome_attributable_eq5d` folder.
- Run the R scripts step by step.
- Tables and figures will be saved in `outcome_attributable_eq5d/output/table/` and `outcome_attributable_eq5d/output/figure/` directories.

---

### Troubleshooting
For any issues with code execution, please contact Xinxin (xx_hao@nus.edu.sg).

--- 


