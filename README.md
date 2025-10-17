<p align="center">
  <img src="./logo/mlearner_logo.png" alt="mlearner Logo" width="500px" />
</p>

# 🧬 Polygenic prediction of treatment efficacy

M-Learner is a framework for detecting and interpreting genetically driven **heterogeneous treatment effects (HTE)** in randomized controlled trials (RCTs). It introduces the concept of **polygenic efficacy**, defined by the aggregate contribution of genome-wide genetic variants to treatment response.

**Manuscript**: https://jiachengmiao.com/paper/M_Learner.pdf

## 🚀 Quick Start

> **⚠️ Prerequisites**: Complete the [installation & setup](#️-installation--setup) below before running M-Learner.

**M-Learner-I** provides causal transfer learning with individual-level data. This repository implements the M-learner framework for estimating treatment effect heterogeneity (HTE) and supports:

- Propensity score estimation (`constant`, `lasso`, `random_forest`, `tree`, `svm`, `logit`)
- Proxy outcome models for baseline control adjustment (BCA) and conditional average treatment effect (CATE)
- Fold-splitting for bias-robust estimation
- Estimated PES-by-treatment interaction effect and average treatment effect in each subgroup
- FIS (Feature Importance Scores, per-covariate association with CATE)
- Optional subgroup effect plots

### Basic Example
Estimate with 5 folds, SVM learner, and 3 subgroups using the example data:

```bash
cd mlearner
Rscript M_Learner_I.R \
  --outcome example_data/Y.txt \
  --treatment example_data/D.txt \
  --prs example_data/prs.txt \
  --num_folds 5 \
  --seed 123 \
  --prop_learner constant \
  --base_learners svm \
  --quantile_cutoffs 0.3333,0.6667 \
  --significance_level 0.05
```

## 🛠️ Installation & Setup

### Requirements
- **R ≥ 4.2** (tested with R 4.2.1)

### Install Required Packages

Run the following command in R to install all required dependencies:

```r
install.packages(c("optparse", "mlr3", "mlr3learners", "sandwich", "lmtest", "ggplot2", "scales", "data.table"))
```

## 📋 Input Format

All inputs must be plain `.txt` files:

- **`Y.txt`** — Outcome vector (numeric). One value per line, no header.
- **`D.txt`** — Binary treatment indicator (0/1). One number per line, no header.
- **`Z.txt`** — PRS matrix. Must include a header row with PRS names.

## 💻 Usage

```bash
Rscript M_Learner_I.R \
  --outcome <file> \
  --treatment <file> \
  --prs <file> \
  --num_folds <int> \
  --seed <int> \
  --prop_learner <method> \
  --base_learners <method> \
  --quantile_cutoffs <values> \
  --significance_level <value> \
  [--plot_file <file>] \
  [--fis_file <file>]
```

### Required Arguments

- `--outcome` — Outcome file path (e.g., `Y.txt`)
- `--treatment` — Treatment indicator file path (e.g., `D.txt`)
- `--prs` — PRS/covariates file path (e.g., `Z.txt`)
- `--num_folds` — Number of folds for cross-fitting (e.g., `5`)
- `--seed` — Random seed for reproducibility (e.g., `123`)
- `--prop_learner` — Propensity score learner: `constant`, `lasso`, `random_forest`, `tree`, `svm`, or `logit`
- `--base_learners` — CATE learner(s), comma-separated (e.g., `svm` or `svm,ranger`)
- `--quantile_cutoffs` — Subgroup cutoffs, comma-separated decimals (e.g., `0.3333,0.6667` for 3 groups)
- `--significance_level` — Significance level for confidence intervals (e.g., `0.05`)

### Optional Arguments

- `--plot_file` — Output path for subgroup treatment effect plot (PNG format)
- `--fis_file` — Output path for Feature Importance Scores (text file)

## 📊 Output

- **Console output** — Prints tables (including estimate, confidence interval and p value) for estimated PES-by-treatment interaction effect and average treatment effect in each subgroup

## 📚 Reference
