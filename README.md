<p align="center">
  <img src="./logo/mlearner_logo.png" alt="mlearner Logo" width="500px" />
</p>

# 🧬 Polygenic prediction of treatment efficacy

``M-Learner`` is a AI/ML framework for detecting and interpreting genetically-driven **heterogeneous treatment effects (HTE)** in randomized controlled trials (RCTs). 

It introduces the concept of **polygenic efficacy**, defined by the aggregate contribution of genome-wide genetic variants to treatment response.

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
  --significance_level 0.05 \
  --output_path /tmp/ 
```

## 🛠️ Installation & Setup

### Requirements
- **R ≥ 4.2** (tested with R 4.2.1)

### Install Required Packages

Run the following command in R to install all required dependencies:

```r
install.packages(c("optparse", "mlr3", "mlr3learners", "sandwich", "lmtest", "ggplot2", "scales", "data.table", "e1071"))
```

## 📋 Input Format

All inputs must be plain `.txt` files:

- **`Y.txt`** — An Nx1 treatment outcome vector (numeric). One value per line with header.
Example:
```
Y
-0.783
1.266
0.573
0.905
0.735
-0.616
1.038
0.0571
-0.591
```
- **`D.txt`** — An Nx1 binary treatment indicator (0/1). One number per line with header.
Example:
```
treatment
0
1
0
0
1
0
0
0
1
```
- **`prs.txt`** — An NxP PRS matrix, where P is the number of PRS. Must include a header row with PRS names.
Example:
```
PRS1	PRS2	PRS3	PRS4	PRS5	PRS6
-0.560	-0.789	0.020	-1.042	0.784	0.620
-0.230	-0.502	0.314	-1.728	2.300	-0.758
1.559	1.496	1.328	0.642	0.157	0.852
0.071	-1.137	0.121	-1.529	0.047	-0.748
0.129	-0.179	0.713	0.002	0.097	0.630
1.715	1.902	0.779	0.250	0.070	1.097
0.461	-0.101	0.915	0.564	-1.848	-0.988
-1.265	-1.360	-0.574	0.189	-1.671	1.108
-0.687	-0.665	1.627	-0.733	-0.078	-0.490
```

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
  [--fis_file <file>] \
  --output_path <value> 
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
- `--output_path` — Path prefix to save PESxT, subgroup ATE, and feature importance results in CSV format

### Optional Arguments

- `--plot_file` — Output path for subgroup treatment effect plot (PNG format)
- `--fis_file` — Output path for Feature Importance Scores (text file)

## 📊 Output

- **Console output** — Prints tables (including estimate, confidence interval and p-value) for estimated PES-by-treatment interaction effect and average treatment effect in each subgroup
```
===== PESxT interaction results =====
        Estimate  CB lower CB upper     P value
beta.2 0.8584659 0.2170776 1.499854 0.008707927

===== ATE in each subgroup =====
        Estimate  CB lower CB upper      P value
gamma.1 0.962239 0.3056268 1.618851 4.075618e-03
gamma.2 1.157403 0.5087509 1.806056 4.701655e-04
gamma.3 2.199665 1.5694140 2.829917 7.889066e-12

CSV files saved:
 - /tmp/PESxT.csv
 - /tmp/ATE_in_each_group.csv
 - /tmp/feature_importance.csv
```
### Interpretation
* We observed a significant PES×T interaction (beta = 0.86, 95% CB: 0.22–1.50, P = 8.7×10⁻³), indicating that genetic variation captured by PES meaningfully modifies treatment response. Individuals with higher PES values show greater benefit from treatment (assuming that a higher treatment response value corresponds to a better clinical outcome).

* Subgroup analyses further support this pattern. The ATE increases across PES-defined strata:
   * Low PES group: ATE = 0.96 (95% CB: 0.31–1.62, P = 4.1×10⁻³)
   * Medium PES group: ATE = 1.16 (95% CB: 0.51–1.81, P = 4.7×10⁻⁴)
   * High PES group: ATE = 2.20 (95% CB: 1.57–2.83, P = 7.9×10⁻¹²)

Treatment efficacy is consistently positive in every group, but the magnitude is larger in high PES group. This provides evidence that the PES can stratify patients by expected treatment benefit and potentially guide more personalized therapeutic decisions.

## 📚 Reference
If you use ``M-Learner``, please cite
```
@article{miao2025polygenic,
  title={Polygenic prediction of treatment efficacy with causal transfer learning},
  author={Miao, Jiacheng and Mu, Jin and Yang, Xiaoyu and Fletcher, Jason and Schmitz, Lauren L and Lu, Qiongshi},
  journal={medRxiv},
  pages={2025--10},
  year={2025},
  publisher={Cold Spring Harbor Laboratory Press}
}
```
