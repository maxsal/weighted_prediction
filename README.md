# **Working title:** Weighty predictions: do selection weights improve generalizability of risk prediction models from non-probabilistic EHR samples?

## Objective

* Should researchers using EHR data from non-probabilistic cohorts consider selection weights to develop more generalizable risk prediction models?
    1.  Develop several one- and two-step risk scores selecting features and tuning hyperparameters in a discovery cohort (MGI) and estimating feature weights in a validation cohort (AOU). 
    2.  Risk scores are estimated in an assessment cohort (UKB) and assessed for several performance metrics and risk stratification capacity. 
    3.  Switch the roles of MGI and AOU and assess the performance of these models to determine whether prediction models built using selection weights are more robust to selection of discovery and validation cohorts.
    4.  Discuss the results and make recommendations regarding the use of selection weights in the development of EHR-based risk scores.

## Setup
Make sure you have the `maxsal/ms` and `maxsal/aimTwo` R packages installed

```
# install.packages("remotes")
remotes::install_github("maxsal/ms")
remotes::install_github("maxsal/aimTwo")
```