# libraries, functions, and options --------------------------------------------
ms::libri(
    ms, data.table, MatchIt, glue, qs, cli, optparse, tidyverse, survey, logistf,
    brglm2
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
    make_option("--outcome",
        type = "character", default = "CA_101.6",
        help = "Outcome phecode [default = %default]"
    ),
    make_option("--mgi_version",
        type = "character", default = "20220822",
        help = "Version of MGI data [default = %default]"
    ),
    make_option("--mgi_cohort",
        type = "character", default = "comb",
        help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"
    ),
    make_option("--time_thresholds",
        type = "character", default = "0,1,2,5",
        help = glue(
            "Time thresholds for the phenome data ",
            "[default = %default]"
        )
    ),
    make_option("--matching_ratio",
        type = "numeric", default = "2",
        help = glue(
            "Number of controls per case [default = %default]"
        )
    ),
    make_option("--weights",
        type = "character", default = "ip_selection",
        help = glue(
            "Weighting variable to use for weighted analyses - ",
            "selection, all, or list of named weight variables [default = %default]"
        )
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])
cooccur_covs <- c("age_at_threshold", "female", "length_followup")

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version)

# read data --------------------------------------------------------------------
## mgi
mgi_tr_pims <- lapply(
    seq_along(time_thresholds),
    \(x) {
        glue(
            "data/private/mgi/{opt$mgi_version}/{opt$outcome}/",
            "time_restricted_phenomes/mgi_{opt$outcome}_t",
            "{time_thresholds[x]}_r{opt$matching_ratio}_{opt$mgi_version}.qs"
        ) |>
            read_qs()
    }
) |> set_names(glue("t{time_thresholds}_threshold"))

mgi_covariates <- read_qs(glue(
    "data/private/mgi/{opt$mgi_version}/{opt$outcome}/",
    "matched_covariates_r{opt$matching_ratio}.qs"
))

mgi_demo <- read_qs(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_{opt$mgi_cohort}.qs"))

mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_{opt$mgi_cohort}.qs"))

## phenome
pheinfo <- ms::pheinfox
if (pheinfo[phecode == opt$outcome, sex] != "Both") {
    if (any(c("sex", "female", "male") %in% tolower(cooccur_covs))) {
        cooccur_covs <- cooccur_covs[-which(tolower(cooccur_covs) %in% c("sex", "female", "male"))]
    }
}

# weights ----------------------------------------------------------------------
if (opt$weights == "all") {
    weight_vars <- names(mgi_weights)[!names(mgi_weights) %in% c("id", "DeID_PatientID")]
} else if (opt$weights == "selection") {
    weight_vars <- grep("selection_c", names(mgi_weights), value = TRUE)
} else {
    weight_vars <- unlist(strsplit(opt$weights, ","))
}

# merge data -------------------------------------------------------------------
weight_id <- c("id", weight_vars)

mgi_tr_merged <- lapply(
    names(mgi_tr_pims),
    \(x) {
        merge_list(list(mgi_tr_pims[[x]], mgi_covariates[, !c("case")], mgi_weights[, ..weight_id], mgi_demo[, .(id, race_eth, smoker, drinker, nhw)]), by_var = "id", join_fn = dplyr::left_join)[, `:=`(
            age_at_threshold = round(get(x) / 365.25, 1)
        )][group == "test", ]
    }
)
names(mgi_tr_merged) <- glue("t{time_thresholds}_threshold")


# data[, `:=` (
#     race_eth = relevel(factor(fcase(
#     race == "Caucasian" & ethn == "Non-Hispanic", "NHW",
#     race == "African America" & ethn == "Non-Hispanic", "NHB",
#     race == "Asian" & ethn == "Non-Hispanic", "NHA",
#     ethn == "Hispanic", "Hisp",
#     race %in% c("Native American") & ethn != "", "Other",
#     default = "Unknown"
#     )), ref = "NHW")
# )]

#### CASCADE
outcome_phecode <- opt$outcome
outcome <- "case"
weight_var <- "ip_selection"
data <- "" # contains id, outcome, weight, covariates, risk factors, symptoms

data <- mgi_tr_merged[[1]]

# prepare covariates
covariates <- c("age_at_threshold", "female", "nhw")

risk_factor_table <- fread("data/public/dig_can_risk_factors.csv") # add to github
risk_factors <- risk_factor_table[outcome_phecode == outcome_phecode, unique(risk_factor_variable)]
risk_factors[risk_factors == "alcohol_ever"] <- "drinker"
risk_factors[risk_factors == "smoke_ever"] <- "smoker"
risk_factors <- risk_factors[risk_factors %in% names(data)]

symptoms_table <- fread("data/public/dig_can_symptoms.csv") # add to github
symptoms <- symptoms_table[outcome_phecode == outcome_phecode, unique(symptom_phecode)]
symptoms <- symptoms[symptoms != ""]
symptoms <- symptoms[symptoms %in% names(data)]

### THESE ARE ONCE PER OUTCOME
# covariates (non-modifiable) ONCE PER OUTCOME
cov_f <- paste0(covariates, collapse = " + ")
## unweighted
cov_un <- glm(
    formula = paste0(outcome, " ~ ", cov_f),
    data = data,
    family = binomial(link = "logit"),
    method = "brglmFit"
)
## weighted
cov_dsn <- svydesign(
    id = ~ 1,
    weights = ~ get(weight_var),
    data = data[!is.na(get(weight_var)), ]
)
cov_w <- svyglm(
    formula = paste0(outcome, " ~ ", cov_f),
    family = "quasibinomial",
    design = cov_dsn
)

# risk factors (modifiable) ONCE PER OUTCOME
risk_f <- paste0(risk_factors, collapse = " + ")
## unweighted
risk_un <- glm(
    formula = paste0(outcome, " ~ ", risk_f),
    data = data,
    family = binomial(link = "logit"),
    method = "brglmFit"
)
## weighted
risk_dsn <- svydesign(
    id = ~ 1,
    weights = ~ get(weight_var),
    data = data[!is.na(get(weight_var)), ]
)
risk_w <- svyglm(
    formula = paste0(outcome, " ~ ", risk_f),
    family = "quasibinomial",
    design = risk_dsn
)

# symptoms ONCE PER OUTCOME
symptoms_f <- paste0(symptoms, collapse = " + ")
## unweighted
symptoms_un <- glm(
    formula = paste0(outcome, " ~ ", symptoms_f),
    data = data,
    family = binomial(link = "logit"),
    method = "brglmFit"
)
## weighted
symptoms_dsn <- svydesign(
    id = ~ 1,
    weights = ~ get(weight_var),
    data = data[!is.na(get(weight_var)), ]
)
symptoms_w <- svyglm(
    formula = paste0(outcome, " ~ ", symptoms_f),
    family = "quasibinomial",
    design = symptoms_dsn
)

# covariates, risk factors, symptoms ONCE PER OUTCOME
crs_f <- paste0(unique(c(cov_f, risk_f, symptoms_f)), collapse = " + ")
crs_f <- paste0(cov_f, " + ", risk_f, " + ", symptoms_f)
## unweighted
crs_un <- glm(
    formula = paste0(outcome, " ~ ", crs_f),
    data = data,
    family = binomial(link = "logit"),
    method = "brglmFit"
)
## weighted
crs_design <- survey::svydesign(
    id = ~ 1,
    weights = ~ get(weight_var),
    data = data[!is.na(get(weight_var)), ]
)
crs_w <- svyglm(
    formula = paste0(outcome, " ~ ", crs_f),
    design = crs_design,
    family = "quasibinomial"
)

tmp <- data.table(
    case = data[!is.na(get(weight_var)), case],
    pred = crs_w$fitted.values
)
pROC::ci.auc(pROC::roc(tmp$case, tmp$pred))

### THESE ARE SEVERAL PER OUTCOME
# phers SEVERAL PER OUTCOME
## unweighted
## weighted

### for EACH phers, also obtain fitted (standardized) predictions for use in
### following section

# covariates, risk factors, symptoms, phers SEVERAL PER OUTCOME
# crsp_f <- paste0(crs_f, " + ", phers_f)
# must predict phers_f and add as covariate
## unweighted
## weighted

