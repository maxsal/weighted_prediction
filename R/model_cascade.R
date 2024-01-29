# libraries, functions, and options --------------------------------------------
ms::libri(
    ms, data.table, MatchIt, glue, qs, cli, optparse, tidyverse, survey, logistf
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
    make_option("--outcome",
        type = "character", default = "CA_101.1",
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

# functions --------------------------------------------------------------------
try_glmnet <- function(
    x,
    y,
    weights = NULL,
    alpha = 1,
    lambda = NULL,
    alt_lambda = NULL,
    family = "binomial",
    maxit = 1e5,
    ...) {
            tryCatch(
                {
                    mod <- glmnet(
                        x         = x,
                        y         = y,
                        weights   = weights,
                        alpha     = alpha,
                        lambda    = lambda,
                        family    = family,
                        maxit = maxit,
                        ...
                    )
                    if (mod$dev.ratio < 0) {
                        warning("dev.ratio < 0")
                    } else {
                    return(mod)
                    }
                },
                error = function(e) {
                    message("Here's the original error message:")
                    message(conditionMessage(e))
                    return(NA)
                },
                warning = function(w) {
                    message("Here's the original warning message:")
                    message(conditionMessage(w))
                    tryCatch({
                        mod <- glmnet(
                            x         = x,
                            y         = y,
                            weights   = weights,
                            alpha     = alpha,
                            lambda    = alt_lambda,
                            family    = family,
                        maxit = maxit,
                            ...
                        )
                        if (mod$dev.ratio < 0) {
                            warning("dev.ratio < 0")
                        } else {
                        return(mod)
                        }
                    },
                    error = function(e) {
                        message("Here's the original error message:")
                        message(conditionMessage(e))
                        return(NA)
                    },
                    warning = function(w) {
                        message("Here's the original warning message:")
                        message(conditionMessage(w))
                        tryCatch({
                            tmp_x <- cbind(x, weight = weights)
                            mod <- glmnet(
                                x = tmp_x,
                                y = y,
                                alpha = alpha,
                                lambda = lambda,
                                family = family,
                                maxit = maxit,
                                penalty.factor = as.numeric(unlist(dimnames(tmp_x)) != "weight"),
                                ...
                        )
                        if (mod$dev.ratio < 0) {
                            warning("dev.ratio < 0")
                        } else {
                        return(mod)
                        }
                        },
                        error = function(e) {
                            message("Here's the original error message:")
                            message(conditionMessage(e))
                            return(NA)
                        },
                        warning = function(w) {
                            message("Here's the original warning message:")
                            message(conditionMessage(w))
                            tmp_x <- cbind(x, weight = weights)
                            mod <- glmnet(
                                x = tmp_x,
                                y = y,
                                alpha = alpha,
                                lambda = alt_lambda,
                                family = family,
                                maxit = maxit,
                                penalty.factor = as.numeric(unlist(dimnames(tmp_x)) != "weight"),
                                ...)
                            if (mod$dev.ratio < 0) {
                                warning("dev.ratio < 0")
                            } else {
                            return(mod)
                            }
                        })
                    })

                }
            )
                    }

more_than_one_unique <- function(dt, var_names) {
    # Check if var_names are in the columns of dt
    if (!all(var_names %in% names(dt))) {
        stop("Some variables are not in the data.table.")
    }

    # Apply the function to each variable and return those with more than one unique value
    return(var_names[sapply(var_names, function(v) length(unique(dt[[v]])) > 1)])
}

keep_top_phecodes <- function(vec) {
    # Identify 'integer' and 'decimal' versions separately
    has_decimal <- grepl("\\.", vec)
    vec_decimal <- vec[has_decimal]
    vec_integer <- vec[!has_decimal]

    # Keep decimal version only when the integer version is missing
    vec_decimal_to_keep <- sapply(vec_decimal, function(x) {
        prefix <- sub("_.*", "", x)
        integer_version <- paste0(prefix, "_", floor(as.numeric(sub(".*_", "", x))))
        if (!(integer_version %in% vec_integer)) {
            return(x)
        }
    })

    # Combine integers and decimals that should be kept
    vec_final <- c(vec_integer, vec_decimal_to_keep[!is.null(vec_decimal_to_keep)])

    # Return cleaned vector
    return(as.character(unlist(vec_final)))
}


# read data --------------------------------------------------------------------
## mgi
mgi_tr_pims <- map(
    seq_along(time_thresholds),
    \(x) {
        glue(
            "data/private/mgi/{opt$mgi_version}/{opt$outcome}/",
            "time_restricted_phenomes/mgi_{opt$mgi_version}_{opt$outcome}_t",
            "{time_thresholds[x]}_pim_r{opt$matching_ratio}.qs"
        ) |>
            qread()
    }
) |> set_names(glue("t{time_thresholds}_threshold"))

mgi_covariates <- qread(glue(
    "data/private/mgi/{opt$mgi_version}/{opt$outcome}/",
    "mgi_{opt$mgi_version}_{opt$outcome}_match_data_r{opt$matching_ratio}.qs"
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
        merge_list(list(mgi_tr_pims[[x]], mgi_covariates[, !c("case")], mgi_weights[, ..weight_id]), by_var = "id", join_fn = dplyr::left_join)[, `:=`(
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

for (i in seq_along(time_thresholds)) {
    cli_progress_step(glue("fitting phers for {opt$outcome} at t = {time_thresholds[i]}"))
    data <- mgi_tr_merged[[i]][group == "test", ]

    # prepare covariates
    covariates <- c("age_at_threshold", "female", "nhw")

    risk_factor_table <- fread("data/public/dig_can_risk_factors.csv") # add to github
    risk_factors <- risk_factor_table[outcome_phecode == outcome_phecode, unique(risk_factor_variable)]
    risk_factors[risk_factors == "alcohol_ever"] <- "drinker"
    risk_factors[risk_factors == "smoke_ever"] <- "smoker"
    risk_factors <- unique(risk_factors[risk_factors %in% names(data)])

    symptoms_table <- fread("data/public/dig_can_symptoms.csv") # add to github
    symptoms <- symptoms_table[outcome_phecode == outcome_phecode, unique(symptom_phecode)]
    symptoms <- symptoms[symptoms != ""]
    symptoms <- unique(symptoms[symptoms %in% names(data)])

    ### THESE ARE ONCE PER OUTCOME
    # covariates (non-modifiable) ONCE PER OUTCOME
    covariates <- more_than_one_unique(data, covariates)
    cov_f <- paste0(covariates, collapse = " + ")
    ## unweighted
    cov_un <- logistf::logistf(
        formula = paste0(outcome, " ~ ", cov_f),
        data = data,
        control = logistf.control(maxit = 100, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 10000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)
    ## weighted
    cov_w <- logistf::logistf(
        formula = paste0(outcome, " ~ ", cov_f),
        data = data[!is.na(get(weight_var)), ],
        weights = data[!is.na(get(weight_var)), ][[weight_var]],
        control = logistf.control(maxit = 100, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 10000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)

    # risk factors (modifiable) ONCE PER OUTCOME
    risk_factors <- more_than_one_unique(data, risk_factors) |>
        keep_top_phecodes()
    risk_f <- paste0(risk_factors, collapse = " + ")
    ## unweighted
    risk_un <- logistf::logistf(
        formula = paste0(outcome, " ~ ", risk_f),
        data = data,
        control = logistf.control(maxit = 100000, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 100000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)
    ## weighted
    risk_w <- logistf::logistf(
        formula = paste0(outcome, " ~ ", risk_f),
        data = data[!is.na(get(weight_var)), ],
        weights = data[!is.na(get(weight_var)), ][[weight_var]],
        control = logistf.control(maxit = 100, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 10000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)

    # symptoms ONCE PER OUTCOME
    symptoms <- more_than_one_unique(data, symptoms) |>
        keep_top_phecodes()
    symptoms_f <- paste0(symptoms, collapse = " + ")
    ## unweighted
    symptoms_un <- logistf::logistf(
        formula = paste0(outcome, " ~ ", symptoms_f),
        data = data,
        control = logistf.control(maxit = 100, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 10000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)
    ## weighted
    symptoms_w <- logistf::logistf(
        formula = paste0(outcome, " ~ ", symptoms_f),
        data = data[!is.na(get(weight_var)), ],
        weights = data[!is.na(get(weight_var)), ][[weight_var]],
        control = logistf.control(maxit = 100, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 10000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)

    # covariates, risk factors, symptoms ONCE PER OUTCOME
    crs_f <- paste0(keep_top_phecodes(unique(c(covariates, risk_factors, symptoms))), collapse = " + ")
    ## unweighted
    crs_un <- logistf::logistf(
        formula = paste0(outcome, " ~ ", crs_f),
        data = data,
        control = logistf.control(maxit = 100, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 10000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)
    ## weighted
    crs_w <- logistf::logistf(
        formula = paste0(outcome, " ~ ", crs_f),
        data = data[!is.na(get(weight_var)), ],
        weights = data[!is.na(get(weight_var)), ][[weight_var]],
        control = logistf.control(maxit = 100, maxstep = 0.5),
        plcontrol = logistf.control(maxit = 10000, maxstep = 0.5)
    ) |> aimTwo::betas_from_mod(intercept = TRUE)

    cascade_models <- list(
        covariates_unweighted = cov_un,
        covariates_weighted = cov_w,
        risk_factors_unweighted = risk_un,
        risk_factors_weighted = risk_w,
        symptoms_unweighted = symptoms_un,
        symptoms_weighted = symptoms_w,
        covariates_risk_factors_symptoms_unweighted = crs_un,
        covariates_risk_factors_symptoms_weighted = crs_w
    )

    qsave(
        cascade_models,
        file = glue(
            "results/mgi/{opt$mgi_version}/{opt$outcome}/",
            "mgi_{opt$mgi_version}_{opt$outcome}_t{time_thresholds[i]}_cascade_models.qs"
        )
    )

    ### THESE ARE SEVERAL PER OUTCOME
    # phers SEVERAL PER OUTCOME
    ex_vars <- names(mgi_tr_merged[[i]])[names(mgi_tr_merged[[i]]) %in% ms::pheinfox[["phecode"]]]
    # hyperparameters from mgi -----------------------------------------------------
    hyperparameters <- fread(paste0(
        "results/mgi/", opt$mgi_version, "/",
        opt$outcome, "/mgi_", opt$mgi_version, "_", opt$outcome,
        "_t", time_thresholds[i], "_train_hyperparameters_ip_selection.csv"
    ))

        for (j in paste0("ridge_", c("alpha", "lambda.min", "lambda.1se"))) {
            if (length(which(hyperparameters[, parameter] == j)) == 1) next
            hyperparameters[
                which(hyperparameters[, parameter] == j)[2],
                parameter := paste0("w", parameter)
            ]
        }

    for (j in paste0("lasso_", c("alpha", "lambda.min", "lambda.1se"))) {
        if (length(which(hyperparameters[, parameter] == j)) == 1) next
        hyperparameters[
            which(hyperparameters[, parameter] == j)[2],
            parameter := paste0("w", parameter)
        ]
    }

    for (j in paste0("enet_", c("alpha", "lambda.min", "lambda.1se"))) {
        if (length(which(hyperparameters[, parameter] == j)) == 1) next
        hyperparameters[
            which(hyperparameters[, parameter] == j)[2],
            parameter := paste0("w", parameter)
        ]
    }

    ## unweighted
    # ridge
    ridge_un <- try_glmnet(
        x = as.matrix(data[, ..ex_vars]),
        y = data[, case],
        alpha = 0,
        lambda = hyperparameters[parameter == "ridge_lambda.min", value],
        alt_lambda = hyperparameters[parameter == "ridge_lambda.1se", value],
        family = "binomial"
    )

    # lasso
    lasso_un <- try_glmnet(
        x = as.matrix(data[, ..ex_vars]),
        y = data[, case],
        alpha = 1,
        lambda = hyperparameters[parameter == "lasso_lambda.min", value],
        alt_lambda = hyperparameters[parameter == "lasso_lambda.1se", value],
        family = "binomial"
    )

    # enet
    enet_un <- try_glmnet(
        x = as.matrix(data[, ..ex_vars]),
        y = data[, case],
        alpha = hyperparameters[parameter == "enet_alpha", value],
        lambda = hyperparameters[parameter == "enet_lambda.min", value],
        alt_lambda = hyperparameters[parameter == "enet_lambda.1se", value],
        family = "binomial"
    )

    # rf
    f <- as.formula(paste0("case ~ ", paste0(ex_vars, collapse = " + ")))
    rf_un <- ranger(
        formula       = f,
        data          = data |> select(any_of(c("case", ex_vars))),
        num.trees     = 500,
        mtry          = hyperparameters[parameter == "rf.mtry", value],
        min.node.size = hyperparameters[parameter == "rf.node_size", value],
        num.threads   = 12,
        importance    = "permutation",
        write.forest  = TRUE
    )

    ## weighted
    # ridge
    ridge_w <- try_glmnet(
        x = as.matrix(data[!is.na(get(weight_var)), ..ex_vars]),
        y = data[!is.na(get(weight_var)), case],
        alpha = 0,
        weights = data[!is.na(get(weight_var)), get(weight_var)],
        lambda = hyperparameters[parameter == "wridge_lambda.min", value],
        alt_lambda = hyperparameters[parameter == "wridge_lambda.1se", value],
        family = "binomial"
    )

    # lasso
    lasso_w <- try_glmnet(
        x = as.matrix(data[!is.na(get(weight_var)), ..ex_vars]),
        y = data[!is.na(get(weight_var)), case],
        alpha = 1,
        weights = data[!is.na(get(weight_var)), get(weight_var)],
        lambda = hyperparameters[parameter == "wlasso_lambda.min", value],
        alt_lambda = hyperparameters[parameter == "wlasso_lambda.1se", value],
        family = "binomial"
    )

    # enet
    enet_w <- try_glmnet(
        x = as.matrix(data[!is.na(get(weight_var)), ..ex_vars]),
        y = data[!is.na(get(weight_var)), case],
        alpha = hyperparameters[parameter == "wenet_alpha", value],
        weights = data[!is.na(get(weight_var)), get(weight_var)],
        lambda = hyperparameters[parameter == "wenet_lambda.1se", value],
        alt_lambda = hyperparameters[parameter == "wenet_lambda.1se", value],
        weight_as_pred = TRUE,
        family = "binomial"
    )

    # rf
    rf_w <- ranger(
        formula       = f,
        data          = data[!is.na(get(weight_var)), ] |> select(any_of(c("case", ex_vars))),
        num.trees     = 500,
        mtry          = hyperparameters[parameter == "wrf_mtry", value],
        min.node.size = hyperparameters[parameter == "wrf_min.node.size", value],
        num.threads   = 12,
        importance    = "permutation",
        write.forest  = TRUE,
        case.weights       = data[!is.na(get(weight_var)), get(weight_var)]
    )


    # ### for EACH phers, also obtain fitted (standardized) predictions for use in
    # ### following section
    # cbind(data, rf_w_pred = rf_w_pred)
    # phers_name <- "rf_w_pred"
    # crsp_f <- paste0(unique(c(cov_f, risk_f, symptoms_f, phers_name)), collapse = " + ")
    # ### need to add models that incorporate PheRSs
    # full_mod <- logistf::logistf(
    #     formula = paste0(outcome, " ~ ", crsp_f),
    #     data = cbind(data, rf_w_pred = rf_w_pred),
    #     control = logistf.control(maxit = 200, maxstep = 0.5)
    # )

    # covariates, risk factors, symptoms, phers SEVERAL PER OUTCOME
    # crsp_f <- paste0(crs_f, " + ", phers_f)
    # must predict phers_f and add as covariate
    ## unweighted
    ## weighted

    ## save models
    phers_models <- list(
        unweighted_ridge = ridge_un,
        unweighted_lasso = lasso_un,
        unweighted_enet = enet_un,
        unweighted_rf = rf_un,
        weighted_ridge = ridge_w,
        weighted_lasso = lasso_w,
        weighted_enet = enet_w,
        weighted_rf = rf_w
    )

    qsave(
        phers_models,
        file = glue(
            "results/mgi/{opt$mgi_version}/{opt$outcome}/",
            "mgi_{opt$mgi_version}_{opt$outcome}_t{time_thresholds[i]}_phers_models.qs"
        )
    )

}

cli_alert_success("Done!")
