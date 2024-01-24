# libraries and such
options(stringsAsFactors = FALSE)
set.seed(61787)

ms::libri(ms, data.table, qs, cli, glue, purrr, ggplot2, dplyr, optparse, aimTwo)

# load functions
walk(list.files("fn/", full.names = TRUE), source)

# optparse list ----------------------------------------------------------------
option_list <- list(
    make_option("--outcome",
        type = "character", default = "CA_101.1",
        help = "Outcome phecode [default = %default]"
    ),
    make_option("--aou_version",
        type = "character", default = "20230309",
        help = "Version of AOU data [default = %default]"
    ),
    make_option("--ukb_version",
        type = "character", default = "20221117",
        help = "Version of UKB data [default = %default]"
    ),
    make_option("--mgi_version",
        type = "character", default = "20220822",
        help = "Version of MGI data [default = %default]"
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
        help = "Number of non-cases to match per case [default = %default]"
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

ukb_version     <- opt$ukb_version
aou_version     <- opt$aou_version
mgi_version     <- opt$mgi_version
outcome         <- opt$outcome
matching_ratio  <- opt$matching_ratio
time_thresholds <- strsplit(opt$time_thresholds, ",")[[1]]

# load data
cli_progress_step("Loading data")
## aou
aou_models <- map(
    seq_along(time_thresholds),
    \(i) {
        qread(glue("results/aou/{aou_version}/{outcome}/aou_{aou_version}_{outcome}_t{time_thresholds[i]}_models.qs"))
    }
) |> set_names(glue("t{time_thresholds}_threshold"))
aou_results <- qread(glue("results/aou/{aou_version}/{outcome}/aou_{outcome}_r{matching_ratio}_results_list.qs"))

## ukb
ukb_tr_pims <- map(
    seq_along(time_thresholds),
    \(i) {
        qread(glue(
            "data/private/ukb/{ukb_version}/",
            "{outcome}/time_restricted_phenomes/ukb_",
            "{ukb_version}_{outcome}_",
            "t{time_thresholds[i]}_pim_",
            "r{matching_ratio}.qs"
        )) |>
        as.data.table()
    }
) |>
    set_names(glue("t{time_thresholds}_threshold"))

ukb_covariates <- qread(
    glue(
        "data/private/ukb/{ukb_version}/",
        "{outcome}/ukb_",
        "{ukb_version}_{outcome}_",
        "match_data_",
        "r{matching_ratio}.qs"
    )
)

ukb_tr_merged <- map(
    names(ukb_tr_pims),
    \(x) {
        merge.data.table(
            ukb_tr_pims[[x]],
            ukb_covariates[, !c("case", "age", "sex", "female", "length_followup")],
            by = "id",
            all.x = TRUE
        ) |> 
            dplyr::mutate(age_at_threshold = round(get(x) / 365.25, 1))
    }
) |>
    set_names(glue("t{time_thresholds}_threshold"))


# construct PheRS
cli_progress_step("Constructing PheRS")

## functions -------------------------------------------------------------------
calculate_phers <- function(
    pim,
    weight_data,
    id_var     = "id",
    trait_var  = "phecode",
    weight_var = "beta",
    name
) {
    codes <- 0
    if (nrow(weight_data) == 0) {
        cli_alert_warning("No results for {name}, skipping")
        return(data.table(phers = NA, name = name, codes = codes))
    }
    out <- data.table::copy(as.data.table(pim))[, pred := 0]
    for(i in weight_data[[trait_var]]) {
        if (!(i %in% names(out))) {
            cli_alert_warning("{i} not in phecode indicator matrix, skipping")
            next
        }
        codes <- codes + 1
        out[, pred := pred + (weight_data[get(trait_var) == i, get(weight_var)] * get(i))]
    }
    out[, phers := scale(pred)]
    out[, .(id, case, pred, phers, name = name, codes = codes)]
}

one_step_phers <- function(
    test_data,
    fitted_models,
    glmnet_models,
    rf_models,
    id_var      = "id",
    outcome_var = "case"
) {
    data <- data.table::copy(test_data)
    # glmnet
    # glmnet_models <- names(fitted_models)[sapply(fitted_models, function(x) "glmnet" %in% class(x))]
    for (i in seq_along(glmnet_models)) {
        if (fitted_models[[i]]$df == 0) {
            next
        }
        vars <- rownames(fitted_models[[i]]$beta)
        tmp2 <- predict(fitted_models[[i]], newx = as.matrix(data0(data, vars)))
        data[, paste0("phers_", glmnet_models[i]) := tmp2]
    }

    # rf
    for (i in seq_along(rf_models)) {
        vars <- fitted_models[[rf_models[i]]]$independent.variable.names
        data[, paste0("phers_", rf_models[i]) := predict(fitted_models[[rf_models[i]]], data = data0(data, vars))$predictions]
    }

    phers_vars    <- names(data)[startsWith(names(data), "phers_")]
    phers_id_vars <- c(id_var, outcome_var, phers_vars)

    data_phers_diagnostics <- map(
        phers_vars,
        \(x) {
            aimTwo::get_bin_diagnostics(data, outcome = "case", exposure = x)
        }
    ) |> set_names(phers_vars)

    return(
        list(
            data = data[, ..phers_id_vars],
            diagnostics = data_phers_diagnostics
        )
    )

}

###


## two-step PheRS --------------------------------------------------------------
ukb_phers_list <- list()
for (i in seq_along(time_thresholds)) {
    tmp_names <- names(aou_results)[startsWith(names(aou_results), glue("t{time_thresholds[i]}_"))]
    for (w in seq_along(tmp_names)) {
        ukb_phers_list[[length(ukb_phers_list) + 1]] <- calculate_phers(
            pim         = ukb_tr_merged[[glue("t{time_thresholds[i]}_threshold")]],
            weight_data = aou_results[[tmp_names[w]]],
            name = tmp_names[w]
        )
    }
}

phers_calculated <- map_dfr(
    ukb_phers_list,
    \(x) {
        data.table(
            name  = x[["name"]][1],
            phers = (\(y) {
                if (all(is.na(x[["phers"]])) | all(is.nan(x[["phers"]]))) {
                    0
                } else {
                    1
                }
            })()
        )
    }
)

diag_results <- map(
    which(phers_calculated[, phers == 1]),
    \(i) {
        aimTwo::get_bin_diagnostics(ukb_phers_list[[i]], outcome = "case", exposure = "phers")
    },
    .progress = TRUE
) |> purrr::set_names(phers_calculated[which(phers_calculated[, phers == 1]), name])

## one-step phers --------------------------------------------------------------
one_step_res <- map(
    seq_along(time_thresholds),
    \(i) {
        one_step_phers(
            test_data     = ukb_tr_merged[[glue("t{time_thresholds[i]}_threshold")]],
            fitted_models = aou_models[[glue("t{time_thresholds[i]}_threshold")]],
            glmnet_models = names(aou_models[[glue("t{time_thresholds[i]}_threshold")]])[sapply(aou_models[[glue("t{time_thresholds[i]}_threshold")]], function(x) "glmnet" %in% class(x))],
            rf_models     = "unweighted_rf"
        )
    },
    .progress = TRUE
) |> set_names(glue("t{time_thresholds}_threshold"))


# save results
cli_progress_step("Saving results")
if (!dir.exists(glue("results/ukb/{ukb_version}/{outcome}/phers/"))) {
    dir.create(glue("results/ukb/{ukb_version}/{outcome}/phers/"), recursive = TRUE)
}
qsave(
    ukb_phers_list,
    glue("results/ukb/{ukb_version}/{outcome}/phers/ukb_{outcome}_two_step_phers.qs")
)
walk(
    seq_along(time_thresholds),
    \(i) {
        qsave(
            one_step_res[[i]]$data,
            glue("results/ukb/{ukb_version}/{outcome}/phers/ukb_{outcome}_t{time_thresholds[i]}_one_step_phers.qs")
        )
    }
)

qsave(
    diag_results,
    glue("results/ukb/{ukb_version}/{outcome}/phers/ukb_{outcome}_two_step_phers_diag_results.qs")
)
walk(
    seq_along(time_thresholds),
    \(i) {
        qsave(
            one_step_res[[i]]$diagnostics,
            glue("results/ukb/{ukb_version}/{outcome}/phers/ukb_{outcome}_t{time_thresholds[i]}_one_step_phers_diag_results.qs")
        )
    }
)

cli_alert_success("Done! 🎉")
