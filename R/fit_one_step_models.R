### SCRIPT NEEDS TO BE GENERALIZED
###  [ ] add in options for mgi version, outcome, time threshold, etc.
###  [ ] add in looping to handle multiple time thresholds

# libraries
ms::libri(ms, data.table, glue, cli, qs, ranger, wlasso, doParallel,
          parallelly, parallel, optparse, glmnet, tidyverse, aimTwo)

set.seed(61787)

walk(list.files("./fn/", full.names = TRUE), source)

# optparse list ----------------------------------------------------------------
option_list <- list(
    make_option("--outcome",
        type = "character", default = "CA_101.41",
        help = "Outcome phecode [default = %default]"
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
    make_option("--weights",
        type = "character", default = "ip_selection",
        help = glue(
            "Weighting variable to use for weighted analyses - ",
            "selection, all, or list of named weight variables [default = %default]"
        )
    ),
    make_option("--matching_ratio",
        type = "numeric", default = "2",
        help = "Number of non-cases to match per case [default = %default]"
    ),
    make_option("--parallel",
        type = "logical", default = "TRUE",
        help = glue(
            "Whether to parallelize [default = %default]"
        )
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

# functions --------------------------------------------------------------------
try_glmnet <- function(
    x,
    y,
    weights    = NULL,
    alpha      = 1,
    lambda     = NULL,
    alt_lambda = NULL,
    family     = "binomial",
    ...) {
    tryCatch(
        {
            glmnet(
                x         = x,
                y         = y,
                weights   = weights,
                alpha     = alpha,
                lambda    = lambda,
                family    = family,
                ...
            )
        },
        error = function(e) {
            message("Here's the original error message:")
            message(conditionMessage(e))
            return(NA)
        },
        warning = function(w) {
            message("Here's the original warning message:")
            message(conditionMessage(w))
            glmnet(
                x         = x,
                y         = y,
                weights   = weights,
                alpha     = alpha,
                lambda    = alt_lambda,
                family    = family,
                ...
            )
        }
    )
}

extract_rf_vip <- function(ranger_mod) {
    data.table(
        variable = names(ranger_mod$variable.importance),
        importance = ranger_mod$variable.importance
    ) |> arrange(desc(importance))
}

# read data --------------------------------------------------------------------
mgi_tr_pims <- map(
  seq_along(time_thresholds),
  \(i) {
    qread(glue(
      "data/private/mgi/{opt$mgi_version}/{opt$outcome}/",
      "time_restricted_phenomes/mgi_{opt$mgi_version}_{opt$outcome}_t",
      "{time_thresholds[i]}_pim_r{opt$matching_ratio}.qs"
    ))
  }
) |>
  set_names(glue("t{time_thresholds}_threshold")) 

mgi_weights <- qread(
  glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs")
)

## phenome
pheinfo <- ms::pheinfox

# weights ----------------------------------------------------------------------
if (opt$weights == "all") {
    weight_vars <- names(mgi_weights)[!names(mgi_weights) %in% c("id", "DeID_PatientID")]
} else if (opt$weights == "selection") {
    weight_vars <- grep("selection_c", names(mgi_weights), value = TRUE)
} else {
    weight_vars <- unlist(strsplit(opt$weights, ","))
}

setnames(mgi_weights, opt$weights, "weight")
weight_id <- c("id", "weight")

mgi_tr_merged <- map(
  names(mgi_tr_pims),
  \(x) {
    list(
      mgi_tr_pims[[x]],
      mgi_weights[, ..weight_id]
    ) |>
    reduce(merge.data.table, by = "id", all.x = TRUE)
  }
) |>
  set_names(glue("t{time_thresholds}_threshold"))

## prep
for (i in seq_along(time_thresholds)) {
    cli_progress_step(paste0("Fitting models at t=", time_thresholds[i]))
    ex_vars <- names(mgi_tr_merged[[i]])[names(mgi_tr_merged[[i]]) %in% ms::pheinfox[["phecode"]]]
    # hyperparameters from mgi -----------------------------------------------------
    hyperparameters <- fread(paste0("results/mgi/", opt$mgi_version, "/",
        opt$outcome, "/mgi_", opt$mgi_version, "_", opt$outcome,
        "_t", time_thresholds[i], "_hyperparameters_ip_selection.csv"))

    # fit models -------------------------------------------------------------------
    # ridge
    ## unweighted
    unweighted_ridge_fit <- try_glmnet(
        x         = as.matrix(mgi_tr_merged[[i]][, ..ex_vars]),
        y         = mgi_tr_merged[[i]]$case,
        alpha     = 0,
        lambda    = hyperparameters[parameter == "ridge_lambda.min", value],
        alt_lambda    = hyperparameters[parameter == "ridge_lambda.1se", value],
        family    = "binomial"
    )
    ## weighted
    tmp_data <- mgi_tr_merged[[i]][!is.na(weight), ]
    weighted_ridge_fit <- try_glmnet(
        x         = as.matrix(tmp_data[, ..ex_vars]),
        y         = tmp_data$case,
        weights   = tmp_data$weight,
        alpha     = 0,
        lambda    = hyperparameters[parameter == "wwridge_lambda.min", value],
        alt_lambda    = hyperparameters[parameter == "wwridge_lambda.1se", value],
        family    = "binomial"
    )

    # lasso
    ## unweighted
    unweighted_lasso_fit <- try_glmnet(
        x         = as.matrix(mgi_tr_merged[[i]][, ..ex_vars]),
        y         = mgi_tr_merged[[i]]$case,
        alpha     = 1,
        lambda    = hyperparameters[parameter == "lasso_lambda.min", value],
        alt_lambda    = hyperparameters[parameter == "lasso_lambda.1se", value],
        family    = "binomial"
    )

    newx_vars <- rownames(unweighted_lasso_fit$beta)
    tmp_lpred <- predict(unweighted_lasso_fit,
        newx = as.matrix(data0(mgi_tr_merged[[i]], newx_vars)),
        type = "response"
    ) |> scale()

    mgi_tr_merged[[i]][, lpred := tmp_lpred]
    glm(case ~ lpred, data = mgi_tr_merged[[i]], family = "binomial")
    print(suppressMessages(pROC::roc(case ~ lpred, data = mgi_tr_merged[[i]])))

    ## weighted
    weighted_lasso_fit <- try_glmnet(
        x         = as.matrix(tmp_data[, ..ex_vars]),
        y         = tmp_data$case,
        weights   = tmp_data$weight,
        alpha     = 1,
        lambda    = hyperparameters[parameter == "wlasso_lambda.min", value],
        alt_lambda    = hyperparameters[parameter == "wlasso_lambda.1se", value],
        family    = "binomial"
    )

    # elastic net
    ## unweighted
    unweighted_enet_fit <- try_glmnet(
        x         = as.matrix(mgi_tr_merged[[i]][, ..ex_vars]),
        y         = mgi_tr_merged[[i]]$case,
        alpha     = hyperparameters[parameter == "enet_alpha", value],
        lambda    = hyperparameters[parameter == "enet_lambda.min", value],
        alt_lambda    = hyperparameters[parameter == "enet_lambda.1se", value],
        family    = "binomial"
    )

    ## weighted
    weighted_enet_fit <- try_glmnet(
        x          = as.matrix(tmp_data[, ..ex_vars]),
        y          = tmp_data$case,
        weights    = tmp_data$weight,
        alpha      = hyperparameters[parameter == "wwenet_alpha", value],
        lambda     = hyperparameters[parameter == "wwenet_lambda.min", value],
        alt_lambda = hyperparameters[parameter == "wwenet_lambda.1se", value],
        family     = "binomial"
    )

    # random forest
    ## unweighted
    f <- as.formula(paste0("case ~ ", paste0(ex_vars, collapse = " + ")))
    unweighted_rf_fit <- ranger(
        formula       = f,
        data          = mgi_tr_merged[[i]] |> select(-any_of(c("weight", "unique "))),
        num.trees     = 500,
        mtry          = hyperparameters[parameter == "rf.mtry", value],
        min.node.size = hyperparameters[parameter == "rf.node_size", value],
        num.threads   = 12,
        importance    = "permutation",
        write.forest  = TRUE
    )

    ## weighted
    ### FORTHCOMING
    # weighted_rf_fit <- ranger(
    #     formula       = f,
    #     data          = tmp_data |> select(-any_of(c("weight", "unique "))),
    #     num.trees     = 500,
    #     mtry          = hyperparameters[parameter == "rf.mtry", value],
    #     min.node.size = hyperparameters[parameter == "rf.node_size", value],
    #     case.weights  = tmp_data[, weight],
    #     num.threads   = 12,
    #     importance    = "permutation",
    #     write.forest  = TRUE
    # )

    # prepare output ---------------------------------------------------------------
    out_models <- list(
        "unweighted_ridge" = unweighted_ridge_fit,
        "weighted_ridge"   = weighted_ridge_fit,
        "unweighted_lasso" = unweighted_lasso_fit,
        "weighted_lasso"   = weighted_lasso_fit,
        "unweighted_enet"  = unweighted_enet_fit,
        "weighted_enet"    = weighted_enet_fit,
        "unweighted_rf"    = unweighted_rf_fit$forest
    )

    # save output ------------------------------------------------------------------
    qsave(out_models, glue("results/mgi/{opt$mgi_version}/{opt$outcome}/mgi_{opt$mgi_version}_{opt$outcome}_t{time_thresholds[i]}_models.qs"))
}
cli_progress_done()

cli::cli_alert_success("Done! 🎉")
