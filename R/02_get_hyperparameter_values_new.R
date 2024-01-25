# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds with weights
# author:   max salvatore

# libraries, functions, and options --------------------------------------------
ms::libri(
  ms, data.table, MatchIt, glue, qs, cli, optparse, tidyverse, maxsal/aimTwo,
  maxsal/wglmnet
)

set.seed(61787)

walk(list.files("fn/", full.names = TRUE), source)

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
    help = "Number of non-cases to
     match per case [default = %default]"
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
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version)

# read data --------------------------------------------------------------------
## mgi
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

for (i in seq_along(time_thresholds)) {
  for (w in seq_along(weight_vars)) {

  data = mgi_tr_merged[[i]][group == "train", ]
  outcome = "case"
  exposures = names(mgi_tr_merged[[i]])[names(mgi_tr_merged[[i]]) %in% ms::pheinfox[, phecode]]
  weight = weight_vars[w]
  parallel      = TRUE
  methods       = c("ridge", "lasso", "enet", "rf")
  n_folds       = 10
  n_cores       = 16
  alpha         = NULL
  mtry_seq      = NULL
  node_size_seq = NULL
  n_trees       = 500
  verbose       = TRUE
  return_mod    = TRUE


  if (parallel) {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }

  if (is.null(weight)) {
    dataset <- data.table::copy(data[stats::complete.cases(data), ])
    if (is.null(exposures)) {
      exposures <- names(dataset)[names(dataset) != outcome]
    }
  } else {
    use_these_vars <- names(data)[!names(data) %in% c("id", weight)]
    dataset <- data.table::copy(
      data[
        stats::complete.cases(data |> dplyr::select(tidyselect::any_of(use_these_vars))),
        ]
    )
    dataset[["id"]] <- NULL
    wdataset <- data.table::copy(dataset[stats::complete.cases(dataset), ])
    if (is.null(exposures)) {
      exposures <- names(dataset)[!names(dataset) %in% c(outcome, weight)]
    }
  }

  out <- list()

  if ("ridge" %in% methods) {
    ridge_mod <- tryCatch({
      tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      .alpha    = 0,
      parallel  = parallel,
      verbose   = verbose,
      return_mod = return_mod
    )
    },
    error = function(e) NULL)
    out[["ridge"]] <- ridge_mod
    if (!is.null(weight)) {
      wridge_mod <- tryCatch({tune_wglmnet(
        data           = wdataset,
        outcome        = outcome,
        exposures      = exposures,
        weight         = weight,
        n_folds        = n_folds,
        alpha          = 0,
        parallel       = parallel,
        verbose        = verbose,
      return_mod = return_mod
      ) },
      error = function(e) NULL)
      out[["wridge"]] <- wridge_mod
    }
  }

  if ("lasso" %in% methods) {
    lasso_mod <- tryCatch({tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      .alpha    = 1,
      parallel  = parallel,
      verbose   = verbose,
      return_mod = return_mod
    )},
    error = function(e) NULL)
    out[["lasso"]] <- lasso_mod
    if (!is.null(weight)) {
      wlasso_mod <- tryCatch({tune_wglmnet(
        data      = wdataset,
        outcome   = outcome,
        exposures = exposures,
        weight    = weight,
        alpha     = 1,
        n_folds   = n_folds,
        verbose   = verbose,
      return_mod = return_mod
      )},
      error = function(e) NULL)
      out[["wlasso"]] <- wlasso_mod
    }
  }

  if ("enet" %in% methods) {
    if (is.null(alpha)) alpha <- seq(0, 1, length.out = 20)
    enet_mod <- tryCatch({tune_glmnet(
      data      = dataset,
      outcome   = outcome,
      exposures = exposures,
      n_folds   = n_folds,
      .alpha    = alpha,
      parallel  = parallel,
      verbose   = verbose,
      return_mod = return_mod
    )},
    error = function(e) NULL)
    out[["enet"]] <- enet_mod
    if (!is.null(weight)) {
      wenet_mod <- tryCatch({tune_wglmnet(
        data           = wdataset,
        outcome        = outcome,
        exposures      = exposures,
        weight         = weight,
        n_folds        = n_folds,
        alpha          = alpha,
        parallel       = parallel,
        verbose        = verbose,
      return_mod = return_mod
      )},
      error = function(e) NULL)
      out[["wenet"]] <- wenet_mod
    }
  }

  if ("rf" %in% methods) {
    rf_mod <- tryCatch({aimTwo::tune_ranger(
      data          = dataset,
      outcome       = outcome,
      exposures     = exposures,
      n_cores       = n_cores,
      mtry_seq      = NULL,
      node_size_seq = NULL,
      n_trees       = n_trees,
      verbose       = verbose,
      return_mod = return_mod
    )},
    error = function(e) NULL)
    out[["rf"]] <- rf_mod
    if (!is.null(weight)) {
      wrf_mod <- tryCatch({wglmnet::wranger(
        data = wdataset,
        col.y = outcome,
        col.x = exposures,
        weights = weight,
        mtry_grid             = NULL,
    min.node.size_grid = NULL,
    num_trees             = n_trees
      ),
      error = function(e) NULL)
      out[["wrf"]] <- wrf_mod
    }
  }

  ### save out
  qsave(
    out,
    file = paste0("results/mgi/", opt$mgi_version, "/", opt$outcome, "/mgi_", opt$mgi_version, "_", opt$outcome, "_t", time_thresholds[i], "_train_hyperparameters_", weight_vars[w], ".qs")
  )

  fwrite(
    map_dfr(out, \(x) x[["param"]]),
    file = paste0("results/mgi/", opt$mgi_version, "/", opt$outcome, "/mgi_", opt$mgi_version, "_", opt$outcome, "_t", time_thresholds[i], "_train_hyperparameters_", weight_vars[w], ".csv")
  )

  }
}

cli_progress_done("done! 🥳")
