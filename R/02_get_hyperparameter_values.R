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

# get hyperparameters ----------------------------------------------------------
out <- map(
  seq_along(weight_vars),
  \(w) {
    map(
      seq_along(time_thresholds),
      \(i) {
        cli_alert_info("t{time_thresholds[i]}")
        tune_models(
          data      = mgi_tr_merged[[i]][group == "train", ],
          outcome   = "case",
          exposures = names(mgi_tr_merged[[i]])[names(mgi_tr_merged[[i]]) %in% ms::pheinfox[, phecode]],
          weight    = weight_vars[w],
          n_cores   = 16
        )
      }
    ) |>
      set_names(glue("t{time_thresholds}_threshold"))
  }
) |> set_names(weight_vars)

# save results -----------------------------------------------------------------
walk(
  seq_along(weight_vars),
  \(w) {
    walk(
      seq_along(time_thresholds),
      \(x) {
        fwrite(
          out[[w]][[x]],
          paste0("results/mgi/", opt$mgi_version, "/", opt$outcome, "/mgi_", opt$mgi_version, "_", opt$outcome, "_t", time_thresholds[x], "_hyperparameters_", weight_vars[w], ".csv")
        )
      }
    )
  }
)

cli_alert_success("done! 🎉")
