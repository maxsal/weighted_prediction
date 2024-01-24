# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds with weights
# author:   max salvatore

# libraries, functions, and options --------------------------------------------
ms::libri(data.table, MatchIt, logistf, glue, qs, cli, optparse, ms)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "CA_105.1",
    help = "Outcome phecode [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Version of UKB data [default = %default]"
  ),
  make_option("--time_thresholds",
    type = "character", default = "0,0.5,1,2,5",
    help = glue(
      "Time thresholds for the phenome data ",
      "[default = %default]"
    )
  ),
  make_option("--mod_type",
    type = "character", default = "glm",
    help = glue(
      "Type of model to use in cooccurrence analysis - ",
      "glm, logistf or SPAtest [default = %default]"
    )
  ),
  make_option("--weights",
    type = "character", default = "weight",
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
cooccur_covs    <- c("age_at_threshold", "female", "length_followup")

## extract file paths
file_paths <- get_files(ukb_version = opt$ukb_version)

# read data --------------------------------------------------------------------
## mgi
ukb_tr_pims <- lapply(
  seq_along(time_thresholds),
  \(x) {
    glue(
      "data/private/ukb/{opt$ukb_version}/{opt$outcome}/",
      "time_restricted_phenomes/ukb_{opt$outcome}_t",
      "{time_thresholds[x]}_{opt$ukb_version}.qs"
    ) |>
      read_qs()
  }
)
names(ukb_tr_pims) <- glue("t{time_thresholds}_threshold")

ukb_covariates <- read_qs(glue(
  "data/private/ukb/{opt$ukb_version}/{opt$outcome}/",
  "matched_covariates.qs"
))

ukb_weights <- fread(file_paths[["ukb"]][["weight_file"]], colClasses = "character")[
    ,
    .(id = f.eid, weight = LassoWeight)
]
ukb_weights[, weight := as.numeric(weight)]

## phenome
pheinfo <- ms::pheinfox
if (pheinfo[phecode == opt$outcome, sex] != "Both") {
  if (any(c("sex", "female", "male") %in% tolower(cooccur_covs))) {
    cooccur_covs <- cooccur_covs[-which(tolower(cooccur_covs) %in% c("sex", "female", "male"))]
  }
}

# weights ----------------------------------------------------------------------
weight_vars <- unlist(strsplit(opt$weights, ","))

weight_id <- c("id", weight_vars)

ukb_tr_merged <- lapply(
  names(ukb_tr_pims),
  \(x) {
    merge_list(list(ukb_tr_pims[[x]], ukb_covariates[, !c("case")], ukb_weights[, ..weight_id]), by_var = "id", join_fn = dplyr::left_join)[, `:=`(
      age_at_threshold = round(get(x) / 365.25, 1)
    )]
  }
)
names(ukb_tr_merged) <- glue("t{time_thresholds}_threshold")


# cooccurrence analysis --------------------------------------------------------
out <- list()
for (w in seq_along(weight_vars)) {
  cli_alert(glue("cooccurrence using {weight_vars[w]} [{w}/{length(weight_vars)}]..."))
  out[[w]] <- lapply(
    seq_along(time_thresholds),
    \(x) {
      cli_alert_info("{x}")
      cooccurrence_analysis(
        data               = ukb_tr_merged[[x]],
        covariates         = cooccur_covs,
        possible_exposures = pheinfo[, phecode],
        weight_var         = weight_vars[w],
        min_overlap_count  = 10,
        n_cores            = parallelly::availableCores() / 4,
        parallel           = TRUE
      )
    }
  )
}
names(out) <- weight_vars

# save results -----------------------------------------------------------------
## mgi
for (w in seq_along(weight_vars)) {
  for (i in seq_along(time_thresholds)) {
    save_qs(
      x = out[[w]][[i]],
      file = glue(
        "results/ukb/{opt$ukb_version}/{opt$outcome}/",
        "ukb_{opt$outcome}_t{time_thresholds[i]}_",
        "{opt$ukb_version}_{weight_vars[w]}_results.qs"
      )
    )
  }
}

cli_alert_success("done! 🥳")
