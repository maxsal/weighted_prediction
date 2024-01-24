# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable
# author:   max salvatore

# libraries, functions, and options --------------------------------------------
ms::libri(
  data.table, ms, MatchIt, logistf, glue, qs, optparse, EValue, cli,
  tidyverse
)

set.seed(61787)

for (i in list.files("./fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "CA_101.8",
    help = "Outcome phecode [default = %default]"
  ),
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Version of MGI data [default = %default]"
  ),
  make_option("--matching_ratio",
    type = "numeric", default = "2",
    help = "Case:control matching ratio [default = %default]"
  ),
  make_option("--mgi_cohort",
    type = "character", default = "comb",
    help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"
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
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])
cooccur_covs    <- c("age_at_threshold", "female", "length_followup")

## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# read data --------------------------------------------------------------------
## mgi
mgi_tr_pims <- map(
  time_thresholds,
  # read in data
  ~qread(
    glue(
      "data/private/mgi/{opt$mgi_version}/{opt$outcome}/",
      "time_restricted_phenomes/mgi_{opt$outcome}_t",
      "{.x}_r{opt$matching_ratio}_{opt$mgi_version}.qs"
    )
  ),
) |>
  # standardize
  map(
    \(x) {
      vars <- names(x)[!names(x) %in% c("id", "case")]
      standardize_variables(x, cols = vars)
    }
  )
names(mgi_tr_pims) <- glue("t{time_thresholds}_threshold")

mgi_covariates <- qread(glue(
  "data/private/mgi/{opt$mgi_version}/{opt$outcome}/",
  "matched_covariates_r{opt$matching_ratio}.qs"
))

mgi_tr_merged <- map(
  names(mgi_tr_pims),
  \(x) {
    reduce(
      list(mgi_tr_pims[[x]], mgi_covariates[, !c("case")]),
      left_join,
      by = "id"
    ) |>
      mutate(
        age_at_threshold = round(get(x) / 365.25, 1)
      )
  }
)
names(mgi_tr_merged) <- glue("t{time_thresholds}_threshold")

## ukb
ukb_tr_pims <- map(
  time_thresholds,
  # read in data
  ~qread(
    glue(
      "data/private/ukb/{opt$ukb_version}/{opt$outcome}/",
      "time_restricted_phenomes/ukb_{opt$outcome}_t",
      "{.x}_r{opt$matching_ratio}_{opt$ukb_version}.qs"
    )
  )
) |>
  map(
    \(x) {
      vars <- names(x)[!names(x) %in% c("id", "case")]
      standardize_variables(x, cols = vars)
    }
  )
names(ukb_tr_pims) <- glue("t{time_thresholds}_threshold")

ukb_covariates <- read_qs(
  glue(
    "data/private/ukb/{opt$ukb_version}/{opt$outcome}/",
    "matched_covariates_r{opt$matching_ratio}.qs"
  )
)

ukb_tr_merged <- map(
  names(ukb_tr_pims),
  \(x) {
    reduce(
      list(ukb_tr_pims[[x]], ukb_covariates[, !c("case")]),
      left_join,
      by = "id"
    ) |>
      mutate(
        age_at_threshold = round(get(x) / 365.25, 1)
      )
  }
)
names(ukb_tr_merged) <- glue("t{time_thresholds}_threshold")

## phenome
pheinfo <- ms::pheinfox
if (pheinfo[phecode == opt$outcome, sex] != "Both") {
  if (any(c("sex", "female", "male") %in% tolower(cooccur_covs))) {
    cooccur_covs <- cooccur_covs[-which(tolower(cooccur_covs) %in% c("sex", "female", "male"))]
  }
}

# cooccurrence analysis --------------------------------------------------------
## mgi

mgi_results <- list()
cli_progress_bar(name = "MGI PheWAS", total = length(time_thresholds))
for (i in seq_along(time_thresholds)) {
  mgi_results[[i]] <- cooccurrence_analysis(
    data               = mgi_tr_merged[[i]],
    covariates         = cooccur_covs,
    possible_exposures = pheinfo[, phecode],
    min_overlap_count  = 10,
    n_cores            = 16,
    parallel           = TRUE,
    detect_separation  = FALSE,
    logistf_pl         = FALSE
  )
  cli_progress_update()
}
cli_progress_done()
names(mgi_results) <- glue("t{time_thresholds}")

## ukb
ukb_results <- list()
cli_progress_bar(name = "UKB PheWAS", total = length(time_thresholds))
for (i in seq_along(time_thresholds)) {
  ukb_results[[i]] <- cooccurrence_analysis(
    data               = ukb_tr_merged[[i]],
    covariates         = cooccur_covs,
    possible_exposures = pheinfo[, phecode],
    min_overlap_count  = 10,
    n_cores            = parallelly::availableCores() / 4,
    parallel           = TRUE,
    detect_separation  = FALSE,
    logistf_pl         = FALSE
  )
  cli_progress_update()
}
cli_progress_done()
names(ukb_results) <- glue("t{time_thresholds}")

# remove correlated phecodes ---------------------------------------------------
## mgi
mgi_pwide_hits <- lapply(mgi_results, get_pwide_hits)
mgi_topn_hits  <- lapply(mgi_results, get_top_hits)

mgi_pwide_results_post_cor <- mapply(
  \(x, y, z) {
    remove_by_correlation(
      pim = x, cooccurrence_results = y, exposures_to_conside = z, 
      exposure_var = "exposure", p_value_var = "p_value", corr_thresh = 0.5
    )
  },
  mgi_tr_pims, mgi_results, mgi_pwide_hits
)

mgi_pwide_results_post_cor <- mgi_pwide_results_post_cor |>
  as.data.table() |>
  lapply(
    \(x) {
      out        <- as.data.table(x)
      names(out) <- rownames(mgi_pwide_results_post_cor)
      out
    }
  )

mgi_topn_results_post_cor <- mapply(
  \(x, y, z) {
    remove_by_correlation(
      pim = x, cooccurrence_results = y, exposures_to_conside = z, 
      exposure_var = "exposure", p_value_var = "p_value", corr_thresh = 0.5
    )
  },
  mgi_tr_pims, mgi_results, mgi_topn_hits
)

mgi_topn_results_post_cor <- mgi_topn_results_post_cor |>
  as.data.table() |>
  lapply(
    \(x) {
      out        <- as.data.table(x)
      names(out) <- rownames(mgi_topn_results_post_cor)
      out
    }
  )

## ukb
ukb_pwide_hits <- lapply(ukb_results, get_pwide_hits)
ukb_topn_hits  <- lapply(ukb_results, get_top_hits)

ukb_pwide_results_post_cor <- mapply(
  \(x, y, z) {
    remove_by_correlation(
      pim = x, cooccurrence_results = y, exposures_to_conside = z, 
      exposure_var = "exposure", p_value_var = "p_value", corr_thresh = 0.5
    )
  },
  ukb_tr_pims, ukb_results, ukb_pwide_hits
)

ukb_pwide_results_post_cor <- ukb_results_post_cor |>
  as.data.table() |>
  lapply(
    \(x) {
      out        <- as.data.table(x)
      names(out) <- rownames(ukb_results_post_cor)
      out
    }
  )

ukb_topn_results_post_cor <- mapply(
  \(x, y, z) {
    remove_by_correlation(
      pim = x, cooccurrence_results = y, exposures_to_conside = z, 
      exposure_var = "exposure", p_value_var = "p_value", corr_thresh = 0.5
    )
  },
  ukb_tr_pims, ukb_results, ukb_topn_hits
)

ukb_topn_results_post_cor <- ukb_topn_results_post_cor |>
  as.data.table() |>
  lapply(
    \(x) {
      out        <- as.data.table(x)
      names(out) <- rownames(ukb_topn_results_post_cor)
      out
    }
  )

# save results -----------------------------------------------------------------
## mgi
for (i in seq_along(time_thresholds)) {
  save_qs(
    x = mgi_results[[i]],
    file = glue(
      "results/mgi/{opt$mgi_version}/{opt$outcome}/",
      "mgi_{opt$outcome}_t{time_thresholds[i]}_",
      "r{opt$matching_ratio}_{opt$mgi_version}_results.qs"
    )
  )
}

## mgi pwide
qsave(
  mgi_pwide_results_post_cor,
  file = glue(
    "results/mgi/{opt$mgi_version}/{opt$outcome}/",
    "mgi_{opt$outcome}_pwide_r{opt$matching_ratio}_{opt$mgi_version}_results_post_cor.qs"
  )
)
## mgi topn
qsave(
  mgi_topn_results_post_cor,
  file = glue(
    "results/mgi/{opt$mgi_version}/{opt$outcome}/",
    "mgi_{opt$outcome}_topn_r{opt$matching_ratio}_{opt$mgi_version}_results_post_cor.qs"
  )
)

## ukb
for (i in seq_along(time_thresholds)) {
  save_qs(
    x = ukb_results[[i]],
    file = glue(
      "results/ukb/{opt$ukb_version}/{opt$outcome}/",
      "ukb_{opt$outcome}_t{time_thresholds[i]}_r{opt$matching_ratio}_",
      "{opt$ukb_version}_results.qs"
    )
  )
}

## ukb pwide
qsave(
  ukb_pwide_results_post_cor,
  file = glue(
    "results/ukb/{opt$ukb_version}/{opt$outcome}/",
    "ukb_{opt$outcome}_pwide_r{opt$matching_ratio}_{opt$ukb_version}_results_post_cor.qs"
  )
)
## ukb topn
qsave(
  ukb_topn_results_post_cor,
  file = glue(
    "results/ukb/{opt$ukb_version}/{opt$outcome}/",
    "ukb_{opt$outcome}_topn_r{opt$matching_ratio}_{opt$ukb_version}_results_post_cor.qs"
  )
)

cli_alert_success("script complete! 🎉🎉🎉")
