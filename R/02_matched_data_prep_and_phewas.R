# prepare mgi and ukb data and conduct unweighted phewas
# author:  max salvatore
# date:    20231201

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)
set.seed(61787)

# load libraries
ms::libri(data.table, qs, MatchIt, optparse, glue, cli, ms, tidyverse)

# load functions
walk(list.files("fn/", full.names = TRUE), source)

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
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Version of UKB data [default = %default]"
  ),
  make_option("--time_thresholds",
    type = "character", default = "0,1,2,5",
    help = glue(
      "Time thresholds for the phenome data ",
      "[default = %default]"
    )
  ),
  make_option("--nearest_matching_var",
    type = "character",
    default = "age,length_followup",
    help = glue(
      "Matching variables by nearest  [default = ",
      "%default]"
    )
  ),
  make_option("--exact_matching_var",
    type = "character", default = "female",
    help = "Matching variables by exact [default = %default]"
  ),
  make_option("--matching_caliper",
    type = "numeric", default = "0.25",
    help = "Matching caliper [default = %default]"
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

# 2. specifications (specifies outcome) --------------------------------------
time_thresholds       <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])
nearest_matching_vars <- strsplit(opt$nearest_matching_var, ",")[[1]]
exact_matching_vars   <- strsplit(opt$exact_matching_var, ",")[[1]]

# process outcome information
outcome_sex <- ms::pheinfox[phecode == opt$outcome, sex]

if (tolower(outcome_sex) %in% c("female", "f", "male", "m")) {
  exact_matching_vars <- exact_matching_vars[!exact_matching_vars %in% c("female", "male", "m", "f")]
  if (length(exact_matching_vars) == 0) {
    exact_matching_vars <- NULL
  }
}

# 3. extra preparations (outcome-specific) -------------------------------------
## pull file paths corresponding to the data version specified
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

## check that outcome exists in both datasets
cli_progress_step("checking that outcome exists in both phenomes...")
common_codes <- fread("data/public/phecodex_20plus.csv")[plus20 == 1, phecode]
if (!opt$outcome %in% common_codes) stop("outcome is not defined in AOU, MGI, and UKB. stopping.")

## confirm file structure for a given outcome exists - if not, create paths
### mgi
check_folder_structure(
  cohort          = "mgi",
  data_version    = opt$mgi_version,
  outcome_phecode = opt$outcome
)

dir.create(glue("data/private/mgi/{opt$mgi_version}/{opt$outcome}/time_restricted_phenomes/"), recursive = TRUE)


### ukb
check_folder_structure(
  cohort          = "ukb",
  data_version    = opt$ukb_version,
  outcome_phecode = opt$outcome
)

dir.create(glue("data/private/ukb/{opt$ukb_version}/{opt$outcome}/time_restricted_phenomes/"), recursive = TRUE)




# 4. read data -----------------------------------------------------------------
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")

## mgi
cli_progress_step("loading mgi data")
### demographics
mgi_demo <- qread(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))


### icd-phecode data
mgi_full_phe  <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs")) |>
  dplyr::filter(id %in% mgi_demo[, id] & phecode %in% common_codes)

### mgi weights
mgi_weights <- qread(glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs"))

## ukb
cli_progress_step("loading ukb data")
### demographics
ukb_demo <- qread(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

### phecode-dsb data
ukb_full_phe  <- qread(glue("data/private/ukb/{opt$ukb_version}/UKB_FULL_PHECODEX_DSB_{opt$ukb_version}.qs")) |>
  dplyr::filter(id %in% ukb_demo[, id] & phecode %in% common_codes)

### ukb weights
ukb_weights <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/UKBSelectionWeights.tab")

source("R/tr_prepare-utils.R")

# prepare matched data
## mgi
mgi_prepped_data <- prepare_matched_phenomes(
  phe_dsb_data          = mgi_full_phe,
  demo_data             = mgi_demo,
  outcome               = opt$outcome,
  case_definition       = "first",
  malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
  specific_phecodes     = cancer_phecodesx[specific == 1, phecode],
  matched               = TRUE,
  match_ratio           = opt$matching_ratio,
  match_caliper         = opt$matching_caliper,
  nearest_vars          = strsplit(opt$nearest_matching_var, ",")[[1]],
  exact_vars            = strsplit(opt$exact_matching_var, ",")[[1]],
  control_definition    = "no_cancer",
  time_thresholds       = as.numeric(strsplit(opt$time_thresholds, ",")[[1]]),
  exclude_other_cancers = TRUE
)

# time restricted pim
walk(
  seq_along(time_thresholds),
  \(x) qsave(mgi_prepped_data$tr_pim[[x]], glue(
    "data/private/mgi/{opt$mgi_version}/",
    "{opt$outcome}/time_restricted_phenomes/mgi_",
    "{opt$mgi_version}_{opt$outcome}_",
    "t{time_thresholds[x]}_pim_",
    "r{opt$matching_ratio}.qs"
  ))
)

# match data
qsave(
  mgi_prepped_data$match_data,
  glue(
    "data/private/mgi/{opt$mgi_version}/",
    "{opt$outcome}/mgi_",
    "{opt$mgi_version}_{opt$outcome}_",
    "match_data_",
    "r{opt$matching_ratio}.qs"
  )
)

# case data
qsave(
  x = mgi_prepped_data$case_data,
  file = glue(
    "data/private/mgi/{opt$mgi_version}/",
    "{opt$outcome}/mgi_",
    "{opt$mgi_version}_{opt$outcome}_",
    "case_data_",
    "r{opt$matching_ratio}.qs"
  )
)


## ukb
ukb_prepped_data <- prepare_matched_phenomes(
  phe_dsb_data          = ukb_full_phe,
  demo_data             = ukb_demo |>
    mutate(
      length_followup = last_dsb - first_dsb,
      female = as.numeric(sex == "Female")
    ) |>
    rename(age = age_at_first_diagnosis),
  outcome               = opt$outcome,
  case_definition       = "first",
  malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
  specific_phecodes     = cancer_phecodesx[specific == 1, phecode],
  matched               = TRUE,
  match_ratio           = opt$matching_ratio,
  match_caliper         = opt$matching_caliper,
  nearest_vars          = strsplit(opt$nearest_matching_var, ",")[[1]],
  exact_vars            = strsplit(opt$exact_matching_var, ",")[[1]],
  control_definition    = "no_cancer",
  time_thresholds       = as.numeric(strsplit(opt$time_thresholds, ",")[[1]]),
  exclude_other_cancers = TRUE
)

# time restricted pim
walk(
  seq_along(time_thresholds),
  \(x) qsave(ukb_prepped_data$tr_pim[[x]], glue(
    "data/private/ukb/{opt$ukb_version}/",
    "{opt$outcome}/time_restricted_phenomes/ukb_",
    "{opt$ukb_version}_{opt$outcome}_",
    "t{time_thresholds[x]}_pim_",
    "r{opt$matching_ratio}.qs"
  ))
)

# match data
qsave(
  ukb_prepped_data$match_data,
  glue(
    "data/private/ukb/{opt$ukb_version}/",
    "{opt$outcome}/ukb_",
    "{opt$ukb_version}_{opt$outcome}_",
    "match_data_",
    "r{opt$matching_ratio}.qs"
  )
)

# case data
qsave(
  x = ukb_prepped_data$case_data,
  file = glue(
    "data/private/ukb/{opt$ukb_version}/",
    "{opt$outcome}/ukb_",
    "{opt$ukb_version}_{opt$outcome}_",
    "case_data_",
    "r{opt$matching_ratio}.qs"
  )
)

# phewas
## mgi
mgi_phewas_res <- map(
  seq_along(mgi_prepped_data$tr_pim),
  \(x) {
    cli_progress_step(paste0("running phewas for time threshold ", time_thresholds[x]))
    mgi_data       <- mgi_prepped_data$tr_pim[[x]]
    mgi_phecodes   <- names(mgi_data)[names(mgi_data) %in% ms::pheinfox[, phecode]]
    mgi_phewas_res <- ms::map_phewas(
      data       = mgi_data,
      outcome    = "case",
      exposures  = mgi_phecodes,
      covariates = c("age", "length_followup"),
      method     = "glm",
      workers    = 16
    )
  }
) |> set_names(names(mgi_prepped_data$tr_pim))

walk(
  seq_along(mgi_phewas_res),
  \(x) qsave(mgi_phewas_res[[x]], glue("results/mgi/{opt$mgi_version}/",
                                     "{opt$outcome}/mgi_",
                                     "{opt$mgi_version}_phewas_{opt$outcome}_",
                                     "t{time_thresholds[x]}_",
                                     "r{opt$matching_ratio}.qs"))
)


mgi_phewas_plots <- map(
  seq_along(mgi_phewas_res),
  \(x) {
    plot_phewasx(
      mgi_phewas_res[[x]],
      phe_var = "phecode",
      title = stringr::str_wrap(glue("PheWAS for {ms::pheinfox[phecode == opt$outcome, description]} ",
                   "[{opt$outcome}] at t = {time_thresholds[x]} in MGI"), width = 50)
      )
  }
)

walk(
  seq_along(mgi_phewas_plots),
  \(x) ggsave(
    filename = glue("results/mgi/{opt$mgi_version}/{opt$outcome}/mgi_",
                    "{opt$mgi_version}_phewas_plot_{opt$outcome}_",
                    "t{time_thresholds[x]}_",
                    "r{opt$matching_ratio}.pdf"),
    plot = mgi_phewas_plots[[x]],
    width = 8,
    height = 6,
    device = cairo_pdf
  )
)

## ukb
ukb_phewas_res <- map(
  seq_along(ukb_prepped_data$tr_pim),
  \(x) {
    cli_progress_step(paste0("running phewas for time threshold ", time_thresholds[x]))
    ukb_data <- ukb_prepped_data$tr_pim[[x]]
    ukb_phecodes <- names(ukb_data)[names(ukb_data) %in% ms::pheinfox[, phecode]]
    ukb_phewas_res <- ms::map_phewas(
      data       = ukb_data,
      outcome    = "case",
      exposures  = ukb_phecodes,
      covariates = c("age", "length_followup"),
      method     = "glm",
      workers    = 16
    )
  }
) |> set_names(names(ukb_prepped_data$tr_pim))

walk(
  seq_along(ukb_phewas_res),
  \(x) qsave(ukb_phewas_res[[x]], glue(
    "results/ukb/{opt$ukb_version}/{opt$outcome}/ukb_",
    "{opt$ukb_version}_phewas_{opt$outcome}_",
    "t{time_thresholds[x]}_",
    "r{opt$matching_ratio}.qs"
  ))
)

ukb_phewas_plots <- map(
  seq_along(ukb_phewas_res),
  \(x) {
    plot_phewasx(
      ukb_phewas_res[[x]],
      phe_var = "phecode",
      title = stringr::str_wrap(glue("PheWAS for {ms::pheinfox[phecode == opt$outcome, description]} ",
                   "[{opt$outcome}] at t = {time_thresholds[x]} in UKB"), width = 50)
    )
  }
)

walk(
  seq_along(ukb_phewas_plots),
  \(x) ggsave(
    filename = glue(
      "results/ukb/{opt$ukb_version}/{opt$outcome}/ukb_",
      "{opt$ukb_version}_phewas_plot_{opt$outcome}_",
      "t{time_thresholds[x]}_",
      "r{opt$matching_ratio}.pdf"
    ),
    plot = ukb_phewas_plots[[x]],
    width = 8,
    height = 6,
    device = cairo_pdf
  )
)

#### correlation thresholding of top 50
mgi_after_cor <- map(
  seq_along(time_thresholds),
  \(i) {
    top_after_correlation <- remove_by_correlation(
      cooccurrence_results = as.data.table(mgi_phewas_res[[i]]),
      pim = mgi_prepped_data$tr_pim[[i]],
      exposure_var = "phecode",
      p_value_var = "p_value",
      top_n = 50,
      corr_thresh = 0.5,
      weights = NULL
    )
  }
) |> set_names(glue("t{time_thresholds}"))

ukb_after_cor <- map(
  seq_along(time_thresholds),
  \(i) {
    top_after_correlation <- remove_by_correlation(
      cooccurrence_results = as.data.table(ukb_phewas_res[[i]]),
      pim = ukb_prepped_data$tr_pim[[i]],
      exposure_var = "phecode",
      p_value_var = "p_value",
      top_n = 50,
      corr_thresh = 0.5,
      weights = NULL
    )
  }
) |> set_names(glue("t{time_thresholds}"))

walk(
  seq_along(time_thresholds),
  \(i) {
    qsave(
      mgi_after_cor[[i]],
      glue(
        "results/mgi/{opt$mgi_version}/",
        "{opt$outcome}/mgi_",
        "{opt$mgi_version}_post_cor_phe_{opt$outcome}_",
        "t{time_thresholds[i]}_",
        "r{opt$matching_ratio}.qs"
      )
    )
    qsave(
      ukb_after_cor[[i]],
      glue(
        "results/ukb/{opt$ukb_version}/",
        "{opt$outcome}/ukb_",
        "{opt$ukb_version}_post_cor_phe_{opt$outcome}_",
        "t{time_thresholds[i]}_",
        "r{opt$matching_ratio}.qs"
      )
    )
  }
)
####

#### WEIGHTED PHEWAS ----------------------------------------------------------
# MERGE IN WEIGHTS
## mgi
mgi_tr_merged <- map(
  seq_along(time_thresholds),
  \(i) {
    left_join(
      mgi_prepped_data$tr_pim[[i]],
      mgi_weights,
      by = "id"
    )
  }
) |> set_names(glue("t{time_thresholds}"))

# ukb
ukb_tr_merged <- map(
  seq_along(time_thresholds),
  \(i) {
    left_join(
      ukb_prepped_data$tr_pim[[i]],
      ukb_weights[, .(id = as.character(f.eid), weight = LassoWeight)],
      by = "id"
    )
  }
) |> set_names(glue("t{time_thresholds}"))

# IP
## mgi
mgi_ip_phewas_res <- map(
  seq_along(mgi_tr_merged),
  \(x) {
    cli_progress_step(paste0("running phewas for time threshold ", time_thresholds[x]))
    mgi_data <- mgi_tr_merged[[x]]
    mgi_phecodes <- names(mgi_data)[names(mgi_data) %in% ms::pheinfox[, phecode]]
    mgi_dsn <- survey::svydesign(
      id      = ~1,
      weights = ~ip_selection,
      data    = mgi_data[!is.na(ip_selection), ]
    )
    mgi_phewas_res <- ms::map_phewas(
      data        = mgi_data,
      design      = mgi_dsn,
      outcome     = "case",
      exposures   = mgi_phecodes,
      covariates  = c("age", "length_followup"),
      method      = "weighted",
      .weight_var = "ip_selection",
      workers     = 16
    )
  }
) |> set_names(names(mgi_tr_merged))

### save
walk(
  seq_along(mgi_ip_phewas_res),
  \(x) qsave(mgi_ip_phewas_res[[x]], glue(
    "results/mgi/{opt$mgi_version}/",
    "{opt$outcome}/mgi_",
    "{opt$mgi_version}_ip_phewas_{opt$outcome}_",
    "t{time_thresholds[x]}_",
    "r{opt$matching_ratio}.qs"
  ))
)

mgi_ip_phewas_plots <- map(
  seq_along(mgi_ip_phewas_res),
  \(x) {
    plot_phewasx(
      mgi_ip_phewas_res[[x]],
      phe_var = "phecode",
      title = stringr::str_wrap(glue(
        "IP-weighted PheWAS for {ms::pheinfox[phecode == opt$outcome, description]} ",
        "[{opt$outcome}] at t = {time_thresholds[x]} in MGI"
      ), width = 50)
    )
  }
)

walk(
  seq_along(mgi_ip_phewas_plots),
  \(x) ggsave(
    filename = glue(
      "results/mgi/{opt$mgi_version}/{opt$outcome}/mgi_",
      "{opt$mgi_version}_ip_phewas_plot_{opt$outcome}_",
      "t{time_thresholds[x]}_",
      "r{opt$matching_ratio}.pdf"
    ),
    plot = mgi_ip_phewas_plots[[x]],
    width = 8,
    height = 6,
    device = cairo_pdf
  )
)

mgi_ip_after_cor <- map(
  seq_along(time_thresholds),
  \(i) {
    top_after_correlation <- remove_by_correlation(
      cooccurrence_results = as.data.table(mgi_ip_phewas_res[[i]]),
      pim = mgi_tr_merged[[i]],
      exposure_var = "phecode",
      p_value_var = "p_value",
      top_n = 50,
      corr_thresh = 0.5,
      weights = NULL
    )
  }
) |> set_names(glue("t{time_thresholds}"))

walk(
  seq_along(time_thresholds),
  \(i) {
    qsave(
      mgi_ip_after_cor[[i]],
      glue(
        "results/mgi/{opt$mgi_version}/",
        "{opt$outcome}/mgi_",
        "{opt$mgi_version}_ip_post_cor_phe_{opt$outcome}_",
        "t{time_thresholds[i]}_",
        "r{opt$matching_ratio}.qs"
      )
    )
  }
)

## ukb
ukb_ip_phewas_res <- map(
  seq_along(ukb_tr_merged),
  \(x) {
    cli_progress_step(paste0("running phewas for time threshold ", time_thresholds[x]))
    ukb_data <- ukb_tr_merged[[x]]
    ukb_phecodes <- names(ukb_data)[names(ukb_data) %in% ms::pheinfox[, phecode]]
    ukb_dsn <- survey::svydesign(
      id      = ~1,
      weights = ~weight,
      data    = ukb_data[!is.na(weight), ]
    )
    ukb_phewas_res <- ms::map_phewas(
      data        = ukb_data,
      design      = ukb_dsn,
      outcome     = "case",
      exposures   = ukb_phecodes,
      covariates  = c("age", "length_followup"),
      method      = "weighted",
      .weight_var = "weight",
      workers     = 16
    )
  }
) |> set_names(names(ukb_tr_merged))

### save
walk(
  seq_along(ukb_ip_phewas_res),
  \(x) qsave(ukb_ip_phewas_res[[x]], glue(
    "results/ukb/{opt$ukb_version}/",
    "{opt$outcome}/ukb_",
    "{opt$ukb_version}_ip_phewas_{opt$outcome}_",
    "t{time_thresholds[x]}_",
    "r{opt$matching_ratio}.qs"
  ))
)

ukb_ip_phewas_plots <- map(
  seq_along(ukb_ip_phewas_res),
  \(x) {
    plot_phewasx(
      ukb_ip_phewas_res[[x]],
      phe_var = "phecode",
      title = stringr::str_wrap(glue(
        "IP-weighted PheWAS for {ms::pheinfox[phecode == opt$outcome, description]} ",
        "[{opt$outcome}] at t = {time_thresholds[x]} in ukb"
      ), width = 50)
    )
  }
)

walk(
  seq_along(ukb_ip_phewas_plots),
  \(x) ggsave(
    filename = glue(
      "results/ukb/{opt$ukb_version}/{opt$outcome}/ukb_",
      "{opt$ukb_version}_ip_phewas_plot_{opt$outcome}_",
      "t{time_thresholds[x]}_",
      "r{opt$matching_ratio}.pdf"
    ),
    plot = ukb_ip_phewas_plots[[x]],
    width = 8,
    height = 6,
    device = cairo_pdf
  )
)

ukb_ip_after_cor <- map(
  seq_along(time_thresholds),
  \(i) {
    top_after_correlation <- remove_by_correlation(
      cooccurrence_results = as.data.table(ukb_ip_phewas_res[[i]]),
      pim = ukb_tr_merged[[i]],
      exposure_var = "phecode",
      p_value_var = "p_value",
      top_n = 50,
      corr_thresh = 0.5,
      weights = NULL
    )
  }
) |> set_names(glue("t{time_thresholds}"))

walk(
  seq_along(time_thresholds),
  \(i) {
    qsave(
      ukb_ip_after_cor[[i]],
      glue(
        "results/ukb/{opt$ukb_version}/",
        "{opt$outcome}/ukb_",
        "{opt$ukb_version}_ip_post_cor_phe_{opt$outcome}_",
        "t{time_thresholds[i]}_",
        "r{opt$matching_ratio}.qs"
      )
    )
  }
)

# PS
## mgi
mgi_ps_phewas_res <- map(
  seq_along(mgi_tr_merged),
  \(x) {
    cli_progress_step(paste0("running phewas for time threshold ", time_thresholds[x]))
    mgi_data <- mgi_tr_merged[[x]]
    mgi_phecodes <- names(mgi_data)[names(mgi_data) %in% ms::pheinfox[, phecode]]
    mgi_dsn <- survey::svydesign(
      id      = ~1,
      weights = ~ps_selection,
      data    = mgi_data[!is.na(ps_selection), ]
    )
    mgi_phewas_res <- ms::map_phewas(
      data        = mgi_data,
      design      = mgi_dsn,
      outcome     = "case",
      exposures   = mgi_phecodes,
      covariates  = c("age", "length_followup"),
      method      = "weighted",
      .weight_var = "ps_selection",
      workers     = 16
    )
  }
) |> set_names(names(mgi_tr_merged))

### save
walk(
  seq_along(mgi_ps_phewas_res),
  \(x) qsave(mgi_ps_phewas_res[[x]], glue(
    "results/mgi/{opt$mgi_version}/",
    "{opt$outcome}/mgi_",
    "{opt$mgi_version}_ps_phewas_{opt$outcome}_",
    "t{time_thresholds[x]}_",
    "r{opt$matching_ratio}.qs"
  ))
)

mgi_ps_phewas_plots <- map(
  seq_along(mgi_ps_phewas_res),
  \(x) {
    plot_phewasx(
      mgi_ps_phewas_res[[x]],
      phe_var = "phecode",
      title = stringr::str_wrap(glue(
        "PS-weighted PheWAS for {ms::pheinfox[phecode == opt$outcome, description]} ",
        "[{opt$outcome}] at t = {time_thresholds[x]} in MGI"
      ), width = 50)
    )
  }
)

walk(
  seq_along(mgi_ps_phewas_plots),
  \(x) ggsave(
    filename = glue(
      "results/mgi/{opt$mgi_version}/{opt$outcome}/mgi_",
      "{opt$mgi_version}_ps_phewas_plot_{opt$outcome}_",
      "t{time_thresholds[x]}_",
      "r{opt$matching_ratio}.pdf"
    ),
    plot = mgi_ps_phewas_plots[[x]],
    width = 8,
    height = 6,
    device = cairo_pdf
  )
)

mgi_ps_after_cor <- map(
  seq_along(time_thresholds),
  \(i) {
    top_after_correlation <- remove_by_correlation(
      cooccurrence_results = as.data.table(mgi_ps_phewas_res[[i]]),
      pim = mgi_tr_merged[[i]],
      exposure_var = "phecode",
      p_value_var = "p_value",
      top_n = 50,
      corr_thresh = 0.5,
      weights = NULL
    )
  }
) |> set_names(glue("t{time_thresholds}"))

walk(
  seq_along(time_thresholds),
  \(i) {
    qsave(
      mgi_ps_after_cor[[i]],
      glue(
        "results/mgi/{opt$mgi_version}/",
        "{opt$outcome}/mgi_",
        "{opt$mgi_version}_ps_post_cor_phe_{opt$outcome}_",
        "t{time_thresholds[i]}_",
        "r{opt$matching_ratio}.qs"
      )
    )
  }
)

## ukb
### NOT APPLICABLE

####

cli_alert_success("script done! 🙌")
