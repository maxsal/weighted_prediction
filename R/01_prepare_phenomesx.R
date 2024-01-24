# prepare mgi and ukb data for time-restricted phecode-phecode phewas
# author:  max salvatore
# date:    20230809

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)
set.seed(61787)

# load libraries
ms::libri(data.table, qs, MatchIt, optparse, glue, cli, ms, tidyverse)

# load functions
walk(list.files("fn/", full.names = TRUE), source)
source("R/match-utils.R")

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
    default = "age_at_first_diagnosis,length_followup",
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
common_codes_tab <- fread("data/public/phecodex_20plus.csv")
common_codes <- common_codes_tab[aou >= 20 & mgi >= 20, phecode]
also_ukb <- opt$outcome %in% common_codes_tab[ukb >= 20, phecode]
if (!also_ukb) {
  cli_alert("outcome is not defined in UKB. skipping UKB analyses...")
} else {
  common_codes_ukb <- common_codes_tab[ukb >= 20, phecode]
}
if (!opt$outcome %in% common_codes) stop("outcome is not defined in MGI and AOU. stopping.")

## confirm file structure for a given outcome exists - if not, create paths
### mgi
check_folder_structure(
  cohort          = "mgi",
  data_version    = opt$mgi_version,
  outcome_phecode = opt$outcome
)

### ukb
if (also_ukb) {
  check_folder_structure(
    cohort          = "ukb",
    data_version    = opt$ukb_version,
    outcome_phecode = opt$outcome
  )
}





# 4. read data -----------------------------------------------------------------
## mgi
cli_progress_step("loading mgi data")
### demographics
mgi_demo <- qread(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))
if (outcome_sex != "Both") {
  mgi_demo <- mgi_demo[sex %in% c(outcome_sex, substr(outcome_sex, 1, 1)), ]
}

### pim - IS PIM DATA NECESSARY?
mgi_pim0 <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_PIM0X_{opt$mgi_version}.qs"))[id %in% mgi_demo[, id], ]
# restrict to common phecodes
pim_vars <- c("id", common_codes)
mgi_pim0 <- mgi_pim0[, ..pim_vars]
setnames(mgi_pim0, old = opt$outcome, new = "outcome")

### icd-phecode data
mgi_full_phe  <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs")) |>
  dplyr::filter(id %in% mgi_demo[, id] & phecode %in% common_codes)
mgi_first_phe <- mgi_full_phe[ mgi_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1 ]

## ukb
if (also_ukb) {
  cli_progress_step("loading ukb data")
  ### demographics
  ukb_demo <- qread(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))
  if (outcome_sex != "Both") {
    ukb_demo <- ukb_demo[sex %in% c(outcome_sex, substr(outcome_sex, 1, 1)), ]
  }
  ukb_demo <- ukb_demo[, .(id, age_at_first_diagnosis = age_at_first_diagnosisx, age_at_last_diagnosisx, ethn = race_eth, sex)]
  ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

  ### pim
  ukb_pim_vars <- c("id", common_codes_ukb)
  ukb_pim0 <- qread(glue("data/private/ukb/{opt$ukb_version}/UKB_PIM0X_{opt$ukb_version}.qs"))[id %in% ukb_demo[, id], ]
  ukb_pim0 <- ukb_pim0[, ..ukb_pim_vars]
  setnames(ukb_pim0, old = opt$outcome, new = "outcome")

  ### phecode-dsb data
  ukb_full_phe  <- qread(glue("data/private/ukb/{opt$ukb_version}/UKB_FULL_PHECODEX_DSB_{opt$ukb_version}.qs")) |>
    dplyr::filter(id %in% ukb_demo[, id] & phecode %in% common_codes_ukb)
  ukb_first_phe <- ukb_full_phe[ ukb_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1 ]
}





# 5. identify cases ------------------------------------------------------------
## mgi
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")

mgi_case_data <- prepare_case_data(
  phe_dsb_data = mgi_full_phe,
  demo_data    = mgi_demo,
  outcome      = opt$outcome,
  malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
  specific_phecodes     = cancer_phecodesx[specific == 1, phecode]
)

if (also_ukb) {
  ukb_case_data <- prepare_case_data(
    phe_dsb_data = ukb_full_phe,
    demo_data    = ukb_demo,
    outcome      = opt$outcome,
  malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
  specific_phecodes     = cancer_phecodesx[specific == 1, phecode]
  )
}





# 6. calculate diagnostic metrics ----------------------------------------------
## mgi
mgi_demo <- mgi_demo[!(id %in% mgi_case_data$exclude_ids), ][, `:=`(
  female = fcase(sex == "Female", 1, sex == "Male", 0),
  case = fifelse(id %in% mgi_case_data$case[, id], 1, 0),
  length_followup = round((last_dsbx - first_dsbx) / 365.25, 1)
)][, .(
  id, female, case, length_followup, age_at_first_diagnosis = age_at_first_diagnosisx,
  first_dsb = first_dsbx
)]

## ukb
if (also_ukb) {
  ukb_diag_metrics <- get_icd_phecode_metrics(
    full_phe_data = ukb_first_phe
  )[, id := as.character(id)]
  ukb_matching_cov <- merge.data.table(
    ukb_demo[!(id %in% ukb_case_data$exclude_ids), .(id, female = fcase(sex == "Female", 1, sex == "Male", 0))],
    ukb_diag_metrics,
    by = "id"
  )[, case := fifelse(id %in% ukb_case_data$case[, id], 1, 0)]
}





# 7. perform matching ----------------------------------------------------------
## mgi
cli_progress_step(glue(
  "performing 1:{opt$matching_ratio} case:non-case ",
  "matching in mgi..."
))

mgi_post_match_cov <- perform_matching(
  matching_data   = mgi_demo,
  case_data       = mgi_case_data$case,
  match_ratio     = opt$matching_ratio,
  nearest_vars    = nearest_matching_vars,
  exact_vars      = exact_matching_vars,
  match_caliper   = opt$matching_caliper,
  time_thresholds = time_thresholds 
)

if (!dir.exists(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "{opt$outcome}/",
  "time_restricted_phenomes/"
))) {
  dir.create(
    glue(
      "data/private/mgi/{opt$mgi_version}/",
      "{opt$outcome}/time_restricted_phenomes/"
    ),
    recursive = TRUE
  )
}
### save mgi matching data
save_qs(
  x = mgi_post_match_cov,
  file = glue(
    "data/private/mgi/{opt$mgi_version}/",
    "{opt$outcome}/matched_covariates_r{opt$matching_ratio}.qs"
  )
)

## ukb
if (also_ukb) {
  cli_progress_step(glue(
    "performing 1:{opt$matching_ratio} case:non-case ",
    "matching in ukb..."
  ))

  ukb_post_match_cov <- perform_matching(
    matching_data   = ukb_matching_cov,
    case_data       = ukb_case_data$case,
    match_ratio     = opt$matching_ratio,
    nearest_vars    = nearest_matching_vars,
    exact_vars      = exact_matching_vars,
    match_caliper   = opt$matching_caliper,
    time_thresholds = time_thresholds
  )

  if (!dir.exists(glue(
    "data/private/ukb/{opt$ukb_version}/",
    "{opt$outcome}/",
    "time_restricted_phenomes/"
  ))) {
    dir.create(
      glue(
        "data/private/ukb/{opt$ukb_version}/",
        "{opt$outcome}/time_restricted_phenomes/"
      ),
      recursive = TRUE
    )
  }
  ## save ukb matching data
  save_qs(
    x = ukb_post_match_cov,
    file = glue(
      "data/private/ukb/{opt$ukb_version}/",
      "{opt$outcome}/matched_covariates_r{opt$matching_ratio}.qs"
    )
  )
}





# 8. create time-restricted phenomes -------------------------------------------
## mgi
cli_alert("constructing time-restricted phenomes in mgi...")
mgi_matched_phe <- merge.data.table(
  mgi_first_phe[id %in% mgi_post_match_cov[, id]],
  mgi_post_match_cov,
  by = "id"
)

mgi_pims <- map(
  time_thresholds,
  \(i) generate_restricted_phenome(
    phe_data    = mgi_matched_phe,
    threshold   = i,
    cases       = mgi_case_data$case[, id],
    outcome_phe = opt$outcome
  ),
  .progress = TRUE
) |> set_names(glue("t{time_thresholds}"))

mgi_pims <- map(
  seq_along(time_thresholds),
  \(i) {
    partition_data(
      data = mgi_pims[[i]],
      prop = c(0.5, 0.5),
      group_names = c("train", "test")
    )
  }
) |> set_names(glue("t{time_thresholds}"))

walk(
  seq_along(mgi_pims),
  \(i) {
    qsave(
      x = mgi_pims[[i]],
      file = glue(
        "data/private/mgi/{opt$mgi_version}/",
        "{opt$outcome}/time_restricted_phenomes/",
        "mgi_{opt$outcome}",
        "_{names(mgi_pims)[i]}_r{opt$matching_ratio}_{opt$mgi_version}.qs"
      )
    )
  }
)

## ukb
if (also_ukb) {
  cli_progress_step("constructing time-restricted phenomes in ukb...")
  ukb_matched_phe <- merge.data.table(
    ukb_first_phe[id %in% ukb_post_match_cov[, id]],
    ukb_post_match_cov,
    by = "id"
  )

  ukb_pims <- map(
    time_thresholds,
    \(i) generate_restricted_phenome(
      phe_data    = ukb_matched_phe,
      threshold   = i,
      cases       = ukb_case_data$case[, id],
      outcome_phe = opt$outcome
    )
  ) |> set_names(glue("t{time_thresholds}"))

  ukb_pims <- map(
    seq_along(time_thresholds),
    \(i) {
      partition_data(
        data = ukb_pims[[i]],
        prop = c(0.5, 0.5),
        group_names = c("train", "test")
      )
    }
  ) |> set_names(glue("t{time_thresholds}"))

  walk(
    seq_along(ukb_pims),
    \(i) {
      qsave(
        x = ukb_pims[[i]],
        file = glue(
          "data/private/ukb/{opt$ukb_version}/",
          "{opt$outcome}/time_restricted_phenomes/",
          "ukb_{opt$outcome}",
          "_{names(ukb_pims)[i]}_r{opt$matching_ratio}_{opt$ukb_version}.qs"
        )
      )
    }
  )
}
cli_progress_done()

# # 9. explore obesity-outcome relationship by threshold -------------------------
# exposure <- "EM_236"

# mgi_case_ex_plot <- mgi_pims |>
#     plot_quick_tv_logor(
#         outcome = opt$outcome,
#         exposure = exposure,
#         cohort = "MGI"
#     )
# ukb_case_ex_plot <- ukb_pims |>
#     plot_quick_tv_logor(
#         outcome = opt$outcome,
#         exposure = exposure,
#         cohort = "UKB"
#     )

# ggsave(
#     filename = glue(
#         "data/private/mgi/{opt$mgi_version}/",
#         "{opt$outcome}/time_restricted_phenomes/",
#         "mgi_{opt$outcome}_{exposure}_logor_r{opt$matching_ratio}_{opt$mgi_version}.pdf"
#     ),
#     plot   = mgi_case_ex_plot,
#     width  = 7,
#     height = 4
# )

# ggsave(
#     filename = glue(
#         "data/private/ukb/{opt$ukb_version}/",
#         "{opt$outcome}/time_restricted_phenomes/",
#         "ukb_{opt$outcome}_{exposure}_logor_r{opt$matching_ratio}_{opt$ukb_version}.pdf"
#     ),
#     plot   = ukb_case_ex_plot,
#     width  = 7,
#     height = 4
# )

cli_alert("script success! 🎉")
