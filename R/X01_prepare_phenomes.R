# prepare mgi and ukb data for time-restricted phecode-phecode phewas
# author:  max salvatore
# date:    20230816

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

ms::libri(ms, data.table, qs, MatchIt, optparse, glue)

set.seed(61787)

lapply(list.files("fn/", full.names = TRUE), source) |> # load functions
  invisible()

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "157",
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
    type = "character", default = "0,0.5,1,2,3,5",
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
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications (specifies outcome) --------------------------------------
time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])
nearest_matching_vars <- strsplit(opt$nearest_matching_var, ",")[[1]]
exact_matching_vars <- strsplit(opt$exact_matching_var, ",")[[1]]

# 3. extra preparations (outcome-specific) -------------------------------------
## pull file paths corresponding to the data version specified
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

## check that outcome exists in both datasets
cli_alert("checking that outcome exists in both phenomes...")
mgi_pim0 <- qread(file_paths[["mgi"]][["pim0_file"]])
ukb_pim0 <- qread(file_paths[["ukb"]][["pim0_file"]])

if (paste0("X", gsub("X", "", opt$outcome)) %in% names(mgi_pim0) &
  paste0("X", gsub("X", "", opt$outcome)) %in% names(ukb_pim0)) {
  cli_alert_success("outcome exists! procceeding...")
} else {
  cli_alert_danger()("outcome does not exist in both datasets! stopping...")
  stop(paste0(
    "outcome phecode ", opt$outcome, " does not exist in ",
    fcase(
      !(paste0("X", gsub("X", "", opt$outcome)) %in% names(mgi_pim0)) &
        !(paste0("X", gsub("X", "", opt$outcome)) %in% names(ukb_pim0)), "mgi and ukb",
      !(paste0("X", gsub("X", "", opt$outcome)) %in% names(mgi_pim0)), "mgi",
      !(paste0("X", gsub("X", "", opt$outcome)) %in% names(ukb_pim0)), "ukb"
    )
  ))
}

setnames(mgi_pim0,
  old = c("IID", paste0("X", opt$outcome)),
  new = c("id", "outcome")
)
setnames(ukb_pim0,
  old = paste0("X", gsub("X", "", opt$outcome)),
  new = "outcome"
)

## confirm file structure for a given outcome exists - if not, create paths
### mgi
check_folder_structure(
  cohort          = "mgi",
  data_version    = opt$mgi_version,
  outcome_phecode = opt$outcome
)

### ukb
check_folder_structure(
  cohort          = "ukb",
  data_version    = opt$ukb_version,
  outcome_phecode = opt$outcome
)

# 4. read data -----------------------------------------------------------------
## mgi
cli_alert("loading mgi data...")
### demographics
mgi_cov <- read_qs(glue("data/private/mgi/{opt$mgi_version}/data_{opt$mgi_version}_comb.qs"))
setnames(mgi_cov,
  old = c(
    "DeID_PatientID", "Age", "AliveYN", "Deceased_DaysSinceBirth",
    "Ethnicity", "MaritalStatusCode", "Sex", "Race", "AgeFirstEntry",
    "AgeLastEntry", "YearsInEHR", "FirstDaySinceBirth", "LastDaySinceBirth"
  ),
  new = c(
    "id", "age", "alive", "dead_dsb", "ethn", "marital",
    "sex", "race", "age_at_first_diagnosis", "age_at_last_diagnosis",
    "length_followup", "first_dsb", "last_dsb"
  )
)

### icd-phecode data
mgi_full_phe <- get(load(file_paths[["mgi"]][["phecode_dsb_file"]]))
if ("IID" %in% names(mgi_full_phe)) setnames(mgi_full_phe, "IID", "id")
if ("DaysSinceBirth" %in% names(mgi_full_phe)) setnames(mgi_full_phe, "DaysSinceBirth", "dsb")

mgi_first_phe <- mgi_full_phe[
  mgi_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
]

## ukb
cli_alert("loading ukb data...")
### demographics
ukb_demo <- qread(file_paths[["ukb"]][["demo_file"]])
ukb_demo <- ukb_demo[, .(id, age = age_at_consent, ethn = race_eth, sex)]
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

### icd-phecode data
ukb_full_phe <- qread(file_paths[["ukb"]][["icd_phecode_file"]])
ukb_first_phe <- ukb_full_phe[
  ukb_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
]

# 5. identify cases ------------------------------------------------------------
## mgi
mgi_case <- unique(mgi_first_phe[phecode == opt$outcome])
mgi_case_ids <- mgi_case[, unique(id)]

## ukb
ukb_case <- unique(ukb_first_phe[phecode == opt$outcome])[, id := as.character(id)]
ukb_case_ids <- ukb_case[, unique(id)]

# 6. calculate diagnostic metrics ----------------------------------------------
## mgi
cli_alert("calculating diagnostic metrics in mgi...")
mgi_cov[, `:=`(
  female = as.numeric(sex == "F"),
  case = fifelse(id %in% mgi_case_ids, 1, 0)
)]

## ukb
cli_alert("calculating diagnostic metrics in ukb...")
ukb_diag_metrics <- get_icd_phecode_metrics(
  full_phe_data = ukb_first_phe
)[, id := as.character(id)]
ukb_matching_cov <- merge.data.table(
  ukb_demo[, .(id, age, female = as.numeric(sex == "Female"))],
  ukb_diag_metrics,
  by = "id"
)[, case := fifelse(id %in% ukb_case_ids, 1, 0)]

# 7. perform matching ----------------------------------------------------------
## mgi
cli_alert(glue(
  "performing 1:{opt$matching_ratio} case:non-case ",
  "matching in mgi..."
))
mgi_match_text <- glue(
  "MatchIt::matchit(case ~ ",
  "{glue_collapse(c(nearest_matching_vars, ",
  "exact_matching_vars), sep = ' + ')}, ",
  "data = mgi_cov, calclosest = TRUE, ",
  "mahvars = c({paste0(sapply(nearest_matching_vars, ",
  "\\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
  "caliper = {opt$matching_caliper}, ",
  "exact = c({paste0(sapply(exact_matching_vars, ",
  "\\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
  "ratio = {opt$matching_ratio})"
)
mgi_match <- eval(parse(text = mgi_match_text))
mgi_matched <- MatchIt::match.data(mgi_match)
mgi_post_match_cov <- merge.data.table(
  mgi_matched,
  merge.data.table(
    mgi_case,
    mgi_matched,
    by = "id",
    all.x = TRUE
  )[, .(subclass, case_dsb = dsb)],
  by = "subclass"
)

lapply(
  time_thresholds,
  \(i) {
    mgi_post_match_cov[, (glue("t{i}_threshold")) :=
      floor(case_dsb - (365.25 * i))]
  }
) |> invisible()

lapply(
  time_thresholds,
  \(i) {
    mgi_post_match_cov[, (glue("t{i}_indicator")) :=
      as.numeric(all(get(glue("t{i}_threshold")) > first_dsb)),
    by = subclass
    ]
  }
) |> invisible()


if (!dir.exists(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "X{gsub('X','', opt$outcome)}/",
  "time_restricted_phenomes/"
))) {
  dir.create(
    glue(
      "data/private/mgi/{opt$mgi_version}/",
      "X{gsub('X','', opt$outcome)}/time_restricted_phenomes/"
    ),
    recursive = TRUE
  )
}
### save mgi matching data
save_qs(
  x = mgi_post_match_cov,
  file = glue(
    "data/private/mgi/{opt$mgi_version}/",
    "X{gsub('X', '', opt$outcome)}/matched_covariates.qs"
  )
)

## ukb
cli_alert(glue(
  "performing 1:{opt$matching_ratio} case:non-case ",
  "matching in ukb..."
))
ukb_match_text <- glue(
  "matchit(case ~ ",
  "{glue_collapse(c(nearest_matching_vars, ",
  "exact_matching_vars), sep = ' + ')}, ",
  "data = ukb_matching_cov, calclosest = TRUE, ",
  "mahvars = c({paste0(sapply(nearest_matching_vars,",
  "\\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
  "caliper = {opt$matching_caliper}, ",
  "exact = c({paste0(sapply(exact_matching_vars,",
  "\\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
  "ratio = {opt$matching_ratio})"
)
ukb_match <- eval(parse(text = ukb_match_text))
ukb_matched <- match.data(ukb_match)
ukb_post_match_cov <- merge.data.table(
  ukb_matched,
  merge.data.table(
    ukb_case,
    ukb_matched,
    by = "id",
    all.x = TRUE
  )[, .(subclass, case_dsb = dsb)],
  by = "subclass"
)[, id := as.character(id)]
for (i in time_thresholds) {
  ukb_post_match_cov[, (glue("t{i}_threshold")) :=
    floor(case_dsb - (365.25 * i))]
}
for (i in time_thresholds) {
  ukb_post_match_cov[, (glue("t{i}_indicator")) :=
    as.numeric(all(get(glue("t{i}_threshold")) > first_dsb)),
  by = subclass
  ]
}
if (!dir.exists(glue(
  "data/private/ukb/{opt$ukb_version}/",
  "X{gsub('X','', opt$outcome)}/",
  "time_restricted_phenomes/"
))) {
  dir.create(
    glue(
      "data/private/ukb/{opt$ukb_version}/",
      "X{gsub('X','', opt$outcome)}/time_restricted_phenomes/"
    ),
    recursive = TRUE
  )
}
### save ukb matching data
save_qs(
  x = ukb_post_match_cov,
  file = glue(
    "data/private/ukb/{opt$ukb_version}/",
    "X{gsub('X', '', opt$outcome)}/matched_covariates.qs"
  )
)

# 8. create time-restricted phenomes -------------------------------------------
## mgi
cli_alert("constructing time-restricted phenomes in mgi...")
mgi_matched_phe <- merge.data.table(
  mgi_first_phe[id %in% mgi_post_match_cov[, id]],
  mgi_post_match_cov,
  by = "id"
)
mgi_pims <- lapply(
  seq_along(time_thresholds),
  \(i) generate_restricted_phenome(
    phe_data = mgi_matched_phe,
    threshold = time_thresholds[i],
    cases = mgi_case_ids,
    outcome_phe = opt$outcome
  )
)
names(mgi_pims) <- glue("t{time_thresholds}")

for (i in 1:length(mgi_pims)) {
  save_qs(
    x = mgi_pims[[i]],
    file = glue(
      "data/private/mgi/{opt$mgi_version}/",
      "X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/",
      "mgi_X{gsub('X', '', opt$outcome)}",
      "_{names(mgi_pims)[i]}_{opt$mgi_version}.qs"
    )
  )
}

## ukb
cli_alert("constructing time-restricted phenomes in ukb...")
ukb_matched_phe <- merge.data.table(
  ukb_first_phe[id %in% ukb_post_match_cov[, id]],
  ukb_post_match_cov,
  by = "id"
)
ukb_pims <- list()
for (i in seq_along(time_thresholds)) {
  ukb_pims[[i]] <- generate_restricted_phenome(
    phe_data = ukb_matched_phe,
    threshold = time_thresholds[i],
    cases = ukb_case_ids,
    outcome_phe = opt$outcome
  )
}
names(ukb_pims) <- glue("t{time_thresholds}")
for (i in 1:length(ukb_pims)) {
  save_qs(
    x = ukb_pims[[i]],
    file = glue(
      "data/private/ukb/{opt$ukb_version}/",
      "X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/",
      "ukb_X{gsub('X', '', opt$outcome)}",
      "_{names(ukb_pims)[i]}_{opt$ukb_version}.qs"
    )
  )
}

cli_alert("script success! 🎉")
