# Calculate simplex regression-based inverse probability weights using NHANES
# data and postratificatiod-based weights using Census, CDC, and SEER data in
# MGI
# author:   max salvatore
# date:     20230809

# libraries --------------------------------------------------------------------
ms::libri(
  maxsal/ms, haven, survey, dplyr, pracma, simplexreg, data.table,
  glue, cli, qs, optparse
)

# optparse list ---
option_list <- list(
  make_option("--cohort_version",
    type = "character", default = "20220822",
    help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = '20220822']"
  ),
  make_option("--mgi_cohort",
    type = "character", default = "comb",
    help = "Specific MGI cohort [default = %default]"
  ),
  make_option("--nhanes_wave_letter",
    type = "character",
    default = "J",
    help = glue(
      "NHANES data prefix corresponding to wave ",
      "[default = %default]"
    )
  ),
  make_option("--nhanes_wave_years",
    type = "character",
    default = "2017-2018",
    help = glue(
      "NHANES wave years corresponding to wave ",
      "[default = %default]"
    )
  ),
  make_option("--nhanes_survey_names",
    type = "character",
    default = "DEMO,BMX,SMQ,DIQ,MCQ,DPQ,BPQ,HIV",
    help = glue(
      "NHANES wave years corresponding to wave ",
      "[default = %default]"
    )
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

data_path <- glue("data/private/mgi/{opt$cohort_version}/")

for (i in list.files("fn", full.names = TRUE)) source(i)

# load data --------------------------------------------------------------------
cli_alert("loading data...")
mgi <- read_qs(glue("{data_path}data_{opt$cohort_version}_{opt$mgi_cohort}.qs"))[, `:=` (
  smoking_current = as.numeric(SmokingStatus == "Current"),
  smoking_former  = as.numeric(SmokingStatus == "Former")
  )]
setnames(mgi, "DeID_PatientID", "id", skip_absent = TRUE)

# nhanes_datasets <- unlist(strsplit(opt$nhanes_survey_names, ","))
# nhanes_merged <- download_nhanes_data(
#   wave_letter = opt$nhanes_wave_letter,
#   wave_years  = opt$nhanes_wave_years,
#   datasets    = nhanes_datasets
# )

# keep_vars <- c(
#   "SEQN", "RIAGENDR", "WTINT2YR", "RIDAGEYR", "RIDRETH1", "RIDRETH3", "MCQ220",
#   "BMXBMI", "SMQ040", "SMQ020", "DIQ010", "MCQ160C", "WTMEC2YR",
#   "SDMVSTRA", "SDMVPSU", paste0("DPQ0", 1:9, "0"), "BPQ020", "LBXHIVC"
# )

# if ("WTMECPRP" %in% names(nhanes_merged)) {
#   setnames(
#     nhanes_merged,
#     "WTMECPRP",
#     "WTMEC2YR"
#   )
# }
# if ("WTINTPRP" %in% names(nhanes_merged)) {
#   setnames(
#     nhanes_merged,
#     "WTINTPRP",
#     "WTINT2YR"
#   )
# }

# nhanes_merged <- nhanes_merged[, ..keep_vars]

prepped_nhanes <- prepare_nhanes_data()
prepped_nhanes[, smoker := fcase(
  smoking_status %in% c("Current", "Former"), 1,
  smoking_status == "Never", 0,
  default = NA
)]

stacked <- rbindlist(list(
  prepped_nhanes,
  mgi[][, dataset := "MGI"]
), use.names = TRUE, fill = TRUE)

# estimate ipw and postratification weights ------------------------------------
cli_alert("estimating ipw weights...")

ip_weights_list <- list(
  "simple" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
               "smoker",
               "bmi_under", "bmi_overweight", "bmi_obese", "nhw"),
  "simple_f" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
                 "smoker", "bmi_under",
                 "bmi_overweight", "bmi_obese", "nhw", "female"),
  "selection" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
                  "cad", "diabetes", "smoker",
                  "bmi_under", "bmi_overweight", "bmi_obese", "nhw"),
  "selection_c" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
                  "cad", "diabetes", "smoker",
                  "bmi_under", "bmi_overweight", "bmi_obese", "nhw",
                  "cancer"),
  "selection_f" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
                  "cad", "diabetes", "smoker",
                  "bmi_under", "bmi_overweight", "bmi_obese", "nhw",
                  "female", "cancer"),
  "cancer" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
               "smoker", "bmi_under",
               "bmi_overweight", "bmi_obese", "nhw", "cancer"),
  "depression" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
                   "smoker", "bmi_under",
                   "bmi_overweight", "bmi_obese", "nhw", "depression"),
  "cad" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
            "smoker", "bmi_under",
            "bmi_overweight", "bmi_obese", "nhw", "cad"),
  "diabetes" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
                 "smoker", "bmi_under",
                 "bmi_overweight", "bmi_obese", "nhw", "diabetes"),
  "hypertension" = c("as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
                     "smoker", "bmi_under",
                     "bmi_overweight", "bmi_obese", "nhw", "hypertension")
)

ip_weights <- list()
cli_progress_bar(name = "estimating ipw weights", total = length(names(ip_weights_list)))
for (i in seq_along(names(ip_weights_list))) {
  tmp_weights <- ipw(stacked_data = stacked,
                     covs = ip_weights_list[[names(ip_weights_list)[i]]])
  setnames(tmp_weights, "ip_weight", paste0("ip_", names(ip_weights_list)[i]))
  ip_weights[[i]] <- tmp_weights
  cli_progress_update()
}
cli_progress_done()

ipws <- ms::merge_list(ip_weights, join_fn = dplyr::full_join)

cli_alert("estimating poststratification weights...")
post_weights_list <- list(
  "selection"    = c("smoker", "cad", "diabetes"),
  "selection_c"  = c("smoker", "cad", "diabetes", "cancer"),
  "selection_f"  = c("smoker", "cad", "diabetes", "cancer", "female"),
  "nhw"          = c("smoker", "cad", "diabetes", "cancer", "nhw"),
  "nhw_f"        = c("smoker", "cad", "diabetes", "cancer", "nhw", "female"),
  "depression"   = c("smoker", "cad", "diabetes", "cancer", "depression"),
  "diabetes"     = c("smoker", "cad", "diabetes", "cancer", "diabetes"),
  "hypertension" = c("smoker", "cad", "diabetes", "cancer", "hypertension")
)

post_weights <- list()
cli_progress_bar(name = "estimating poststratification weights", total = length(names(post_weights_list)))
for (i in seq_along(names(post_weights_list))) {
  tmp_weights <- poststrat_nhanes(
    int_data       = mgi,
    nhanes_data    = prepped_nhanes,
    chop           = TRUE,
    age_bin        = TRUE,
    age_bin_breaks = c(seq(0, 80, 10), 150),
    covs           = post_weights_list[[names(post_weights_list)[i]]]
  )
  setnames(tmp_weights, "ps_weight", paste0("ps_", names(post_weights_list)[i]))
  post_weights[[i]] <- tmp_weights
  cli_progress_update()
}
cli_progress_done()

posts <- merge_list(post_weights, join_fn = dplyr::full_join)

weights <- merge.data.table(ipws, posts, by = "id")
merged  <- merge_list(list(mgi, weights), by = "id", join_fn = dplyr::full_join)

# estimate cancer~female log(OR) -----------------------------------------------
cli_alert("estimating cancer~female log(OR) as sanity check...")

extract_estimates <- function(
  outcome,
  exposure,
  data,
  weights = NULL
) {
  if (is.null(weights)) {
    mod <- glm(
      as.formula(paste0(outcome, " ~ ", exposure)),
      family  = "quasibinomial",
      data    = data)
  } else {
    tmp_design <- svydesign(
      id      = ~1,
      weights = ~get(weights),
      data    = data[!is.na(get(weights)), ]
    )
    mod <- svyglm(
      as.formula(paste0(outcome, " ~ ", exposure)),
      family = "quasibinomial",
      design = tmp_design
    )
  }

  est     <- summary(mod)$coef[exposure, 1:2]
  beta_lo <- est[[1]] - qnorm(0.975) * est[[2]]
  beta_hi <- est[[1]] + qnorm(0.975) * est[[2]]
  rare <- ifelse(sum(mod$data[[outcome]], na.rm = TRUE) /
    length(na.omit(mod$data[[outcome]])) >= 0.15, FALSE, TRUE)
  evals <- suppressMessages(evalues.OR(
    est  = exp(est[[1]]),
    lo   = exp(beta_lo),
    hi   = exp(beta_hi),
    rare = rare
  ))
  data.table(
    outcome  = outcome,
    exposure = exposure,
    weights  = ifelse(is.null(weights), "None", weights),
    beta     = est[[1]],
    beta_lo  = beta_lo,
    beta_hi  = beta_hi,
    or       = exp(est[[1]]),
    or_lo    = exp(beta_lo),
    or_hi    = exp(beta_hi),
    eval     = evals[2, 1],
    eval_lo  = evals[2, 2],
    eval_hi  = evals[2, 3]
  )
}

est_out <- list()
weight_vars <- names(weights)[names(weights) != "id"]
cli_progress_bar(name = "estimating cancer~female log(OR)", total = length(weight_vars))
for (i in seq_along(weight_vars)) {
  est_out[[i]] <- extract_estimates(
    outcome  = "cancer",
    exposure = 'as.numeric(Sex == "F")',
    data     = merged,
    weights  = weight_vars[i]
  )
  cli_progress_update()
}
cli_progress_done()

est_out[[i+1]] <- extract_estimates(
        outcome  = "cancer",
        exposure = 'as.numeric(Sex == "F")',
        data     = merged
      )

est_out <- rbindlist(est_out)

fwrite(
  x    = est_out,
  file = glue("{data_path}cancer_female_logor_est_{opt$cohort_version}_{opt$mgi_cohort}.txt"),
  sep  = "\t"
)

# save -------------------------------------------------------------------------
save_qs(
  x    = weights,
  file = glue("{data_path}weights_{opt$cohort_version}_{opt$mgi_cohort}.qs")
)

cli_alert_success("script success! see {.path {data_path}} and suffix {opt$mgi_cohort} 🥳")
