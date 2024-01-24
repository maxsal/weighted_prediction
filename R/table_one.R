# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable

# libraries, functions, and options --------------------------------------------
ms::libri(
    data.table, tidyverse, glue, qs, optparse, ms, cli
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--mgi_version", type = "character", default = "20220822",
              help = "Version of MGI data [default = %default]"),
  make_option("--mgi_cohort", type = "character", default = "comb",
              help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# read data --------------------------------------------------------------------
cli_alert("loading mgi data...")
## mgi
### demographics
mgi_cov <- read_qs(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))
setnames(mgi_cov,
         old = c("DeID_PatientID", "Age", "AliveYN",
                 "Ethnicity", "MaritalStatusCode", "Sex", "Race",
                 "YearsInEHR", "FirstDaySinceBirth", "LastDaySinceBirth"),
         new = c("id", "age", "alive", "ethn", "marital",
                 "sex", "race", "length_followup", "first_dsb", "last_dsb"),
         skip_absent = TRUE
        )

### icd-phecode data
mgi_full_phe <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs"))
if ("IID" %in% names(mgi_full_phe)) setnames(mgi_full_phe, "IID", "id")
if ("DaysSinceBirth" %in% names(mgi_full_phe)) setnames(mgi_full_phe, "DaysSinceBirth", "dsb")
mgi_full_phe <- mgi_full_phe[id %in% mgi_cov[, unique(id)], ]

cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv", showProgress = FALSE)

mgi_cancer       <- mgi_full_phe[phecode %in% cancer_phecodesx[keep == 1 | specific == 1, phecode], ]
mgi_first_cancer <- mgi_cancer[ mgi_cancer[, .I[dsb == min(dsb)], by = c("id")][["V1"]], ]

cancer_outcomes <- c("CA_101.1", "CA_101.41", "CA_101.6", "CA_101.8")
mgi_cancer_outcome_counts <- map(
    cancer_outcomes,
    function(x) {
        mgi_first_cancer[phecode == x, .N]
    }
) |> set_names(cancer_outcomes)

### weights data
mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs"))

mgi <- merge.data.table(
    mgi_cov,
    mgi_weights,
    by = "id",
    all.x = TRUE
)

## ukb
cli_alert("loading ukb data...")
### demographics
ukb_demo <- read_qs(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))
ukb_demo <- ukb_demo[, .(
  id,
  dob  = birth_date,
  age  = age_at_consent,
  race_eth,
  race_eth2 = fcase(race_eth == "White", "White", race_eth == "Black", "Black", race_eth == "Asian", "Asian", race_eth %in% c("Other", "Unknown", NA), "Other/Unknown"),
  sex,
  smoker = fcase(smk_status %in% c("Previous", "Current"), 1, smk_status == "Never", 0, default = NA),
  bmi = bmi_med,
  bmi_cat = bmi_med_cat,
  drinker = as.numeric(alc_ev == "Ever"),
  hypertension = hypertensionx,
  cancer, diabetes, cad, anxiety, depression, in_phenome)]

### icd-phecode data
ukb_full_phe <- read_qs(glue("data/private/ukb/{opt$ukb_version}/UKB_FULL_PHECODEX_DSB_{opt$ukb_version}.qs"))
if ("IID" %in% names(ukb_full_phe)) setnames(ukb_full_phe, "IID", "id")
if ("DaysSinceBirth" %in% names(ukb_full_phe)) setnames(ukb_full_phe, "DaysSinceBirth", "dsb")
ukb_full_phe[, dsb := as.numeric(dsb)]

ukb_cancer <- ukb_full_phe[phecode %in% cancer_phecodesx[keep == 1 | specific == 1, phecode], ]
ukb_first_cancer <- ukb_cancer[ukb_cancer[, .I[dsb == min(dsb)], by = c("id")][["V1"]], ]

ukb_cancer_outcome_counts <- map(
    cancer_outcomes,
    function(x) {
        ukb_first_cancer[phecode == x, .N]
    }
) |> set_names(cancer_outcomes)

### weights data
ukb_weights <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/UKBSelectionWeights.tab",
    colClasses = "character")[, .(id = f.eid, ip_weight = as.numeric(LassoWeight))]

ukb <- merge.data.table(
    ukb_demo[id %in% unique(ukb_full_phe[, id]), ],
    ukb_weights,
    by = "id",
    all.x = TRUE
)

## other
cancer_phecodes <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv",
                         colClasses = "character")[keep == 1, phecode]

# prep data --------------------------------------------------------------------
## age categories
age_cats    <- seq(0, 80, 10)
tmp_age_cat <- age_grp_table(
    lower_ages   = age_cats,
    num_vec      = rep(NA, length(age_cats)),
    num_var_name = "counts"
)
## mgi
mgi[, `:=` (
    race_eth = fcase(
        race == "Caucasian" & ethn == "Non-Hispanic", "NH White",
        race == "African American" & ethn == "Non-Hispanic", "NH Black",
        race == "Asian" & ethn == "Non-Hispanic", "NH Asian",
        ethn == "Hispanic", "Hispanic",
        default = "Other/Unknown"
    ),
    SmokingStatus = relevel(factor(SmokingStatus, levels = c("Never", "Past", "Current", "Unknown")), ref = "Never"),
    bmi_verbose   = factor(fcase(
        bmi_cat == 1, "Underweight (<18.5)",
        bmi_cat == 2, "Healthy [18.5, 25)",
        bmi_cat == 3, "Overweight [25, 30)",
        bmi_cat == 4, "Obese [30+)"
    ), levels = c("Underweight (<18.5)", "Healthy [18.5, 25)", "Overweight [25, 30)", "Obese [30+)")),
    age_verbose = factor(cut(age_at_last_diagnosisx,
                              breaks = c(0, tmp_age_cat[["upper"]]),
                              labels = tmp_age_cat[["group"]],
                              right  = FALSE),
                          levels = tmp_age_cat[["group"]])
)]
## ukb
ukb[, `:=` (
    age_verbose = factor(cut(age,
                            breaks = c(0, tmp_age_cat[["upper"]]),
                            labels = tmp_age_cat[["group"]],
                            right  = FALSE),
                        levels = tmp_age_cat[["group"]])
                          )
    ]

# unweighted demographics summary ----------------------------------------------
mgi_demo_summary <- mgi[id %in% unique(mgi_full_phe[, id]), ] |>
    summarizer(col_names = c("age_at_last_diagnosisx", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
                             "cancerx", "diabetesx", "anxietyx", "depressionx", "hypertensionx"))

mgi_cancer_outcome_summary <- map(
    seq_along(mgi_cancer_outcome_counts),
    \(x) {
        data.table(
            val   = names(mgi_cancer_outcome_counts)[x],
            print = paste0(
                format(round(mgi_cancer_outcome_counts[[x]] * 100 / mgi_demo_summary[variable == "N", as.numeric(gsub(",", "", print))], 1), nsmall = 1),
                " (", format(mgi_cancer_outcome_counts[[x]], big.mark = ","), ")"
            ),
            print_unit = "% (n)"
        )
    }
) |> bind_rows() |>
    left_join(
        ms::pheinfox[, .(val = phecode, variable = description)],
        by = "val"
    ) |>
    select(variable, everything())

mgi_demo_summary <- bind_rows(
    mgi_demo_summary,
    mgi_cancer_outcome_summary
)

ukb_demo_summary <- ukb[id %in% unique(ukb_full_phe[, id]), ] |>
    summarizer(col_names = c("age", "age_verbose", "sex", "race_eth2", "bmi", "bmi_cat",
                             "cancer", "diabetes", "anxiety", "depression", "hypertension"))

ukb_cancer_outcome_summary <- map(
    seq_along(ukb_cancer_outcome_counts),
    \(x) {
        data.table(
            val = names(ukb_cancer_outcome_counts)[x],
            print = paste0(
                format(round(ukb_cancer_outcome_counts[[x]] * 100 / ukb_demo_summary[variable == "N", as.numeric(gsub(",", "", print))], 1), nsmall = 1),
                " (", format(ukb_cancer_outcome_counts[[x]], big.mark = ","), ")"
            ),
            print_unit = "% (n)"
        )
    }
) |>
    bind_rows() |>
    left_join(
        ms::pheinfox[, .(val = phecode, variable = description)],
        by = "val"
    ) |>
    select(variable, everything())

ukb_demo_summary <- bind_rows(
    ukb_demo_summary,
    ukb_cancer_outcome_summary
)

# weighted demographics summary ------------------------------------------------
(mgi_demo_summary_ip <- weighted_summary_wrapper(
    mgi,
    weight = "ip_selection",
    vars = c(
        "age_at_last_diagnosisx", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
        "cancerx", "diabetesx", "anxietyx", "depressionx", "hypertensionx"
    )
))
(mgi_demo_summary_ps <- weighted_summary_wrapper(
    mgi,
    weight = "ps_selection",
    vars = c(
        "age_at_last_diagnosisx", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
        "cancerx", "diabetesx", "anxietyx", "depressionx", "hypertensionx"
    )
))
(ukb_demo_summary_w <- weighted_summary_wrapper(
    data = ukb,
    weight = "ip_weight",
    vars = c(
        "age", "age_verbose", "sex", "race_eth2", "bmi", "bmi_cat",
        "cancer", "diabetes", "cad", "anxiety", "depression", "hypertension"
    )
))

# summarize ehr data -----------------------------------------------------------
(mgi_ehr_summary <- phecode_dsb_summarizer(mgi_full_phe))

(ukb_ehr_summary <- phecode_dsb_summarizer(ukb_full_phe))

# stack summaries --------------------------------------------------------------
(mgi_stacked_summary <- rbindlist(list(
    mgi_demo_summary,
    mgi_ehr_summary
), use.names = TRUE, fill = TRUE))
(ukb_stacked_summary <- rbindlist(list(
    ukb_demo_summary,
    ukb_ehr_summary
), use.names = TRUE, fill = TRUE))

# save -------------------------------------------------------------------------
fwrite(
    x    = mgi_stacked_summary,
    file = glue("data/private/mgi/{opt$mgi_version}/mgi_demox_ehr_summary.txt"),
    sep  = "\t"
)
fwrite(
    x    = mgi_demo_summary_ip,
    file = glue("data/private/mgi/{opt$mgi_version}/mgi_demox_summary_ip.txt"),
    sep  = "\t"
)
fwrite(
    x    = mgi_demo_summary_ps,
    file = glue("data/private/mgi/{opt$mgi_version}/mgi_demox_summary_ps.txt"),
    sep  = "\t"
)
fwrite(
    x    = ukb_stacked_summary,
    file = glue("data/private/ukb/{opt$ukb_version}/ukb_demox_ehr_summary.txt"),
    sep  = "\t"
)
fwrite(
    x    = ukb_demo_summary_w,
    file = glue("data/private/ukb/{opt$ukb_version}/ukb_demox_summary_w.txt"),
    sep  = "\t"
)

cli_alert_success("done! 🎉")
