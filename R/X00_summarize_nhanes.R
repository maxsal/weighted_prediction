
# Calculate simplex regression-based inverse probability weights using NHANES
# data and postratificatiod-based weights using Census, CDC, and SEER data in
# MGI
# author:   max salvatore
# date:     20230418

# libraries --------------------------------------------------------------------
ms::libri(
    ms, haven, survey, data.table, glue, optparse, cli, qs
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
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

data_path <- glue("data/private/mgi/{opt$cohort_version}/")

for (i in list.files("fn", full.names = TRUE)) source(i)

nhanes_datasets <- unlist(strsplit(opt$nhanes_survey_names, ","))
nhanes_merged <- download_nhanes_data(
  wave_letter = opt$nhanes_wave_letter,
  wave_years  = opt$nhanes_wave_years,
  datasets    = nhanes_datasets
)

keep_vars <- c(
  "SEQN", "RIAGENDR", "WTINT2YR", "RIDAGEYR", "RIDRETH1", "RIDRETH3", "MCQ220",
  "BMXBMI", "SMQ040", "SMQ020", "DIQ010", "MCQ160C", "WTMEC2YR",
  "SDMVSTRA", "SDMVPSU", paste0("DPQ0", 1:9, "0"), "BPQ020", "LBXHIVC"
)

if ("WTMECPRP" %in% names(nhanes_merged)) {
  setnames(
    nhanes_merged,
    "WTMECPRP",
    "WTMEC2YR"
  )
}
if ("WTINTPRP" %in% names(nhanes_merged)) {
  setnames(
    nhanes_merged,
    "WTINTPRP",
    "WTINT2YR"
  )
}

nhanes_merged <- nhanes_merged[, ..keep_vars]

prepped_nhanes <- prepare_nhanes_data(
  nhanes_data = nhanes_merged,
  mec_wt_var = "WTMEC2YR"
)

## age categories
age_cats <- seq(0, 80, 10)
tmp_age_cat <- age_grp_table(
    lower_ages   = age_cats,
    num_vec      = rep(NA, length(age_cats)),
    num_var_name = "counts"
)
prepped_nhanes[, `:=`(
    age_verbose = factor(cut(age,
        breaks = c(0, tmp_age_cat[["upper"]]),
        labels = tmp_age_cat[["group"]],
        right  = FALSE
    ),
    levels = tmp_age_cat[["group"]]
    )
)]

# estimate weighted NHANES
dsn <- svydesign(
    id = ~psu_nhanes,
    strata = ~strata_nhanes,
    weights = ~weight_nhanes,
    data = prepped_nhanes,
    nest = TRUE
)

(nhanes_demo_summary <- weighted_summary_wrapper(
    prepped_nhanes,
    dsn = dsn,
    vars = c(
        "age", "age_verbose", "sex", "race_eth", "bmi", "bmi_cat",
        "cancer", "diabetes", "cad", "depression", "smoking_status"
    )
))

fwrite(
    x    = nhanes_demo_summary,
    file = glue("data/public/nhanes_weighted_demo_summary.txt"),
    sep  = "\t"
)
