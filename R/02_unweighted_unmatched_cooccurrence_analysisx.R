# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable
# author:   max salvatore

# libraries, functions, and options --------------------------------------------
ms::libri(
  data.table, ms, MatchIt, logistf, glue, qs, optparse, EValue, cli,
  tidyverse, gtsummary
)

set.seed(61787)

for (i in list.files("./fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "CA_105.1",
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
cooccur_covs    <- c("age", "female", "length_followup")

## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)





# read data --------------------------------------------------------------------

## helper data
### codes with at least 20 cases in AOU, MGI, and UKB
common_codes     <- fread("data/public/phecodex_20plus.csv")[plus20 == 1, phecode]
### qualifying cancer phecodes (malignant)
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")

## mgi
cli_progress_step("Reading MGI data")
### demographics
mgi_demo <- qread(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))
ukb_demo <- qread(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))




### time-stamped phenome data
mgi_full_phe <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs"))[phecode %in% common_codes, ]
ukb_full_phe <- qread(glue("data/private/ukb/{opt$ukb_version}/UKB_FULL_PHECODEX_DSB_{opt$ukb_version}.qs"))[phecode %in% common_codes, ]


mgi_first_phe <- mgi_full_phe[
    mgi_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
]

# identify cases and exclusions
mgi_first_cancer_dsb <- mgi_first_phe[phecode %in% cancer_phecodesx, ] # first occurrence of each cancer
mgi_case             <- mgi_first_cancer_dsb[
    mgi_first_cancer_dsb[, .I[dsb == min(dsb, na.rm = TRUE)], id]$V1
    ][phecode == opt$outcome, ] # first occurrence of any cancer

# mgi_case     <- unique(mgi_first_cancer[get(paste0(opt$outcome, "_first")) == 1, ])
mgi_case         <- unique(mgi_first_phe[phecode == opt$outcome, ])
mgi_case_ids     <- mgi_case[, unique(id)]
mgi_other_cancer <- setdiff(mgi_first_cancer_dsb[, id], mgi_case_ids)

# truncate phenome
## identify index dsb in cases
mgi_case_index <- mgi_first_cancer_dsb[phecode == opt$outcome, ][, .(id, index_dsb = dsb)]

mgi_case_phe <- merge.data.table(
    mgi_first_phe[id %in% mgi_case_ids, ],
    mgi_case_index,
    by = "id"
)

for (i in seq_along(time_thresholds)) {
    mgi_case_phe[[paste0("t", time_thresholds[i], "_threshold")]] <- mgi_case_phe$index_dsb - (time_thresholds[i] * 365.25)
    mgi_case_phe[[paste0("t", time_thresholds[i], "_indicator")]] <- as.numeric(
        mgi_case_phe$dsb <= 
        mgi_case_phe[[paste0("t", time_thresholds[i], "_threshold")]]
        )
}

mgi_case_pims <- lapply(
    seq_along(time_thresholds),
    \(i) generate_restricted_phenome(
        phe_data    = mgi_case_phe,
        threshold   = time_thresholds[i],
        cases       = mgi_case_ids,
        outcome_phe = opt$outcome
    )
)

# identify whether there is sex-restriction for controls based on outcome
outcome_female_only <- ms::phecodex_sex[phecode == opt$outcome, female_only]
outcome_male_only   <- ms::phecodex_sex[phecode == opt$outcome, male_only]

if (outcome_female_only) {
    control_ids <- mgi_cov[!(id %in% c(mgi_case_ids, mgi_other_cancer)) & sex == "F", id]
} else if (outcome_male_only) {
    control_ids <- mgi_cov[!(id %in% c(mgi_case_ids, mgi_other_cancer)) & sex == "M", id]
} else {
    control_ids <- mgi_cov[!(id %in% c(mgi_case_ids, mgi_other_cancer)), id]
}


mgi_control_pim <- dcast(
    mgi_first_phe[id %in% control_ids, ],
    id ~ phecode,
    value.var = "phecode",
    fun.aggregate = length,
    fill = 0
)

mgi_tr_pims <- lapply(
    seq_along(time_thresholds),
    \(i) {
        tmp <- rbindlist(
            list(
                mgi_case_pims[[i]],
                mgi_control_pim[, `:=`(
                    case = 0
                )]
            ),
            use.names = TRUE,
            fill      = TRUE
        )
        tmp[is.na(tmp)] <- 0
        tmp
    }
)

# standardize
mgi_tr_pims <- lapply(
  mgi_tr_pims,
  \(x) {
    vars <- names(x)[!names(x) %in% c("id", "case")]
    standardize_variables(x, cols = vars)
  }
)

mgi_case_sum <- mgi_case_phe[ mgi_case_phe[, .I[which.min(dsb)], by = c("id")]$V1 ]

mgi_t_covs <- lapply(
    seq_along(time_thresholds),
    \(i) {
        mgi_case_sum[get(paste0("t", time_thresholds[i], "_indicator")) == 1, ][,
        `:=` (
                age             = round(dsb / 365.25, 1),
                length_followup = round((get(paste0("t", time_thresholds[i], "_threshold")) - dsb) / 365.25, 3)
        )][, .(id, age, length_followup)] |>
        merge.data.table(
            mgi_cov[, .(id, sex)],
            by    = "id",
            all.x = TRUE
        )
    }
)

mgi_tr_merged <- lapply(
    seq_along(time_thresholds),
    \(i) {
        merge.data.table(
            mgi_tr_pims[[i]],
            rbindlist(list(
                mgi_t_covs[[i]],
                mgi_cov[id %in% control_ids, .(id, age, sex, length_followup)]
            ), use.names = TRUE, fill = TRUE)
        )
    }
)

names(mgi_tr_merged) <- glue("t{time_thresholds}_threshold")

## phenome
pheinfo <- ms::pheinfox
if (pheinfo[phecode == opt$outcome, sex] != "Both") {
  if (any(c("sex", "female", "male") %in% tolower(cooccur_covs))) {
    cooccur_covs <- cooccur_covs[-which(tolower(cooccur_covs) %in% c("sex", "female", "male"))]
  }
}


logistf(case ~ BI_160 + age + length_followup, data = mgi_tr_merged[[1]]) |> summary()
glm(case ~ BI_160 + age + length_followup,
    data = mgi_tr_merged[[1]], method = brglm_fit, family = "binomial",
    type = "AS_mean") |> tidy()

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
    n_cores            = parallelly::availableCores() / 4,
    parallel           = TRUE,
    verbose            = TRUE
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
    parallel           = TRUE
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
      "{opt$mgi_version}_results.qs"
    )
  )
}

## mgi pwide
qsave(
  mgi_pwide_results_post_cor,
  file = glue(
    "results/mgi/{opt$mgi_version}/{opt$outcome}/",
    "mgi_{opt$outcome}_pwide_{opt$mgi_version}_results_post_cor.qs"
  )
)
## mgi topn
qsave(
  mgi_topn_results_post_cor,
  file = glue(
    "results/mgi/{opt$mgi_version}/{opt$outcome}/",
    "mgi_{opt$outcome}_topn_{opt$mgi_version}_results_post_cor.qs"
  )
)

## ukb
for (i in seq_along(time_thresholds)) {
  save_qs(
    x = ukb_results[[i]],
    file = glue(
      "results/ukb/{opt$ukb_version}/{opt$outcome}/",
      "ukb_{opt$outcome}_t{time_thresholds[i]}_",
      "{opt$ukb_version}_results.qs"
    )
  )
}

## ukb pwide
qsave(
  ukb_pwide_results_post_cor,
  file = glue(
    "results/ukb/{opt$ukb_version}/{opt$outcome}/",
    "ukb_{opt$outcome}_pwide_{opt$ukb_version}_results_post_cor.qs"
  )
)
## ukb topn
qsave(
  ukb_topn_results_post_cor,
  file = glue(
    "results/ukb/{opt$ukb_version}/{opt$outcome}/",
    "ukb_{opt$outcome}_topn_{opt$ukb_version}_results_post_cor.qs"
  )
)

cli_alert_success("script complete! 🎉🎉🎉")
