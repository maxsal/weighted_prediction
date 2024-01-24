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
  make_option("--matching_ratio",
    type = "numeric", default = "2",
    help = "Number of non-cases to match per case [default = %default]"
  ),
  make_option("--time_thresholds",
    type = "character", default = "0,1,2,5",
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
    type = "character", default = "ip_selection,ip_crossfit,ps_selection",
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

# 2. specifications (specifies outcome) --------------------------------------
time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])
cooccur_covs    <- c("age", "female", "length_followup")

# process outcome information
outcome_sex <- ms::pheinfox[phecode == opt$outcome, sex]

if (tolower(outcome_sex) %in% c("female", "f", "male", "m")) {
  cooccur_covs <- cooccur_covs[!cooccur_covs %in% c("female", "male", "m", "f")]
  if (length(cooccur_covs) == 0) {
    cooccur_covs <- NULL
  }
}

# 3. extra preparations (outcome-specific) -------------------------------------
## pull file paths corresponding to the data version specified
file_paths <- get_files(
  mgi_version = opt$mgi_version
)


# 4. read data -----------------------------------------------------------------
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")

## mgi
cli_progress_step("loading mgi data")

mgi_tr_pims <- map(
    seq_along(time_thresholds),
    \(x) {
        qread(
            glue(
                "data/private/mgi/{opt$mgi_version}/",
                "{opt$outcome}/time_restricted_phenomes/mgi_",
                "{opt$mgi_version}_{opt$outcome}_",
                "t{time_thresholds[x]}_pim_",
                "r{opt$matching_ratio}.qs"
            )
        )
    }
)
names(mgi_tr_pims) <- glue("t{time_thresholds}_threshold")

# standardize
mgi_tr_pims_std <- map(
    mgi_tr_pims,
    \(x) {
        vars <- names(x)[!names(x) %in% c("id", "case")]
        standardize_variables(x, cols = vars)
    }
)
names(mgi_tr_pims_std) <- glue("t{time_thresholds}_threshold")

mgi_covariates <- qread(
      glue(
          "data/private/mgi/{opt$mgi_version}/",
          "{opt$outcome}/mgi_",
          "{opt$mgi_version}_{opt$outcome}_",
          "match_data_",
          "r{opt$matching_ratio}.qs"
      )
)

mgi_weights <- qread(
    glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs")
)

## phenome
pheinfo <- ms::pheinfox
if (pheinfo[phecode == opt$outcome, sex] != "Both") {
    if (any(c("sex", "female", "male") %in% tolower(cooccur_covs))) {
        cooccur_covs <- cooccur_covs[-which(tolower(cooccur_covs) %in% c("sex", "female", "male"))]
    }
}

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
    names(mgi_tr_pims_std),
    \(x) {
        list(mgi_tr_pims_std[[x]], mgi_covariates[, !c("age", "length_followup", "sex", "case")], mgi_weights[, ..weight_id]) |>
            reduce(left_join, by = "id") |>
            (\(y) copy(y)[, `:=`(age_at_threshold = round(get(x) / 365.25, 1))])()
    }
)
names(mgi_tr_merged) <- glue("t{time_thresholds}_threshold")

# phewas
mgi_phewas_res <- map(
    seq_along(weight_vars),
    \(w) {
        cli_progress_step(paste0("running phewas for weight ", weight_vars[w]))
        map(
            seq_along(time_thresholds),
            \(i) {
                print(paste0("t = ", time_thresholds[i]))
                mgi_data <- mgi_tr_merged[[i]]
                mgi_phecodes <- names(mgi_data)[names(mgi_data) %in% ms::pheinfox[, phecode]]
                mgi_design <- survey::svydesign(
                    ids     = ~1,
                    data    = mgi_data[!is.na(mgi_data[[weight_vars[w]]]), ],
                    weights = ~get(weight_vars[w])
                )
                ms::map_phewas(
                    data        = mgi_data,
                    design      = mgi_design,
                    outcome     = "case",
                    exposures   = mgi_phecodes,
                    covariates  = c("age", "length_followup"),
                    method      = "weighted",
                    .weight_var = weight_vars[w],
                    workers     = 16
                )
            }
        )
    }
)
names(mgi_phewas_res) <- weight_vars
for (i in seq_along(mgi_phewas_res)) {
    names(mgi_phewas_res[[i]]) <- glue("t{time_thresholds}_threshold")
}

mgi_pwide_hits <- list()
mgi_topn_hits <- list()
for (i in seq_along(weight_vars)) {
    mgi_pwide_hits[[i]] <- lapply(mgi_phewas_res[[i]], \(x) get_pwide_hits(data = as.data.table(x), exposure_var = "phecode"))
    mgi_topn_hits[[i]]  <- lapply(mgi_phewas_res[[i]],\(x) get_top_hits(data = as.data.table(x), exposure_var = "phecode", top_n = 50))
}
names(mgi_pwide_hits) <- weight_vars
names(mgi_topn_hits)  <- weight_vars

pwide_post_cor <- list()
for (i in seq_along(weight_vars)) {
    pwide_post_cor[[i]] <- mapply(
        \(x, y, z) {
            x <- x[!is.na(get(weight_vars[i])), ]
            remove_by_correlation(
                pim = x, cooccurrence_results = as.data.table(y), exposures_to_consider = z,
                exposure_var = "phecode", p_value_var = "p_value", corr_thresh = 0.5,
                weights = x[, get(weight_vars[i])]
            )
        },
        x = mgi_tr_merged, y = mgi_phewas_res[[i]], z = mgi_pwide_hits[[i]]
    )
}

topn_post_cor <- list()
for (i in seq_along(weight_vars)) {
    topn_post_cor[[i]] <- mapply(
        \(x, y, z) {
            x <- x[!is.na(get(weight_vars[i])), ]
            remove_by_correlation(
                pim = x, cooccurrence_results = as.data.table(y), exposures_to_consider = z,
                exposure_var = "phecode", p_value_var = "p_value", corr_thresh = 0.5,
                weights = x[, get(weight_vars[i])]
            )
        },
        x = mgi_tr_merged, y = mgi_phewas_res[[i]], z = mgi_topn_hits[[i]]
    )
}

for (i in seq_along(weight_vars)) {
    pwide_post_cor[[i]] <- pwide_post_cor[[i]] |>
        as.data.table() |>
        lapply(
            \(x) {
                out        <- as.data.table(x)
                names(out) <- rownames(pwide_post_cor[[i]])
                out
            }
        )
}
names(pwide_post_cor) <- weight_vars

for (i in seq_along(weight_vars)) {
    topn_post_cor[[i]] <- topn_post_cor[[i]] |>
        as.data.table() |>
        lapply(
            \(x) {
                out <- as.data.table(x)
                names(out) <- rownames(topn_post_cor[[i]])
                out
            }
        )
}
names(topn_post_cor) <- weight_vars

# save results -----------------------------------------------------------------
for (w in seq_along(weight_vars)) {
    for (i in seq_along(time_thresholds)) {
        save_qs(
            x = mgi_phewas_res[[w]][[i]],
            file = glue(
                "results/mgi/{opt$mgi_version}/{opt$outcome}/",
                "mgi_{opt$outcome}_t{time_thresholds[i]}_",
                "{opt$mgi_version}_r{opt$matching_ratio}_{weight_vars[w]}_results.qs"
            )
        )
    }
}

mgi_weighted_phewas_plots <- map(
    seq_along(weight_vars),
    \(w) {
        map(seq_along(time_thresholds),
            \(i) {
                plot_phewasx(
                    mgi_phewas_res[[w]][[i]],
                    phe_var = "phecode",
                    title = stringr::str_wrap(glue(
                        "PheWAS for {ms::pheinfox[phecode == opt$outcome, description]} ",
                        "[{opt$outcome}] at t = {time_thresholds[i]} in MGI ",
                        "with {weight_vars[w]} weighting"
                    ), width = 50)
                )
            }
        )
    }
)

walk(
    seq_along(weight_vars),
    \(w) {
        map(
            seq_along(time_thresholds),
            \(i) {
                ggsave(
                    filename = glue(
                        "results/mgi/{opt$mgi_version}/{opt$outcome}/mgi_",
                        "{opt$outcome}_t{time_thresholds[i]}_",
                        "{opt$mgi_version}_r{opt$matching_ratio}_{weight_vars[w]}_phewas_plot.pdf"
                    ),
                    plot = mgi_weighted_phewas_plots[[w]][[i]],
                    width = 8,
                    height = 6,
                    device = cairo_pdf
                )
            }
        )
    }
)

qsave(
    pwide_post_cor,
    file = glue(
        "results/mgi/{opt$mgi_version}/{opt$outcome}/",
        "mgi_{opt$outcome}_pwide_{opt$mgi_version}_r{opt$matching_ratio}_",
        "weighted_post_cor_results.qs"
    )
)

qsave(
    topn_post_cor,
    file = glue(
        "results/mgi/{opt$mgi_version}/{opt$outcome}/",
        "mgi_{opt$outcome}_top50_{opt$mgi_version}_r{opt$matching_ratio}_",
        "weighted_post_cor_results.qs"
    )
)

cli_alert_success("script done! 🙌")
