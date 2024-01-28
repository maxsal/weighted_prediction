# libraries and such
ms::libri(ms, data.table, qs, glue, cli, optparse, survey, tidyverse)

# optparse list ----------------------------------------------------------------
option_list <- list(
    make_option("--outcome",
        type = "character", default = "CA_101.1",
        help = "Outcome phecode [default = %default]"
    ),
    make_option("--mgi_version",
        type = "character", default = "20220822",
        help = "Version of AOU data [default = %default]"
    ),
    make_option("--time_thresholds",
        type = "character", default = "0,1,2,5",
        help = glue(
            "Time thresholds for the phenome data ",
            "[default = %default]"
        )
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

outcome         <- opt$outcome
mgi_version     <- opt$mgi_version
time_thresholds <- as.character(strsplit(opt$time_thresholds, ",")[[1]])
cooccur_covs    <- c("age", "female", "length_followup")
matching_ratio  <- opt$matching_ratio

walk(
    list.files("fn", full.names = TRUE),
    source
)

# load data
cli_progress_step("Loading data")
## mgi results
### unweighted
mgi_data <- map(
    seq_along(time_thresholds),
    \(i) { as.data.table(qread(glue("results/mgi/{mgi_version}/{outcome}/mgi_{mgi_version}_train_phewas_{outcome}_t{time_thresholds[i]}_r{matching_ratio}.qs"))) }
) |> set_names(glue("t{time_thresholds}"))
mgi_post_cor <- map(
    seq_along(time_thresholds),
    \(i) {
        as.data.table(qread(glue("results/mgi/{mgi_version}/{outcome}/mgi_{mgi_version}_train_post_cor_phe_{outcome}_t{time_thresholds[i]}_r{matching_ratio}.qs")))
    }
)

### ip weighted
mgi_ip_data <- map(
    seq_along(time_thresholds),
    \(i) {
        as.data.table(qread(glue("results/mgi/{mgi_version}/{outcome}/mgi_{mgi_version}_train_ip_phewas_{outcome}_t{time_thresholds[i]}_r{matching_ratio}.qs")))
    }
) |> set_names(glue("t{time_thresholds}"))
mgi_ip_post_cor <- map(
    seq_along(time_thresholds),
    \(i) {
        as.data.table(qread(glue("results/mgi/{mgi_version}/{outcome}/mgi_{mgi_version}_train_ip_post_cor_phe_{outcome}_t{time_thresholds[i]}_r{matching_ratio}.qs")))
    }
)

## mgi data
### tr pims
mgi_tr_pims <- map(
    seq_along(time_thresholds),
    \(i) {
        qread(glue(
            "data/private/mgi/{opt$mgi_version}/",
            "{opt$outcome}/time_restricted_phenomes/mgi_",
            "{opt$mgi_version}_{opt$outcome}_",
            "t{time_thresholds[i]}_pim_",
            "r{opt$matching_ratio}.qs"
        ))
    }
) |>
    set_names(paste0("t", time_thresholds, "_threshold"))

### demographics
mgi_covariates <- read_qs(glue(
    "data/private/mgi/{opt$mgi_version}/",
    "{opt$outcome}/mgi_",
    "{opt$mgi_version}_{opt$outcome}_",
    "match_data_",
    "r{opt$matching_ratio}.qs"
))

### weights
mgi_weights <- qread(
    glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs")
)

mgi_tr_merged <- map(
    names(mgi_tr_pims),
    \(x) {
        list(
            mgi_tr_pims[[x]],
            mgi_covariates[, !c("case", "age", "sex", "length_followup")],
            mgi_weights
        ) |>
        reduce(merge.data.table, by = "id", all.x = TRUE) |>
        mutate(
            age_at_threshold = round(get(x) / 365.25, 1)
        )
    }
) |> set_names(glue("t{time_thresholds}"))

## pheinfo
pheinfo <- ms::pheinfox
if (pheinfo[phecode == outcome, sex] != "Both") {
    if (any(c("sex", "female", "male") %in% tolower(cooccur_covs))) {
        cooccur_covs <- cooccur_covs[-which(tolower(cooccur_covs) %in% c("sex", "female", "male"))]
    }
}

# estimate disease weights
## univariable
cli_progress_step("Estimating univariable disease weights")

mgi_topn_hits <- map(
    mgi_post_cor,
    \(x) {
        x[, phecode]
    }
) |> set_names(glue("t{time_thresholds}"))
mgi_ip_topn_hits <- map(
    mgi_ip_post_cor,
    \(x) {
        x[, phecode]
    }
) |> set_names(glue("t{time_thresholds}"))

topn_univariable_results <- map(
    seq_along(mgi_ip_topn_hits),
    \(i) {
        print(glue("t{time_thresholds[i]}"))
        tmp_hits <- mgi_ip_topn_hits[[i]][mgi_ip_topn_hits[[i]] %in% names(mgi_tr_merged[[names(mgi_ip_topn_hits)[i]]])]
        if (length(tmp_hits) == 0) {
            return(data.table())
        }
        ms::map_phewas(
            data               = mgi_tr_merged[[names(mgi_ip_topn_hits)[i]]],
            outcome            = "case",
            exposures          = tmp_hits,
            covariates         = cooccur_covs,
            method             = "logistf",
            workers            = 12
        ) |> as.data.table()
    }
) |> set_names(glue("t{time_thresholds}"))

topn_ip_univariable_results <- map(
    seq_along(mgi_ip_topn_hits),
    \(i) {
        print(glue("t{time_thresholds[i]}"))
        tmp_hits <- mgi_ip_topn_hits[[i]][mgi_ip_topn_hits[[i]] %in% names(mgi_tr_merged[[names(mgi_ip_topn_hits)[i]]])]
        if (length(tmp_hits) == 0) {
            return(data.table())
        }
        tmp_dsn <- survey::svydesign(
            ids     = ~1,
            weights = ~ip_selection,
            data    = mgi_tr_merged[[names(mgi_ip_topn_hits)[i]]][!is.na(ip_selection), ]
        )
        ms::map_phewas(
            data               = mgi_tr_merged[[names(mgi_ip_topn_hits)[i]]],
            design             = tmp_dsn,
            outcome            = "case",
            exposures          = tmp_hits,
            covariates         = cooccur_covs,
            method             = "weighted",
            .weight_var        = "ip_selection",
            workers            = 12
        ) |> as.data.table()
    }
) |> set_names(glue("t{time_thresholds}"))

## multivariable
### helper function
glm_coefficient_summary <- function(
    data,
    covs,
    hits) {
    hits <- hits[hits %in% names(data)]
    if (length(hits) == 0) return(data.table())
    multi_mod_covs <- paste0(c(covs, hits), collapse = " + ")
    multi_mod <- glm(
        formula = as.formula(paste0("case ~ ", multi_mod_covs)),
        data    = data,
        family  = "binomial"
    )
    cbind(
        data.table(phecode = rownames(summary(multi_mod)$coefficients)),
        as.data.table(summary(multi_mod)$coefficients)
    )[
        !phecode %in% c("(Intercept)", covs),
        .(phecode, beta = Estimate, se_beta = `Std. Error`, p_value = `Pr(>|z|)`)
    ][, `:=` (log10p = log10(p_value))][]
}
svyglm_coefficient_summary <- function(
    data,
    covs,
    hits,
    weight_var
) {
    hits <- hits[hits %in% names(data)]
    if (length(hits) == 0) return(data.table())
    multi_mod_covs <- paste0(c(covs, hits), collapse = " + ")
    multi_mod_dsn <- svydesign(
        id = ~1,
        weights = ~get(weight_var),
        data = data[!is.na(get(weight_var)), ]
    )
    multi_mod <- svyglm(
        formula = as.formula(paste0("case ~ ", multi_mod_covs)),
        design  = multi_mod_dsn,
        family  = "quasibinomial"
    )
    cbind(
        data.table(phecode = rownames(summary(multi_mod)$coefficients)),
        as.data.table(summary(multi_mod)$coefficients)
    )[
        !phecode %in% c("(Intercept)", covs),
        .(phecode, beta = Estimate, se_beta = `Std. Error`, p_value = `Pr(>|t|)`, weight = weight_var)
    ][, `:=` (log10p = log10(p_value))][]
}
###

topn_multivariable_results <- map(
    seq_along(mgi_topn_hits),
    \(i) {
        glm_coefficient_summary(
            data = mgi_tr_merged[[names(mgi_topn_hits)[i]]],
            covs = cooccur_covs,
            hits = mgi_topn_hits[[i]]
        )
    }
) |> set_names(glue("t{time_thresholds}"))

ip_topn_multivariable_results <- map(
    seq_along(mgi_ip_topn_hits),
    \(i) {
        svyglm_coefficient_summary(
            data = mgi_tr_merged[[names(mgi_ip_topn_hits)[i]]],
            covs = cooccur_covs,
            hits = mgi_ip_topn_hits[[i]],
            weight_var = "ip_selection"
        )
    }
) |> set_names(glue("t{time_thresholds}"))


# save disease weights
cli_progress_step("Saving results")
## unweighted
results_list <- list()
for (i in seq_along(time_thresholds)) {
    # unweighted
    ## univariable
    results_list[[length(results_list) + 1]] <- topn_univariable_results[[i]]
    names(results_list)[length(results_list)] <- glue("t{time_thresholds[i]}_top50_univariable_results")

    ## multivariable
    results_list[[length(results_list) + 1]] <- topn_multivariable_results[[i]]
    names(results_list)[length(results_list)] <- glue("t{time_thresholds[i]}_top50_multivariable_results")

    # ip-weighted
    ## univariable
    results_list[[length(results_list) + 1]] <- topn_ip_univariable_results[[i]]
    names(results_list)[length(results_list)] <- glue("t{time_thresholds[i]}_ip_top50_univariable_results")

    ## multivariable
    results_list[[length(results_list) + 1]] <- ip_topn_multivariable_results[[i]]
    names(results_list)[length(results_list)] <- glue("t{time_thresholds[i]}_ip_top50_multivariable_results")
}

qsave(
    x = results_list,
    file = glue(
        "results/mgi/{mgi_version}/{outcome}/",
        "mgi_{outcome}_r{matching_ratio}_results_list.qs",
    )
)

cli_progress_done()

cli_alert_success("Done! 🎉")
