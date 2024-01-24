# prepare mgi and ukb data and conduct unweighted phewas
# author:  max salvatore
# date:    20231201

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)
set.seed(61787)

# load libraries
ms::libri(data.table, qs, optparse, glue, cli, ms, tidyverse)

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
    make_option("--matching_ratio",
        type = "numeric", default = "2",
        help = "Number of non-cases to match per case [default = %default]"
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

time_thresholds <- strsplit(opt$time_thresholds, ",")[[1]]

# load mgi data ----------------------------------------------------------------
mgi_phe_res <- map(
    time_thresholds,
    \(i) {
        qread(
            glue(
                "results/mgi/{opt$mgi_version}/{opt$outcome}/mgi_{opt$mgi_version}_phewas_{opt$outcome}_t{i}_r{opt$matching_ratio}.qs"
            )
        )
    }
)|> set_names(glue("t{time_thresholds}"))

mgi_tr_pim <- map(
    time_thresholds,
    \(i) {
        qread(
            glue(
                "data/private/mgi/{opt$mgi_version}/{opt$outcome}/time_restricted_phenomes/mgi_{opt$mgi_version}_{opt$outcome}_t{i}_pim_r{opt$matching_ratio}.qs"
            )
        )
    }
) |> set_names(glue("t{time_thresholds}"))

mgi_after_cor <- map(
    seq_along(time_thresholds),
    \(i) {
        top_after_correlation <- remove_by_correlation(
            cooccurrence_results = as.data.table(mgi_phe_res[[i]]),
            pim = mgi_tr_pim[[i]],
            exposure_var = "phecode",
            p_value_var = "p_value",
            top_n = 50,
            corr_thresh = 0.5,
            weights = NULL
        )
    }
) |> set_names(glue("t{time_thresholds}"))

# load ukb data ----------------------------------------------------------------
ukb_phe_res <- map(
    time_thresholds,
    \(i) {
        qread(
            glue(
                "results/ukb/{opt$ukb_version}/{opt$outcome}/ukb_{opt$ukb_version}_phewas_{opt$outcome}_t{i}_r{opt$matching_ratio}.qs"
            )
        )
    }
)|> set_names(glue("t{time_thresholds}"))

ukb_tr_pim <- map(
    time_thresholds,
    \(i) {
        qread(
            glue(
                "data/private/ukb/{opt$ukb_version}/{opt$outcome}/time_restricted_phenomes/ukb_{opt$ukb_version}_{opt$outcome}_t{i}_pim_r{opt$matching_ratio}.qs"
            )
        )
    }
)|> set_names(glue("t{time_thresholds}"))

ukb_after_cor <- map(
    seq_along(time_thresholds),
    \(i) {
        top_after_correlation <- remove_by_correlation(
            cooccurrence_results = as.data.table(ukb_phe_res[[i]]),
            pim = ukb_tr_pim[[i]],
            exposure_var = "phecode",
            p_value_var = "p_value",
            top_n = 50,
            corr_thresh = 0.5,
            weights = NULL
        )
    }
)  |> set_names(glue("t{time_thresholds}"))

# save results -----------------------------------------------------------------
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

cli_alert_success("Done! 🎉")