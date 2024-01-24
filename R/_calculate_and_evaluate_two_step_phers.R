# libraries and such
ms::libri(ms, data.table, qs, cli, glue, purrr, ggplot2, dplyr)

ukb_version     <- "20221117"
aou_version     <- "20230309"
outcome         <- "CA_101.8"
matching_ratio  <- 2
time_thresholds <- as.character(c(0, 1, 2, 5))

# load data
cli_progress_step("Loading data")
## aou
aou_results <- qread(glue("results/aou/{aou_version}/{outcome}/aou_{outcome}_r{matching_ratio}_results_list.qs"))

## ukb
ukb_tr_pims <- map(
    seq_along(time_thresholds),
    \(i) {
        qread(glue(
            "data/private/ukb/{ukb_version}/",
            "{outcome}/time_restricted_phenomes/ukb_",
            "{ukb_version}_{outcome}_",
            "t{time_thresholds[i]}_pim_",
            "r{matching_ratio}.qs"
        )) |>
        as.data.table()
    }
) |>
    set_names(glue("t{time_thresholds}_threshold"))

ukb_covariates <- qread(
    glue(
        "data/private/ukb/{ukb_version}/",
        "{outcome}/ukb_",
        "{ukb_version}_{outcome}_",
        "match_data_",
        "r{matching_ratio}.qs"
    )
)

ukb_tr_merged <- map(
    names(ukb_tr_pims),
    \(x) {
        merge.data.table(
            ukb_tr_pims[[x]],
            ukb_covariates[, !c("case", "age", "sex", "female", "length_followup")],
            by = "id",
            all.x = TRUE
        ) |> 
            dplyr::mutate(age_at_threshold = round(get(x) / 365.25, 1))
    }
) |>
    set_names(glue("t{time_thresholds}_threshold"))

# construct PheRS
cli_progress_step("Constructing PheRS")

calculate_phers <- function(
    pim,
    weight_data,
    id_var     = "id",
    trait_var  = "phecode",
    weight_var = "beta",
    name
) {
    codes <- 0
    if (nrow(weight_data) == 0) {
        cli_alert_warning("No results for {name}, skipping")
        return(data.table(phers = NA, name = name, codes = codes))
    }
    out <- data.table::copy(as.data.table(pim))[, pred := 0]
    for(i in weight_data[[trait_var]]) {
        if (!(i %in% names(out))) {
            cli_alert_warning("{i} not in phecode indicator matrix, skipping")
            next
        }
        codes <- codes + 1
        out[, pred := pred + (weight_data[get(trait_var) == i, get(weight_var)] * get(i))]
    }
    out[, phers := scale(pred)]
    out[, .(id, case, pred, phers, name = name, codes = codes)]
}

ukb_phers_list <- list()
for (i in seq_along(time_thresholds)) {
    tmp_names <- names(aou_results)[startsWith(names(aou_results), glue("t{time_thresholds[i]}_"))]
    for (w in seq_along(tmp_names)) {
        ukb_phers_list[[length(ukb_phers_list) + 1]] <- calculate_phers(
            pim         = ukb_tr_merged[[glue("t{time_thresholds[i]}_threshold")]],
            weight_data = aou_results[[tmp_names[w]]],
            name = tmp_names[w]
        )
    }
}

phers_calculated <- map_dfr(
    ukb_phers_list,
    \(x) {
        data.table(
            name  = x[["name"]][1],
            phers = (\(y) {
                if (all(is.na(x[["phers"]])) | all(is.nan(x[["phers"]]))) {
                    0
                } else {
                    1
                }
            })()
        )
    }
)

auc_results <- map_dfr(
    which(phers_calculated[, phers == 1]),
    \(i) {
        tmp <- suppressMessages(pROC::roc(
            response  = ukb_phers_list[[i]][["case"]],
            predictor = ukb_phers_list[[i]][["phers"]],
            family    = binomial(),
            ci        = TRUE
        ))
        data.table(
            data   = "UKB",
            type   = unique(ukb_phers_list[[i]][["name"]]),
            time   = as.numeric(gsub("t", "", stringr::str_extract(ukb_phers_list[[i]][["name"]][1], "t[0-9].[0-9]|t[0-9]"))),
            auc    = tmp$auc[1],
            auc_lo = tmp$ci[1],
            auc_hi = tmp$ci[3]
        )
    }
) |>
    dplyr::mutate(
        variable = ifelse(grepl("univariable", type), "Univariable", "Multivariable"),
        weighted = fcase(
            grepl("ip", type), "IP",
            grepl("ps", type), "PS",
            default = "None"
        )
    )

quick_pwide_auc_plot <- auc_results[!grepl("top50", type), ] |>
    ggplot(aes(x = reorder(time, -as.numeric(time)), y = auc, color = variable)) +
    geom_hline(yintercept = 0.5, linewidth = 1) +
    geom_pointrange(aes(shape = weighted, ymin = auc_lo, ymax = auc_hi), position = position_dodge(0.4)) +
    scale_color_ms() +
    coord_flip(ylim = c(0, 1)) +
    labs(
        x     = "Time threshold (years)",
        y     = "AUC (95% CI)",
        title = glue("UKB {outcome} pwide PheRS performance"),
        subtitle = glue("{ms::pheinfox[phecode == outcome, description]}")
    ) +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )

quick_topn_auc_plot <- auc_results[!grepl("pwide", type), ] |>
    ggplot(aes(x = reorder(time, -as.numeric(time)), y = auc, color = variable)) +
    geom_hline(yintercept = 0.5, linewidth = 1) +
    geom_pointrange(aes(shape = weighted, ymin = auc_lo, ymax = auc_hi), position = position_dodge(0.4)) +
    scale_color_ms() +
    coord_flip(ylim = c(0, 1)) +
    labs(
        x = "Time threshold (years)",
        y = "AUC (95% CI)",
        title = glue("UKB {outcome} top N PheRS performance"),
        subtitle = glue("{ms::pheinfox[phecode == outcome, description]}")
    ) +
    theme_ms() +
    theme(
        legend.title = element_blank()
    )

# save results
cli_progress_step("Saving results")
if (!dir.exists(glue("results/ukb/{ukb_version}/{outcome}/phers/"))) {
    dir.create(glue("results/ukb/{ukb_version}/{outcome}/phers/"), recursive = TRUE)
}
qsave(
    ukb_phers_list,
    glue("results/ukb/{ukb_version}/{outcome}/phers/ukb_{outcome}_phers.qs")
)


ggsave(
    filename = glue("results/ukb/{ukb_version}/{outcome}/phers/ukb_{outcome}_pwide_phers_auc.pdf"),
    plot     = quick_pwide_auc_plot,
    width    = 6,
    height   = 6,
    device   = cairo_pdf
)
ggsave(
    filename = glue("results/ukb/{ukb_version}/{outcome}/phers/ukb_{outcome}_pwide_phers_auc.pdf"),
    plot     = quick_topn_auc_plot,
    width    = 6,
    height   = 6,
    device   = cairo_pdf
)

cli_alert_success("Done! 🎉")
