cancer_timing <- function(phe_dsb_data) {
    first_phe_dsb <- phe_dsb_data[ phe_dsb_data[, .I[dsb == min(dsb)], by = c("id", "phecode")][["V1"]] ]
    first_dx <- first_phe_dsb |>
        select(id, dsb) |>
        filter(dsb == min(dsb), .by = "id") |> 
        distinct()
    first_cancer_dx <- first_phe_dsb |>
        filter(phecode %in% cancer_phecodesx) |>
        select(id, dsb) |>
        filter(dsb == min(dsb), .by = "id") |>
        distinct()
    timing <- full_join(
        first_dx |> rename(first_dsb = dsb),
        first_cancer_dx |> rename(cancer_dsb = dsb),
        by = "id"
    ) |>
        mutate(precancer_followup = cancer_dsb - first_dsb) |>
        as_tibble()
    return(timing)
}

mgi_timing <- cancer_timing(mgi_full_phe)
ukb_timing <- cancer_timing(ukb_full_phe)

plot_timing <- function(timing_data) {
    timing_data |>
        # filter(precancer_followup > 0) |>
        ggplot(aes(x = precancer_followup)) +
        geom_histogram(aes(y = after_stat(density)), bins = 100) +
        geom_vline(xintercept = median(timing_data$precancer_followup, na.rm = TRUE), color = "red") +
        geom_vline(xintercept = mean(timing_data$precancer_followup, na.rm = TRUE), color = "blue") +
        geom_text(
            x = 2500,
            y = 0.001,
            label = paste0("Median: ", round(median(timing_data$precancer_followup, na.rm = TRUE), 2)),
            color = "red",
            size = 5,
            hjust = 0,
            vjust = 0
        ) +
        geom_text(
            x = 2500,
            y = 0.0005,
            label = paste0("Mean: ", round(mean(timing_data$precancer_followup, na.rm = TRUE), 2)),
            color = "blue",
            size = 5,
            hjust = 0,
            vjust = 0
        ) +
        geom_density()
}

plot_timing(mgi_timing)
plot_timing(ukb_timing)

phe_dsb_data <- mgi_full_phe

specific_cancer_codes <- cancer_phecodesx[!(cancer_phecodesx %in% c("CA_130", "CA_105", "CA_112", "CA_103"))]

prior_cancers <- function(phe_dsb_data, index_diagnosis = "CA_105.1") {
    first_phe_dsb <- phe_dsb_data[phe_dsb_data[, .I[dsb == min(dsb)], by = c("id", "phecode")][["V1"]]]
    first_index <- first_phe_dsb |>
        filter(phecode == index_diagnosis) |>
        distinct()
    specific_cancer_codes <- cancer_phecodesx[!(cancer_phecodesx %in% c("CA_130", "CA_112", "CA_130", "CA_103", "CA_116", remove_decimal_places(index_diagnosis)))]
    first_other_cancer <- first_phe_dsb |>
        filter(phecode %in% specific_cancer_codes & phecode != index_diagnosis) |>
        inner_join(first_index |> select(id, index_dsb = dsb), by = "id") |>
        filter(dsb < index_dsb) |>
        distinct()
    first_index |>
        mutate(
            other_first = ifelse(id %in% first_other_cancer$id, 1, 0)
        )
    # phe_dsb_to_pim(first_other_cancer)
}

prior_cancers(mgi_full_phe, "CA_101.8") |>
    count(other_first) |>
    mutate(
        prop = round(n * 100 / sum(n), 1)
    )

    first_index |>
        mutate(
            other_first = ifelse(id %in% first_other_cancer$id, 1, 0)
        ) |>
    count(other_first) |>
    mutate(
        prop = round(n * 100 / sum(n), 1)
    )

sub("\\..*", "", "CA_105.1")
nchar("CA_105.1")

outcome <- "CA_107.2"

specific_cancer_codes <- cancer_phecodesx[!(cancer_phecodesx %in% c("CA_130", "CA_112", "CA_130", "CA_103", "CA_116", remove_decimal_places(outcome)))]


first_other_cancer |>
    count(phecode) |>
    arrange(desc(n)) |>
    left_join(ms::pheinfox[, .(phecode, description)]) |>
    as_tibble()

ms::phecodex_icdcm_map |>
    filter(phecode %in% c("CA_101.4")) |>
    select(-one_of("phecode")) |>
    distinct() |>
    ms::icd_desc() |>
        as_tibble() |>
        print(n = 48)

tibble(
    "phecode" = cancer_phecodesx
) |>
    left_join(
        ms::pheinfox[, .(phecode, description)],
        by = "phecode"
    ) |>
    print(n = 151)


remove_decimal_places <- function(x) {
    # split the string into two parts: base and decimal part.
    split_string <- strsplit(x, ".", fixed = TRUE)[[1]]

    # If the string doesn't contain a decimal part, return the original string
    if (length(split_string) < 2) {
        return(x)
    }

    # otherwise, iterate to remove each decimal digit progressively
    base_part    <- split_string[1]
    decimal_part <- split_string[2]

    results <- vector() # container for results

    for (i in 1:nchar(decimal_part)) {
       new_string <- paste0(base_part, ".", substring(decimal_part, 1, nchar(decimal_part) - i))
       # Remove the string if it ends with decimal
       if (!grepl("\\.$", new_string)) {
           results <- c(results, new_string)
       }
    }

    results <- c(results, base_part) # add the base part without any decimal

    return(results)
}

remove_decimal_places("CA_105.12345")

#######################################
outcome <- "CA_105.1"
mgi_brca <- prepare_case_data(
    phe_dsb_data       = mgi_full_phe,
    demo_data          = mgi_demo,
    outcome            = outcome,
    malignant_phecodes = cancer_phecodesx[keep == 1, phecode],
    specific_phecodes  = cancer_phecodesx[specific == 1, phecode],
)
ukb_brca <- prepare_case_data(
    phe_dsb_data       = ukb_full_phe,
    demo_data          = ukb_demo,
    outcome            = outcome,
    malignant_phecodes = cancer_phecodesx[keep == 1, phecode],
    specific_phecodes  = cancer_phecodesx[specific == 1, phecode],
)

unmatched_analytic_data <- function(case_data) {
    bind_rows(
        case_data$case_pim0 |> mutate(case = 1),
        case_data$unmatched_pim |> mutate(case = 0)
    ) |>
    left_join(
        bind_rows(
            case_data$case_demo,
            case_data$unmatched_demo
        ),
        by = "id"
    )
}

mgi_brca_analytic <- unmatched_analytic_data(mgi_brca)
mgi_brca_phecodes <- names(mgi_brca_analytic)[names(mgi_brca_analytic) %in% ms::pheinfox[, phecode]]

ukb_brca_analytic <- unmatched_analytic_data(ukb_brca)
ukb_brca_phecodes <- names(ukb_brca_analytic)[names(ukb_brca_analytic) %in% ms::pheinfox[, phecode]]

options(future.globals.maxSize = (object.size(mgi_brca_analytic) + 50 * 1024^2))
plan(multicore, workers = 16)
mgi_brca_res <- mgi_brca_phecodes |>
    future_map(
        \(x) {
            phewas_analysis(
                data       = mgi_brca_analytic,
                outcome    = "case",
                exposure   = x,
                covariates = c("age", "length_followup"),
                method     = "SPAtest"
            )
        },
        .progress = TRUE
    ) |>
    bind_rows() |>
    arrange(p_value)

options(future.globals.maxSize = (object.size(ukb_brca_analytic)+ 50 * 1024^2))
ukb_brca_res <- ukb_brca_phecodes |>
    future_map(
        \(x) {
            phewas_analysis(
                data       = ukb_brca_analytic,
                outcome    = "case",
                exposure   = x,
                covariates = c("age", "length_followup"),
                method     = "SPAtest"
            )
        },
        .progress = TRUE
    ) |>
    bind_rows() |>
    arrange(p_value)

mgi_brca_plot <- plot_phewasx(mgi_brca_res |> filter(is.finite(log10p)), title = glue("Unmatched {outcome} PheWAS in MGI at 0 years"), phe_var = "phecode")
ukb_brca_plot <- plot_phewasx(ukb_brca_res |> filter(is.finite(log10p)), title = glue("Unmatched {outcome} PheWAS in UKB at 0 years"), phe_var = "phecode")

ggsave(filename = glue("bin/mgi_{outcome}_plot.pdf"), plot = mgi_brca_plot,
        width = 8, height = 6, device = cairo_pdf)
ggsave(
    filename = glue("bin/ukb_{outcome}_plot.pdf"), plot = ukb_brca_plot,
    width = 8, height = 6, device = cairo_pdf
)

quick_summary_table <- function(
    data,
    outcome,
    cohort,
    case_var  = "case",
    desc_vars = c("age", "sex", "length_followup", "unique_phecodes", "unique_encounters", "EM_236", "SS_809", "CV_401")
) {
    phecode_vars <- names(data)[names(data) %in% ms::pheinfox[, phecode]]
    for (i in phecode_vars) {
        labelled::var_label(data[[i]]) <- ms::pheinfox[ms::pheinfox[, phecode] == i, description]
    }
    data |>
        mutate(case = factor(ifelse(get(case_var) == 0, 0, 1), levels = c(0, 1), labels = c("Control", "Case"))) |>
        select(all_of(c(case_var, desc_vars))) |>
        tbl_summary(
            by = case,
            label = list(
                age ~ "Age",
                sex ~ "Sex",
                length_followup ~ "Length of Follow-up (days)",
                unique_phecodes ~ "Unique phecodes",
                unique_encounters ~ "Unique encounters"
            )
        ) |>
        add_overall() |>
        add_n() |>
        add_p() |>
        bold_labels() |>
        modify_spanning_header(c("stat_1", "stat_2") ~ paste0("**", ms::pheinfox[phecode == outcome, description],"**")) |>
        modify_caption(paste0("**Table 1.** Summary of ", ms::pheinfox[phecode == outcome, description], " [", outcome, "] in ", cohort, " at 0 years"))
}

mgi_brca_table <- quick_summary_table(mgi_brca_analytic, outcome = outcome, cohort = "MGI")
ukb_brca_table <- quick_summary_table(ukb_brca_analytic, outcome = outcome, cohort = "UKB")

mgi_brca_table |>
    as_gt() |>
    gt::gtsave(glue("bin/mgi_{outcome}_table.docx"))
ukb_brca_table |>
    as_gt() |>
    gt::gtsave(glue("bin/ukb_{outcome}_table.docx"))
