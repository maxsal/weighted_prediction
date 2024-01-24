outcome <- "CA_101.6"
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")

mgi_test <- prepare_matched_phenomes(
    phe_dsb_data          = mgi_full_phe,
    demo_data             = mgi_demo,
    outcome               = opt$outcome,
    case_definition       = "first", # "first" | "first_occurrence" or "all" | "all_occurrences" | "any" | "anytime" | "any_occurrence"
    malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
    specific_phecodes     = cancer_phecodesx[specific == 1, phecode],
    matched               = TRUE,
    match_ratio           = 2,
    match_caliper         = 0.25,
    nearest_vars          = c("age", "length_followup"),
    exact_vars            = c("female"),
    control_definition    = "no_cancer",
    time_thresholds       = c(0, 0.5, 1, 2, 3, 5),
    exclude_other_cancers = TRUE
)

ukb_test <- prepare_matched_phenomes(
    phe_dsb_data          = ukb_full_phe,
    demo_data             = ukb_demo |>
        rename(age = age_at_first_diagnosisx) |>
        mutate(
            length_followup = last_dsbx - first_dsbx,
            female = as.numeric(sex == "Female")
        ),
    outcome               = outcome,
    case_definition       = "first", # "first" | "first_occurrence" or "all" | "all_occurrences" | "any" | "anytime" | "any_occurrence"
    malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
    specific_phecodes     = cancer_phecodesx[specific == 1, phecode],
    matched               = TRUE,
    match_ratio           = 2,
    match_caliper         = 0.25,
    nearest_vars          = c("age", "length_followup"),
    exact_vars            = c("female"),
    control_definition    = "no_cancer",
    time_thresholds       = c(0, 0.5, 1, 2, 3, 5),
    exclude_other_cancers = TRUE
)

time_threshold <- 1
time_threshold_index <- which(time_thresholds == time_threshold)

mgi_test_analytic <- mgi_test$tr_pim[[time_threshold_index]]
ukb_test_analytic <- ukb_test$tr_pim[[time_threshold_index]]

mgi_test_phecodes <- names(mgi_test_analytic)[names(mgi_test_analytic) %in% ms::pheinfox[, phecode]]
ukb_test_phecodes <- names(ukb_test_analytic)[names(ukb_test_analytic) %in% ms::pheinfox[, phecode]]

options(future.globals.maxSize = (object.size(mgi_test_analytic) + 50 * 1024^2))
plan(multicore, workers = 16)
mgi_test_res2 <- mgi_test_phecodes |>
    future_map(
        \(x) {
            phewas_analysis(
                data       = mgi_test_analytic,
                outcome    = "case",
                exposure   = x,
                covariates = c("age", "length_followup"),
                method     = "glm"
            )
        },
        .progress = TRUE
    ) |>
    bind_rows() |>
    arrange(p_value)

options(future.globals.maxSize = (object.size(ukb_test_analytic) + 50 * 1024^2))
ukb_test_res <- ukb_test_phecodes |>
    future_map(
        \(x) {
            phewas_analysis(
                data       = ukb_test_analytic,
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

mgi_test_plot2 <- plot_phewasx(mgi_test_res2 |> filter(is.finite(log10p)), title = glue("Matched {outcome} PheWAS in MGI at {time_threshold} years"), phe_var = "phecode")
ukb_test_plot <- plot_phewasx(ukb_test_res |> filter(is.finite(log10p)), title = glue("Matched {outcome} PheWAS in UKB at {time_threshold} years"), phe_var = "phecode")

ggsave(
    filename = glue("bin/mgi_matched_{outcome}_t{time_threshold}_plot.pdf"), plot = mgi_test_plot,
    width = 8, height = 6, device = cairo_pdf
)
ggsave(
    filename = glue("bin/ukb_matched_{outcome}_t{time_threshold}_plot.pdf"), plot = ukb_test_plot,
    width = 8, height = 6, device = cairo_pdf
)

quick_summary_table <- function(
    data,
    outcome,
    cohort,
    case_var = "case",
    desc_vars = c("age", "sex", "length_followup", "unique_phecodes", "unique_encounters", "EM_236", "SS_809", "CV_401")) {
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
        modify_spanning_header(c("stat_1", "stat_2") ~ paste0("**", ms::pheinfox[phecode == outcome, description], "**")) |>
        modify_caption(paste0("**Table 1.** Summary of ", ms::pheinfox[phecode == outcome, description], " [", outcome, "] in ", cohort, " at 0 years"))
}

mgi_test_table <- quick_summary_table(mgi_test_analytic, outcome = outcome, cohort = "MGI")
ukb_test_table <- quick_summary_table(ukb_test_analytic, outcome = outcome, cohort = "UKB")

mgi_test_table |>
    as_gt() |>
    gt::gtsave(glue("bin/mgi_matched_{outcome}_t{time_threshold}_table.docx"))
ukb_test_table |>
    as_gt() |>
    gt::gtsave(glue("bin/ukb_matched_{outcome}_t{time_threshold}_table.docx"))


