case_data <- prepare_case_data(
    phe_dsb_data          = mgi_full_phe,
    demo_data             = mgi_demo,
    outcome               = opt$outcome,
    cancer_codes          = cancer_phecodesx,
    exclude_other_cancers = TRUE
)

tmp <- bind_rows(
    case_data$case_pim0,
    case_data$unmatched_pim
) |>
    left_join(
        bind_rows(
            case_data$case_demo,
            case_data$unmatched_demo
        )
    )

names(case_data)
case_data$case
case_data$case_demo

tmp <- left_join(
    phe_dsb_data[id %in% case_data$case[, id], ],
    case_data$case,
    by = "id"
) |>
    filter(dsb <= index_dsb) |>
    summarize(
        length_followup   = max(dsb) - min(dsb),
        unique_encounters = length(unique(dsb)),
        unique_phecodes   = length(unique(phecode)),
        .by = "id"
    ) |>
    as_tibble()

tmp |>
    ggplot(aes(x = "", y = length_followup)) +
    geom_violin(
        draw_quantiles = c(0.25, 0.5, 0.75)
    ) +
    geom_jitter(alpha = 0.25)

stacked_pim <- bind_rows(
    case_data$case_pim0 |> mutate(case = 1),
    case_data$unmatched_pim |> mutate(case = 0)
)
stacked_pim[is.na(stacked_pim)] <- 0

stacked_demo <- bind_rows(
    case_data$case_demo,
    case_data$unmatched_demo
)

merged <- left_join(
    stacked_pim,
    stacked_demo,
    by = "id"
)

merged <- mgi_test$tr_pim[[3]]
table(
    merged[, case],
    merged[, sex],
    useNA = "always"
)

phecodes <- names(merged)[names(merged) %in% ms::pheinfox[, phecode]]

for (i in phecodes) {
    labelled::var_label(merged[[i]]) <- ms::pheinfox[phecode == i, description]
}

sjlabelled::get_label(merged$CV_401)


glm(case ~ EM_236 + age + length_followup, data = merged, family = quasibinomial())

merged 
mgi_test$tr_pim[[1]]|>
    mutate(case = factor(case, levels = c(0, 1), labels = c("Control", "Case"))) |>
    select(case, age, sex, length_followup, unique_phecodes, unique_encounters, SS_809, CV_401, EM_236) |>
    tbl_summary(
        by = case,
        label = list(
            age ~ "Age",
            length_followup ~ "Length of follow-up (days)",
            unique_phecodes ~ "Unique phecodes",
            unique_encounters ~ "Unique encounters"
        )) |>
    add_p() |>
    add_overall() |>
    bold_labels() |>
    modify_spanning_header(c("stat_1", "stat_2") ~ paste0("**", ms::pheinfox[phecode == "CA_107.2", description], "**")) |>
    modify_caption("**Table 1.** Demographic and clinical characteristics of cases and controls")

phecodes <- names(merged)[names(merged) %in% ms::pheinfox[, phecode]]
options(future.globals.maxSize = (object.size(merged) + 50 * 1024^2))
plan(multicore, workers = 16)

print((object.size(merged) + 50 * 1024^2), units = "auto")

quick_res <- phecodes |>
    future_map(
        \(x) phewas_analysis(
            data       = merged,
            outcome    = "case",
            exposure   = x,
            covariates = NULL,
            logistf_pl = TRUE,
            method     = "SPAtest"
        ),
        .progress = TRUE
    ) |>
    bind_rows() |>
    arrange(p_value)


quick_res_dt <- as.data.table(quick_res)
for (i in quick_res$phecode) {
    quick_res_dt[phecode == i, `:=` (
        n_case = merged[case == 1, ][[i]] |> sum(),
        n_ctrl = merged[case == 0, ][[i]] |> sum()
    )]
}

plot_phewasx(quick_res |> filter(is.finite(log10p)), phe_var = "phecode",
    title = paste0("PheWAS of cancer (", opt$outcome, ")"))

###
tmp <- phe_dsb_to_pim(tr_full_phe[[x]] |> filter(!(phecode %in% cancer_codes))) |>
    (\(y) {
        if (time_thresholds[x] == 0) {
            bind_rows(
                y,
                case_data$case[!(id %in% y$id), .(id)]
            ) |>
                (\(z) {
                    z[is.na(z)] <- 0
                    z
                })()
        } else {
            y
        }
    })() |>
    (\(a) {
        if (time_thresholds[x] == 0) {
            bind_rows(
                a,
                match_data[!(id %in% a$id), .(id)]
            ) |>
                (\(b) {
                    b[is.na(b)] <- 0
                    b
                })()
        } else {
            a
        }
    })() |>
    mutate(case = ifelse(id %in% case_data$case[, id], 1, 0))

match_data[, id %in% tmp$id] |> table()
tmp[, id %in% match_data$id] |> table()
