prepare_matched_phenomes <- function(
    phe_dsb_data,
    demo_data,
    outcome,
    case_definition       = "first_occurrence", # "first" | "first_occurrence" or "all" | "all_occurrences" | "any" | "anytime" | "any_occurrence"
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
) {
    
    # PREPARE CASE DATA
    cli_progress_step("Preparing case data")
    case_data <- prepare_case_data(
        phe_dsb_data          = phe_dsb_data,
        demo_data             = demo_data,
        outcome               = outcome,
        case_definition       = case_definition,
        malignant_phecodes    = malignant_phecodes,
        specific_phecodes     = specific_phecodes,
        exclude_other_cancers = exclude_other_cancers
    )

    # MATCH
    if (matched) {
        cli_progress_step("Matching")
        match_data <- perform_matching(
            matching_data   = demo_data[!id %in% case_data$exclude_ids, ],
            case_data       = case_data$case,
            match_ratio     = match_ratio,
            nearest_vars    = nearest_vars,
            exact_vars      = exact_vars,
            match_caliper   = match_caliper,
            time_thresholds = time_thresholds
        )
        full_phe_sub <- phe_dsb_data[id %in% match_data[, id], ]
        cli_progress_step("Preparing time-restricted phenomes")
        merged_full_phe <- merge.data.table(
            full_phe_sub,
            match_data[, .(id, index_dsb)],
            by = "id",
            all.x = TRUE
        )
        tr_full_phe <- map(
            seq_along(time_thresholds),
            \(x) {
                merged_full_phe |>
                    dplyr::filter(dsb <= index_dsb - (365.25 * time_thresholds[x])) |>
                    dplyr::select(-one_of("index_dsb")) 
            }
        )
        tr_demo <- map(
            seq_along(time_thresholds),
            \(x) {
                tr_full_phe[[x]] |>
                    summarize(
                        length_followup   = max(dsb) - min(dsb),
                        unique_encounters = length(unique(dsb)),
                        unique_phecodes   = length(unique(phecode)),
                        .by = "id"
                    ) |>
                    (\(y) {
                        if (time_thresholds[x] == 0) {
                            bind_rows(
                                y,
                                match_data[!(id %in% y$id), .(id)]
                            ) |>
                                (\(z) {
                                    z[is.na(z)] <- 0
                                    z
                                })()
                        } else {
                            y
                        }
                    })() |>
                    full_join(
                        match_data[, .(id, age = round((index_dsb - (365.25 * time_thresholds[x])) / 365.25, 1))],
                        by = "id"
                    ) |>
                    full_join(
                        demo_data[, .(id, sex)],
                        by = "id"
                    )
            } 
        )
        tr_pim <- map(
            seq_along(time_thresholds),
            \(x) {
                phe_dsb_to_pim(tr_full_phe[[x]] |> filter(!(phecode %in% unique(c(malignant_phecodes, specific_phecodes))))) |>
                    (\(y) {
                        if (time_thresholds[x] == 0) {
                            bind_rows(
                                y,
                                match_data[!(id %in% y$id), .(id)]
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
            }
        )
        names(tr_pim) <- paste0("t", time_thresholds)

        merged_tr <- map(
            seq_along(time_thresholds),
            \(x) {
                tr_pim[[x]] |>
                    left_join(
                        tr_demo[[x]],
                        by = join_by(id)
                    )
            }
        )

        out <- list(
            case_data  = case_data,
            match_data = match_data,
            tr_pim     = merged_tr
        )

    } else {
        cli_progress_step("Preparing time-restricted phenomes")
        full_phe_sub  <- phe_dsb_data[!(id %in% case_data$exclude_ids), ]
        case_full_phe <- full_phe_sub[id %in% case_data$case[, id], ]
        case_full_phe <- merge.data.table(case_full_phe, case_data$case[, .(id, case_dsb = index_dsb)], by = "id", all.x = TRUE)
        tr_case_phe <- map(
            seq_along(time_thresholds),
            \(x) {
                case_full_phe[dsb <= case_dsb - (365.25 * time_thresholds[x]), ][, !c("case_dsb")]
            }
        )
        tr_case_pim <- map(
            seq_along(time_thresholds),
            \(x) {
                phe_dsb_to_pim(tr_case_phe[[x]] |>
                    (\(q) if (exclude_other_cancers) {
                        q |>
                        filter(!(phecode %in% unique(c(malignant_phecodes, specific_phecodes))))
                    } else {
                        q
                    })()) |>
                (\(y) {
                    if (time_thresholds[x] == 0) {
                        bind_rows(
                            y,
                            case_data$case[!(id %in% y$id), ]
                        ) |>
                            (\(z) {
                                z[is.na(z)] <- 0
                                z
                            })()
                    } else {
                        y
                    }
                })()
            }
        )

        tr_case_demo <- map(
            seq_along(time_thresholds),
            \(x) {
                tr_case_phe[[x]] |>
                    summarize(
                        length_followup   = max(dsb) - min(dsb),
                        unique_encounters = length(unique(dsb)),
                        unique_phecodes   = length(unique(phecode)),
                        .by = "id"
                    ) |>
                    left_join(
                        case_data$case[, .(id, age = round((index_dsb - (365.25 * time_thresholds[x])) / 365.25, 1))],
                        by = "id"
                    ) |>
                    left_join(
                        demo_data[, .(id, sex)],
                        by = "id"
                    )
            }
        )
        control_demo <- full_phe_sub[!(id %in% case_data$case[, id]), ] |>
            summarize(
                length_followup   = max(dsb) - min(dsb),
                unique_encounters = length(unique(dsb)),
                unique_phecodes   = length(unique(phecode)),
                age               = round(max(dsb) / 365.25, 1),
                .by = "id"
            ) |>
            left_join(
                demo_data[, .(id, sex)],
                by = "id"
            )
        control_pim <- full_phe_sub[!(id %in% case_data$case[, id]), ] |>
            phe_dsb_to_pim()
        
        tr_pim <- map(
            seq_along(time_thresholds), 
            \(x) bind_rows(
                tr_case_pim[[x]],
                control_pim
            ) |>
            left_join(
                bind_rows(
                    tr_case_demo[[x]],
                    control_demo
                ),
                by = join_by(id)
            ) |>
            mutate(
                case = ifelse(id %in% case_data$case[, id], 1, 0)
            )
        )
        names(tr_pim) <- paste0("t", time_thresholds)

        out <- list(
            case_data  = case_data,
            tr_pim     = tr_pim
        )
    }

    return(out)

}

test <- prepare_matched_phenomes(
    phe_dsb_data          = mgi_full_phe,
    demo_data             = mgi_demo,
    outcome               = opt$outcome,
    case_definition       = "first_occurrence", # "first" | "first_occurrence" or "all" | "all_occurrences" | "any" | "anytime" | "any_occurrence"
    cancer_codes          = cancer_phecodesx,
    matched               = TRUE,
    match_ratio           = 2,
    match_caliper         = 0.25,
    nearest_vars          = c("age", "length_followup"),
    exact_vars            = c("female"),
    control_definition    = "no_cancer",
    time_thresholds       = c(0, 0.5, 1, 2, 3, 5),
    exclude_other_cancers = TRUE
)

mgi_test <- prepare_matched_phenomes(
    phe_dsb_data          = mgi_full_phe,
    demo_data             = mgi_demo,
    outcome               = opt$outcome,
    case_definition       = "all", # "first" | "first_occurrence" or "all" | "all_occurrences" | "any" | "anytime" | "any_occurrence"
    cancer_codes          = cancer_phecodesx,
    matched               = TRUE,
    match_ratio           = 2,
    match_caliper         = 0.25,
    nearest_vars          = c("age", "length_followup"),
    exact_vars            = c("female"),
    control_definition    = "no_cancer",
    time_thresholds       = c(0, 0.5, 1, 2, 3, 5),
    exclude_other_cancers = FALSE
)

ukb_test <- prepare_matched_phenomes(
    phe_dsb_data          = ukb_full_phe,
    demo_data             = ukb_demo |>
        rename(age = age_at_consent) |>
        mutate(
            length_followup = last_dsbx - first_dsbx,
            female = as.numeric(sex == "Female")
        ) |>
        drop_na(age, length_followup, female),
    outcome               = opt$outcome,
    case_definition       = "first_occurrence", # "first" | "first_occurrence" or "all" | "all_occurrences" | "any" | "anytime" | "any_occurrence"
    cancer_codes          = cancer_phecodesx,
    matched               = FALSE,
    match_ratio           = 2,
    match_caliper         = 0.25,
    nearest_vars          = c("age", "length_followup"),
    exact_vars            = c("female"),
    control_definition    = "no_cancer",
    time_thresholds       = c(0, 0.5, 1, 2, 3, 5),
    exclude_other_cancers = TRUE
)

ukb_demo |>
    rename(age = age_at_consent) |>
    mutate(
        length_followup = last_dsbx - first_dsbx,
        female = as.numeric(sex == "Female")
    ) |>
    select(age, length_followup, female) |>
    map_dbl(\(x) sum(is.na(x)))

map(
    seq_along(test$tr_pim),
    \(x) {
        glm(case ~ EM_236, data = test$tr_pim[[x]], family = quasibinomial()) |>
            broom::tidy() |>
            filter(term == "EM_236") |>
            mutate(time = names(test$tr_pim)[x])
    }
) |> bind_rows()

table(
    case = test$tr_pim[[1]][, case],
    test = test$tr_pim[[1]][, EM_236]
)
