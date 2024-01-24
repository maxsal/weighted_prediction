### LIBRARIES
suppressPackageStartupMessages({
    library(data.table)
    library(glue)
    library(cli)
    library(tidyverse)
    library(broom)
    library(MatchIt)
})

### HELPER FUNCTIONS
### REMOVE DECIMAL PLACES
.remove_decimal_places <- function(x) {
    # split the string into two parts: base and decimal part.
    split_string <- strsplit(x, ".", fixed = TRUE)[[1]]

    # If the string doesn't contain a decimal part, return the original string
    if (length(split_string) < 2) {
        return(x)
    }

    # otherwise, iterate to remove each decimal digit progressively
    base_part <- split_string[1]
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
### PHE DSB TO PIM
phe_dsb_to_pim <- function(phe_dsb_data) {
    pim <- dcast(
        as.data.table(phe_dsb_data)[, .(id, phecode)],
        id ~ phecode,
        value.var     = "phecode",
        fun.aggregate = length,
        fill          = 0
    )
    int_vars <- names(pim)[map_lgl(pim, is.numeric)]
    walk(int_vars, \(x) set(pim, j = x, value = as.numeric(pim[[x]] > 0)))
    pim[]
}


### PREPARE CASE DATA
.prepare_case_data <- function(
    phe_dsb_data,
    demo_data,
    case_definition       = "first_occurrence",
    outcome               = opt$outcome,
    malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
    specific_phecodes     = cancer_phecodesx[specific == 1, phecode],
    exclude_other_cancers = TRUE
) {
    # identify whether outcome phecode is sex-specific and, if so, exclude
    # individuals with discordant sex
    outcome_sex <- ms::pheinfox[phecode == outcome, sex]
    exclude_ids <- switch(
        outcome_sex,
        "Female" = unique(demo_data[sex %in% c("M", "Male"), id]),
        "Male"   = unique(demo_data[sex %in% c("F", "Female"), id]),
        "Both"   = vector(mode = "character", length = 0L)
    )

    # restrict to first occurrence of outcome
    first_phe_dsb <- phe_dsb_data[
        phe_dsb_data[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
    ]
    if (tolower(outcome_sex) != "both") {
        first_phe_dsb <- first_phe_dsb[!(id %in% exclude_ids), ]
    }

    # restrict to cancers
    first_cancer_dsb <- first_phe_dsb[
        phecode %in% intersect(malignant_phecodes, specific_phecodes),
    ]

    # restrict to cases
    ## first occurrence of outcome
    if (tolower(case_definition) %in%
            c("first", "first_occurrence", "first occurrence")
        ) {
            case_dsb <- first_cancer_dsb[
                first_cancer_dsb[, .I[dsb == min(dsb, na.rm = TRUE)], id]$V1, 
            ][phecode == outcome, .(id, index_dsb = dsb)]
    }
    ## any occurrence of outcome
    if (tolower(case_definition) %in%
            c("any", "any_occurrence", "any occurrence",
              "all", "all_occurrences", "all occurrences")
        ) {
        case_dsb <- first_cancer_dsb[phecode == outcome, .(id, index_dsb = dsb)]
    }

    # exclude other cancers?
    if (exclude_other_cancers) {
        exclude_ids <- unique(
            c(exclude_ids, setdiff(first_cancer_dsb[, id], case_dsb[, id]))
        )
    }

    # make case pim
    case_pim0 <- phe_dsb_to_pim(
        left_join(
            first_phe_dsb |>
                dplyr::filter(id %in% case_dsb[, id]) |> # restrict to cases
                dplyr::filter(!(phecode %in% 
                                unique(
                                    c(
                                        malignant_phecodes,
                                        specific_phecodes
                                        )
                                    )
                                )
                            ), # exclude other cancers diagnosed at t=0
            case_dsb,
            by = "id"
        ) |>
            dplyr::filter(dsb <= index_dsb)
    )

    case_pim0 <- left_join(
        case_dsb[, !c("index_dsb")],  # add individuals with no phecodes at t=0
        case_pim0,
        by = "id"
    ) |>
        (\(x) {
            x[is.na(x)] <- 0 # replace missing with 0
            x
        })() 

    # length follow up
    case_follow_up <- left_join(
        phe_dsb_data[id %in% case_dsb[, id], ],
        case_dsb,
        by = "id"
    ) |>
        filter(dsb <= index_dsb) |>
        summarize(
            length_followup   = max(dsb) - min(dsb),
            unique_encounters = length(unique(dsb)),
            unique_phecodes   = length(unique(phecode)),
            .by = "id"
        )

    case_demo <- demo_data[id %in% case_dsb[, id], .(id, sex)] |>
        left_join(
            case_follow_up,
            by = "id"
        ) |>
        left_join(
            case_dsb[, .(id, age = round(index_dsb / 365.25, 1))],
            by = "id"
        )

    # unmatched, unrestricted pim
    unmatched_pim <- first_phe_dsb |>
        dplyr::filter(!(id %in% c(case_dsb[, id], exclude_ids))) |>
        phe_dsb_to_pim()

    unmatched_follow_up <- phe_dsb_data |>
        dplyr::filter(!(id %in% c(case_dsb[, id], exclude_ids))) |>
        summarize(
            length_followup   = max(dsb) - min(dsb),
            unique_encounters = length(unique(dsb)),
            unique_phecodes   = length(unique(phecode)),
            .by = "id"
        )

    # unmatched demo
    unmatched_demo <- demo_data |>
        dplyr::filter(!(id %in% c(case_dsb[, id], exclude_ids))) |>
        mutate(age = age_at_last_diagnosisx) |>
        select(id, age, sex) |>
        left_join(
            unmatched_follow_up,
            by = "id"
        )

    return(list(
        case           = case_dsb,
        case_pim0      = case_pim0,
        case_demo      = case_demo,
        exclude_ids    = exclude_ids,
        unmatched_pim  = unmatched_pim,
        unmatched_demo = unmatched_demo
    ))
}

### PERFORM MATCHING
.perform_matching <- function(
    matching_data,
    case_data,
    match_ratio     = 2,
    nearest_vars    = c("age", "length_followup"),
    exact_vars      = "sex",
    match_caliper   = 0.25,
    time_thresholds = NULL
) {
    matching_data[, case := ifelse(id %in% case_data[, id], 1, 0)]
    match_text <- paste0(
        "MatchIt::matchit(case ~ ",
        paste0(c(nearest_vars, exact_vars), collapse = " + "),
        ", data = matching_data, calclosest = TRUE, ",
        "mahvars = c(",
        paste0(sapply(nearest_vars, \(x) paste0("'", x, "'")), collapse = ", "),
        "), caliper = ", match_caliper, ", exact = c(",
        paste0(sapply(exact_vars, \(x) paste0("'", x, "'")), collapse = ", "),
        "), ratio = ", match_ratio, ")"
    )

    match   <- eval(parse(text = match_text))
    matched <- MatchIt::match.data(match)

    post_match_cov <- merge.data.table(
        matched,
        merge.data.table(
            case_data[, .(id, index_dsb)],
            matched,
            by    = "id",
            all.x = TRUE
        )[, .(subclass, index_dsb)],
        by = "subclass"
    )

    if (!is.null(time_thresholds)) {
        walk(
            time_thresholds,
            \(i) {
                post_match_cov[, (glue("t{i}_threshold")) :=
                    floor(index_dsb - (365.25 * i))]
            }
        )

        walk(
            time_thresholds,
            \(i) {
                post_match_cov[, (glue("t{i}_indicator")) :=
                    as.numeric(all(get(glue("t{i}_threshold")) > first_dsb)),
                by = subclass
                ]
            }
        )
    }

    return(post_match_cov[])
}


### PRIMARY FUNCTION
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
    case_data <- .prepare_case_data(
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
        match_data <- .perform_matching(
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
        names(merged_tr) <- paste0("t", time_thresholds)

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
