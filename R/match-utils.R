suppressPackageStartupMessages({
    library(data.table)
    library(glue)
    library(cli)
    library(tidyverse)
    library(broom)
    library(MatchIt)
})

plot_quick_tv_logor <- function(.data_list, outcome, exposure, cohort) {
    imap(.data_list,
        ~ {
            f <- as.formula(glue("case ~ {exposure}"))
            glm(f, family = quasibinomial(), data = .data_list[[.y]]) |>
                broom::tidy() |>
                filter(term != "(Intercept)") |>
                mutate(
                    time = .y,
                    lower = estimate - 1.96 * std.error,
                    upper = estimate + 1.96 * std.error
                )
        }) |>
        bind_rows() |>
            ggplot(aes(x = time, y = estimate, ymin = lower, ymax = upper)) +
            geom_pointrange() +
            geom_hline(yintercept = 0, linetype = 2) +
            labs(
                title = glue("{exposure} log(OR) (95% CI) for {outcome} status by time_threshold in {cohort}"),
                subtitle = glue("
                    {outcome}: {ms::pheinfox[phecode == outcome, description]}
                    {exposure}: {ms::pheinfox[phecode == exposure, description]}
                "),
                x = "Time threshold",
                y = "Log(OR) estimate",
                caption = glue("Matching ratio: 1:{opt$matching_ratio}")
            ) +
            theme_ms()

}

perform_matching <- function(
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
        paste0(c(nearest_vars, exact_vars), collapse = ' + '),
        ", data = matching_data, calclosest = TRUE, ",
        "mahvars = c(", paste0(sapply(nearest_vars, \(x) paste0('\'', x, '\'')), collapse = ', '), 
        "), caliper = ", match_caliper, ", exact = c(",
        paste0(sapply(exact_vars, \(x) paste0('\'', x, '\'')), collapse = ', '),
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

phe_dsb_to_pim <- function(phe_dsb_data) {
    pim <- dcast(
        as.data.table(phe_dsb_data)[, .(id, phecode)],
        id ~ phecode,
        value.var = "phecode",
        fun.aggregate = length,
        fill = 0
    )
    int_vars <- names(pim)[map_lgl(pim, is.numeric)]
    walk(int_vars, \(x) set(pim, j = x, value = as.numeric(pim[[x]] > 0)))
    pim
}

prepare_case_data <- function(
    phe_dsb_data,
    demo_data,
    case_definition       = "first_occurrence",
    outcome               = opt$outcome,
    malignant_phecodes    = cancer_phecodesx[keep == 1, phecode],
    specific_phecodes     = cancer_phecodesx[specific == 1, phecode],
    exclude_other_cancers = TRUE
) {

    outcome_sex <- ms::pheinfox[phecode == outcome, sex]
    if (outcome_sex %in% c("Female", "Male")) {
        if (outcome_sex == "Female") {
            exclude_ids <- unique(demo_data[sex %in% c("M", "Male"), id])
        } else {
            exclude_ids <- unique(demo_data[sex %in% c("F", "Female"), id])
        }
    } else {
        exclude_ids <- vector(mode = "character", length = 0L)
    }

    # restrict to first occurrence of outcome
    first_phe_dsb <- phe_dsb_data[ phe_dsb_data[, .I[which.min(dsb)], by = c("id", "phecode")]$V1 ]
    if (tolower(outcome_sex) != "both") {
        first_phe_dsb <- first_phe_dsb[id %in% demo_data[tolower(sex) %in% tolower(c(outcome_sex, substr(outcome_sex, 1, 1))), unique(id)]]
    }

    # restrict to cancers
    first_cancer_dsb <- first_phe_dsb[phecode %in% intersect(malignant_phecodes, specific_phecodes), ]

    # restrict to cases
    ## first occurrence of outcome
    if (tolower(case_definition) %in% c("first", "first_occurrence", "first occurrence")) {
        case_dsb <- first_cancer_dsb[ first_cancer_dsb[, .I[dsb == min(dsb, na.rm = TRUE)], id]$V1 ][phecode == outcome, .(id, index_dsb = dsb)] # first occurrence of any cancer
    }
    ## any occurrence of outcome
    if (tolower(case_definition) %in% c("any", "any_occurrence", "any occurrence", "all", "all_occurrences", "all occurrences")) {
        case_dsb <- first_cancer_dsb[phecode == outcome, .(id, index_dsb = dsb)]
    }

    # exclude other cancers?
    if (exclude_other_cancers) {
        exclude_ids <- unique(c(exclude_ids, setdiff(first_cancer_dsb[, id], case_dsb[, id])))
    }
    
    # make case pim
    case_pim0 <- phe_dsb_to_pim(
        left_join(
            first_phe_dsb[id %in% case_dsb[, id] & !(phecode %in% unique(c(malignant_phecodes, specific_phecodes))), ],
            case_dsb,
            by = "id"
        ) |>
            dplyr::filter(dsb <= index_dsb)
    )
    case_pim0 <- left_join(
        case_dsb,
        case_pim0,
        by = "id"
    ) |>
    (\(x) {
        x[is.na(x)] <- 0
        x
        })() |>
    dplyr::select(-one_of("index_dsb"))

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
        merge.data.table(case_follow_up, by = "id") |>
        merge.data.table(case_dsb[, .(id, age = round(index_dsb / 365.25, 1))], by = "id")
    
    # unmatched, unrestricted pim
    unmatched_pim <- first_phe_dsb[!(id %in% case_dsb[, id]) & !(id %in% exclude_ids), ] |>
        phe_dsb_to_pim()
    
    unmatched_follow_up <- phe_dsb_data[!(id %in% case_dsb[, id]) & !(id %in% exclude_ids), ] |>
        summarize(
            length_followup = max(dsb) - min(dsb),
            unique_encounters = length(unique(dsb)),
            unique_phecodes = length(unique(phecode)),
            .by = "id"
        )
    
    # unmatched demo
    unmatched_demo <- demo_data[!(id %in% case_dsb[, id]) & !(id %in% exclude_ids),
        .(id, age = age_at_last_diagnosisx, sex)] |>
        merge.data.table(unmatched_follow_up, by = "id")

    return(list(
        case           = case_dsb,
        case_pim0      = case_pim0,
        case_demo      = case_demo,
        exclude_ids    = exclude_ids,
        unmatched_pim  = unmatched_pim,
        unmatched_demo = unmatched_demo
    ))
}
