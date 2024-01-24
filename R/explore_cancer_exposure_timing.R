# libraries, functions, and options --------------------------------------------
ms::libri(
    data.table, ms, MatchIt, logistf, glue, qs, optparse, EValue, cli,
    tidyverse, ggstance, patchwork
)

set.seed(61787)

for (i in list.files("./fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
    make_option("--outcome",
        type = "character", default = "CA_105.1",
        help = "Outcome phecode [default = %default]"
    ),
    make_option("--exposure",
        type = "character", default = "EM_236",
        help = "Exposure phecode [default = %default]"
    ),
    make_option("--mgi_version",
        type = "character", default = "20220822",
        help = "Version of MGI data [default = %default]"
    ),
    make_option("--mgi_cohort",
        type = "character", default = "comb",
        help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"
    ),
    make_option("--ukb_version",
        type = "character", default = "20221117",
        help = "Version of UKB data [default = %default]"
    ),
    make_option("--time_thresholds",
        type = "character", default = "0,0.5,1,2,5",
        help = glue(
            "Time thresholds for the phenome data ",
            "[default = %default]"
        )
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])
cooccur_covs <- c("age", "female", "length_followup")

## extract file paths
file_paths <- get_files(
    mgi_version = opt$mgi_version,
    ukb_version = opt$ukb_version
)


# read data --------------------------------------------------------------------

## helper data
### codes with at least 20 cases in AOU, MGI, and UKB
common_codes <- fread("data/public/phecodex_20plus.csv")[plus20 == 1, phecode]
### qualifying cancer phecodes (malignant)
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")[
    keep == 1, phecode
]

## mgi
cli_progress_step("Reading MGI data")
### demographics
mgi_cov <- qread(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))

### time-stamped phenome data
mgi_full_phe <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs"))[phecode %in% common_codes, ]

mgi_first_phe <- mgi_full_phe[
    mgi_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
]

# identify cases and exclusions
mgi_first_cancer_dsb <- mgi_first_phe[phecode %in% cancer_phecodesx, ] # first occurrence of each cancer
mgi_case <- mgi_first_cancer_dsb[
    mgi_first_cancer_dsb[, .I[dsb == min(dsb, na.rm = TRUE)], id]$V1
][phecode == opt$outcome, ] # first occurrence of any cancer

mgi_case_index <- mgi_case[phecode == opt$outcome, ][, .(id, index_dsb = dsb)]

mgi_case_phe <- merge.data.table(
    mgi_first_phe[id %in% mgi_case[, id], ],
    mgi_case_index,
    by = "id"
)

mgi_first_phe <- mgi_first_phe[id %in% mgi_case_index[, id] | id %in% setdiff(mgi_first_phe[, id], mgi_first_cancer_dsb[, id]), ]

if (ms::pheinfox[phecode == opt$outcome, sex] %in% c("Female", "Male")) {
   mgi_first_phe <-  mgi_first_phe[id %in% mgi_cov[sex == ms::pheinfox[phecode == opt$outcome, substr(sex, 1, 1)], id], ]
}



### UKB
ukb_cov <- qread(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))

### time-stamped phenome data
ukb_full_phe <- qread(glue("data/private/ukb/{opt$ukb_version}/UKB_FULL_PHECODEX_DSB_{opt$ukb_version}.qs"))[phecode %in% common_codes, ]

ukb_first_phe <- ukb_full_phe[
    ukb_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
]

# identify cases and exclusions
ukb_first_cancer_dsb <- ukb_first_phe[phecode %in% cancer_phecodesx, ] # first occurrence of each cancer
ukb_case <- ukb_first_cancer_dsb[
    ukb_first_cancer_dsb[, .I[dsb == min(dsb, na.rm = TRUE)], id]$V1
][phecode == opt$outcome, ] # first occurrence of any cancer

ukb_case_index <- ukb_case[phecode == opt$outcome, ][, .(id, index_dsb = dsb)]

ukb_case_phe <- merge.data.table(
    ukb_first_phe[id %in% ukb_case[, id], ],
    ukb_case_index,
    by = "id"
)

ukb_first_phe <- ukb_first_phe[id %in% ukb_case_index[, id] | id %in% setdiff(ukb_first_phe[, id], ukb_first_cancer_dsb[, id]), ]

if (ms::pheinfox[phecode == opt$outcome, sex] %in% c("Female", "Male")) {
    ukb_first_phe <- ukb_first_phe[id %in% ukb_cov[sex == ms::pheinfox[phecode == opt$outcome, sex], id], ]
}

###

time_between <- function(.data, x, y) {
    first_x <- .data[phecode == x, ][, .(id, x_dsb = dsb)]
    first_y <- .data[phecode == y, ][, .(id, y_dsb = dsb)]
    inner_join(first_x, first_y, by = "id")[
        , .(id, time_between = x_dsb - y_dsb)
    ]
}

plot_time_between <- function(.data, case, exposure, save = FALSE, path = NULL) {
    d <- .data |>
        time_between(exposure, case) |>
        arrange(desc(time_between)) |>
        mutate(
            order = 1:n(),
            time = case_when(
                time_between >= 31 ~ "After",
                time_between <= -31 ~ "Before",
                TRUE ~ "Within 1 month"
            ),
            time = factor(time, levels = c("Before", "Within 1 month", "After"))
        )

    cols <- c(
        "Before"         = "#009E73",
        "Within 1 month" = "#1a1a1a",
        "After"          = "#D55E00"
    )

    d_prop <- d[, .N, time][, prop := round(N / sum(N) * 100, 1)][, label := paste0(time, ": ", prop, "%")]

    time_between_plot <- d |>
        ggplot(aes(x = order, y = time_between)) +
        geom_segment(
            aes(x = order, xend = order, y = 0, yend = time_between, color = time)
        ) +
        coord_flip() +
        geom_text(
            aes(
                x = 0, y = min(d[, time_between]),
                label = paste0(d_prop[, label], collapse = "\n")
            ),
            hjust = 0, vjust = 0
        ) +
        scale_color_manual(values = cols) +
        scale_y_continuous(labels = scales::comma) +
        labs(
            title    = glue("Time between first occurrence of {case} and {exposure}"),
            subtitle = glue("{case}: {ms::pheinfox[phecode == case, description]}\n{exposure}: {ms::pheinfox[phecode == exposure, description]}"),
            x        = "Participant",
            y        = glue("Days ({case} - {exposure})"),
            caption  = glue("
        {round((nrow(mgi_case_index) - nrow(d))/nrow(mgi_case_index) * 100, 1)}% of {case} cases never had {exposure} diagnosis.
        ")
        ) +
        theme_ms() +
        theme(
            legend.title = element_blank()
        )

    if (save) {
        if (is.null(path)) {
            cli_alert_info("No path specified. Saving to bin/")
            path <- "bin/"
        }
        ggsave(
            plot = time_between_plot,
            filename = glue("{path}{case}_{exposure}_between.pdf"),
            width = 6, height = 6, device = cairo_pdf)
    }

    return(time_between_plot)
}

plot_case_distribution <- function(.data, case, save = TRUE, path = NULL, subtitle = FALSE) {
    case_plot <- .data |>
        mutate(age_at_dx = round(index_dsb / 365.25, 1)) |>
        ggplot(aes(x = age_at_dx, y = -0.1)) +
        geom_boxploth(width = 0.1, fill = "#D55E00") +
        geom_density(aes(x = age_at_dx, y = stat(scaled)),
            linewidth = 1, inherit.aes = FALSE, fill = "#D55E00", alpha = 0.6
        ) +
        labs(
            title = glue("Distribution of age at {case} diagnosis"),
            x = glue("Age at {case} diagnosis"),
            y = "Density"
        ) +
        (\() if (subtitle) labs(subtitle = glue("{case}: {ms::pheinfox[phecode == case, description]}")))() +
        theme_ms()

     if (save) {
         if (is.null(path)) {
             cli_alert_info("No path specified. Saving to bin/")
             path <- "bin/"
         }
         ggsave(
             plot = case_plot,
             filename = glue("{path}{case}_dist_plot.pdf"),
             width = 6, height = 6, device = cairo_pdf
         )
     }

     return(case_plot)
}

plot_exposure_distribution <- function(.data, case, exposure, case_index_data, save = FALSE, path = NULL, vline = TRUE, subtitle = FALSE) {          
    exp_dist_data <- .data |>
        filter(phecode == exposure) |>
        mutate(
            age_at_dx = round(dsb / 365.25, 1),
            case = ifelse(id %in% case_index_data[, id], "Case", "Control")
        )
    cols2 <- c(
        "Case" = "#D55E00",
        "Control" = "#009E73"
    )
    exposure_plot <- exp_dist_data |>
        ggplot(aes(x = age_at_dx, y = -0.2, group = case, fill = case)) +
        geom_boxploth(width = 0.2, show.legend = FALSE) +
        geom_density(aes(x = age_at_dx, y = stat(scaled), fill = case),
            alpha = 0.6,
            linewidth = 1, inherit.aes = FALSE
        ) +
        (\() if (vline) geom_vline(
            aes(xintercept = median(case_index_data[, index_dsb / 365.25])),
            color = "#1a1a1a", linetype = 2, linewidth = 1
        ))() +
        scale_fill_manual(values = cols2) +
        labs(
            title    = glue("Distribution of age at {exposure} diagnosis"),
            x        = glue("Age at {exposure} diagnosis"),
            y        = "Density",
            caption  = glue("Vertical dashed line represents median age at {case} diagnosis")
        ) +
        (\() if (subtitle) labs(subtitle = glue("{exposure}: {ms::pheinfox[phecode == exposure, description]}")))() +
        theme_ms() +
        theme(
            legend.title = element_blank()
        )
    
    if (save) {
        if (is.null(path)) {
            cli_alert_info("No path specified. Saving to bin/")
            path <- "bin/"
        }
        ggsave(plot = exposure_plot,
            filename = glue("{path}{exposure}_dist_plot.pdf"),
            width = 6, height = 6, device = cairo_pdf
        )
    }

    return(exposure_plot)
}

plot_age_at_first_dx <- function(.data, case, case_index_data, vline = TRUE, save = FALSE, path = NULL) {
    age_at_first_dx_data <- .data |>
        mutate(
            case = ifelse(id %in% case_index_data[, id], "Case", "Control")
        )
    cols2 <- c(
        "Case" = "#D55E00",
        "Control" = "#009E73"
    )
    age_at_first_dx_plot <- age_at_first_dx_data |>
        ggplot(aes(x = age_at_first_diagnosisx, y = -0.2, group = case, fill = case)) +
        geom_boxploth(width = 0.2, show.legend = FALSE) +
        geom_density(aes(x = age_at_first_diagnosisx, y = stat(scaled), fill = case),
            alpha = 0.6,
            linewidth = 1, inherit.aes = FALSE
        ) +
        (\() if (vline) {
            geom_vline(
                aes(xintercept = median(case_index_data[, index_dsb / 365.25])),
                color = "#1a1a1a", linetype = 2, linewidth = 1
            )
        })() +
        scale_fill_manual(values = cols2) +
        labs(
            title    = glue("Distribution of age at first diagnosis"),
            x        = glue("Age at first diagnosis"),
            y        = "Density",
            caption  = glue("Vertical dashed line represents median age at {case} diagnosis")
        ) +
        theme_ms() +
        theme(
            legend.title = element_blank()
        )

    if (save) {
        if (is.null(path)) {
            cli_alert_info("No path specified. Saving to bin/")
            path <- "bin/"
        }
        ggsave(
            plot = age_at_first_dx_plot,
            filename = glue("{path}{case}_age_at_first_dx_dist_plot.pdf"),
            width = 6, height = 6, device = cairo_pdf
        )
    }

    return(age_at_first_dx_plot)
}

# MGI --------------------------------------------------------------------------
time_between_plot <- mgi_case_phe |>
    plot_time_between(
        case = opt$outcome,
        exposure = opt$exposure,
        save = FALSE
    )

case_dist_plot <- mgi_case_index |>
    plot_case_distribution(case = opt$outcome, save = FALSE)

exposure_dist_plot <- mgi_first_phe |>
    plot_exposure_distribution(
        case            = opt$outcome,
        exposure        = opt$exposure,
        case_index_data = mgi_case_index,
        save            = FALSE,
        vline           = TRUE
    )

age_at_first_dx_plot <- mgi_cov |>
    plot_age_at_first_dx(case = opt$outcome, case_index_data = mgi_case_index)

patched <- (time_between_plot / case_dist_plot / exposure_dist_plot / age_at_first_dx_plot) +
    plot_layout(height = c(3, 1, 1, 1))

ggsave(
    plot = patched,
    filename = glue("bin/{opt$outcome}_{opt$exposure}_time_and_dist_plot.pdf"),
    width = 6, height = 10, device = cairo_pdf
)

# UKB --------------------------------------------------------------------------
ukb_time_between_plot <- ukb_case_phe |>
    plot_time_between(
        case = opt$outcome,
        exposure = opt$exposure,
        save = FALSE
    )

ukb_case_dist_plot <- ukb_case_index |>
    plot_case_distribution(case = opt$outcome, save = FALSE)

ukb_exposure_dist_plot <- ukb_first_phe |>
    plot_exposure_distribution(
        case            = opt$outcome,
        exposure        = opt$exposure,
        case_index_data = ukb_case_index,
        save            = FALSE,
        vline           = TRUE
    )

ukb_age_at_first_dx_plot <- ukb_cov |>
    plot_age_at_first_dx(case = opt$outcome, case_index_data = ukb_case_index)

ukb_patched <- (ukb_time_between_plot / ukb_case_dist_plot / ukb_exposure_dist_plot / ukb_age_at_first_dx_plot) +
    plot_layout(height = c(3, 1, 1, 1))

ggsave(
    plot = ukb_patched,
    filename = glue("bin/ukb_{opt$outcome}_{opt$exposure}_time_and_dist_plot.pdf"),
    width = 6, height = 10, device = cairo_pdf
)

