suppressPackageStartupMessages({
    library(data.table)
    library(survey)
    library(ms)
    library(cli)
    library(scales)
})

# poststrat function -----------------------------------------------------------
poststrat_nhanes <- function(
    int_data,
    nhanes_data    = NULL,
    id_var         = "id",
    age_var        = "AgeLastEntry",
    age_bin        = FALSE,
    age_bin_breaks = c(0, 18, 35, 65, 80, 150),
    covs           = c("age_bin", "cad", "smoker", "diabetes", "female"),
    not_num_vars   = NULL,
    chop = TRUE) {

    if (age_bin) {
        cli_alert("constructing categorical age variable using breaks: {age_bin_breaks}")
        covs <- unique(c(covs, "age_bin"))
    }
    cli_alert_info("estimating poststratification weights for covariates: {.field {covs}}")
    if (chop) cli_alert_info("truncating weights to 2.75 and 97.5 percentiles")

    # 1. load and prep nhanes data ----------------------------------------------
    if (is.null(nhanes_data)) {
        nhanes <- download_nhanes_data(datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ", "BPQ", "DPQ"))
        phanes <- prepare_nhanes_data(nhanes)
    } else {
        phanes <- copy(nhanes_data)
    }
    phanes[, `:=`(
        smoker = ifelse(smoking_former == 1 | smoking_current == 1, 1, 0),
        age_bin = cut(age, c(seq(0, 80, by = 10), 150), right = FALSE)
    )]
    setnames(phanes, "nhanes_nhw", "nhw", skip_absent = TRUE)

    if (age_bin == TRUE) {
        phanes[, age_bin := cut(age, age_bin_breaks, right = FALSE)]
    }

    nhanes_design <- svydesign(
        id      = ~psu_nhanes,
        strata  = ~strata_nhanes,
        weights = ~weight_nhanes,
        nest    = TRUE,
        data    = phanes
    )

    # 2. get nhanes proportions -------------------------------------------------
    population_proportions <- svytable(
        formula = as.formula(paste0("~", paste(covs, collapse = " + "))),
        design  = nhanes_design
    ) |>
    prop.table() |>
    as.data.table()
    setnames(population_proportions, "N", "pop_prob")

    if (!is.null(not_num_vars)) {
        num_vars <- setdiff(covs, c(not_num_vars, "age_bin"))
    } else {
        num_vars <- setdiff(covs, "age_bin")
    }
    
    population_proportions[, (num_vars) := lapply(.SD, as.numeric), .SDcols = num_vars]

    # 3. get internal proportions -----------------------------------------------
    if (age_bin == TRUE) {
        int_data[, age_bin := cut(get(age_var), age_bin_breaks, right = FALSE)]
    }

    internal_probabilities <- int_data[, ..covs] |>
        table() |>
        prop.table() |>
        as.data.table()
    setnames(internal_probabilities, "N", "int_prob")
    internal_probabilities[, (num_vars) := lapply(.SD, as.numeric), .SDcols = num_vars]

    # 4. merge probabilities into internal data ---------------------------------
    sub_vars <- c(id_var, covs)
    merged <- int_data[, ..sub_vars] |>
        merge.data.table(population_proportions, by = covs) |>
        merge.data.table(internal_probabilities, by = covs)
    
    merged[, ps_weight := pop_prob / int_prob]

    # 5. process ----------------------------------------------------------------
    if (chop == TRUE)  merged[, ps_weight := chopr(ps_weight)]

    merged[, ps_weight := (.N * ps_weight) / sum(ps_weight, na.rm = TRUE)]

    return_vars <- c(id_var, "ps_weight")
    return(
        merged[, ..return_vars]
    )
    
}
