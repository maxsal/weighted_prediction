# replace NA with value --------------------------------------------------------
## inspiration from: https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
replace_missing <- function(data, cols = NULL, new_value = 0) {
  if (!is.null(cols)) {
    if (!all(cols %in% names(data))) {
      missing_cols <- col[!(cols %in% names(data))]
      message(paste0("Some names in cols argument not present in data set: ", paste0(missing_cols, collapse = ", ")))
    }
    for (i in cols) {
      data[is.na(get(i)), (i) := new_value]
    }
  } else {
    message(paste0("Replacing missing in *all* columns with: ", new_value))
    for (i in names(data)) {
      data[is.na(get(i)), (i) := new_value]
    }
  }
}

# replace an existing value with NA --------------------------------------------
make_missing <- function(data, cols = NULL, old_value = "Unknown") {
  if (!is.null(cols)) {
    if (!all(cols %in% names(data))) {
      missing_cols <- col[!(cols %in% names(data))]
      message(paste0("Some names in cols argument not present in data set: ", paste0(missing_cols, collapse = ", ")))
    }
    for (i in cols) {
      data[get(i) == old_value, (i) := NA]
    }
  } else {
    message(paste0("Replacing '", old_value, "' in *all* columns with missing (NA)"))
    for (i in names(data)) {
      data[get(i) == old_value, (i) := NA]
    }
  }
}

# quickly generate an age group prevalence table
age_grp_table <- function(lower_ages,
                          upper_char = "+",
                          upper_offset = 0,
                          upper_val = 150,
                          num_vec,
                          num_var_name = "prevalences") {
  out <- data.table(
    group      = paste0(lower_ages, c(paste0("-", lower_ages[-1] - 1), upper_char)),
    lower      = lower_ages,
    upper      = c(lower_ages[-1] - upper_offset, upper_val),
    num_var    = num_vec
  )
  setnames(out, "num_var", num_var_name)
  out
}

## helper for truncating probabilities
chopr <- function(x, probs = c(0.025, 0.975)) {
  quant <- quantile(x, probs = probs, na.rm = TRUE)
  x[x < quant[1L]] <- quant[1L]
  x[x > quant[2L]] <- quant[2L]
  x
}

## help summarize demographic data
demo_summarizr <- function(x,
                           age_var = "AgeLastEntry",
                           age_cats = seq(0, 80, 10),
                           female_var = "female",
                           nhw_var = "nhanes_nhw",
                           cancer_var = "cancer") {
  # continuous age (pretty_print() is in eval-utils.R)--------------------------
  age_cont <- x[, .(mean = mean(get(age_var), na.rm = TRUE), sd = sd(get(age_var), na.rm = TRUE))][, `:=`(
    print    = paste0(pretty_round(mean, 1), " (", pretty_round(sd, 1), ")"),
    var_name = age_var,
    variable = "Age (continuous)"
  )]

  # categorical age (age_grp_table() is in cleaning-utils.R) -------------------
  tmp_age_cat <- age_grp_table(
    lower_ages = age_cats,
    num_vec = rep(NA, length(age_cats)),
    num_var_name = "counts"
  )

  age_cat <- data.table(
    sub_variable = factor(
      cut(x[[age_var]],
        breaks = c(0, tmp_age_cat[["upper"]]),
        labels = tmp_age_cat[["group"]],
        right  = FALSE
      ),
      levels = tmp_age_cat[["group"]]
    )
  )[, .N, sub_variable][, `:=`(
    prop = round(N * 100 / x[, .N], 1)
  )][, `:=`(
    print    = paste0(prop, " (", trimws(format(N, big.mark = ",")), ")"),
    variable = "Age (categorical)"
  )][order(sub_variable), ][]

  # female ---------------------------------------------------------------------
  fem_tab <- x[, .N, female_var][, `:=`(
    prop         = round(N * 100 / x[, .N], 1),
    sub_variable = get(female_var)
  )][, `:=`(
    print    = paste0(prop, " (", trimws(format(N, big.mark = ",")), ")"),
    variable = "Female",
    var_name = female_var
  )]

  # Non-Hispanic White ---------------------------------------------------------
  nhw_tab <- x[, .N, nhw_var][, `:=`(
    prop         = round(N * 100 / x[, .N], 1),
    sub_variable = get(nhw_var)
  )][, `:=`(
    print    = paste0(prop, " (", trimws(format(N, big.mark = ",")), ")"),
    variable = "Non-Hispanic White",
    var_name = nhw_var
  )]

  # Cancer ---------------------------------------------------------------------
  cancer_tab <- x[, .N, cancer_var][, `:=`(
    prop         = round(N * 100 / x[, .N], 1),
    sub_variable = get(cancer_var)
  )][, `:=`(
    print    = paste0(prop, " (", trimws(format(N, big.mark = ",")), ")"),
    variable = "Cancer",
    var_name = cancer_var
  )]

  rbindlist(list(
    age_cont,
    age_cat,
    fem_tab,
    nhw_tab,
    cancer_tab
  ), fill = TRUE, use.names = TRUE)[, .(var_name, variable, sub_variable, print)]
}
