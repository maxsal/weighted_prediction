# functions for obtaining quick summaries from datasets
suppressPackageStartupMessages({
  library(data.table)
  library(survey)
  library(cli)
})

summarizer_bin <- function(x, var_name = NULL) {
  if (!is.numeric(x)) stop("variable must be numeric to use summarizer_bin()")
  if (is.numeric(x) & all(x %in% c(0:1, NA_real_))) {
    # dichotomous
    out <- data.table(
      val = as.character(c(0, 1)),
      N   = c(sum(x == 0, na.rm = TRUE), sum(x == 1, na.rm = TRUE)),
      nas = sum(is.na(x))
    )[,
      prop := round(N * 100 / length(na.omit(x)), 1)
    ][, `:=` (
      print      = paste0(format(round(prop, 1), nsmall = 1), " (", trimws(format(N, big.mark = ",")), ")"),
      print_unit = "% (n)"
    )
    ][]
  } else if (is.numeric(x)) {
    # numeric
    out <- data.table(
      mean   = mean(x, na.rm = TRUE),
      sd     = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      N      = length(na.omit(x)),
      val    = "Continuous",
      nas    = sum(is.na(x))
    )[, `:=` (
      print = paste0(
        format(round(mean, 2), big.mark = ",", nsmall = 2),
        " (",
        format(round(sd, 2), big.mark = ",", nsmall = 2),
        ")"),
      print_unit = "mean (sd)"
    )
    ][]
  } else {
    out <- data.table()
  }
  if (!is.null(var_name)) {
    out[, variable := var_name]
  }
  first_cols <- c("variable", "val", "print", "print_unit")
  setcolorder(out, c(first_cols, setdiff(names(out), first_cols)))
  return(out)
}

summarizer_ch <- function(x, var_name = NULL) {
  if (is.factor(x)) {
    out <- data.table(
      val = levels(x),
      N   = map_dbl(levels(x), \(level) sum(x == level, na.rm = TRUE)),
      nas = sum(is.na(x))
    )[,
      prop := round(N * 100 / length(na.omit(x)), 1)
    ][, `:=` (
      print      = paste0(format(round(prop, 1), nsmall = 1), " (", trimws(format(N, big.mark = ",")), ")"),
      print_unit = "% (n)"
    )
    ][]
  } else if (is.character(x)) {
    out <- data.table(
      val = unique(x),
      N   = map_dbl(unique(x), \(level) sum(x == level, na.rm = TRUE)),
      nas = sum(is.na(x))
    )[,
      prop := round(N * 100 / length(na.omit(x)), 1)
    ][, `:=` (
      print      = paste0(format(round(prop, 1), nsmall = 1), " (", trimws(format(N, big.mark = ",")), ")"),
      print_unit = "% (n)"
    )
    ][]
  } else {
    out <- data.table()
  }
  if (!is.null(var_name)) {
    out[, variable := var_name]
  }
  first_cols <- c("variable", "val", "print", "print_unit")
  setcolorder(out, c(first_cols, setdiff(names(out), first_cols)))
  return(out)
}

summarizer_func <- function(x, var_name = NULL) {
  if (is.numeric(x)) {
    out <- summarizer_bin(x, var_name)
  } else if (is.factor(x) | is.character(x)) {
    out <- summarizer_ch(x, var_name)
  } else {
    out <- data.table()
  }
  out
}

summarizer <- function(x, col_names = NULL) {
  if (is.null(col_names)) {
    out <- rbindlist(list(
        data.table(variable = "N", val = "Count", print = format(x[, .N], big.mark = ",")),
        rbindlist(mapply(summarizer_func, x, names(x)), use.names = TRUE, fill = TRUE)
    ), use.names = TRUE, fill = TRUE)
  } else {
    out <- rbindlist(list(
        data.table(variable = "N", val = "Count", print = format(x[, .N], big.mark = ",")),        
        rbindlist(mapply(summarizer_func, x[, ..col_names], col_names), use.names = TRUE, fill = TRUE)
    ), use.names = TRUE, fill = TRUE)
  }
  out
}

# summarize EHR data
phecode_dsb_summarizer <- function(x) {

    first_phe <- x[ x[, .I[which.min(dsb)], id][, V1] ]
    last_phe  <- x[ x[, .I[which.max(dsb)], id][, V1] ]
    ehr_cov   <- merge.data.table(
        first_phe[, .(id, first = dsb)],
        last_phe[, .(id, last = dsb)],
        by = "id"
    )[,
    ehr_days := as.numeric(last - first)
    ][,
    ehr_years := round(ehr_days/365.25, 1)]

    rbindlist(list(
    # number of people
    data.table(variable = "Number of people",
               print    = format(unique(x[, .(id)])[, .N], big.mark = ","),
               val      = "count"),

    # encounters per person
    v_sum(unique(x[, .(id, dsb)])[, .N, id][, N],
          var_name = "Encounters per person, unique"),
          
    # total phecodes per person
    v_sum(summary(x[, .N, id][, N]),
          var_name = "Phecodes per person, total"),

    # unique phecodes per person
    v_sum(unique(x[, .(id, phecode)])[, .N, id][, N],
          var_name = "Phecodes per person, unique")
    ,

    # length follow-up in days
    v_sum(ehr_cov[, ehr_days],
          var_name = "Length of EHR follow-up, days"),
    
    # length of follow-up in years
    v_sum(ehr_cov[, ehr_years],
          var_name = "Length of EHR follow-up, years", other_round = 1)
    ), use.names = TRUE, fill = TRUE)
}

v_sum <- function(x, var_name, mean_round = 1, other_round = 0) {
    data.table(
        variable = var_name,
        val      = c("min", "median", "mean", "max"),
        print    = c(
            format(round(min(x, na.rm = TRUE), other_round), big.mark = ",", nsmall = other_round),
            format(round(median(x, na.rm = TRUE), other_round), big.mark = ",", nsmall = other_round),
            format(round(mean(x, na.rm = TRUE), mean_round), big.mark = ",", nsmall = mean_round),
            format(round(max(x, na.rm = TRUE), other_round), big.mark = ",", nsmall = other_round)
            )
    )
}

# count number of occurrences of phecodes in a phecode indicator matrix
pim_counter <- function(pim, id_var) {
  keep <- setdiff(names(pim), id_var)
  short_pim <- pim[, ..keep]
  counter <- colSums(short_pim, na.rm = TRUE)
  data.table(
    phecode = names(counter),
    N = counter
  )
}

# weighted summary statistics --------------------------------------------------
weighted_summary <- function(var, dsn) {
  tmp <- svymean(as.formula(paste0("~", var)), dsn, na.rm = TRUE)
  out <- as.data.table(tmp)
  if (is.character(model.frame(dsn)[[var]]) | is.factor(model.frame(dsn)[[var]])) {
    names <- gsub(var, "", names(tmp))
    out[, `:=`(
      variable = var,
      val      = names,
      print    = paste0(format(round(mean * 100, 1), nsmall = 1), " (", trimws(format(round(SE * 100, 1), big.mark = ",", nsmall = 1)), ")")
    )]
  } else {
    out[, `:=`(
      variable = var,
      val      = var
    )]
    setnames(out, var, "SE")
    out[, `:=`(
      print = paste0(format(round(mean, 1), nsmall = 1), " (", trimws(format(round(SE, 1), big.mark = ",", nsmall = 1)), ")")
    )]
  }
  out[]
}

weighted_summary_wrapper <- function(data,
                                     weight,
                                     vars,
                                     dsn = NULL) {
  if (is.null(dsn)) {
    pre <- nrow(data)
    data2 <- data[!is.na(get(weight)), ]
    post <- nrow(data2)
    if (pre != post) {
      cli_alert_warning(paste0("Dropped ", pre - post, "(", round((pre - post)*100/pre, 1), "%) rows due to missingness in ", weight))
    }
    dsn <- svydesign(
      id = ~1,
      weights = ~ get(weight),
      data = data2
    )
  }
  res <- lapply(
    vars,
    \(x) weighted_summary(var = x, dsn = dsn)
  )
  rbindlist(res)
}
