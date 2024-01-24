# functions to assist in evaluation of cooccurence and phers results
# author: max salvatore
# date:   20221212

suppressPackageStartupMessages({
  library(data.table)
  library(progress)
  library(glue)
  library(logistf)
  library(ResourceSelection)
  library(DescTools)
  library(pROC)
  library(pwr)
})

# log10toP: converts log10 p-values to decimal format --------------------------
log10toP <- function(log10P) {
  log10P <- abs(as.numeric(log10P))
  if (is.na(log10P)) {
    return(NA)
  }
  if (log10P > 300) {
    part1 <- log10P %/% 100 * 100
    part2 <- log10P - part1
    if (part2 != 0) {
      P <- format(signif(10^-part2, 3), scientific = T)
      P <- paste(as.numeric(gsub("e-.+", "", P)), "e-",
        as.numeric(gsub(".+-", "", P), sep = "") + part1,
        sep = ""
      )
    } else {
      P <- paste("1e-", part1, sep = "")
    }
  } else {
    P <- signif(10^-log10P, 3)
  }
  as.character(P)
}

# getTopEffects (from Lars) ----------------------------------------------------
getTopEffects <- function(prob, pred, cov, dat) {
  add_cov <- ifelse(length(cov) == 0, "",
    paste0(" + ", paste0(cov, collapse = " + "))
  )
  riskBin <- glue("Top_{prob}")
  dat[[riskBin]] <- ifelse(
    dat[[pred]] >= quantile(dat[case == 0, ][[pred]], probs = 1 - prob),
    1,
    0
  )
  tmp_vars <- c("case", riskBin, cov)
  fitAsso2 <- logistf(
    as.formula(paste0("case ~ ", riskBin, add_cov)),
    data = na.omit(dat[, ..tmp_vars]),
    control = logistf.control(maxit = 1000)
  )
  getValues(fitAsso2, riskBin)
}

# getTopPower (from Lars) ------------------------------------------------------
getTopPower <- function(prob, dat, pred) {
  riskBin <- glue("Top_{prob}")
  dat[[riskBin]] <- ifelse(
    dat[[pred]] >= quantile(dat[case == 0, ][[pred]], probs = 1 - prob),
    1,
    0
  )
  ptable <- table("Top" = dat[[riskBin]], "trait" = dat[["case"]])
  if (nrow(ptable) != 2) {
    return(c("power" = NA, "h" = NA))
  }
  ntop <- sum(ptable[, 2])
  nrest <- sum(ptable[, 1])
  ptop <- ptable[2, 2] / ntop
  prest <- ptable[2, 1] / nrest
  h <- ES.h(ptop, prest)
  unlist(pwr.2p2n.test(
    h = h,
    n1 = ntop,
    n2 = nrest,
    sig.level = 0.05,
    power = NULL
  )[c("power", "h")])
}

# getValues (from Lars) --------------------------------------------------------
getValues <- function(fitAsso, predictor) {
  sebeta <- sqrt(diag(vcov(fitAsso)))
  names(sebeta) <- fitAsso$terms
  BETA <- round(fitAsso$coefficient[predictor], 4)
  SEBETA <- round(sebeta[predictor], 4)

  vars <- diag(fitAsso$var)
  names(vars) <- fitAsso$terms
  LOG10P <- -pchisq((fitAsso$coefficient[predictor]^2 / vars[predictor]),
    1,
    lower.tail = F, log = T
  ) / log(10)
  P <- log10toP(LOG10P)

  OR <- round(exp(fitAsso$coefficient[predictor]), 4)
  CI1 <- round(exp(fitAsso$ci.lower[predictor]), 4)
  CI2 <- round(exp(fitAsso$ci.upper[predictor]), 4)
  or_print <- glue("{format(round(OR, 2), nsmall = 2)} ({format(round(CI1, 2), nsmall = 2)}, {format(round(CI2, 2), nsmall = 2)})")
  out <- c(BETA, SEBETA, P, OR, CI1, CI2, or_print, signif(LOG10P, 4))
  names(out) <- paste0(
    predictor, "_",
    c("beta", "sebeta", "p_val", "or", "or_lo", "or_hi", "or_print", "log10p")
  )
  return(out)
}

# calculate_phers --------------------------------------------------------------
calculate_phers <- function(pim, # phecode indicator matrix
                            res, # results matrix - phecode, beta, p-value
                            method, # tophits or pwide_sig
                            tophits_n = 50, # n of hits to select based on p-value
                            bonf_tests = NULL, # denominator of bonferroni correction
                            reverse_code = FALSE, # reverse code negatives to keep phers above 0?
                            corr_remove = 0.25) {
  # check that columns exist in res data
  res_data_cols <- c("phecode", "beta", "p_value")
  if (!any(res_data_cols %in% names(res))) {
    stop(glue(
      "Missing columns in 'res': ",
      "{glue_collapse(res_data_cols[!(res_data_cols %in% names(res))]",
      ", sep = ', ')}"
    ))
  }

  ## select significant hits
  if (method == "tophits") {
    phers_hits <- res[order(p_value)][1:min(tophits_n, nrow(res))]
  }
  if (method == "pwide_sig") {
    if (is.null(bonf_tests)) {
      phers_hits <- res[p_value < 0.05 / nrow(pim)]
    } else {
      phers_hits <- res[p_value < 0.05 / bonf_tests]
    }
  }

  if (!is.null(corr_remove)) {
    phers_hits <- remove_by_correlation(
      pim         = pim,
      co_res      = phers_hits,
      phecodes    = phers_hits[, phecode],
      corr_thresh = corr_remove
    )
    message(glue("{nrow(phers_hits)} phecodes remain after correlation thresholding (r2 < {corr_remove})"))
  }

  out <- copy(pim)

  out[, pred := 0]
  message("calculating phers...")
  pb <- txtProgressBar(max = length(phers_hits[, phecode]), width = 50, style = 3)
  if (reverse_code == FALSE) {
    for (i in phers_hits[, phecode]) {
      if (!(i %in% names(pim))) {
        cli_alert_warning("{i} not in phecode indicator matrix, skipping")
        next
      }
      out[, pred := pred + (phers_hits[phecode == i, beta] * get(i))]
      setTxtProgressBar(pb, i)
    }
  } else {
    for (i in phers_hits[, phecode]) {
      if (!(i %in% names(pim))) {
        cli_alert_warning("{i} not in phecode indicator matrix, skipping")
        next
      }
      out[, pred := pred + (
        phers_hits[phecode == i, beta] * get(i) *
          as.numeric(phers_hits[phecode == i, beta] > 0)) - (
        phers_hits[phecode == i, beta] * (1 - get(i)) *
          as.numeric(phers_hits[phecode == i, beta] < 0)
      )]
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)

  out[, phers := scale(pred)]

  list(
    n_phecodes = length(phers_hits[, phecode]),
    method = method,
    phecodes = phers_hits,
    ata = out[, .(id, case, pred, phers)]
  )
}

# quick AUCs -------------------------------------------------------------------
quick_naive_aucs <- function(x) {
  out <- data.table()
  message("AUC")
  pb <- txtProgressBar(max = length(names(x)[grepl("phers", names(x))]))
  for (i in names(x)[grepl("phers", names(x))]) {
    out <- rbindlist(list(
      out,
      data.table(
        phers = i,
        auc = suppressMessages(pROC::roc(
          response = x[["case"]],
          predictor = x[[i]],
          family = binomial(),
          ci = TRUE
        )$auc)
      )
    ), fill = TRUE)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  out[order(-auc)]
}

# remove phecodes by correlation -----------------------------------------------
remove_by_correlation <- function(
    pim,
    cooccurrence_results,
    exposures_to_consider = NULL,
    exposure_var          = "phecode",
    p_value_var           = "p_value",
    top_n                 = 50,
    corr_thresh           = 0.5,
    weights               = NULL
) {
  if (is.null(exposures_to_consider)) {
    sub_co <- cooccurrence_results[order(p_value), ][1:top_n, ]
  } else {
    sub_co <- cooccurrence_results[get(exposure_var) %in% exposures_to_consider, ]
  }
  if (dim(sub_co)[1] <= 1) {
    cli_alert_info("No exposures to remove at correlation >= {corr_thresh}")
    return(sub_co)
  } else {

    cli_alert_info("{nrow(sub_co)} exposures before correlation thresholding")

    exposure_vars <- sub_co[, get(exposure_var)]
    if (!is.null(weights)) {
      corr_results <- ms::fcor(pim[, ..exposure_vars], w = weights)
    } else {
      corr_results <- ms::fcor(pim[, ..exposure_vars])
    }
    setnames(corr_results, old = c("V1", "V2", "N"), new = c("exposure1", "exposure2", "correlation"))
    corr_results <- merge.data.table(
      corr_results,
      cooccurrence_results[, .(exposure1 = get(exposure_var), exposure1_p = get(p_value_var))],
      by = "exposure1"
    )
    corr_results <- merge.data.table(
      corr_results,
      cooccurrence_results[, .(exposure2 = get(exposure_var), exposure2_p = get(p_value_var))],
      by = "exposure2"
    )
    corr_results[, exclude := NA]

      sub_corr <- corr_results[order(-correlation), ][correlation >= corr_thresh]

      if (dim(sub_corr)[1] != 0) {
        for (k in 1:dim(sub_corr)[1]) {
          if (sub_corr[k, exposure1] %in% sub_corr[, exclude] | sub_corr[k, exposure2] %in% sub_corr[, exclude]) {
            sub_corr$exclude[k] <- "Previously Removed"
          } else {
            sub_corr$exclude[k] <- ifelse(sub_corr[k, exposure1_p] < sub_corr[k, exposure2_p],
              sub_corr[k, exposure2], sub_corr[k, exposure2]
            )
          }
        }
        exclude_corr <- sub_corr$exclude[!sub_corr$exclude %in% "Previously Removed"]

      cli_alert_info("{length(exclude_corr)} exposures removed at correlation >= {corr_thresh}")

      sub_co[!(get(exposure_var) %in% exclude_corr), ]
    } else {
      cli_alert_info("No exposures to remove at correlation >= {corr_thresh}")
    }
  }
}

# performs simple ridge selecting minimum lambda -------------------------------
glmnet_ridge <- function(.x, .y, .alpha = 0, .family = "binomial",
                         .nfolds = 10, .weights = NULL) {
  message("identifying optimal lambda...")
  cv_ridge <- cv.glmnet(
    x       = .x,
    y       = .y,
    alpha   = .alpha,
    family  = .family,
    nfolds  = .nfolds,
    weights = .weights
  )
  lambda_opt <- cv_ridge$lambda.min
  message(glue("optimal lambda: {format(lambda_opt, scientific = TRUE)}. fitting final model..."))
  ridge <- glmnet(
    x       = .x,
    y       = .y,
    alpha   = .alpha,
    family  = .family,
    lambda  = lambda_opt,
    weights = .weights
  )
  return(ridge)
}

# check for missing predictors and add 0 variables if missing ------------------
predictor_checker <- function(data, predictors) {
  out <- data.table::copy(data)

  # see which predictors are not in data
  missing_predictors <- predictors[!(predictors %in% names(out))]
  if (length(missing_predictors) > 0) {
    message(glue("Found {length(missing_predictors)} missing predictor{?s}: {stringr::str_trunc(paste0(missing_predictors, collapse = ', '), width = 300)}"))
    for (i in seq_along(missing_predictors)) {
      out[[missing_predictors[i]]] <- 0
    }
  } else {
    message("No missing predictors found")
  }

  return(out[, ..predictors])
}

# get phenomewide hits ---------------------------------------------------------
get_pwide_hits <- function(
    data,
    pval_var = "p_value",
    exposure_var = "exposure",
    exposure_only = TRUE) {
  out <- data[get(pval_var) < 0.05 / .N, ]
  if (exposure_only) {
    return(out[[exposure_var]])
  } else {
    return(out)
  }
}
# get top N hits by p-value ----------------------------------------------------
get_top_hits <- function(
    data,
    top_n = NULL,
    pval_var = "p_value",
    exposure_var = "exposure",
    exposure_only = TRUE) {
  if (is.null(top_n)) {
    cli::cli_alert_info("{.field top_n} argument not specified, defaulting to {.val 50}")
    top_n <- 50
  }
  out <- data[order(get(pval_var)), ][1:top_n, ]
  if (exposure_only) {
    return(out[[exposure_var]])
  } else {
    return(out)
  }
}

# extract number after t in string ---------------------------------------------
get_time_threshold <- function(string, numeric = FALSE) {
  match <- stringr::str_extract(string, "t([0-9]+(\\.[0-9]+)?)")
  if (!is.na(match)) {
    if (numeric) {
      as.numeric(sub("t", "", match))
    } else {
      sub("t", "", match)
    }
  } else {
    NA
  }
}

# helper functions (see 03_calculated_multivariable_phers.R) -------------------
## round to specific number of decimal places with big mark
pretty_round <- function(x, r) {
  format(round(x, r), big.mark = ",", nsmall = r)
}
## print an estimate (CI) from vector (lower, estimate, upper)
pretty_print <- function(x, r = 3) {
  paste0(
    pretty_round(x[2], r), " (", pretty_round(x[1], r), ", ",
    pretty_round(x[3], r), ")"
  )
}
## extract OR from a model
extractr_or <- function(x, r = 2) {
  suppressMessages(y <- confint(x))
  data.table(
    or_est = exp(coef(x)[["phers"]]),
    or_lo = exp(y["phers", 1]),
    or_hi = exp(y["phers", 2])
  )[, or_print := paste0(
    format(round(or_est, r), big.mark = ",", nsmall = r), " (",
    format(round(or_lo, r), big.mark = ",", nsmall = r), ", ",
    format(round(or_hi, r), big.mark = ",", nsmall = r), ")"
  )][]
}
