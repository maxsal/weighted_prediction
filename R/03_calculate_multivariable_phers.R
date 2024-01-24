# construct multivariable naive phers
# author: max salvatore
# date:   20230816

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)
ms::libri(
  ms, data.table, pROC, glue, logistf, qs, optparse, colorblindr,
  glmnet, cli, parallelly
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "157",
    help = "Outcome phecode [default = %default]"
  ),
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Version of MGI data [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Version of UKB data [default = %default]"
  ),
  make_option("--time_threshold",
    type = "numeric", default = "0",
    help = glue(
      "Time threshold for the phenome data ",
      "[default = %default]"
    )
  ),
  make_option("--discovery_cohort",
    type = "character", default = "mgi",
    help = glue(
      "Use co-occurrence results from discovery cohort in phers ",
      "[default = %default]"
    )
  ),
  make_option("--method",
    type = "character", default = "pwide_sig",
    help = glue(
      "Method for determining phecodes for PheRS ('pwide_sig' or 'tophits') ",
      "[default = %default]"
    )
  ),
  make_option("--tophits_n",
    type = "numeric", default = "50",
    help = glue(
      "Number of top hits to use in top hits PheRS ",
      "[default = %default]"
    )
  ),
  make_option("--weights",
    type = "character", default = "cancer_ipw",
    help = glue(
      "Name of weight suffix for weights used in multivariable regression model (see --cooccur_weights for weights used in cooccurrence results) ",
      "[default = %default]"
    )
  ),
  make_option("--corr_remove",
    type = "numeric", default = NULL,
    help = glue(
      "Correlation threshold for phecode removal (set NULL if no thresholding) ",
      "[default = %default]"
    )
  ),
  make_option("--exclude",
    type = "logical", default = FALSE,
    help = glue(
      "Remove phecodes in exclusion range ",
      "[default = %default]"
    )
  ),
  make_option("--model",
    type = "character", default = "ridge",
    help = glue(
      "Type of logistic regression model - logistic (glm) or ridge ",
      "[default = %default]"
    )
  ),
  make_option("--nfolds",
    type = "numeric", default = 10,
    help = glue(
      "Number of folds for ridge regression cross-validation ",
      "[default = %default]"
    )
  ),
  make_option("--cooccur_weights",
    type = "character", default = NULL,
    help = glue(
      "Weights used in cooccurrence analysis, if desired (see --weights for weights for use in multivariable regression model) ",
      "[default = %default]"
    )
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

if (opt$discovery_cohort == "ukb" & (!is.null(opt$weights) | !is.null(opt$weighted_cooccur))) {
  stop("Script cannot use weights when discovery cohort is UKB. Set --weights=NULL and --weighted_cooccur=NULL or change discovery cohort.")
}
if (!(opt$model %in% c("glm", "ridge"))) {
  stop("Script can only take 'glm' or 'ridge' settings for '--model'")
}

# check output folder exists ---------------------------------------------------
out_path <- glue("results/{coh}/{coh_version}/X{outc}/multivariable/",
  coh = opt$discovery_cohort,
  coh_version = ifelse(opt$discovery_cohort == "mgi",
    opt$mgi_version,
    ifelse(opt$discovery_cohort == "ukb",
      opt$ukb_version, NA
    )
  ),
  outc = gsub("X", "", opt$outcome)
)
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

# 2. specifications (specifies outcome) ----------------------------------------
external_cohort <- ifelse(opt$discovery_cohort == "mgi", "ukb", "mgi")
w               <- ifelse(is.null(opt$weights), "naive", opt$weights) # weights used in multivariable regression
w_co            <- ifelse(is.null(opt$cooccur_weights), "", paste0(opt$cooccur_weights, "_"))

mgi_out_prefix  <- glue("{opt$discovery_cohort}d{w_co}_mgi_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$method}{ifelse(opt$method == 'tophits', opt$tophits_n, '')}_{w}_{opt$model}_")
ukb_out_prefix  <- glue("{opt$discovery_cohort}d{w_co}_ukb_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$method}{ifelse(opt$method == 'tophits', opt$tophits_n, '')}_{w}_{opt$model}_")
comb_out_prefix <- glue("{opt$discovery_cohort}d{w_co}_{external_cohort}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$method}{ifelse(opt$method == 'tophits', opt$tophits_n, '')}_{w}_{opt$model}_")

## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# 3. read data -----------------------------------------------------------------
cli_alert("loading data...")
## mgi
mgi_pim0 <- fread(file_paths[["mgi"]]$pim0_file)
mgi_pim <- read_qs(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/",
  "mgi_X{gsub('X', '', opt$outcome)}_",
  "t{opt$time_threshold}_{opt$mgi_version}.qs"
))
mgi_pim[is.na(mgi_pim)] <- 0
mgi_covariates <- read_qs(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "X{gsub('X', '', opt$outcome)}/",
  "matched_covariates.qs"
))
if (w != "naive") {
  mgi_weights <- read_qs(glue(
    "data/private/mgi/{opt$mgi_version}/",
    "weights_{opt$mgi_version}_comb.qs"
  ))
}

## ukb
ukb_pim0 <- fread(file_paths[["ukb"]]$pim0_file)
ukb_pim <- read_qs(glue(
  "data/private/ukb/{opt$ukb_version}/",
  "X{gsub('X', '', opt$outcome)}/",
  "time_restricted_phenomes/",
  "ukb_X{gsub('X', '', opt$outcome)}_",
  "t{opt$time_threshold}_{opt$ukb_version}.qs"
))
ukb_covariates <- read_qs(glue(
  "data/private/ukb/{opt$ukb_version}/",
  "X{gsub('X', '', opt$outcome)}/",
  "matched_covariates.qs"
))

## cooccur
if (opt$discovery_cohort == "mgi") {
  cooccur <- read_qs(glue(
    "results/mgi/{opt$mgi_version}/",
    "X{gsub('X', '', opt$outcome)}/",
    "mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_{opt$mgi_version}_",
    "{ifelse(is.null(opt$cooccur_weights), '', paste0(opt$cooccur_weights, '_'))}results.qs"
  ))
} else {
  cooccur <- read_qs(glue(
    "results/ukb/{opt$ukb_version}/",
    "X{gsub('X', '', opt$outcome)}/",
    "ukb_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_{opt$ukb_version}",
    "{ifelse(is.null(opt$weights), '', paste0(opt$weights, '_'))}_results.qs"
  ))
}

## phecode info
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character"
)

# 4. phecode exclusions ------------------------------------------------------
if (opt$exclude == TRUE) {
  ## exclusion range from PhewasCatalog
  exclusionRange <- pheinfo[phecode == opt$outcome, phecode_exclude_range]
  exclusions1 <- pheinfo[phecode %in% unlist(unname(sapply(
    strsplit(exclusionRange, ", {0,1}")[[1]],
    expandPhecodes
  ))), phecode]

  ## phecodes not defined in both cohorts
  exclusions2 <- pheinfo[
    !(phecode %in% gsub(
      "X", "",
      intersect(
        names(mgi_pim0),
        names(ukb_pim0)
      )
    )),
    phecode
  ]

  exclusionsX <- glue("X{union(exclusions1, exclusions2)}")

  cooccur <- cooccur[!phecode %in% exclusionsX, ]
}

# 5. calculate naive phers -----------------------------------------------------
cli_alert("calculating naive phers for phecode {gsub('X', '', opt$outcome)}...")
## naive
if (opt$method == "pwide_sig") {
  phes <- cooccur[p_value < 0.05 / .N, length(phecode)]
  if (length(phes) == 0) {
    stop("No phenomewide significant phecodes, stopping...")
  }
}

## select significant hits
if (opt$discovery_cohort == "mgi") {
  pim <- mgi_pim
} else {
  pim <- ukb_pim
}

if (opt$method == "tophits") {
  phers_hits <- cooccur[order(p_value)][1:min(opt$tophits_n, nrow(cooccur))]
}
if (opt$method == "pwide_sig") {
  phers_hits <- cooccur[p_value < 0.05 / nrow(pim)]
}

if (!is.null(opt$corr_remove)) {
  phers_hits <- remove_by_correlation(
    pim         = pim,
    co_res      = phers_hits,
    phecodes    = phers_hits[, phecode],
    corr_thresh = opt$corr_remove
  )
  cli_alert("{nrow(phers_hits)} phecodes remain after correlation thresholding (r2 < {opt$corr_remove})")
}

if (!is.null(opt$weights)) {
  keep_wgts <- c("id", opt$weights)
  wgts <- na.omit(merge.data.table(pim[, .(id)], mgi_weights, by = "id")[, ..keep_wgts])
  pim <- pim[id %in% wgts[, id]]
  wgts <- wgts[[opt$weights]]
} else {
  wgts <- NULL
}

preds <- phers_hits[, phecode]
y <- pim[, case]
xs <- as.matrix(pim[, ..preds])

if (opt$model == "ridge") {
  mod <- glmnet_ridge(.x = xs, .y = y, .nfolds = opt$nfolds, .weights = wgts)
} else {
  mod <- glm(paste0("case ~ ", paste0(preds, collapse = " + ")),
    family = binomial(),
    data = pim, weights = wgts
  )
}

mod_sum <- data.table(
  "phecode" = rownames(coef(mod)),
  "beta"    = coef(mod)[, 1]
)[!phecode %in% c("(Intercept)", "Intercept"), ]

save_qs(
  x = mod,
  file = glue("{out_path}{opt$discovery_cohort}d_{ifelse({opt$discovery_cohort} == 'mgi', 'ukb', 'mgi')}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$model}.qs")
)

mgi_phers <- data.table(
  id   = mgi_pim[, id],
  case = mgi_pim[, case],
  pred = predict(mod, newx = as.matrix(predictor_checker(mgi_pim, preds)), type = "response") |> as.numeric()
)[, phers := scale(pred)][]

ukb_phers <- data.table(
  id   = ukb_pim[, id],
  case = ukb_pim[, case],
  pred = predict(mod, newx = as.matrix(predictor_checker(ukb_pim, preds)), type = "response") |> as.numeric()
)[, phers := scale(pred)][]

cli_alert("generating outputs...")
suppressMessages({
  mgi_roc <- roc(mgi_phers[, case], mgi_phers[, phers], smooth = TRUE, se = FALSE)
  mgi_auc <- ci.auc(mgi_phers[, case], mgi_phers[, phers])

  ukb_roc <- roc(ukb_phers[, case], ukb_phers[, phers], smooth = TRUE, se = FALSE)
  ukb_auc <- ci.auc(ukb_phers[, case], ukb_phers[, phers])
})


mgi_stuff <- data.table(
  sensitivity = mgi_roc$sensitivities,
  specificity = mgi_roc$specificities,
  test_data   = glue("MGI ({ifelse(opt$discovery_cohort == 'mgi', 'discovery', 'external')})")
)
ukb_stuff <- data.table(
  sensitivity = ukb_roc$sensitivities,
  specificity = ukb_roc$specificities,
  test_data   = glue("UKB ({ifelse(opt$discovery_cohort == 'mgi', 'external', 'discovery')})")
)

auc_sum <- data.table(
  test_data = c(glue("MGI ({ifelse(opt$discovery_cohort == 'mgi', 'discovery', 'external')})"), glue("UKB ({ifelse(opt$discovery_cohort == 'mgi', 'external', 'discovery')})")),
  cohort    = c("mgi", "ukb"),
  auc_est   = c(mgi_auc[2], ukb_auc[2]),
  auc_lo    = c(mgi_auc[1], ukb_auc[1]),
  auc_hi    = c(mgi_auc[3], ukb_auc[3]),
  auc_print = c(pretty_print(mgi_auc), pretty_print(ukb_auc))
)

auc_plot <- rbindlist(list(mgi_stuff, ukb_stuff)) |>
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = test_data)) +
  geom_abline(lty = 3) +
  geom_path(linewidth = 2) +
  scale_color_OkabeIto() +
  annotate(
    geom = "label", x = rep(0.75, 2), y = c(0.275, 0.225), label.size = NA,
    label = auc_sum[, auc_print], color = palette_OkabeIto[1:nrow(auc_sum)]
  ) +
  labs(
    title   = glue("AUC for X{gsub('X', '', opt$outcome)} at t{opt$time_threshold}"),
    caption = str_wrap(glue("Discovery cohort: {toupper(opt$discovery_cohort)}; external cohort: {toupper(external_cohort)}; N_phecodes: {nrow(mod_sum)}; method: {opt$method}{ifelse(opt$method == 'tophits', paste0('; N_tophits = ', opt$tophits_n), '')}; exclude codes = {opt$exclude}{ifelse(!is.null(opt$corr_remove), paste0('; correlation threshold = ', opt$corr_remove), '')}"), 80)
  ) +
  coord_equal() +
  cowplot::theme_minimal_grid() +
  theme(
    plot.caption    = element_text(hjust = 0),
    legend.position = "top",
    legend.title    = element_blank()
  )
ggsave(
  plot = auc_plot,
  filename = glue("{out_path}{comb_out_prefix}auc.pdf"),
  width = 7, height = 7, device = cairo_pdf
)

mgi_mod <- glm(case ~ phers, data = mgi_phers, family = "binomial")
ukb_mod <- glm(case ~ phers, data = ukb_phers, family = "binomial")

or_sum <- rbindlist(list(
  cbind(data.table(cohort = "mgi"), extractr_or(mgi_mod)),
  cbind(data.table(cohort = "ukb"), extractr_or(ukb_mod))
))

total_sum <- merge.data.table(
  auc_sum,
  or_sum,
  by = "cohort"
)

total_sum

mgi_phers_dist_plot <- top_or_plotr(
  phers_data = mgi_phers,
  .title = glue(
    "X{gsub('X', '', opt$outcome)}",
    " PheRS distribution by case status at t{opt$time_threshold}"
  ),
  .subtitle = glue("In MGI discovery cohort"),
  .caption = str_wrap(glue("Discovery cohort = {toupper(opt$discovery_cohort)}; N_phecodes: {nrow(mod_sum)}"),
    width = 100
  )
)
ggsave(
  plot = mgi_phers_dist_plot,
  filename = glue("{out_path}{mgi_out_prefix}phers_dist.pdf"),
  width = 8, height = 6, device = cairo_pdf
)


ukb_phers_dist_plot <- top_or_plotr(
  phers_data = ukb_phers,
  .title = glue(
    "X{gsub('X', '', opt$outcome)}",
    " PheRS distribution by case status at t{opt$time_threshold}"
  ),
  .subtitle = glue("In UKB external cohort"),
  .caption = str_wrap(glue("Discovery cohort = {toupper(opt$discovery_cohort)}; N_phecodes: {nrow(mod_sum)}"),
    width = 100
  )
)
ggsave(
  plot = ukb_phers_dist_plot,
  filename = glue("{out_path}{ukb_out_prefix}phers_dist.pdf"),
  width = 8, height = 6, device = cairo_pdf
)

## phers and summary output
save_qs(
  x = mgi_phers,
  file = glue("{out_path}{mgi_out_prefix}phers.qs")
)
save_qs(
  x = ukb_phers,
  file = glue("{out_path}{ukb_out_prefix}phers.qs")
)

list(
  "discovery_cohort" = opt$discovery_cohort,
  "external_cohort"  = external_cohort,
  "total_summary"    = total_sum,
  "model_summary"    = mod_sum,
  "model_call"       = mod$call,
  "n_phecodes"       = nrow(mod_sum),
  "mgi_sense_spec"   = mgi_stuff,
  "ukb_sense_spec"   = mgi_stuff,
  "out_path"         = out_path,
  "optparse_list"    = opt
) |> saveRDS(file = glue("{out_path}{opt$discovery_cohort}d_{ifelse({opt$discovery_cohort} == 'mgi', 'ukb', 'mgi')}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_summary.rds"))

cli_alert_success("script success! 🎉 see output in {.path {out_path}}"))
