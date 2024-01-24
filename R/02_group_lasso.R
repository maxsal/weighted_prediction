# perform group lasso analysis for given threshold and discovery cohort
# author:   max salvatore
# date:     20230418

# libraries --------------------------------------------------------------------
ms::libri(
  ms, qs, data.table, snakecase, stringr, cowplot, optparse,
  glue, colorblindr, rsample, gglasso, prettyunit, cli, parallelly
)

# optparse list ---
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
  make_option("--split_prop",
    type = "numeric", default = "0.7",
    help = glue(
      "Proportion of data in training set ",
      "[default = %default]"
    )
  ),
  make_option("--strata",
    type = "character", default = "case",
    help = glue(
      "Strata (outcome variable) for distributing between train/test sets ",
      "[default = %default]"
    )
  ),
  make_option("--seed",
    type = "numeric", default = "1234",
    help = glue(
      "Set seed for reproducibility ",
      "[default = %default]"
    )
  ),
  make_option("--n_core_prop",
    type = "numeric", default = "0.75",
    help = glue(
      "Proportion of available cores to use ",
      "[default = %default]"
    )
  ),
  make_option("--discovery_cohort",
    type = "character", default = "mgi",
    help = glue(
      "Cohort to use as discovery cohort (mgi / ukb) ",
      "[default = %default]"
    )
  ),
  make_option("--folds",
    type = "numeric", default = 5,
    help = glue(
      "Number of folds in SuperLearner cross validation ",
      "[default = %default]"
    )
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

set.seed(opt$seed)

if (parallelly::availableCores() == 1) {
  cli_alert_info("Only 1 core detected - this could be a while!")
  n_cores <- 1
} else {
  n_cores <- parallelly::availableCores() * opt$n_core_prop
}

options(mc.cores = n_cores)
test_cohort     <- opt$discovery_cohort
external_cohort <- ifelse(opt$discovery_cohort == "mgi", "ukb", "mgi")

for (i in c(
  "expandPhecodes.R", "files-utils.R",
  "super_learner-utils.R", "top_or_plotr.R"
)) {
  source(paste0("fn/", i))
}

# check output folder exists ---------------------------------------------------
out_path <- glue("results/{coh}/{coh_version}/X{outc}/grp_lasso/",
  coh = opt$discovery_cohort,
  coh_version = ifelse(opt$discovery_cohort == "mgi",
    opt$mgi_version,
    ifelse(opt$discovery_cohort == "ukb", opt$ukb_version, NA)
  ),
  outc = gsub("X", "", opt$outcome)
)
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

# data -------------------------------------------------------------------------
cli_alert("reading data...")
d <- read_qs(glue("data/private/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/mgi_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$mgi_version}.qs"))
d_ids <- d[, .(id, case)]


u <- read_qs(glue("data/private/ukb/{opt$ukb_version}/X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/ukb_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$ukb_version}.qs"))
u_ids <- u[, .(id, case)]

p <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character"
)

## exclusion range from PhewasCatalog
exclusionRange <- p[phecode == gsub("X", "", opt$outcome), phecode_exclude_range]
exclusions1 <- p[phecode %in% unlist(unname(sapply(
  strsplit(exclusionRange, ", {0,1}")[[1]],
  expandPhecodes
))), phecode]
exclusions2 <- p[leaf == 0, phecode]

p[, `:=`(
  phecode = paste0("X", phecode),
  group = snakecase::to_sentence_case(group)
)]

keep_these_vars <- intersect(names(d), names(u))
keep_these_vars <- keep_these_vars[!(keep_these_vars %in% paste0("X", union(exclusions1, exclusions2)))]

cli_alert("using {format(length(keep_these_vars) - 1, big.mark = ',')} predictors (only 'leaf' phecodes) available in both cohorts...")
d <- d[, ..keep_these_vars]
u <- u[, ..keep_these_vars]

if (opt$discovery_cohort == "mgi") {
  data <- d
  external <- u
} else if (opt$discovery_cohort == "ukb") {
  data <- u
  external <- d
} else {
  stop("opt$discovery_cohort must be one of 'mgi' or 'ukb'")
}

# below is function or script
cli_alert("splitting data...")

train_obs <- sample(seq_len(nrow(data)), size = round(nrow(data) * opt$split_prop))

data_train <- data[train_obs, ]
data_test <- data[!train_obs, ]

train_ids <- data_train[, .(id, case)]
test_ids <- data_test[, .(id, case)]
external_ids <- external[, .(id, case)]

data_train_covs <- data_train[, !c("id", "case")]
data_train_y <- data_train[, case]

data_test_covs <- data_test[, !c("id", "case")]
data_test_y <- data_test[, case]

external_covs <- external[, !c("id", "case")]
external_y <- external[, case]

# fitting grp lasso ------------------------------------------------------------
## group lasso
v_group <- merge.data.table(
  data.table(phecode = names(data_train_covs)),
  p[, .(phecode, group)],
  by = "phecode"
)[order(group), ][, group_num := .GRP, by = group]
setcolorder(data_train_covs, v_group[, phecode])
data_train_y[data_train_y == 0] <- -1

# gr <- gglasso(x = as.matrix(data_train_covs), y = data_train_y, group = v_group[, group_num],
#               loss = "logit", intercept = FALSE)

# cross validation -------------------------------------------------------------
cli_alert("performing {opt$folds}-fold cross-validation to select lambda...")
gr_cv <- cv.gglasso(
  x         = as.matrix(data_train_covs),
  y         = data_train_y,
  loss      = "logit",
  pred.loss = "L2",
  group     = v_group[, group_num],
  intercept = FALSE, nfolds = opt$folds
)

cli_alert("optimal lambda (within 1sd): {round(gr_cv$lambda.1se, 5)}")

gr_cv_min <- gglasso(
  x         = as.matrix(data_train_covs),
  y         = data_train_y,
  lambda    = gr_cv$lambda.1se,
  group     = v_group[, group_num], loss = "logit",
  intercept = FALSE
)
summary(gr_cv_min)
rownames(gr_cv_min$beta)

data.table(
  phecode = rownames(gr_cv_min$beta),
  beta    = gr_cv_min$beta[, 1]
)[order(beta), ][, or := exp(beta)][]



# predicting -------------------------------------------------------------------
suppressMessages({
  pred_test <- predict(gr_cv_min, data_test_covs, type = "link")[, 1]
  pred_test[pred_test == -1] <- 0
  pred_test_roc <- pROC::roc(data_test_y, pred_test)
  pred_test_auc <- pROC::ci.auc(data_test_y, pred_test)

  pred_ext <- predict(gr_cv_min, external_covs, type = "link")[, 1]
  pred_ext[pred_ext == -1] <- 0
  pred_ext_roc <- pROC::roc(external_y, pred_ext)
  pred_ext_auc <- pROC::ci.auc(external_y, pred_ext)
})

pretty_print <- function(x, r = 3) {
  paste0(format(round(x[2], r), nsmall = r), " (", format(round(x[1], r), nsmall = r), ", ", format(round(x[3], r), nsmall = r), ")")
}

# outputs ----------------------------------------------------------------------

## group lasso object
save_qs(
  x    = gr_cv_min,
  file = glue("{out_path}{opt$discovery_cohort}d_{ifelse({opt$discovery_cohort} == 'mgi', 'ukb', 'mgi')}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_grp_lasso.qs")
)

## dataset containin ID, case status, phers_raw, and phers (mean-standardized)
phers_from_pred <- function(case_data, predicted) {
  out <- cbind(case_data, pred = predicted)
  out[, phers := scale(pred)][]
}
test_phers <- phers_from_pred(case_data = test_ids, pred = pred_test)[, `:=`(cohort = opt$discovery_cohort, split = "test")]
external_phers <- phers_from_pred(case_data = external_ids, pred = pred_ext)[, `:=`(cohort = external_cohort, split = "external")]

save_qs(
  x    = test_phers,
  file = glue("{out_path}{opt$discovery_cohort}d_{opt$discovery_cohort}e_t{opt$time_threshold}_test_phers.qs")
)
save_qs(
  x    = external_phers,
  file = glue("{out_path}{opt$discovery_cohort}d_{external_cohort}e_t{opt$time_threshold}_external_phers.qs")
)

## auc and or summary table
in_stuff <- data.table(
  sensitivity = pred_test_roc$sensitivities,
  specificity = pred_test_roc$specificities,
  test_data   = paste0(toupper(test_cohort), " (test)")
)
out_stuff <- data.table(
  sensitivity = pred_ext_roc$sensitivities,
  specificity = pred_ext_roc$specificities,
  test_data   = paste0(toupper(external_cohort), " (external)")
)

auc_sum <- data.table(
  test_data = c("Hold-out test data", "External data"),
  cohort    = c(opt$discovery_cohort, ifelse(opt$discovery_cohort == "mgi", "ukb", "mgi")),
  auc_est   = c(pred_test_auc[2], pred_ext_auc[2]),
  auc_lo    = c(pred_test_auc[1], pred_ext_auc[1]),
  auc_hi    = c(pred_test_auc[3], pred_ext_auc[3]),
  auc_print = c(pretty_print(pred_test_auc), pretty_print(pred_ext_auc))
)

extractr <- function(x) {
  mod <- glm(case ~ phers, data = x, family = binomial())
  suppressMessages(y <- confint(mod))
  data.table(
    or_est = exp(coef(mod)[["phers"]]),
    or_lo = exp(y["phers", 1]),
    or_hi = exp(y["phers", 2])
  )[, print := paste0(round(or_est, 3), " (", round(or_lo, 3), ", ", round(or_hi, 3), ")")][]
}

or_sum <- rbindlist(list(
  extractr(test_phers)[, cohort := opt$discovery_cohort],
  extractr(external_phers)[, cohort := external_cohort]
))

total_sum <- merge.data.table(
  auc_sum,
  or_sum,
  by = "cohort"
)

total_sum

fwrite(
  x    = total_sum,
  file = glue("{out_path}{opt$discovery_cohort}d_{external_cohort}e_t{opt$time_threshold}_summary.txt")
)

## auc plots --
auc_plot <- rbindlist(list(in_stuff, out_stuff)) |>
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = test_data)) +
  geom_abline(lty = 3) +
  geom_path(linewidth = 2, alpha = 0.2) +
  geom_smooth(
    method = "loess", formula = "y ~ x",
    span = 0.5, se = FALSE, linewidth = 1
  ) +
  scale_color_OkabeIto() +
  annotate(
    geom = "label", x = rep(0.75, 2), y = c(0.275, 0.225), label.size = NA,
    label = auc_sum[, auc_print], color = rev(palette_OkabeIto[1:nrow(auc_sum)])
  ) +
  labs(
    title    = glue("AUC for X{gsub('X', '', opt$outcome)} at t{opt$time_threshold}"),
    subtitle = glue("discovery = {opt$discovery_cohort}, external = {ifelse(opt$discovery_cohort == 'mgi', 'ukb', 'mgi')}"),
    caption  = str_wrap(glue("CV folds = {opt$folds}, train/test prop = {opt$split_prop}"), width = 100)
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
  filename = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_rf_auc.pdf"),
  width = 6, height = 6, device = cairo_pdf
)

# phers distribution by case status ---
test_phers_dist_plot <- top_or_plotr(
  phers_data = test_phers,
  .title = glue("X{gsub('X', '', opt$outcome)} PheRS distribution by case status at t{opt$time_threshold}"),
  .subtitle = glue("In {toupper(test_cohort)} hold out test sample"),
  .caption = str_wrap(glue("Discovery cohort = {opt$discovery_cohort}; CV folds = {opt$folds}, train/test prop = {opt$split_prop}"), width = 100)
)

ggsave(
  plot = test_phers_dist_plot,
  filename = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_test_phers_dist.pdf"),
  width = 6, height = 6, device = cairo_pdf
)

external_phers_dist_plot <- top_or_plotr(
  phers_data = external_phers,
  .title = glue("X{gsub('X', '', opt$outcome)} PheRS distribution by case status at t{opt$time_threshold}"),
  .subtitle = glue("In {toupper(external_cohort)} external sample"),
  .caption = str_wrap(glue("Discovery cohort = {opt$discovery_cohort}; CV folds = {opt$folds}, train/test prop = {opt$split_prop}"), width = 100)
)

ggsave(
  plot = external_phers_dist_plot,
  filename = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_external_phers_dist.pdf"),
  width = 6, height = 6, device = cairo_pdf
)

## summary list object
list(
  "discovery_cohort"     = opt$discovery_cohort,
  "external_cohort"      = external_cohort,
  "total_summary"        = total_sum,
  "group_lasso_summary"  = data.table(
    phecode = rownames(gr_cv_min$beta),
    beta    = gr_cv_min$beta[, 1]
  )[order(beta), ][, or := exp(beta)][],
  "group_lasso_lambda"   = gr_cv_min$lambda,
  "discovery_sense_spec" = in_stuff,
  "external_sense_spec"  = out_stuff,
  "out_path"             = out_path,
  "optparse_list"        = opt
) |> saveRDS(file = glue("{out_path}{opt$discovery_cohort}d_{ifelse({opt$discovery_cohort} == 'mgi', 'ukb', 'mgi')}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_summary.rds"))

message("script success! see output in {.path {out_path}}")
