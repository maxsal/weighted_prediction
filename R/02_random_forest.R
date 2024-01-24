# perform random forest analysis for given threshold and discovery cohort
# requires: time-threshold phecode indicator matrices must already exist
# author:   max salvatore
# date:     20230816

# libraries --------------------------------------------------------------------
ms::libri(
  ms, qs, vip, data.table, snakecase, stringr, cowplot, optparse, glue,
  clauswilke/colorblindr, ranger, caret, ggnewscale, ggpubr, scales, pROC,
  logistf, parallelly, cli
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
  make_option("--n_trees",
    type = "numeric", default = "500",
    help = glue(
      "Number of trees to grow in random forest ",
      "[default = %default]"
    )
  ),
  make_option("--n_vip",
    type = "numeric", default = "20",
    help = glue(
      "Number of variables to include in variable importance plot ",
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
    type = "numeric", default = "0.5",
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
  make_option("--n_mtry",
    type = "numeric", default = 6,
    help = glue(
      "Number of values of mtry to consider ",
      "[default = %default]"
    )
  ),
  make_option("--n_node_size",
    type = "numeric", default = 5,
    help = glue(
      "Number of values of node_size to consider ",
      "[default = %default]"
    )
  ),
  make_option("--n_samp_frac",
    type = "numeric", default = 4,
    help = glue(
      "Number of values of sample to consider ",
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
  message("Only 1 core detected - this could be a while!")
  n_cores <- 1
} else {
  n_cores <- parallelly::availableCores() * opt$n_core_prop
}

test_cohort     <- opt$discovery_cohort
external_cohort <- ifelse(test_cohort == "mgi", "ukb", "mgi")

source("fn/expandPhecodes.R")
source("fn/files-utils.R")
source("fn/top_or_plotr.R")

# check output folder exists ---------------------------------------------------
out_path <- glue("results/{coh}/{coh_version}/X{outc}/random_forest/",
  coh         = opt$discovery_cohort,
  coh_version = ifelse(opt$discovery_cohort == "mgi", opt$mgi_version, ifelse(opt$discovery_cohort == "ukb", opt$ukb_version, NA)),
  outc        = gsub("X", "", opt$outcome)
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

data_train <- data_train[, !c("id")]
data_test <- data_test[, !c("id")]
external <- external[, !c("id")]

data_train_covs <- data_train[, !c("case")]
data_train_y <- data_train[, case]

data_test_covs <- data_test[, !c("case")]
data_test_y <- data_test[, case]

external_covs <- external[, !c("case")]
external_y <- external[, case]

# full grid search
hyper_grid <- expand.grid(
  mtry        = round(seq(5, round(ncol(data_train) / 2), length.out = opt$n_mtry)),
  node_size   = round(seq(3, round(nrow(data_train) / 100), length.out = opt$n_node_size)),
  sample_frac = seq(0.5, 0.8, length.out = opt$n_samp_frac),
  oob_rmse    = 0
) |> as.data.table()

cli_alert("Trying {opt$n_mtry} values of mtry, {opt$n_node_size} of node size, {opt$n_samp_frac} values of sample fraction ({nrow(hyper_grid)} combinations)")

cli_progress_bar(name = "grid search...", total = nrow(hyper_grid))
for (i in 1:nrow(hyper_grid)) {
  model <- ranger(
    formula         = case ~ .,
    data            = data_train,
    num.trees       = opt$n_trees,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_frac[i],
    num.threads     = n_cores
  )
  hyper_grid$oob_rmse[i] <- sqrt(model$prediction.error)
  cli_progress_update()
}
cli_progress_done()

hyper_grid[order(oob_rmse), ][1, ]

oob_rmse <- vector(mode = "numeric", length = 100)

cli_progress_bar(
  name = "fitting random forest with optimal hyperparameters",
  total = length(oob_rmse)
)
for (i in seq_along(oob_rmse)) {
  optimal_ranger <- ranger(
    formula         = case ~ .,
    data            = data_train,
    num.trees       = opt$n_trees,
    mtry            = hyper_grid[order(oob_rmse), ][1, mtry],
    min.node.size   = hyper_grid[order(oob_rmse), ][1, node_size],
    sample.fraction = hyper_grid[order(oob_rmse), ][1, sample_frac],
    importance      = "impurity",
    num.threads     = n_cores
  )
  oob_rmse[i] <- sqrt(optimal_ranger$prediction.error)
  cli_progress_update()
}
cli_progress_done()

# predicting -------------------------------------------------------------------
phers_from_pred <- function(case_data, predicted) {
  out <- cbind(case_data, pred = predicted)
  out[, phers := scale(pred)][]
}

suppressMessages({
  pred_opt     <- predict(optimal_ranger, data_test, type = "response")
  test_phers   <- phers_from_pred(case_data = test_ids, predicted = pred_opt[["predictions"]])
  pred_opt_roc <- roc(test_phers[, case], test_phers[, phers])
  pred_opt_auc <- ci.auc(test_phers[, case], test_phers[, phers])

  pred_oth       <- predict(optimal_ranger, external, type = "response")
  external_phers <- phers_from_pred(case_data = external_ids, predicted = pred_oth[["predictions"]])
  pred_oth_roc   <- roc(external_phers[, case], external_phers[, phers])
  pred_oth_auc   <- ci.auc(external_phers[, case], external_phers[, phers])
})

save_qs(
  x    = test_phers,
  file = glue("{out_path}{test_cohort}d_{test_cohort}e_t{opt$time_threshold}_test_phers.qs")
)
save_qs(
  x    = external_phers,
  file = glue("{out_path}{test_cohort}d_{external_cohort}e_t{opt$time_threshold}_external_phers.qs")
)

pretty_print <- function(x, r = 3) {
  paste0(round(x[2], r), " (", round(x[1], r), ", ", round(x[3], r), ")")
}

in_stuff <- data.table(
  sensitivity = pred_opt_roc$sensitivities,
  specificity = pred_opt_roc$specificities,
  test_data   = paste0(toupper(test_cohort), " (test)")
)
out_stuff <- data.table(
  sensitivity = pred_oth_roc$sensitivities,
  specificity = pred_oth_roc$specificities,
  test_data   = paste0(toupper(external_cohort), " (external)")
)

auc_sum <- data.table(
  test_data = c("Hold-out test data", "External data"),
  cohort    = c(test_cohort, external_cohort),
  auc_est   = c(pred_opt_auc[2], pred_oth_auc[2]),
  auc_lo    = c(pred_opt_auc[1], pred_oth_auc[1]),
  auc_hi    = c(pred_opt_auc[3], pred_oth_auc[3]),
  auc_print = c(pretty_print(pred_opt_auc), pretty_print(pred_oth_auc))
)

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
    caption  = str_wrap(glue("N_trees = {opt$n_trees}, mtry = {hyper_grid[order(oob_rmse), ][1, mtry]}, node size = {hyper_grid[order(oob_rmse), ][1, node_size]}, sample fraction = {hyper_grid[order(oob_rmse), ][1, sample_frac]}"), width = 100)
  ) +
  coord_equal() +
  theme_minimal_grid() +
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

## VIP plot
vip_data <- data.table(
  variable = names(optimal_ranger$variable.importance),
  importance = optimal_ranger$variable.importance
)[order(-importance), ][1:opt$n_vip, ] |>
  merge.data.table(p[, .(phecode, description, group, color)], by.x = "variable", by.y = "phecode", all.x = TRUE)
vip_data[, full_name := ifelse(is.na(description), variable, paste0(variable, ": ", description))]

cols <- unique(vip_data[, .(group, color)])
colors <- cols[, color]
names(colors) <- cols[, group]

vip_plot <- vip_data |>
  ggplot(aes(x = reorder(full_name, importance), y = importance, fill = group)) +
  geom_bar(stat = "identity") +
  labs(
    title = glue("VIP for X{gsub('X', '', opt$outcome)} at t{opt$time_threshold}"),
    subtitle = glue("discovery = {opt$discovery_cohort}, external = {ifelse(opt$discovery_cohort == 'mgi', 'ukb', 'mgi')}"),
    caption = str_wrap(glue("N_trees = {opt$n_trees}, mtry = {hyper_grid[order(oob_rmse), ][1, mtry]}, node size = {hyper_grid[order(oob_rmse), ][1, node_size]}, sample fraction = {hyper_grid[order(oob_rmse), ][1, sample_frac]}"), width = 100),
    x = "",
    y = "Importance"
  ) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = \(x) stringr::str_wrap(x, width = 50)) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE)) +
  coord_flip() +
  cowplot::theme_minimal_grid() +
  theme(
    plot.caption    = element_text(hjust = 0),
    legend.position = "bottom",
    legend.title    = element_blank()
  )
ggsave(
  plot = vip_plot,
  filename = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_vip.pdf"),
  width = 10, height = 8, device = cairo_pdf
)

saveRDS(optimal_ranger, file = glue("{out_path}{opt$discovery_cohort}d_{ifelse({opt$discovery_cohort} == 'mgi', 'ukb', 'mgi')}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_optimal_rf.rds"))


if (opt$discovery_cohort == "mgi") {
  mgi <- test_phers[, type := "test set"]
} else {
  mgi <- external_phers[, type := "external"]
}

if (opt$discovery_cohort == "ukb") {
  ukb <- test_phers[, type := "test set"]
} else {
  ukb <- external_phers[, type := "external"]
}

mgi_mod <- glm(case ~ phers, data = mgi, family = binomial())
ukb_mod <- glm(case ~ phers, data = ukb, family = binomial())
extractr <- function(x) {
  suppressMessages(y <- confint(x))
  data.table(
    or_est = exp(coef(x)[["phers"]]),
    or_lo = exp(y["phers", 1]),
    or_hi = exp(y["phers", 2])
  )[, print := paste0(round(or_est, 3), " (", round(or_lo, 3), ", ", round(or_hi, 3), ")")][]
}
or_sum <- rbindlist(list(
  cbind(data.table(cohort = "mgi"), extractr(mgi_mod)),
  cbind(data.table(cohort = "ukb"), extractr(ukb_mod))
))

total_sum <- merge.data.table(
  auc_sum,
  or_sum,
  by = "cohort"
)

total_sum

# phers distribution by case status ---
test_phers_dist_plot <- top_or_plotr(
  phers_data = test_phers,
  .title = glue("X{gsub('X', '', opt$outcome)} PheRS distribution by case status at t{opt$time_threshold}"),
  .subtitle = glue("In {toupper(test_cohort)} hold out test sample"),
  .caption = str_wrap(glue("Discovery cohort = {test_cohort}; CV folds = {opt$folds}, train/test prop = {opt$split_prop}"), width = 100)
)

ggsave(
  plot = test_phers_dist_plot,
  filename = glue("{out_path}{test_cohort}d__{external_cohort}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_test_phers_dist.pdf"),
  width = 8, height = 6, device = cairo_pdf
)

external_phers_dist_plot <- top_or_plotr(
  phers_data = external_phers,
  .title = glue("X{gsub('X', '', opt$outcome)} PheRS distribution by case status at t{opt$time_threshold}"),
  .subtitle = glue("In {toupper(external_cohort)} external sample"),
  .caption = str_wrap(glue("Discovery cohort = {test_cohort}; CV folds = {opt$folds}, train/test prop = {opt$split_prop}"), width = 100)
)

ggsave(
  plot = external_phers_dist_plot,
  filename = glue("{out_path}{test_cohort}d_{external_cohort}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_external_phers_dist.pdf"),
  width = 8, height = 6, device = cairo_pdf
)

# optimal ranger object
saveRDS(optimal_ranger, file = glue("{out_path}{test_cohort}d_{external_cohort}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_ranger.rds"))

list(
  "discovery_cohort" = opt$discovery_cohort,
  "external_cohort" = ifelse(opt$discovery_cohort == "mgi", "ukb", "mgi"),
  "total_summary" = total_sum,
  "discovery_sense_spec" = in_stuff,
  "external_sense_spec" = out_stuff,
  "vip_table" = data.table(
    variable = names(optimal_ranger$variable.importance),
    importance = optimal_ranger$variable.importance
  )[order(-importance), ],
  "hyper_grid" = hyper_grid,
  "auc_plot" = auc_plot,
  "vip_plot" = vip_plot,
  "out_path" = out_path,
  "optparse_list" = opt
) |> saveRDS(file = glue("{out_path}{opt$discovery_cohort}d_{ifelse({opt$discovery_cohort} == 'mgi', 'ukb', 'mgi')}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_summary.rds"))

cli_alert_success("script success! 🎉 see output in {.path {out_path}}")
