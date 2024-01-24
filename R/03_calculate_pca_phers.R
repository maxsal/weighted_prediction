# phers construction code
# an adaptaion from Lauren Beesley's D1_Construct PheRS PanCan.R script
# author: max salvatore
# date:   20230314

# libraries, functions, and options --------------------------------------------
options(stringsAsFactors = FALSE)
ms::libri(
  ms, data.table, ppcor, ggplot2, gridExtra, MatchIt, glue, SPAtest, ggrepel,
  qs, colorblindr, optparse, stringr, pROC, cowplot, gridExtra, parallelly,
  cli
)

set.seed(61787)

for (i in c(
  "files-utils.R",
  "expandPheocdes.R",
  "eval-utils.R",
  "top_or_plotr.R"
)) {
  source(paste0("fn/", i))
}

mytheme <- ttheme_default(
  core    = list(fg_params = list(cex = 0.5)),
  colhead = list(fg_params = list(cex = 0.5)),
  rowhead = list(fg_params = list(cex = 0.5))
)

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
  make_option("--pc_var_explain",
    type = "numeric", default = 0.95,
    help = glue(
      "Cumulative variation explained by PCs threshold ",
      "[default = %default]"
    )
  ),
  make_option("--discovery_cohort",
    type = "character", default = "mgi",
    help = glue(
      "Which cohort to conduct PCA in ",
      "[default = %default]"
    )
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# specifications ---------------------------------------------------------------
mgi_data_path <- glue(
  "data/private/mgi/{opt$mgi_version}/",
  "X{gsub('X', '', opt$outcome)}"
)
ukb_data_path <- glue(
  "data/private/ukb/{opt$ukb_version}/",
  "X{gsub('X', '', opt$outcome)}"
)
out_path <- glue(
  "results/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
  "pca"
)
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

## pull file paths corresponding to the data version specified
file_paths <- get_files(mgi_version = opt$mgi_version, ukb_version = opt$ukb_version)

# read in data -----------------------------------------------------------------
cli_alert("reading data...")
## mgi covariate data
mgi_cov <- read_qs(glue("{mgi_data_path}/matched_covariates.qs"))[
  get(glue("t{opt$time_threshold}_indicator")) == 1, .(
    id, female,
    age_at_threshold = round(get(glue("t{opt$time_threshold}_threshold")) /
      365.25, 1),
    length_followup, case
  )
]
mgi_cov <- mgi_cov[complete.cases(mgi_cov), ]

# mgi time-restricted pim
mgi_pim <- read_qs(glue(
  "{mgi_data_path}/time_restricted_phenomes/",
  "mgi_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}",
  "_{opt$mgi_version}.qs"
))
mgi_pim <- merge(mgi_pim, data.table(id = mgi_cov$id),
  by = "id", all.x = FALSE
)
short_mgi_pim <- mgi_pim[, !c("id")]
short_mgi_pim[is.na(short_mgi_pim)] <- 0

# ukb time_restricted pim
ukb_cov <- read_qs(glue("{ukb_data_path}/matched_covariates.qst"))[
  get(glue("t{opt$time_threshold}_indicator")) == 1, .(
    id, female,
    age_at_threshold = round(get(glue("t{opt$time_threshold}_threshold")) /
      365.25, 1),
    length_followup, case
  )
]
ukb_cov <- ukb_cov[complete.cases(ukb_cov), ]

ukb_pim <- read_qs(glue(
  "{ukb_data_path}/time_restricted_phenomes/",
  "ukb_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}",
  "_{opt$ukb_version}.qs"
))
ukb_pim <- merge(ukb_pim, data.table(id = ukb_cov$id),
  by = "id", all.x = FALSE
)
ukb_pim[is.na(ukb_pim)] <- 0
ukb_cases <- ukb_pim[, .(id, case)]
short_ukb_pim <- ukb_pim[, !c("id", "case")]

## mgi master PIM - read 1 row for variable names only
mgi_pim0 <- fread(file_paths[["mgi"]]$pim0_file, nrow = 1)

## ukb master PIM - read 1 row for variable names only
ukb_pim0 <- fread(file_paths[["ukb"]]$pim0_file, nrow = 1)

## phecode info
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character"
)

# exclusions -------------------------------------------------------------------
outcome_vector <- short_mgi_pim[["case"]]

## exclusion range phecodes (from PhewasCatalog)
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

exclusionsX <- c("case", glue(
  "X{c(gsub('X', '', opt$outcome), ",
  "union(exclusions1, exclusions2))}"
))

# subset phenotype data --------------------------------------------------------
included <- names(short_mgi_pim)[!(names(short_mgi_pim) %in% exclusionsX)]
exclude_mgi_pim <- short_mgi_pim[, ..included]

# obtain principal components --------------------------------------------------
cli_alert("calculating principal components...")
mgi_pca <- prcomp(exclude_mgi_pim, center = FALSE, scale. = FALSE)

pcs <- mgi_pca$x[
  , 1:which.min(abs(summary(mgi_pca)$importance[3, ] - opt$pc_var_explain))
]

maxs <- apply(abs(pcs), 2, max)
pcs_mod <- sweep(pcs, MARGIN = 2, maxs, "/")
pcs_mod <- pcs_mod + 1

rotations <- data.table(mgi_pca$rotation)
rotations <- rotations[
  , 1:which.min(abs(summary(mgi_pca)$importance[3, ] - opt$pc_var_explain))
]

# run SPAtest ------------------------------------------------------------------
cli_alert("running SPAtest")
## method 2: pca
test2 <- ScoreTest_SPA(
  genos = t(pcs_mod),
  pheno = outcome_vector,
  cov = mgi_cov[id %in% mgi_pim[, id], !c("id")],
  method = "fastSPA",
  beta.out = TRUE, beta.Cutoff = 1
)

# organize results -------------------------------------------------------------
cli_alert("organizing results and generating plots...")
results_pc <- data.frame(
  phecode      = glue("PC{c(1:length(pcs_mod[1, ]))}"),
  pvals_m2     = test2$p.value,
  alphas_m2    = test2$beta,
  alphasmod_m2 = sweep(as.matrix(test2$beta), MARGIN = 1, maxs, "/")
) |> as.data.table()

results_pc[, alphasmod_m2sig := fifelse(
  pvals_m2 > 0.05 / length(results_pc[, 1]) | is.na(pvals_m2),
  rep(0, length(results_pc[, 1])),
  alphasmod_m2
)]

results <- data.table(
  phecode = names(exclude_mgi_pim)
)

results$betas_m2 <- as.matrix(rotations) %*%
  as.matrix(results_pc[, alphasmod_m2])
results$betas_m2sig <- as.matrix(rotations) %*%
  as.matrix(results_pc[, alphasmod_m2sig])

### save rotations, results_pc, and results
save_qs(
  x = rotations,
  file = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_rotations.qs"
  )
)
save_qs(
  x = results_pc,
  file = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pcs.qs"
  )
)
save_qs(
  x = results,
  path = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_results.qs"
  )
)

# manhattan plot ---------------------------------------------------------------
results_short <- data.table(
  log10pvals = -log10(results_pc[, pvals_m2]),
  betas = results_pc[, alphasmod_m2]
)
results_short[, pch := fifelse(
  log10pvals > -log10(0.05 / length(log10pvals)),
  fifelse(
    betas > 0,
    24,
    25
  ),
  21
)]

p <- results_short |>
  ggplot() +
  geom_hline(
    yintercept = -log10(0.05 / length(results_short[, log10pvals])),
    linetype = 1, color = "red"
  ) +
  geom_point(
    aes(
      x = 1:nrow(results_short),
      y = log10pvals
    ),
    shape = results_short[, pch]
  ) +
  labs(
    x = "Principal components",
    y = "-log10(p-value)",
    title = glue(
      "Associations with phecode {opt$outcome}: ",
      "{pheinfo[phecode == gsub('X', '', opt$outcome), description]}"
    ),
    caption = glue(
      "Cumulative variation explained by PCs: {opt$pc_var_explain}; ",
      "n PCs = {format(length(pcs_mod[1, ]), big.mark = ',')}; ",
      "t = {opt$time_threshold}; ",
      "derivation cohort: {toupper(opt$discovery_cohort)}; ",
      "application cohort: MGI"
    )
  ) +
  theme_minimal() +
  theme(
    legend.position  = "",
    legend.text      = element_text(size = 10),
    legend.title     = element_blank(),
    axis.text.x      = element_text(angle = 60, hjust = 1, vjust = 1),
    plot.caption     = element_text(hjust = 0),
    text             = element_text(size = 12)
  )

ggsave(
  filename = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_mahanttan.pdf"
  ),
  plot = p,
  width = 8, height = 6,
  device = cairo_pdf
)

# beta plot --------------------------------------------------------------------
results <- merge.data.table(
  results,
  pheinfo[, .(
    phecode = glue_data(.SD, "X{phecode}"),
    color,
    desc = description
  )],
  by = "phecode"
)

p2 <- results |>
  ggplot() +
  geom_point(aes(
    x = 1:nrow(results),
    y = betas_m2sig,
    col = color
  ), shape = 16) +
  scale_colour_manual(values = unique(as.character(results[, color]))) +
  geom_label_repel(aes(x = c(1:nrow(results)), y = betas_m2sig, label = desc),
    label.size = 0.1, force = 2, size = 2,
    label.padding = 0.1, point.padding = unit(0.2, "lines"),
    segment.alpha = 0.3
  ) +
  labs(
    x = "Phenotypes",
    y = "Beta",
    title = glue(
      "Associations with phecode {opt$outcome}: ",
      "{pheinfo[phecode == gsub('X', '', opt$outcome), description]}"
    ),
    caption = glue(
      "Cumulative variation explained by PCs: {opt$pc_var_explain}; ",
      "n PCs = {format(length(pcs_mod[1, ]), big.mark = ',')}; ",
      "t = {opt$time_threshold}; ",
      "derivation cohort: {toupper(opt$discovery_cohort)}; ",
      "application cohort: MGI"
    )
  ) +
  theme_minimal() +
  theme(
    plot.caption = element_text(hjust = 0),
    legend.position = "",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    text = element_text(size = 12)
  )

ggsave(
  filename = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_betas.pdf"
  ),
  plot = p2,
  width = 8, height = 6,
  device = cairo_pdf
)

# calculate phers --------------------------------------------------------------
message("calculating phers...")
mgi_phers_vec <- as.matrix(exclude_mgi_pim) %*%
  as.matrix(results[, betas_m2sig])
mgi_phers_vec_std <- scale(mgi_phers_vec)

sub_ukb <- predictor_checker(short_ukb_pim, results[, phecode])

ukb_phers_vec <- as.matrix(sub_ukb) %*%
  as.matrix(results[, betas_m2sig])
ukb_phers_vec_std <- scale(ukb_phers_vec)

mgi_phers <- data.frame(
  id    = mgi_pim[, id],
  case  = mgi_pim[, case],
  pred  = mgi_phers_vec,
  phers = mgi_phers_vec_std
) |> as.data.table()

ukb_phers <- data.frame(
  id    = ukb_pim[, id],
  case  = ukb_pim[, case],
  pred  = ukb_phers_vec,
  phers = ukb_phers_vec_std
) |> as.data.table()

save_qs(
  x = ukb_phers,
  file = glue(
    "{out_path}/ukb_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_phers.qs"
  )
)
save_qs(
  x = mgi_phers,
  path = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_phers.qs"
  )
)

cli_alert("generating outputs...")
suppressMessages({
  mgi_roc <- pROC::roc(mgi_phers[, case], mgi_phers[, phers], smooth = TRUE, se = FALSE)
  mgi_auc <- pROC::ci.auc(mgi_phers[, case], mgi_phers[, phers])

  ukb_roc <- pROC::roc(ukb_phers[, case], ukb_phers[, phers], smooth = TRUE, se = FALSE)
  ukb_auc <- pROC::ci.auc(ukb_phers[, case], ukb_phers[, phers])
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
    title    = glue("AUC for X{gsub('X', '', opt$outcome)} at t{opt$time_threshold}"),
    caption = str_wrap(glue("Discovery cohort: {toupper(opt$discovery_cohort)}"), 80)
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
  filename = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_auc.pdf"
  ),
  width = 7, height = 7, device = cairo_pdf
)

mgi_mod <- glm(case ~ phers, data = mgi_phers)
ukb_mod <- glm(case ~ phers, data = ukb_phers)

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
  .caption = str_wrap(glue("Discovery cohort = {toupper(opt$discovery_cohort)}"),
    width = 100
  )
)
ggsave(
  plot = mgi_phers_dist_plot,
  filename = glue(
    "{out_path}/mgi_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_phers_dist.pdf"
  ),
  width = 8, height = 6, device = cairo_pdf
)

ukb_phers_dist_plot <- top_or_plotr(
  phers_data = ukb_phers,
  .title = glue(
    "X{gsub('X', '', opt$outcome)}",
    " PheRS distribution by case status at t{opt$time_threshold}"
  ),
  .subtitle = glue("In UKB external cohort"),
  .caption = str_wrap(glue("Discovery cohort = {toupper(opt$discovery_cohort)}"),
    width = 100
  )
)
ggsave(
  plot = ukb_phers_dist_plot,
  filename = glue(
    "{out_path}/ukb_X{gsub('X', '', opt$outcome)}_",
    "t{opt$time_threshold}_pve{opt$pc_var_explain}_",
    "pc_phers_dist.pdf"
  ),
  width = 8, height = 6, device = cairo_pdf
)

cli_alert("script success! 🎉 see output in {out_path}")
