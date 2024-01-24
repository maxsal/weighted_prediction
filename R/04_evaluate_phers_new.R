# evaluate a phers
# author: max salvatore
# date:   20230418

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)
ms::libri(
  ms, data.table, caret, purrr, progress, pROC, glue, logistf,
  optparse, cli, parallelly
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "157",
    help = "Outcome phecode [default = 157]"
  ),
  make_option("--mgi_version",
    type = "character", default = "20210318",
    help = "Version of MGI data [default = 20210318]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Version of UKB data [default = 20221117]"
  ),
  make_option("--time_threshold",
    type = "numeric", default = "0",
    help = glue(
      "Time threshold for the phenome data ",
      "[default = 0]"
    )
  ),
  make_option("--tophits_n",
    type = "numeric", default = "50",
    help = glue(
      "Number of top hits to use in top hits PheRS ",
      "[default = 50]"
    )
  ),
  make_option("--pctile_or",
    type = "logical", default = "TRUE",
    help = glue(
      "Perform percentile-based OR diagnostics ",
      "which is relatively time consuming ",
      "[default = TRUE]"
    )
  )
)

#### !!! ADD arguments specifying covariates

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

# 2. specifications (specifies outcome) ----------------------------------------
## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# read data --------------------------------------------------------------------
## mgi
### phers
#### naive
mgi_naive_phers <- fread(glue(
  "results/mgi/{opt$mgi_version}/",
  "X{gsub('X', '', opt$outcome)}/phers/",
  "mgi_naive_phers_t{opt$time_threshold}.txt"
))

mgi_naive_phers_info <- list()
for (i in grep(
  list.files(
    glue(
      "results/mgi/{opt$mgi_version}/",
      "X{gsub('X', '', opt$outcome)}/phers"
    ),
    full.names = TRUE
  ),
  pattern = glue("mgi_phers_t{opt$time_threshold}"),
  value = TRUE
)) {
  mgi_naive_phers_info[[i]] <- readRDS(i)
}

### covariates
mgi_covariates <- fread(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "X{gsub('X', '', opt$outcome)}/",
  "matched_covariates.txt"
))[
  id %in% mgi_naive_phers[, id]
][
  , age_at_threshold := round(get(glue("t{opt$time_threshold}_threshold")) / 365.25, 1)
]

mgi_merged <- merge.data.table(
  mgi_naive_phers,
  mgi_covariates[, !c("case")],
  by = "id"
)

## ukb
### phers
ukb_naive_phers <- fread(glue(
  "results/ukb/{opt$ukb_version}/",
  "X{gsub('X', '', opt$outcome)}/phers/",
  "ukb_naive_phers_t{opt$time_threshold}.txt"
))

ukb_naive_phers_info <- list()
for (i in grep(
  list.files(
    glue(
      "results/ukb/{opt$ukb_version}/",
      "X{gsub('X', '', opt$outcome)}/phers"
    ),
    full.names = TRUE
  ),
  pattern = glue("ukb_phers_t{opt$time_threshold}"),
  value = TRUE
)) {
  ukb_naive_phers_info[[i]] <- readRDS(i)
}

### covariates
ukb_covariates <- fread(glue(
  "data/private/ukb/{opt$ukb_version}/",
  "X{gsub('X', '', opt$outcome)}/",
  "matched_covariates.txt"
))[
  id %in% ukb_naive_phers[, id]
][
  , age_at_threshold := round(get(glue("t{opt$time_threshold}_threshold")) / 365.25, 1)
]

ukb_merged <- merge.data.table(
  ukb_naive_phers,
  ukb_covariates[, !c("case")],
  by = "id"
)

## phenome
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character"
)

# scale phers ------------------------------------------------------------------
# remove trailing number from phers - FIX THIS IN NAIVE PHERS SCRIPT
setnames(
  mgi_merged,
  old = grep(x = names(mgi_merged), pattern = "phers", value = TRUE),
  new = gsub("phers[0-9]", "phers", grep(x = names(mgi_merged), pattern = "phers", value = TRUE)),
  skip_absent = TRUE
)
for (i in grep(x = names(mgi_merged), pattern = "phers", value = TRUE)) {
  mgi_merged[, glue("{i}_scaled") := scale(get(i))]
}

setnames(
  ukb_merged,
  old = grep(x = names(ukb_merged), pattern = "phers", value = TRUE),
  new = gsub("phers[0-9]", "phers", grep(x = names(ukb_merged), pattern = "phers", value = TRUE)),
  skip_absent = TRUE
)
for (i in grep(x = names(ukb_merged), pattern = "phers", value = TRUE)) {
  ukb_merged[, glue("{i}_scaled") := scale(get(i))]
}

# evaluate ---------------------------------------------------------------------
# mgi_evals <- list()
cli_progress_bar(
  name = "evaluating mgi",
  total = length(grep("_scaled", names(mgi_merged)))
)
for (i in grep(pattern = "_scaled", names(mgi_merged), value = TRUE)) {
  evaluate_phers(
    merged_data = mgi_merged,
    covars = c("age_at_threshold", "female"),
    outcome = opt$outcome,
    n_phecodes = mgi_naive_phers_info[[grep(
      x = names(mgi_naive_phers_info),
      pattern = gsub("_scaled", "", i),
      value = TRUE
    )]][["n_phecodes"]],
    phers_name = i,
    pctile_or = opt$pctile_or
  ) |>
    fwrite(glue(
      "results/mgi/{opt$mgi_version}/",
      "X{gsub('X', '', opt$outcome)}/phers/",
      "{i}_eval.txt"
    ), sep = "\t")
  cli_progress_update()
}
cli_progress_done()

# ukb_evals <- list()
cli_progress_bar(
  name = "evaluating ukb",
  total = length(grep("_scale", names(ukb_merged)))
)
for (i in grep(pattern = "_scaled", names(ukb_merged), value = TRUE)) {
  evaluate_phers(
    merged_data = ukb_merged,
    covars = c("age_at_threshold", "female"),
    outcome = opt$outcome,
    n_phecodes = ukb_naive_phers_info[[grep(
      x = names(ukb_naive_phers_info),
      pattern = gsub("_scaled", "", i),
      value = TRUE
    )]][["n_phecodes"]],
    phers_name = i,
    pctile_or = opt$pctile_or
  ) |>
    fwrite(glue(
      "results/ukb/{opt$ukb_version}/",
      "X{gsub('X', '', opt$outcome)}/phers/",
      "{i}_eval.txt"
    ), sep = "\t")
  cli_progress_update()
}
cli_progress_done()

cli_alert_success("script complete! 🎉")
