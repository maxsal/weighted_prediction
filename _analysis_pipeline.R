# aim one analysis pipeline
ms::libri(maxsal/ms, glue, cli)

mgi_version     <- "20220822"
ukb_version     <- "20221117"
outcome         <- "157"                 # pancreatic cancer
time_thresholds <- c(0, 0.5, 1, 2, 3, 5)
use_geno_pcs    <- FALSE

# phase 0 scripts --------------------------------------------------------------
## prepare mgi data
system(glue(
  "/usr/bin/time -v -o logs/00_prepare_mgi_data.txt Rscript R/",
  "00_prepare_mgi_data.R --mgi_version={mgi_version}"
))

## prepare ukb data
system(glue(
  "/usr/bin/time -v -o logs/00_create_ukb_datasets.txt Rscript R/",
  "00_create_ukb_datasets.R"
))

# phase 1 scripts --------------------------------------------------------------
## prepare time-restricted phenomes in MGI and UKB #
system(glue("/usr/bin/time -v -o logs/01_prepare_phenomes.txt Rscript R/",
            "01_prepare_phenomes.R --mgi_version={mgi_version} --ukb_version={ukb_version} ",
            "--time_thresholds={paste0(time_thresholds, collapse = ',')} --outcome={outcome}"))

## estimate ipw and poststratification weights in MGI
system(glue("/usr/bin/time -v -o logs/01_estimate_weights.txt Rscript R/",
            "01_estimate_weights.R --cohort_version={mgi_version}"))

## phase 2 scripts -------------------------------------------------------------
cli_alert("running unweighted cooccurrence")
system(glue("/usr/bin/time -v -o logs/02_unweighted_coccurrence_analysis.txt Rscript R/",
            "02_unweighted_cooccurrence_analysis.R --mgi_version={mgi_version} --ukb_version={ukb_version} ",
            "--time_thresholds={paste0(time_thresholds, collapse = ',')} --outcome={outcome}"))

for (w in c("ip_selection_f", "ps_nhw_f")) {
  cli_alert_info("running {w} weighted cooccurrence")
  system(glue("/usr/bin/time -v -o logs/02_weighted_cooccurrence_analysis.txt Rscript R/",
              "02_weighted_cooccurrence_analysis.R --mgi_version={mgi_version} ",
              "--time_thresholds={paste0(time_thresholds, collapse = ',')} --outcome={outcome} ",
              "--weights={w}"))
}
system(glue("/usr/bin/time -v -o logs/02_ukb_weighted_cooccurrence_analysis.txt Rscript R/",
            "02_ukb_weighted_cooccurrence_analysis.R --ukb_version={ukb_version} ",
            "--time_thresholds={paste0(time_thresholds, collapse = ',')} --outcome={outcome}"))

for (i in seq_along(time_thresholds)) {
  cli_alert_info("running random forest at {time_thresholds[i]}")
  system(glue("/usr/bin/time -v -o logs/02_random_forest.txt Rscript R/",
              "02_random_forest.R --mgi_version={mgi_version} --ukb_version={ukb_version} ",
              "--time_threshold={time_thresholds[i]} --outcome={outcome}"))
}

for (i in seq_along(time_thresholds)) {
  cli_alert("running SuperLearner at {time_thresholds[i]}")
  system(glue("/usr/bin/time -v -o logs/02_super_learner.txt Rscript R/",
              "02_super_learner.R --mgi_version={mgi_version} --ukb_version={ukb_version} ",
              "--time_threshold={time_thresholds[i]} --outcome={outcome}"))
}

## phase 3 scripts -------------------------------------------------------------
for (i in time_thresholds) {
  cli_alert("calculating unweighted naive phers at {i}")
  system(glue("/usr/bin/time -v -o logs/03_calculate_naive_phers.txt Rscript R/",
              "03_calculate_naive_phers.R --outcome={outcome} --time_threshold={i}"))
}

for (i in time_thresholds) {
  cli_alert("calculating unweighted naive phers at {i}")
  system(glue("/usr/bin/time -v -o logs/03_calculate_naive_phers.txt Rscript R/",
              "03_calculate_naive_phers.R --outcome={outcome} --time_threshold={i} ",
              "--method=tophits"))
}


for (i in c("cancer_ipw", "cancer_postw")) {
  for (j in time_thresholds) {
    cli_alert("calculating {i} weighted naive pwide sig phers at {j}")
    system(glue("/usr/bin/time -v -o logs/03_calculate_naive_phers.txt Rscript R/",
                "03_calculate_naive_phers.R --outcome={outcome} --time_threshold={j} ",
                "--weights={i}"))
  }
}
for (i in c("cancer_ipw", "cancer_postw")) {
  for (j in time_thresholds) {
    cli_alert("calculating {i} weighted naive top hits phers at {j}")
    system(glue("/usr/bin/time -v -o logs/03_calculate_naive_phers.txt Rscript R/",
                "03_calculate_naive_phers.R --outcome={outcome} --time_threshold={j} ",
                "--weights={i} --method=tophits"))
  }
}

## phase 4 scripts -------------------------------------------------------------
for (i in time_threshold) {
  system(glue("/usr/bin/time -v -o logs/04_evaluate_phers_new.txt Rscript R/",
              "04_evaluate_phers_new.R --outcome={outcome} --time_threshold={i}"))
}


