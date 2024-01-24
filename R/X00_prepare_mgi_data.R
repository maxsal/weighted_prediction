# Prepare MGI data including deriving variables for comorbidity status and
# descriptive variables derived from demographic data
# author:   max salvatore
# date:     20230809

# libraries, paths, and such ---------------------------------------------------
ms::libri(
  data.table, glue, qs, parallel, optparse, PheWAS/PheWAS, cli, ms
)

# optparse list ----
option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

cli_alert(glue("using mgi cohort version {opt$mgi_version}"))

out_path <- glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/private/mgi/{opt$mgi_version}/")
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

source("fn/files-utils.R")
file_paths <- get_files(mgi_version = opt$mgi_version)

# load data --------------------------------------------------------------------
cli_alert("loading data...")
study <- fread("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/MGI_Study_FirstEnrollment_20221102.txt")
MGIcohort <- fread(file_paths[["mgi"]][["cov_file"]])
cancer_phecodes <- fread("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/public/cancer_phecodes.txt",
  colClasses = "character"
)[[1]]

### subset
MGIcohort <- merge.data.table(
  MGIcohort,
  study,
  by = "DeID_PatientID"
)[StudyName != "AOS", ]

### MAP ICD9/ICD10 codes
mgi_icd9  <- get(load(file_paths[["mgi"]][["icd9_file"]]))
mgi_icd10 <- get(load(file_paths[["mgi"]][["icd10_file"]]))

mgi_icd <- rbindlist(list(
  mgi_icd9[, .(IID, vocabulary_id = "ICD9CM", code = DiagnosisCode, DaysSinceBirth)],
  mgi_icd10[, .(IID, vocabulary_id = "ICD10CM", code = DiagnosisCode, DaysSinceBirth)]
))

in_demo_and_phenome <- intersect(unique(MGIcohort[, DeID_PatientID]), unique(mgi_icd[, IID]))

rm(mgi_icd9, mgi_icd10)

save_qs(
  mgi_icd[, .(IID, lexicon = vocabulary_id, DiagnosisCode = code, DaysSinceBirth)],
  file = glue("{out_path}MGI_ICD_{opt$mgi_version}.qs")
)

mapped <- mapCodesToPhecodes(
  input = mgi_icd, vocabulary.map = PheWAS::phecode_map,
  rollup.map = PheWAS::phecode_rollup_map, make.distinct = TRUE
)

id.sex <- MGIcohort[, .(
  IID = DeID_PatientID, sex = Sex
)]

# remove sex-specific phecodes that are discordant with individual's reported sex at birth
restrictPhecodesBySex_mod <- function(phenotypes, id.sex, by_var = "person_id", sex_var = "sex") {
  data <- merge.data.table(
    as.data.table(phenotypes),
    as.data.table(id.sex),
    by = by_var, all.x = TRUE
  )
  # Get the restrictions found in the phenotypes data frame
  current_sex_restriction <- PheWAS::sex_restriction[PheWAS::sex_restriction$phecode %in% unique(data[, phecode]), ] |>
    as.data.table()
  # Get male and female-only phenotypes
  male_only <- current_sex_restriction[current_sex_restriction$male_only, phecode]
  female_only <- current_sex_restriction[current_sex_restriction$female_only, phecode]
  # Set row column matches to NA where inds of a sex meet restricted phenotypes
  data[phecode %in% male_only & sex == "F", phecode := NA]
  data[phecode %in% female_only & sex == "M", phecode := NA]

  na.omit(data)[, (sex_var) := NULL][]
}

restricted <- restrictPhecodesBySex_mod(mapped, id.sex, by_var = "IID")

in_demo_and_phenome <- intersect(unique(MGIcohort[, DeID_PatientID]), unique(restricted[, IID]))

restricted <- restricted[IID %in% in_demo_and_phenome, ]

save_qs(restricted, file = glue("{out_path}MGI_FULL_PHECODE_DSB_{opt$mgi_version}.qs"))
first_restricted <- restricted[ restricted[, .I[which.min(DaysSinceBirth)], by = c("IID", "phecode")][["V1"]] ]
save_qs(first_restricted, file = glue("{out_path}MGI_FIRST_PHECODE_DSB_{opt$mgi_version}.qs"))

pim <- dcast(
  first_restricted[, .(IID, phecode = paste0("X", phecode))],
  IID ~ phecode,
  value.var = "phecode",
  fun.aggregate = length,
  fill = 0
)
# collect integer variable names from data.table
int_vars <- names(pim)[which(sapply(pim, is.integer))]
pim[, (int_vars) := lapply(.SD, \(x) as.numeric(x > 0)), .SDcols = int_vars]

save_qs(pim, file = glue("{out_path}MGI_PIM0_{opt$mgi_version}.qs"))

### replace FirstDaySinceBirth and LastDaySinceBirth
# following a Slack conversation with Lars on 2/21/23, the original values for
# these variables include ICD codes that cannot be mapped to a phecode and
# non-ICD code records (e.g., OrderDate_DaysSinceBirth from LabResults)
MGIcohort <- MGIcohort[, !c("AgeFirstEntry", "AgeLastEntry", "FirstDaySinceBirth", "LastDaySinceBirth")]

first_dsb <- unique(restricted[restricted[
  ,
  .I[which.min(DaysSinceBirth)],
  "IID"
][["V1"]]][
  ,
  .(
    DeID_PatientID         = IID,
    FirstDaySinceBirth     = DaysSinceBirth,
    age_at_first_diagnosis = round(DaysSinceBirth / 365.25, 1)
  )
])
last_dsb <- unique(restricted[restricted[, .I[which.max(DaysSinceBirth)], "IID"][["V1"]]][, .(
  DeID_PatientID        = IID,
  LastDaySinceBirth     = DaysSinceBirth,
  age_at_last_diagnosis = round(DaysSinceBirth / 365.25, 1)
)])

MGIcohort <- Reduce(
  \(x, y) merge.data.table(x, y, by = "DeID_PatientID"),
  list(MGIcohort, first_dsb, last_dsb)
)
###

comorbid <- list(
  "cad"                = list("phecodes" = c("411.4")),
  "diabetes"           = list("phecodes" = c("250")),
  "hypertension"       = list("phecodes" = c("272.12")),
  "mixed_hypertension" = list("phecodes" = c("272.13")),
  "vitamin_d"          = list("phecodes" = c("261.4")),
  "depression"         = list("phecodes" = c("296.2")),
  "anxiety"            = list("phecodes" = c("300")),
  "bipolar"            = list("phecodes" = c("296.1")),
  "cancer"             = list("phecodes" = cancer_phecodes)
)

# identify cases and create indicator variables --------------------------------
cli_alert("identifying cases and creating indicator variables...")
cli_progress_bar("deriving comorbidities", total = length(names(comorbid)))
for (i in names(comorbid)) {
  comorbid[[i]][["ids"]] <- restricted[phecode %in% comorbid[[i]][["phecodes"]], IID] |>
    unique()
  cli_progress_update()
}
cli_progress_done()

cli_progress_bar("deriving comorbidity indicator variables", total = length(names(comorbid)))
for (i in names(comorbid)) {
  set(MGIcohort, j = i, value = fifelse(MGIcohort[["DeID_PatientID"]] %in% comorbid[[i]][["ids"]], 1, 0))
  cli_progress_update()
}
cli_progress_done()

MGIcohort[, triglycerides := fifelse(hypertension == 0 & mixed_hypertension == 0, 0, 1)]

MGIcohort[, `:=` (
  triglycerides = fifelse(hypertension == 0 & mixed_hypertension == 0, 0, 1),
  nhanes_nhw = fifelse(Ethnicity != "Hispanic" & Race == "Caucasian", 1, 0),
  AgeFirstEntry = age_at_first_diagnosis,
  AgeLastEntry = age_at_last_diagnosis
)][, nhw := nhanes_nhw]

MGIcohort[, age_cat := between(AgeLastEntry, 0, 5.99) +
  2 * between(AgeLastEntry, 6, 11.99) +
  3 * between(AgeLastEntry, 12, 19.99) +
  4 * between(AgeLastEntry, 20, 39.99) +
  5 * between(AgeLastEntry, 40, 59.99) +
  6 * between(AgeLastEntry, 60, 150.99)]

setnames(MGIcohort, "BMI", "bmi")
MGIcohort[
  , bmi_cat := fcase(
    between(bmi, 0, 18.499), 1, # underweight
    between(bmi, 18.5, 24.999), 2, # "normal"
    between(bmi, 25.0, 29.999), 3, # overweight
    between(bmi, 30, 120), 4
  ) # obese
][, `:=`(
  bmi_under       = as.numeric(bmi_cat == 1),
  bmi_overweight  = as.numeric(bmi_cat == 3),
  bmi_obese       = as.numeric(bmi_cat == 4),
  smoking_current = as.numeric(SmokingStatus == "Current"),
  smoking_former  = as.numeric(SmokingStatus == "Past")
)]

MGIcohort[, female := as.numeric(Sex == "F")]

# saving files -----------------------------------------------------------------
cli_alert("saving processed files...")

save_qs(MGIcohort, file = glue("{out_path}data_{opt$mgi_version}_comb.qs"))
for (i in c("MGI", "MGI-MEND", "MIPACT", "MHB2")) {
  save_qs(MGIcohort[StudyName == i, ], file = glue("{out_path}data_{opt$mgi_version}_{tolower(i)}.qs"))
}

cli_alert_success("script success! see {.path {out_path}} for output files")
