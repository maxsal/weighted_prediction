# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, PheWAS/PheWAS, glue, optparse, cli,
    labelled, comorbidity
)

# optparse list ---
option_list <- list(
    make_option("--mgi_version",
        type = "character", default = "20220822",
        help = "MGI cohort version [default = '20220822']"
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

# load data --------------------------------------------------------------------
## internal
mgi_icd  <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_ICD_{opt$mgi_version}.qs"))
mgi_demo <- qread(glue("data/private/mgi/{opt$mgi_version}/data_{opt$mgi_version}_comb.qs"))
setnames(mgi_icd, c("IID", "DiagnosisCode", "lexicon", "DaysSinceBirth"), c("id", "code", "vocabulary_id", "dsb"))

## external
### mapping table
phecodex_map <- ms::phecodex_icdcm_map
### rollup map
phecodex_rollup <- ms::phecodex_rollup
### sex restriction map
phecodex_sex <- ms::phecodex_sex
### phecodex info
phecodex_info <- ms::phecodex_labels

# map phenotypes ---------------------------------------------------------------
mapped <- mapCodesToPhecodes(
    input = mgi_icd, vocabulary.map = phecodex_map,
    rollup.map = phecodex_rollup, make.distinct = TRUE
)

id_sex <- mgi_demo[, .(id = DeID_PatientID, sex = Sex)]

# remove sex-specific phecodes that are discordant with individual's reported sex at birth
restrictPhecodesBySex_mod <- function(phenotypes, id.sex, by_var = "person_id", sex_var = "sex") {
    data <- merge.data.table(
        as.data.table(phenotypes),
        as.data.table(id.sex),
        by = by_var, all.x = TRUE
    )
    # Get the restrictions found in the phenotypes data frame
    current_sex_restriction <- phecodex_sex[phecodex_sex$phecode %in% unique(data[, phecode]), ] |>
        as.data.table()
    # Get male and female-only phenotypes
    male_only <- current_sex_restriction[current_sex_restriction$male_only, phecode]
    female_only <- current_sex_restriction[current_sex_restriction$female_only, phecode]
    # Set row column matches to NA where inds of a sex meet restricted phenotypes
    data[phecode %in% male_only & sex == "Female", phecode := NA]
    data[phecode %in% female_only & sex == "Male", phecode := NA]

    na.omit(data)[, (sex_var) := NULL][]
}

## remove sex person-phenotype discordant pairs
restricted <- restrictPhecodesBySex_mod(mapped, id.sex = id_sex, by_var = "id")

in_demo_and_phenome <- intersect(unique(mgi_demo[, DeID_PatientID]), unique(restricted[, id]))

restricted <- restricted[id %in% in_demo_and_phenome, ]

## save
qsave(
    restricted,
    file = glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs")
)

# get first phecode per person -------------------------------------------------
first_restricted <- restricted[ restricted[, .I[which.min(dsb)], by = c("id", "phecode")][["V1"]], ]

## save
qsave(
    first_restricted,
    file = glue("data/private/mgi/{opt$mgi_version}/MGI_FIRST_PHECODEX_DSB_{opt$mgi_version}.qs")
)

# get phecode info -------------------------------------------------------------
phecodex_n <- first_restricted[, .N, phecode][order(-N)]
phecodex_n <- merge.data.table(
    phecodex_n, phecodex_info, by.x = "phecode", by.y = "phenotype",
    all.x = TRUE
)[order(-N), ]
phecodex_n[, prev := N / length(unique(first_restricted$id))]

## save
qsave(
    phecodex_n,
    file = glue("data/private/mgi/{opt$mgi_version}/MGI_PHECODEX_N_{opt$mgi_version}.qs")
)

# pim --------------------------------------------------------------------------
pim <- dcast(
    first_restricted[, .(id, phecode = phecode)],
    id ~ phecode,
    value.var = "phecode",
    fun.aggregate = length,
    fill = 0
)

# collect integer variable names from data.table
int_vars <- names(pim)[which(sapply(pim, is.integer))]
pim[, (int_vars) := lapply(.SD, \(x) as.numeric(x > 0)), .SDcols = int_vars]

for (i in seq_along(int_vars)) {
    var_label(pim[[int_vars[i]]]) <- ms::pheinfox[phecode == int_vars[i], description]
}

qsave(
    pim,
    file = glue("data/private/mgi/{opt$mgi_version}/MGI_PIM0X_{opt$mgi_version}.qs")
)

### first and last dsb
first_dsb <- unique(restricted[restricted[
  ,
  .I[which.min(dsb)],
  "id"
][["V1"]]][
  ,
  .(
    DeID_PatientID = id,
    first_dsbx     = dsb,
    age_at_first_diagnosisx = round(dsb / 365.25, 1)
  )
])
last_dsb <- unique(restricted[restricted[, .I[which.max(dsb)], "id"][["V1"]]][, .(
    DeID_PatientID = id,
    last_dsbx     = dsb,
    age_at_last_diagnosisx = round(dsb / 365.25, 1)
)])

mgi_demo <- Reduce(
    \(x, y) merge.data.table(x, y, by = "DeID_PatientID"),
    list(mgi_demo, first_dsb, last_dsb)
)

###
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")[
    keep == 1, phecode
]

comorbid <- list(
    "cadx"                = list("phecodes" = c("CV_404.2")),
    "diabetesx"           = list("phecodes" = c("EM_202")),
    "hypertensionx"       = list("phecodes" = c("CV_401")),
    "mixed_hypertensionx" = list("phecodes" = c("EM_239.3")),
    "vitamin_dx"          = list("phecodes" = c("EM_232.4")),
    "depressionx"         = list("phecodes" = c("MB_286.2")),
    "anxietyx"            = list("phecodes" = c("MB_288")),
    "bipolarx"            = list("phecodes" = c("MB_286.1")),
    "cancerx"             = list("phecodes" = cancer_phecodesx)
)

# identify cases and create indicator variables --------------------------------
cli_alert("identifying cases and creating indicator variables...")
cli_progress_bar("deriving comorbidities", total = length(names(comorbid)))
for (i in names(comorbid)) {
    comorbid[[i]][["ids"]] <- restricted[phecode %in% comorbid[[i]][["phecodes"]], id] |>
        unique()
    cli_progress_update()
}
cli_progress_done()

cli_progress_bar("deriving comorbidity indicator variables", total = length(names(comorbid)))
for (i in names(comorbid)) {
    set(mgi_demo, j = i, value = fifelse(mgi_demo[["DeID_PatientID"]] %in% comorbid[[i]][["ids"]], 1, 0))
    cli_progress_update()
}
cli_progress_done()

mgi_demo[, triglyceridesx := fifelse(hypertensionx == 0 & mixed_hypertensionx == 0, 0, 1)]

mgi_demo[, age_catx := between(age_at_last_diagnosisx, 0, 5.99) +
    2 * between(age_at_last_diagnosisx, 6, 11.99) +
    3 * between(age_at_last_diagnosisx, 12, 19.99) +
    4 * between(age_at_last_diagnosisx, 20, 39.99) +
    5 * between(age_at_last_diagnosisx, 40, 59.99) +
    6 * between(age_at_last_diagnosisx, 60, 150.99)]

###
# ADD INDICATOR FOR FIRST CANCER DIAGNOSIS
first_cancer <- first_restricted[phecode %in% cancer_phecodesx, ]
first_cancer <- first_cancer[
    first_cancer[, .I[dsb == min(dsb, na.rm = TRUE)], id][["V1"]]
]

mgi_first_cancer <- data.table(
    id = first_cancer[, unique(id)]
)

first_cancer_ids <- list()
for (i in seq_along(cancer_phecodesx)) {
    first_cancer_ids[[i]] <- first_cancer[phecode == cancer_phecodesx[i], unique(id)]
    names(first_cancer_ids)[i] <- cancer_phecodesx[i]
}

for (i in seq_along(cancer_phecodesx)) {
    set(
        mgi_first_cancer,
        j     = paste0(cancer_phecodesx[i], "_first"),
        value = fifelse(mgi_first_cancer[["id"]] %in% first_cancer_ids[[i]], 1, 0)
    )
}

### comorbidity scores
# charlson comorbidity index
charlson10 <- comorbidity(mgi_icd[vocabulary_id == "ICD10CM", ], id = "id", code = "code", map = "charlson_icd10_quan", assign0 = FALSE)
charlson9 <- comorbidity(mgi_icd[vocabulary_id == "ICD9CM", ], id = "id", code = "code", map = "charlson_icd9_quan", assign0 = FALSE)

charlson <- merge.data.table(
    charlson10,
    charlson9,
    by = "id",
    all = TRUE
) |> as.data.table()

charl_comorbid_names <- gsub(".x", "", grep(".x", names(charlson), value = TRUE))
for (i in seq_along(charl_comorbid_names)) {
    charlson[, (charl_comorbid_names[i]) := as.numeric(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c(paste0(charl_comorbid_names[i], ".x"), paste0(charl_comorbid_names[i], ".y"))]
}
keep_charl <- c("id", charl_comorbid_names)
charlson <- charlson[, ..keep_charl]
attr(charlson, "class") <- attributes(charlson10)$class
attr(charlson, "variable.labels") <- attributes(charlson10)$variable.labels
attr(charlson, "map") <- attributes(charlson10)$map
un_charl_scores <- score(charlson, weights = NULL, assign0 = FALSE)
we_charl_scores <- score(charlson, weights = "quan", assign0 = FALSE)

charl_scores <- data.table(
    id = charlson$id,
    charlson_unweighted = un_charl_scores,
    charlson_quan_weighted = we_charl_scores
)

# elixhauser comorbidity index
elixhauser10 <- comorbidity(mgi_icd[vocabulary_id == "ICD10CM", ], id = "id", code = "code", map = "elixhauser_icd10_quan", assign0 = FALSE)
elixhauser9 <- comorbidity(mgi_icd[vocabulary_id == "ICD9CM", ], id = "id", code = "code", map = "elixhauser_icd9_quan", assign0 = FALSE)

elixhauser <- merge.data.table(
    elixhauser10,
    elixhauser9,
    by = "id",
    all = TRUE
) |> as.data.table()

elix_comorbid_names <- gsub(".x", "", grep(".x", names(elixhauser), value = TRUE))
for (i in seq_along(elix_comorbid_names)) {
    elixhauser[, (elix_comorbid_names[i]) := as.numeric(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = c(paste0(elix_comorbid_names[i], ".x"), paste0(elix_comorbid_names[i], ".y"))]
}
keep_elix <- c("id", elix_comorbid_names)
elixhauser <- elixhauser[, ..keep_elix]
attr(elixhauser, "class") <- attributes(elixhauser10)$class
attr(elixhauser, "variable.labels") <- attributes(elixhauser10)$variable.labels
attr(elixhauser, "map") <- attributes(elixhauser10)$map
un_elix_scores <- score(elixhauser, weights = NULL, assign0 = FALSE)
we_elix_scores <- score(elixhauser, weights = "vw", assign0 = FALSE)

elix_scores <- data.table(
    id = elixhauser$id,
    elixhauser_unweighted = un_elix_scores,
    elixhauser_vw_weighted = we_elix_scores
)

comorbid_scores <- merge.data.table(
    charl_scores,
    elix_scores,
    by = "id",
    all = TRUE
) |> as.data.table()

mgi_demo <- merge.data.table(
    mgi_demo,
    comorbid_scores,
    by.x = "DeID_PatientID",
    by.y = "id",
    all = TRUE
)
###

setnames(mgi_demo,
    old = c(
        "DeID_PatientID", "Age", "AliveYN", "DeceasedDaysSinceBirth",
        "Ethnicity", "MaritalStatusCode", "Sex", "Race", "YearsInEHR",
        "FirstDaySinceBirth", "LastDaySinceBirth", "DSB_BMI", "BMI_class",
        "BMI_category", "StudyName", "DaysWithEntries", "DaysInEHR",
        "Enrollment_DaysSinceBirth", "Zip3"
    ),
    new = c(
        "id", "age", "alive", "dead_dsb", "ethn", "marital",
        "sex", "race",
        "length_followup", "first_dsb", "last_dsb",
        "bmi_dsb", "bmi_class", "bmi_category", "study_name",
        "days_with_entries", "days_in_ehr", "enrollment_dsb",
        "zip3"
    ),
    skip_absent = TRUE
)

### add race_eth, nhw, and female
mgi_demo[, `:=` (
    race_eth = relevel(factor(fcase(
        race == "Caucasian" & ethn == "Non-Hispanic", "NHW",
        race == "African American" & ethn == "Non-Hispanic", "NHB",
        ethn == "Hispanic", "Hisp",
        !is.na(race) & !is.na(ethn), "Other",
        default = "Unknown"
    )), ref = "NHW")
)]
mgi_demo[, `:=` (
    nhw = as.numeric(race_eth == "NHW"),
    female = as.numeric(sex == "Female")
)]
###

### drop unused variables
drop_vars <- c(
    "Deceased_DaysSinceBirth", "RaceEthnicity", "RaceEthnicity4",
    "Adult", "AgeCategory", "Drinker", "SmokingStatus", "Smoker",
    "AgeFirstEntry", "AgeLastEntry", "AgeAtDeath")
mgi_demo[, (drop_vars) := NULL]

###
qsave(
    mgi_first_cancer,
    file = glue("data/private/mgi/{opt$mgi_version}/first_cancerx_{opt$mgi_version}.qs")
)
###

qsave(
    mgi_demo,
    file = glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs")
)

cli_alert_success("Finished! 🎉")
