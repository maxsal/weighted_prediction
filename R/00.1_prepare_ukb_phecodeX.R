# libraries --------------------------------------------------------------------
ms::libri(ms, qs, data.table, PheWAS/PheWAS, glue, cli)

# load data --------------------------------------------------------------------
ukb_icd <- qread("data/private/ukb/20221117/UKB_PHENOME_ICD_DSB_20221117.qs")
setnames(ukb_icd, c("lexicon", "diagnosis_code"), c("vocabulary_id", "code"))
ukb_icd[, vocabulary_id := toupper(vocabulary_id)]

ukb_icd <- ukb_icd[vocabulary_id == "ICD10", ]
ukb_icd[, ICD10category := gsub("([A-Z][0-9]{2}).+", "\\1", code)]
ukb_icd[, ICD10suffix := gsub("[A-Z].+", "", gsub("^[A-Z][0-9]{2}", "", code))]
ukb_icd[, code := paste0(ICD10category, ifelse(ICD10suffix == "", "", "."), ICD10suffix)]
ukb_icd[, c("ICD10category", "ICD10suffix") := NULL]
ukb_icd <- unique(ukb_icd)

ukb_demo <- qread("data/private/ukb/20221117/UKB_PHENOME_COV_20221117.qs")

# external
phecodex_map <- ms::phecodex_icd_map
setnames(phecodex_map, "icd", "code")
### rollup map
phecodex_rollup <- ms::phecodex_rollup
### sex restriction map
phecodex_sex <- ms::phecodex_sex
### phecodex info
phecodex_info <- ms::phecodex_labels

# map phenotypes ---------------------------------------------------------------
mapped <- mapCodesToPhecodes(
    input = ukb_icd, vocabulary.map = phecodex_map,
    rollup.map = phecodex_rollup, make.distinct = TRUE
)

id_sex <- ukb_demo[, .(id, sex)]

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

in_demo_and_phenome <- intersect(unique(ukb_demo[, id]), unique(restricted[, id]))

restricted <- restricted[id %in% in_demo_and_phenome, ]
restricted <- unique(restricted[, .(id, phecode, dsb)])

## save
qsave(
    restricted,
    file = glue("data/private/ukb/20221117/UKB_FULL_PHECODEX_DSB_20221117.qs")
)

# get first phecode per person -------------------------------------------------
first_restricted <- restricted[restricted[, .I[which.min(dsb)], by = c("id", "phecode")][["V1"]], ]

## save
qsave(
    first_restricted,
    file = glue("data/private/ukb/20221117/UKB_FIRST_PHECODEX_DSB_20221117.qs")
)

# get phecode info -------------------------------------------------------------
phecodex_n <- first_restricted[, .N, phecode][order(-N)]
phecodex_n <- merge.data.table(
    phecodex_n, phecodex_info,
    by.x = "phecode", by.y = "phenotype",
    all.x = TRUE
)[order(-N), ]
phecodex_n[, prev := N / length(unique(restricted$id))]

## save
qsave(
    phecodex_n,
    file = glue("data/private/ukb/20221117/UKB_PHECODEX_N_20221117.qs")
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

qsave(
    pim,
    file = glue("data/private/ukb/20221117/UKB_PIM0X_20221117.qs")
)

### first and last dsb
first_dsb <- unique(restricted[restricted[
    ,
    .I[which.min(dsb)],
    "id"
][["V1"]]][
    ,
    .(
        id,
        first_dsbx = dsb,
        age_at_first_diagnosisx = round(dsb / 365.25, 1)
    )
])
last_dsb <- unique(restricted[restricted[, .I[which.max(dsb)], "id"][["V1"]]][, .(
    id,
    last_dsbx = dsb,
    age_at_last_diagnosisx = round(dsb / 365.25, 1)
)])

ukb_demo <- Reduce(
    \(x, y) merge.data.table(x, y, by = "id"),
    list(ukb_demo, first_dsb, last_dsb)
)

###
cancer_phecodesx <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv")[
    keep == 1, phecode
]

comorbid <- list(
    "cadx" = list("phecodes" = c("CV_404.2")),
    "diabetesx" = list("phecodes" = c("EM_202")),
    "hypertensionx" = list("phecodes" = c("CV_401")),
    "mixed_hypertensionx" = list("phecodes" = c("EM_239.3")),
    "vitamin_dx" = list("phecodes" = c("EM_232.4")),
    "depressionx" = list("phecodes" = c("MB_286.2")),
    "anxietyx" = list("phecodes" = c("MB_288")),
    "bipolarx" = list("phecodes" = c("MB_286.1")),
    "cancerx" = list("phecodes" = cancer_phecodesx)
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
    set(ukb_demo, j = i, value = fifelse(ukb_demo[["id"]] %in% comorbid[[i]][["ids"]], 1, 0))
    cli_progress_update()
}
cli_progress_done()

ukb_demo[, triglyceridesx := fifelse(hypertensionx == 0 & mixed_hypertensionx == 0, 0, 1)]

ukb_demo[, age_catx := between(age_at_last_diagnosisx, 0, 5.99) +
    2 * between(age_at_last_diagnosisx, 6, 11.99) +
    3 * between(age_at_last_diagnosisx, 12, 19.99) +
    4 * between(age_at_last_diagnosisx, 20, 39.99) +
    5 * between(age_at_last_diagnosisx, 40, 59.99) +
    6 * between(age_at_last_diagnosisx, 60, 150.99)]

qsave(
    ukb_demo,
    file = glue("data/private/ukb/20221117/datax_20221117_comb.qs")
)

cli_alert_success("Finished! 🎉")
