# read and prepare UKB data
# based on lars' scripts from: https://github.com/umich-cphds/createUKBphenome
# author: max salvatore
# date:   20230809

# libraries --------------------------------------------------------------------
# install.packages("remotes")
# remotes::install_github("maxsal/ms")
ms::libri(
  data.table, optparse, janitor, PheWAS/PheWAS, cli, qs
)

# optparse list ----
option_list <- list(
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "UKB cohort version [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

source("/net/junglebook/home/mmsalva/projects/ukb_data_explorations/scripts/function.extractUKBdata.r")
source("/net/junglebook/home/mmsalva/projects/ukb_data_explorations/scripts/function.harmonizeICD9.r")
source("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/fn/files-utils.R")

out_path <- paste0("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/private/ukb/", opt$ukb_version, "/")

# pull data --------------------------------------------------------------------
## icd9 data
cli_alert("Pulling ICD9 data...")
icd9 <- reformat_ukb(
  fields = c(41271, 41281),
  data_coding = FALSE, only_info = FALSE, keep_instances = TRUE, recode = FALSE
)
icd9 <- merge.data.table(
  icd9$data[[1]],
  icd9$data[[2]],
  by.x = c("f.eid", "i.41271"),
  by.y = c("f.eid", "i.41281")
)
setnames(icd9, old = c("f.41271", "f.41281"), new = c("diagnosis_code", "date"))
icd9[, date := as.Date(date)]
icd9 <- icd9[, .(id = as.character(f.eid), diagnosis_code, date)]

## icd10 data
cli_alert("Pulling ICD10 data...")
icd10 <- reformat_ukb(
  fields = c(41270, 41280),
  data_coding = FALSE, only_info = FALSE, keep_instances = TRUE, recode = FALSE
)
icd10 <- merge.data.table(
  icd10$data[[1]],
  icd10$data[[2]],
  by.x = c("f.eid", "i.41270"),
  by.y = c("f.eid", "i.41280")
)
setnames(icd10, old = c("f.41270", "f.41280"), new = c("diagnosis_code", "date"))
icd10[, date := as.Date(date)]
icd10 <- icd10[, .(id = as.character(f.eid), diagnosis_code, date)]

# # pull UKB covariate data by field
# dob <- reformat_ukb(
#   fields = c(31, 34, 52, 22001, 21000, 20116, 20117, 21001, 200, 6142, 680, 728, 709, 2178, 6138),
#   data_coding    = TRUE,
#   only_info      = FALSE,
#   keep_instances = TRUE,
#   recode         = TRUE
# )

# # ### save raw pull
# save_qs(
#   dob,
#   paste0("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/private/ukb/", opt$ukb_version, "/dob_raw.qs")
#   )

### read raw pull
dob <- read_qs(
  paste0("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/private/ukb/", opt$ukb_version, "/dob_raw.qs")
)

###
## smoking [f.20116] -----------------------------------------------------------
smk <- dob$data[["f.20116"]]
# make c.20116 a character variable
smk[, c.20116 := as.character(c.20116)]

# inital assessment
smk_init <- smk[i.20116 == 0, ]

# create ever/never variable
smk[, smk_ev := NA_character_]

# Individuals who responded "Current" or "Previous" will be recoded to Ever
smk2 <- unique(smk[c.20116 %in% c("Current", "Previous"), f.eid])
smk[f.eid %in% smk2, smk_ev := "Ever"]

# Individuals who responded "Never" and have no other smoking status will be recoded to Never
smk3 <- unique(smk[c.20116 %in% c("Never") & !(f.eid %in% smk2), f.eid])
smk[f.eid %in% smk3, smk_ev := "Never"]

# Individuals with missing smk_ev will be recoded to Unknown
smk[is.na(smk_ev), smk_ev := "Unknown"]

# create cleaned smoking data
smk <- merge.data.table(
    smk_init,
    unique(smk[, .(f.eid, smk_ev)]),
)[, .(f.eid, smk_status = c.20116, smk_ev)]

## drinking [f.20117] ----------------------------------------------------------
alc <- dob$data[["f.20117"]]
# make c.20117 a character variable
alc[, c.20117 := as.character(c.20117)]

# inital assessment
alc_init <- alc[i.20117 == 0, ]

# create ever/never variable
alc[, alc_ev := NA_character_]

# Individuals who responded "Current" or "Previous" will be recoded to Ever
alc2 <- unique(alc[c.20117 %in% c("Current", "Previous"), f.eid])
alc[f.eid %in% alc2, alc_ev := "Ever"]

# Individuals who responded "Never" and have no other drinking status will be recoded to Never
alc3 <- unique(alc[c.20117 %in% c("Never") & !(f.eid %in% alc2), f.eid])
alc[f.eid %in% alc3, alc_ev := "Never"]

# Individuals with missing alc_ev will be recoded to Unknown
alc[is.na(alc_ev), alc_ev := "Unknown"]

# create cleaned drinking data
alc <- merge.data.table(
    alc_init,
    unique(alc[, .(f.eid, alc_ev)]),
)[, .(f.eid, alc_status = c.20117, alc_ev)]

## ethnic [f.21000] ------------------------------------------------------------
ethn_coding <- fread("/net/junglebook/home/mmsalva/projects/ukb_data_explorations/data/coding/coding1001.tsv")
ethnicity <- list(
    White = c("White", "British", "Irish", "Any other white background"),
    Other = c(
        "Other ethnic group", "Mixed", "White and Black Caribbean",
        "White and Black African", "White and Asian", "Any other mixed background"
    ),
    Asian = c(
        "Asian or Asian British", "Indian", "Pakistani", "Bangladeshi",
        "Any other Asian background", "Chinese"
    ),
    Black = c("Black or Black British", "Caribbean", "African", "Any other Black background"),
    Unknown = c("Do not know", "Prefer not to answer")
)
ethn <- dob$data$f.21000[, .(id = f.eid, f.21000)]

# recode the data
ethn[, ethnic_group := as.character(factor(f.21000, levels = ethn_coding$coding, labels = ethn_coding$meaning))]

# add column with race_eth
ethn[, race_eth := as.character(NA)]
for (i in 1:5) {
    ethn[ethnic_group %in% ethnicity[[i]], race_eth := names(ethnicity)[i]]
}

# Individuals who responded multiple time with different specific ethnic_group will be recoded to "Other"
ethn2 <- unique(ethn[!race_eth %in% c("Unknown", "Other"), .(id, race_eth)])
dups <- ethn2[duplicated(id), id]
for (i in seq_along(dups)) {
    print(ethn[id == dups[i], race_eth])
    ethn[id == dups[i], `:=`(race_eth, "Other")]
}

# Individuals who responded "Other" and a specific ethnic_group will also be recoded to Other
ethn2 <- unique(ethn[!race_eth %in% c("Unknown"), .(id, race_eth)])
dups <- ethn2[duplicated(id), id]
for (i in seq_along(dups)) {
    print(ethn[id == dups[i], race_eth])
    ethn[id == dups[i], `:=`(race_eth, "Other")]
}

# Remove redundant rows
ethn <- unique(ethn[, .(f.eid = id, race_eth)])

# Individuals who responded "Unknown" but also had a specific answer will be recoded to the specific answer
ethn1 <- ethn[race_eth == "Unknown", ]
ethn2 <- ethn[race_eth != "Unknown", ]
ethn <- rbind(ethn2, ethn1[!f.eid %in% ethn2[, f.eid]])

# bmi [f.21001] ----------------------------------------------------------------
# filter bmi data to only include initial assessment
bmi_init <- dob$data$f.21001[i.21001 == 0, .(f.eid, f.21001)]
bmi_init[, bmi_init_cat := fcase(
    f.21001 < 18.5, "Underweight",
    f.21001 >= 18.5 & f.21001 < 25, "Healthy weight",
    f.21001 >= 25 & f.21001 < 30, "Overweight",
    f.21001 >= 30, "Obese",
    TRUE, NA_character_
    )
]

# calculate median bmi
bmi_med <- dob$data$f.21001[, .(f.eid, f.21001)][, median(f.21001, na.rm = TRUE), by = f.eid]
bmi_med[, bmi_med_cat := fcase(
    V1 < 18.5, "Underweight",
    V1 >= 18.5 & V1 < 25, "Healthy weight",
    V1 >= 25 & V1 < 30, "Overweight",
    V1 >= 30, "Obese",
    TRUE, NA_character_
    )
]

# merge initial assessment and median bmi values
bmi <- merge.data.table(
    bmi_init[, .(f.eid, bmi_init = f.21001, bmi_init_cat)],
    bmi_med[, .(f.eid, bmi_med = V1, bmi_med_cat)],
    by = "f.eid"
    )

# year of birth [f.34] ---------------------------------------------------------
yob <- dob$data$f.34

# make 5-year year of birth bins consistent with Van Alten et al. (2022)
yob[, year_cat := cut(f.34, breaks = seq(1935.5, 1970.5, by = 5), include.lowest = TRUE)]

yob <- yob[, .(f.eid, birth_year = f.34, year_cat)]

# economic status [f.6142] -----------------------------------------------------
econ <- dob$data$f.6142[i.6142 == 0, ]

# “Employed” (1-6, 8), “Retired” (10), “Stay-at-home” (12), “Incapacitated” (13), “Unemployed” (7), and “Student” (9, 11)
econ[, econ_status := as.character(c.6142)]
econ[c.6142 %in% c("Doing unpaid or voluntary work", "Looking after home and/or family"), econ_status := "Stay-at-home"]
econ[c.6142 %in% c("Unable to work because of sickness or disability"), econ_status := "Incapacitated"]
econ[c.6142 %in% c("Full or part-time student"), econ_status := "Student"]
econ[c.6142 %in% c("In paid employment or self-employed"), econ_status := "Employed"]
econ[c.6142 %in% c("None of the above", "Prefer not to answer"), econ_status := NA_character_]

econ <- econ[, .(f.eid, econ_status)]

# tenure of household [f.680] --------------------------------------------------
tenure <- dob$data$f.680[i.680 == 0, ]

tenure[, house_tenure := as.character(c.680)]
tenure[house_tenure %in% c("None of the above", "Prefer not to answer"), house_tenure := NA_character_]

tenure <- tenure[, .(f.eid, house_tenure)]

# vehicle ownership [f.728] ----------------------------------------------------
veh_own <- dob$data$f.728[i.728 == 0, ]
veh_own <- veh_own[, veh_own := as.character(c.728)]
veh_own[c.728 %in% c("Prefers not to answer", "None of the above"), veh_own := NA_character_]

veh_own <- veh_own[, .(f.eid, veh_own)]

# living with [f.709] ----------------------------------------------------------
live_with <- dob$data$f.709[i.709 == 0, ]
live_with[, one_person := as.numeric(f.709 == 1)]
live_with[c.709 %in% c("Prefer not to answer", "None of the above"), one_person := NA_integer_]

live_with <- live_with[, .(f.eid, one_person)]

# self-reported health [f.2178] ------------------------------------------------
srh <- dob$data$f.2178[i.2178 == 0, ]

srh[c.2178 %in% c("Poor"), health := "Bad"]
srh[c.2178 %in% c("Fair"), health := "Fair"]
srh[c.2178 %in% c("Good", "Excellent"), health := "Good/Excellent"]
srh[c.2178 %in% c("Prefer not to answer", "None of the above"), health := NA_character_]

srh <- srh[, .(f.eid, health)]

# education [f.6138] -----------------------------------------------------------
edu <- dob$data$f.6138[i.6138 == 0, ]
edu <- edu[, .(f.eid, edu = c.6138)]

# clean up data ----------------------------------------------------------------
# remove the following datasets from the dob object
multiple <- c(
    "f.20116", "f.20117", "f.21000", "f.21001", "f.34", "f.6142",
    "f.680", "f.728", "f.709", "f.2178", "f.6138"
)
dob$data <- dob$data[!names(dob$data) %in% multiple]

# using dob$info, rename the last columns in the dob$data object
dob$data <- lapply(
    seq_along(dob$data),
    \(x) {
        setnames(
            dob$data[[x]],
            tail(names(dob$data[[x]]), n = 1L),
            make_clean_names(dob$info[`Field ID` == as.numeric(gsub("f.", "", names(dob$data)[x])), Description])
        )
    }
)

# keep first and last column of each dataset in dob$data
dob$data <- lapply(
    dob$data,
    \(x) {
        x[, .SD, .SDcols = c(1L, ncol(x))]
    }
)

# add cleaned data to dob$data
dob$data <- c(
    dob$data,
    list(smk),
    list(alc),
    list(ethn),
    list(bmi),
    list(yob),
    list(econ),
    list(tenure),
    list(veh_own),
    list(live_with),
    list(srh)
)

# combine pulled data ----------------------------------------------------------
dob <- dob$data |>
  purrr::reduce(dplyr::full_join, by = "f.eid")
setnames(dob, c("f.eid", "date_of_consenting_to_join_uk_biobank"), c("id", "consent_date"))
dob[, consent_date := as.Date(consent_date)]

# calculate age at consent
dob[, birth_date := as.Date(paste0(birth_year, "-", month_of_birth, "-15"), "%Y-%B-%d")]
dob[, age_at_consent := round(as.numeric(consent_date - birth_date) / 365.25, 1)]

dob[, id := as.character(id)]
dob[, consent_date := as.Date(consent_date)]

# process data -----------------------------------------------------------------
cli_alert("Calculating DSB for ICD9 data")
icd9 <- merge.data.table(icd9, dob[, .(id, birth_date)], by = "id", all.x = TRUE)
icd9[, dsb := round(as.numeric(((date - birth_date))))]
icd9[, c("date", "birth_date") := NULL][, lexicon := "icd9"]
print(paste(nrow(icd9), "lines of ICD9 codes for", length(unique(icd9[, id])), "people"))

cli_alert("Calculating DSB for ICD10 data")
icd10 <- merge.data.table(icd10, dob[, .(id, birth_date)], by = "id", all.x = TRUE)
icd10[, dsb := round(as.numeric(((date - birth_date))))]
icd10[, c("date", "birth_date") := NULL][, lexicon := "icd10"]
print(paste(nrow(icd10), "lines of ICD10 codes for", length(unique(icd10[, id])), "people"))

## stack icd9 and icd10 data
ukb_icd_dsb <- rbindlist(list(
  icd9, icd10
))

## remove unncessary quotes
ukb_icd_dsb[, diagnosis_code := gsub("\"", "", diagnosis_code)]
print(paste(nrow(ukb_icd_dsb), "lines of ICD codes for", length(unique(ukb_icd_dsb[, id])), "people"))

## save
fwrite(
  x = ukb_icd_dsb,
  file = paste0(out_path, "UKB_PHENOME_ICD_DSB_", opt$ukb_version, ".txt")
)
save_qs(
  x = ukb_icd_dsb,
  file = paste0(out_path, "UKB_PHENOME_ICD_DSB_", opt$ukb_version, ".qs")
)
cli_alert(paste0("File saved as '", out_path, "UKB_PHENOME_ICD_DSB_", opt$ukb_version, ".txt/.qs'!"))

# load necessary information for phecode mapping -------------------------------
# Phecode information
pheinfo <- as.data.table(PheWAS::pheinfo)
pheinfo <- pheinfo[, .(phecode, description, groupnum, group, color)]
phecat_alt_col <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/phecat_alt_colors.txt")
pheinfo <- merge.data.table(
  pheinfo[, !c("color")],
  phecat_alt_col,
  by = "group",
  all.x = TRUE
)


# Sex rules
sex_restriction <- data.table(PheWAS::sex_restriction)
sex_restriction <- rbind(
  data.table("phecode" = sex_restriction[male_only == FALSE & female_only == FALSE, phecode], "sex" = "Both"),
  data.table("phecode" = sex_restriction[male_only == TRUE, phecode], "sex" = "Male"),
  data.table("phecode" = sex_restriction[female_only == TRUE, phecode], "sex" = "Female")
)

# Exclusion / roll up rules
pheinfo2 <- unique(fread("/net/junglebook/home/mmsalva/createUKBphenome/data/phecode_icd9_rolled.csv",
  colClasses = "character",
  select = c("PheCode", "Excl. Phecodes", "Excl. Phenotypes", "Rollup", "Leaf"),
  col.names = c("phecode", "phecode_exclude_range", "phecode_exclude_phenotypes", "rollup", "leaf")
))

# create one data.table with all criteria
pheinfo <- merge(pheinfo, pheinfo2, by = "phecode")
pheinfo <- merge(pheinfo, sex_restriction, by = "phecode")
pheinfo <- pheinfo[, c("phecode", "description", "sex", "rollup", "leaf", "groupnum", "group", "color", "phecode_exclude_range", "phecode_exclude_phenotypes")]

# add phecode information that's missing (collected form earlier versions)
pheinfoOLD <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/Phecode_Definitions_FullTable_Modified.txt", colClasses = "character")
pheinfo <- rbind(pheinfo, pheinfoOLD[!phecode %in% pheinfo$phecode, ])

# Phecode that should not be rolled up
norollup <- pheinfo$phecode[which(pheinfo$rollup == 0)]
# add manual no rollup rules:
norollup <- c(norollup, fread("/net/junglebook/home/mmsalva/createUKBphenome/data/no_rollup_extra.txt", colClasses = "character")$phecode)

# map ICD9 codes to phecodes ---------------------------------------------------
cli_alert("Mapping ICD9 codes to phecodes")

# read PheWAS map (downloaded from https://phewascatalog.org/phecodes; selected "export all" top right corner)
icd9map <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/phecode_icd9_rolled.csv", colClasses = "character")

ICD9codes <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/coding87.tsv")
ICD9codes[!grepl("Block", coding), ICD9 := sapply(coding, harmonizeICD9)]
codeICD9 <- ICD9codes[, sort(unique(ICD9))]

mappedICD9Codes <- NULL
icd9map_new <- list()
cli_progress_bar("mapping ICD9 codes", total = nrow(icd9map))
for (i in seq_len(nrow(icd9map))) {
  cli_progress_update()
  mapped9 <- grep(icd9map$ICD9[i], codeICD9)
  if (length(mapped9) > 0) {
    mappedICD9Codes <- unique(c(mappedICD9Codes, codeICD9[mapped9]))
    icd9map_new[[icd9map$ICD9[i]]] <- data.table(
      "phecode" = icd9map$PheCode[i],
      "ICD9"    = codeICD9[mapped9]
    )
  }
}
cli_progress_done()
icd9key <- rbindlist(icd9map_new)
icd9unmapped <- codeICD9[!codeICD9 %in% mappedICD9Codes]

# roll up PheWAS codes
icd9key$added <- "original"
pcodes <- unique(c(icd9key$phecode, gsub("\\..+", "", icd9key$phecode), gsub("(\\..).", "\\1", icd9key$phecode)))
pcodes <- sort(pcodes)

for (p in seq_along(pcodes)) {
  if (grepl("\\.", pcodes[p])) {
    iSub <- which(icd9key$phecode %in% pcodes[which(grepl(paste0("^", pcodes[p]), pcodes) &
      nchar(pcodes) > nchar(pcodes[p]))] &
      !icd9key$phecode %in% norollup)
  } else {
    iSub <- which(icd9key$phecode %in% pcodes[grep(paste0("^", pcodes[p], "\\."), pcodes)] &
      !icd9key$phecode %in% norollup)
  }
  if (length(iSub) == 0) next
  iTop <- icd9key$ICD9[which(icd9key$phecode == pcodes[p])]
  addTop <- which(icd9key$ICD9 %in% unique(icd9key$ICD9[iSub]) & !icd9key$ICD9 %in% iTop)
  if (length(addTop) == 0) next

  addKey <- icd9key[addTop, ]
  addKey$phecode <- pcodes[p]
  addKey$added <- "rolled up PheWAS code"
  icd9key <- rbind(icd9key, addKey)
}

## Add ICD code description and phecode description
icd9key <- unique(icd9key)
print("Rollup of phewas codes (ICD9 code)")
print(table(icd9key$added))

icd9key <- merge(icd9key, ICD9codes[, .(ICD9, meaning, node_id, parent_id, selectable)], by = "ICD9")
icd9key <- merge(icd9key, pheinfo, by = "phecode")
icd9key <- icd9key[, c(
  "ICD9", "meaning", "node_id", "parent_id", "selectable", "phecode", "description",
  "group", "groupnum", "added", "sex", "rollup", "leaf"
)]

# remove any issues you might have found
icd9map_remove <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/remove_icd9_map.txt", colClasses = "character")
icd9key <- icd9key[!paste(ICD9, phecode, sep = ":") %in% icd9map_remove[, paste(ICD9, phecode, sep = ":")], ]

## map ICD9 data
icd9 <- ukb_icd_dsb[lexicon == "icd9"]
icd9[, diagnosis_code := sapply(diagnosis_code, harmonizeICD9)]
icd9 <- unique(icd9)

phecode1 <- merge(icd9, icd9key[, .(ICD9, phecode)],
  by.x = "diagnosis_code", by.y = "ICD9",
  allow.cartesian = TRUE
)

# map ICD10 codes to phecodes --------------------------------------------------
cli_alert("Mapping ICD10 codes to phecodes")

# read PheWAS map (downloaded from https://phewascatalog.org/phecodes_icd10; selected "export all" top right corner)
icd10map <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/phecode_icd10.csv", colClasses = "character")

# read UKB coding
ICD10codes <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/coding19.tsv")
ICD10codes <- ICD10codes[!grepl("Block", coding), ]
ICD10codes[, ICD10category := gsub("([A-Z][0-9]{2}).+", "\\1", coding)]
ICD10codes[, ICD10suffix := gsub("[A-Z].+", "", gsub("^[A-Z][0-9]{2}", "", coding))]
ICD10codes[, ICD10 := paste0(ICD10category, ifelse(ICD10suffix == "", "", "."), ICD10suffix)]
ICD10codes <- ICD10codes[, c("ICD10category", "ICD10suffix") := NULL]

codeICD10 <- ICD10codes[, sort(unique(ICD10))]
mappedICD10Codes <- NULL
icd10map_new <- list()
cli_progress_bar("mapping ICD10 codes", total = nrow(icd10map))
for (i in seq_len(nrow(icd10map))) {
  cli_progress_update()
  mapped10 <- grep(icd10map$ICD10[i], codeICD10)
  if (length(mapped10) > 0) {
    mappedICD10Codes <- unique(c(mappedICD10Codes, codeICD10[mapped10]))
    icd10map_new[[icd10map$ICD10[i]]] <- data.table(
      "phecode" = icd10map$PheCode[i],
      "ICD10" = codeICD10[mapped10]
    )
  }
}
cli_progress_done()
icd10key <- rbindlist(icd10map_new)
icd10unmapped <- codeICD10[!codeICD10 %in% mappedICD10Codes]

#### roll up phewas codes
icd10key$added <- "original"
pcodes <- unique(c(icd10key$phecode, gsub("\\..+", "", icd10key$phecode), gsub("(\\..).", "\\1", icd10key$phecode)))
pcodes <- sort(pcodes)
for (p in seq_along(pcodes)) {
  if (grepl("\\.", pcodes[p])) {
    iSub <- which(icd10key$phecode %in% pcodes[which(grepl(paste0("^", pcodes[p]), pcodes) &
      nchar(pcodes) > nchar(pcodes[p]))] &
      !icd10key$phecode %in% norollup)
  } else {
    iSub <- which(icd10key$phecode %in% pcodes[grep(paste0("^", pcodes[p], "\\."), pcodes)] &
      !icd10key$phecode %in% norollup)
  }
  if (length(iSub) == 0) next
  iTop <- icd10key$ICD10[which(icd10key$phecode == pcodes[p])]
  addTop <- which(icd10key$ICD10 %in% unique(icd10key$ICD10[iSub]) & !icd10key$ICD10 %in% iTop)
  if (length(addTop) == 0) next

  addKey <- icd10key[addTop, ]
  addKey$phecode <- pcodes[p]
  addKey$added <- "rolled up PheWAS code"
  icd10key <- rbind(icd10key, addKey)
}

icd10key <- unique(icd10key)

print("Rollup of phewas codes (ICD10 code)")
print(table(icd10key$added))

icd10key <- merge(icd10key, ICD10codes[, .(ICD10, meaning, node_id, parent_id, selectable)], by = "ICD10")
icd10key <- merge(icd10key, pheinfo, by = "phecode")
icd10key <- icd10key[, c(
  "ICD10", "meaning", "node_id", "parent_id", "selectable", "phecode", "description",
  "group", "groupnum", "added", "sex", "rollup", "leaf"
)]


## remove any issues you might have found
icd10map_remove <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/remove_icd10_map.txt", colClasses = "character")
icd10key <- icd10key[!paste(ICD10, phecode, sep = ":") %in% icd10map_remove[, paste(ICD10, phecode, sep = ":")], ]

## map ICD10 data
icd10 <- ukb_icd_dsb[lexicon == "icd10", ]
icd10[, ICD10category := gsub("([A-Z][0-9]{2}).+", "\\1", diagnosis_code)]
icd10[, ICD10suffix := gsub("[A-Z].+", "", gsub("^[A-Z][0-9]{2}", "", diagnosis_code))]
icd10[, diagnosis_code := paste0(ICD10category, ifelse(ICD10suffix == "", "", "."), ICD10suffix)]
icd10[, c("ICD10category", "ICD10suffix") := NULL]
icd10 <- unique(icd10)

phecode2 <- merge(icd10, icd10key[, .(ICD10, phecode)],
  by.x = "diagnosis_code", by.y = "ICD10",
  allow.cartesian = TRUE
)

## stack icd-phecode mapped dsb data
ukb_phecode <- unique(rbindlist(list(
  phecode1,
  phecode2
)))

# save mapped data -------------------------------------------------------------
fwrite(
  x = ukb_phecode,
  file = paste0(out_path, "UKB_PHECODE_DSB_MAPPED_", opt$ukb_version, ".txt")
)
save_qs(
  x = ukb_phecode,
  file = paste0(out_path, "UKB_PHECODE_DSB_MAPPED_", opt$ukb_version, ".qs")
)
cli_alert_info(paste0("File saved as '", out_path, "UKB_PHECODE_DSB_MAPPED_", opt$ukb_version, ".txt/.qs'!"))

# save sex data ----------------------------------------------------------------
first_phe <- ukb_phecode[ukb_phecode[, .I[which.min(dsb)], id]$V1]
last_phe  <- ukb_phecode[ukb_phecode[, .I[which.max(dsb)], id]$V1]
ehr_followup <- merge.data.table(
  first_phe[, .(id, first_dsb = dsb, age_at_first_diagnosis = round(dsb / 365.25, 1))],
  last_phe[, .(id, last_dsb = dsb, age_at_last_diagnosis = round(dsb / 365.25, 1))],
  by = "id"
)[, ehr_days := as.numeric(last_dsb - first_dsb)][, ehr_years := round(ehr_days / 365.25, 1)]
dob <- merge.data.table(
  dob,
  ehr_followup,
  by = "id",
  all.x = TRUE
)

## add comorbidity indicators
cancer_phecodes <- fread("/net/junglebook/home/mmsalva/projects/dissertation/aim_two/data/public/cancer_phecodes.txt",
                         colClasses = "character"
)[[1]]

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

cli_alert("identifying cases and creating indicator variables...")
for (i in names(comorbid)) {
  comorbid[[i]][["ids"]] <- ukb_phecode[phecode %in% comorbid[[i]][["phecodes"]], id] |>
    unique()
}

for (i in names(comorbid)) {
  set(dob, j = i, value = fifelse(dob[["id"]] %in% comorbid[[i]][["ids"]], 1, 0))
}

dob[, triglycerides := fifelse(hypertension == 0 & mixed_hypertension == 0, 0, 1)]

# keep individuals whose sex == genetic_sex as in_phenome
dob[, `:=` (
    in_phenome  = fifelse(id %in% ukb_phecode[, unique(id)], 1, 0),
    in_phenome2 = fifelse(id %in% ukb_phecode[, unique(id)] | genetic_sex == sex, 1, 0)
  )]

comorbids <- c(names(comorbid), "triglycerides")
for (i in comorbids) {
  dob[in_phenome == 0, (i) := NA_real_]
}

dob[, `:=` (
  smoker = fcase(
    smk_ev == "Ever", 1,
    smk_ev == "Never", 0,
    default = NA
  ),
  drinker = fcase(
    alc_ev == "Ever", 1,
    alc_ev == "Never", 0,
    default = NA
  )
)]

## save
fwrite(
  x = dob,
  file = paste0(out_path, "UKB_PHENOME_COV_", opt$ukb_version, ".txt")
)
save_qs(
  x = dob,
  file = paste0(out_path, "UKB_PHENOME_COV_", opt$ukb_version, ".qs")
)
cli_alert_info(paste0("File saved as '", out_path, "UKB_PHENOME_COV_", opt$ukv_version, ".txt/.qs'!"))

# save phecode indicator matrix ------------------------------------------------
pim <- dcast(
    ukb_phecode[, .(id, phecode = paste0("X", phecode))],
    id ~ phecode,
    value.var     = "phecode",
    fun.aggregate = length,
    fill          = 0
  )

# collect integer variable names from data.table
int_vars <- names(pim)[which(sapply(pim, is.integer))]
pim[, (int_vars) := lapply(.SD, \(x) as.numeric(x > 0)), .SDcols = int_vars]

save_qs(
  x = pim,
  file = paste0(out_path, "UKB_PHENOME_PIM0_", opt$ukb_version, ".qs")
)

pim2 <- rbindlist(
  list(
    pim,
    dob[!(id %in% pim[, id]) & in_phenome2 == 1, .(id)]
  ),
  fill = TRUE,
  use.names = TRUE
)
pim2[is.na(pim2)] <- 0 
save_qs(
  x = pim2,
  file = paste0(out_path, "UKB_PHENOME_PIM_w_no_codes_", opt$ukb_version, ".qs")
)

cli_alert_success("Done! 🎉")
