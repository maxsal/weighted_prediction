get_full_icd_dsb_phecode <- function(
    icd_file,
    more_icd_file = NULL
) {
  
  if (tolower(tools::file_ext(icd_file)) %in% c("txt", "csv", "tsv")) {
    icd <- data.table::fread(icd_file, colClasses = "character")
  }
  if (tolower(tools::file_ext(icd_file)) == "rsav") {
    icd <- get(load(icd_file))
  }
  
  if (!is.null(more_icd_file)) {
    if (tolower(tools::file_ext(more_icd_file)) %in% c("txt", "csv", "tsv")) {
      more_icd <- data.table::fread(more_icd_file, colClasses = "character")
    }
    if (tolower(tools::file_ext(more_icd_file)) == "rsav") {
      more_icd <- get(load(more_icd_file))
    }
    icd <- data.table::rbindlist(list(icd, more_icd))
  }
  
  if ("IID" %in% names(icd)) {data.table::setnames(icd, "IID", "id")}
  if ("DaysSinceBirth" %in% names(icd)) {data.table::setnames(icd, "DaysSinceBirth", "dsb")}
  icd <- icd[, .(id, dsb, phecode)][, dsb := as.numeric(dsb)]
  
  return(icd)
  
}
