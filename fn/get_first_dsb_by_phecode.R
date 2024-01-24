get_first_dsb_by_phecode <- function(
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
    icd <- rbindlist(list(icd, more_icd))
  }
  
  icd <- icd[, .(id = IID, dsb = DaysSinceBirth, phecode)][, dsb := as.numeric(dsb)]
  
  out <- out[ out[, .I[which.min(dsb)], by = c("id", "phecode")]$V1 ]
  return(out)
  
}