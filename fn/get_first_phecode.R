# function to extract the first occurrence of each phecode by individual
get_first_phecode <- function(cohort, icd_data = NULL, mgi_version = "20210318", ukb_file, icd9_file, icd10_file, force = FALSE) {
  
  if (tolower(cohort) == "mgi") {
    if (!is.null(icd_data)) {
      out <- icd_data[, .(id = IID, dsb = DaysSinceBirth, phecode)]
      out <- out[ out[, .I[which.min(dsb)], by = c("id", "phecode")]$V1 ]
    } else {
      # first, check if file already exists
      if( file.exists( paste0("data/", mgi_version, "/processed/first_phecode.txt") ) ){
        cli::cli_alert_info("congrats! first_phecode already exists")
        out <- data.table::fread(paste0("data/", mgi_version, "/processed/first_phecode.txt"), colClasses = "character")
      } else {
        icd9 <- get(load(icd9_file))
        icd10 <- get(load(icd10_file))
        out <- rbindlist(list(
          icd9, icd10
        ))[, .(id = IID, dsb = DaysSinceBirth, phecode)]
        out <- out[ out[, .I[which.min(dsb)], by = c("id", "phecode")]$V1 ]
        fwrite(out, paste0("data/", mgi_version, "/processed/first_phecode.txt"))
      }
    }
  }
  
  if (tolower(cohort) == "ukb") {
    out <- data.table::fread(ukb_file, na.strings = c(""), colClasses = "character")
  }
  
  return(out)
  
}

