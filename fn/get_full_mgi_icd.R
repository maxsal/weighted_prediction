# function to get icd-phecode 9 and 10 data on all individuals in MGI
get_full_mgi_icd <- function(mgi_version, icd9_file, icd10_file, force = FALSE) {
  
  # first, check if file already exists
  if( file.exists( paste0("data/", data_version, "/processed/full_icd.fst") ) & force == FALSE ){
    cli::cli_alert_info("congrats! full_icd already exists")
    out <- fst::read_fst(paste0("data/", data_version, "/processed/full_icd.fst"), as.data.table = TRUE)
  } else {
    icd9 <- get(load(icd9_file))
    icd10 <- get(load(icd10_file))
    out <- rbindlist(list(
      icd9, icd10
    ))
    fst::write_fst(x = out, path = paste0("data/", data_version, "/processed/full_icd.fst"), compress = 100)
  }
  
  return(out)
  
}
