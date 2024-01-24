get_icd_phecode_metrics <- function(full_phe_data) {
  
  first_phe    <- full_phe_data[ full_phe_data[, .I[which.min(dsb)], by = "id"]$V1 ][, .(id, first_dsb = dsb)]
  last_phe     <- full_phe_data[ full_phe_data[, .I[which.max(dsb)], by = "id"]$V1 ][, .(id, last_dsb = dsb)]
  n_encounters <- unique(full_phe_data[, .(id, dsb)])[, .N, by = "id"][, .(id, n_encounters = N)]
  
  out <- Reduce(merge.data.table, list(first_phe, last_phe, n_encounters))
  out[, `:=` (
    age_at_first_diagnosis = round(first_dsb / 365.25, 1),
    age_at_last_diagnosis  = round(last_dsb / 365.25, 1)
  )][, `:=` (
    length_followup     = round((age_at_last_diagnosis - age_at_first_diagnosis)),
    encounters_per_year = round(n_encounters / (age_at_last_diagnosis - age_at_first_diagnosis), 1)
  )]
  out[is.infinite(encounters_per_year), encounters_per_year := NA]
  
  return(out[])
  
}