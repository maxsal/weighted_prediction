suppressPackageStartupMessages({
  require(data.table)
  require(glue)
})
generate_restricted_phenome <- function(phe_data, threshold, cases, outcome_phe) {
  
  if (is.na(suppressWarnings({as.numeric(substr(outcome_phe, 1, 1))}))) {
    outcome_phe <- glue::glue("X{outcome_phe}")
  }
  
  out <- data.table::dcast(
    unique(phe_data[get(glue::glue("t{threshold}_indicator")) == 1, ][dsb < get(glue::glue("t{threshold}_threshold")), ][, .(id, phecode = fifelse(
      is.na(suppressWarnings({
        as.numeric(substr(phecode, 1, 1))
      })), phecode, paste0("X", phecode)
    ))]),
    id ~ phecode,
    value.var     = "phecode",
    fun.aggregate = length,
    fill          = 0
  )
  
  out[, case := data.table::fifelse(id %in% cases, 1, 0)]
  if (outcome_phe %in% names(out)) {
    out[, c(outcome_phe) := NULL]
  }
  
  return(out)
  
}