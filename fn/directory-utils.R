## function that checks for presence of folders corresponding to data date
## and, optionally, outcome phecode and makes missing folders if absent
suppressPackageStartupMessages({
  require(glue)
})

is_character_numeric <- function(char) {
  # Try to convert the character to numeric
  num <- suppressWarnings({as.numeric(substr(char, 1, 1))})

  # Check if the conversion is successful and not NA
  if (!is.na(num)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

check_folder_structure <- function(cohort = "mgi", data_version, outcome_phecode) {
  created <- 0

  if (is_character_numeric(outcome_phecode)) {
    outcome_phecode <- paste0("X", outcome_phecode)
  }

  # data -----------------------------------------------------------------------
  if (!dir.exists(glue("data/private/{cohort}/{data_version}/{outcome_phecode}"))) {
    dir.create(glue("data/private/{cohort}/{data_version}/{outcome_phecode}"), recursive = TRUE)
    message(glue("data/private/{cohort}/{data_version}/{outcome_phecode}/ folder created"))
    created <- created + 1
  }

  if (!dir.exists(glue::glue("results/{cohort}/{data_version}/"))) {
    dir.create(glue::glue("results/{cohort}/{data_version}/"), recursive = TRUE)
    message(glue("results/{cohort}/{data_version}/ folder created"))
    created <- created + 1
  }

  # results --------------------------------------------------------------------
  if (!is.null(outcome_phecode)) {
    if (!dir.exists(glue("results/{cohort}/{data_version}/{outcome_phecode}/"))) {
      dir.create(glue("results/{cohort}/{data_version}/{outcome_phecode}/"), recursive = TRUE)
      message(glue("results/{cohort}/{data_version}/{outcome_phecode}/ folder created"))
      created <- created + 1
    }
  }

  # summary --------------------------------------------------------------------
  if (created == 0) {
    if (is.null(outcome_phecode)) {
      message(glue("file structure already exists for {cohort} version {data_version}"))
    } else {
      message(glue("file structure already exists for {cohort} version {data_version} and outcome phecode {outcome_phecode}"))
    }
  }
}
