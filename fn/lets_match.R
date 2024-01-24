# function for matching cases to noncases
lets_match <- function(dataset, nearest_vars, exact_vars, calip = 0.25, r = matching_ratio) {
  
  nearest_vars = c("age_at_first_diagnosis", "length_followup")
  exact_vars = c("female")
  calip = 0.25
  r = matching_ratio
  dataset = deparse(substitute(matching_cov))
  match_text <- paste0("MatchIt::matchit(case ~ ",
                       paste0(c(nearest_vars, exact_vars), collapse = " + "),
                       ", data = ", dataset,
                       ", calclosest = TRUE, mahvars = c(", paste0(purrr::map_chr(nearest_vars, ~paste0("\'", .x,"\'")), collapse = ", "),
                       "), caliper = ", calip,
                       ", exact = c(", paste0(purrr::map_chr(exact_vars, ~paste0("\'", .x,"\'")), collapse = ", "), "), ratio = ", r, ")")
  
  return(eval(parse(text = match_text)))
  
}
