# suppressPackageStartupMessages({
#   library(simplexreg)
#   library(data.table)
#   library(survey)
#   library(ms)
#   library(cli)
#   library(scales)
# })
# # calculate IPW from stacked data
# # ADAPTED FROM: /net/junglebook/home/kundur/EHR/Processed Code/Weighted_using_lauren_code_bb.R
# ipw <- function(
#     stacked_data,
#     weight_outcome_var = "WTFA_A",
#     samp_var = "samp_WTFA_A",
#     external_dataset = "NHIS",
#     dataset_name = "MGI",
#     id_var = "id",
#     cancer_factor = FALSE,
#     cancer_factor_var = "cancer",
#     covs = c(
#       "age_50", "female", "nhw", "hypertension",
#       "diabetes", "cancer", "anxiety", "depression",
#       "bmi_cat"
#     ),
#     chop = TRUE) {
#   stacked_data[dataset == external_dataset, (weight_outcome_var) := .N * get(weight_outcome_var) /
#     sum(get(weight_outcome_var), na.rm = TRUE)]

#   if (cancer_factor) {
#     if (cancer_factor_var %in% covs) {
#       cancer_TF <- TRUE
#       covs <- covs[covs != cancer_factor_var]
#     } else {
#       cancer_TF <- FALSE
#     }
#   }
#   select_mod_covs <- paste0(covs, collapse = " + ")

#   # modeling nhanes sampling weights
#   nhanes_select_mod <- simplexreg(
#     as.formula(paste0(
#       samp_var, " ~ ",
#       select_mod_covs
#     )),
#     data = stacked_data[dataset == external_dataset, ]
#   )


#   # selection model into internal data
#   internal_select_mod <- glm(
#     paste0(
#       "as.numeric(dataset == '",
#       dataset_name, "') ~ ",
#       select_mod_covs
#     ),
#     data = stacked_data, family = quasibinomial()
#   )

#   # obtain fitted values from nhanes and internal models
#   p_nhanes <- predict(nhanes_select_mod,
#     newdata = stacked_data[dataset == dataset_name, ],
#     type = "response"
#   )[, 1]
#   p_internal <- predict(internal_select_mod,
#     newdata = stacked_data[dataset == dataset_name, ],
#     type = "response"
#   )

#   ###
#   temp <- rep(0, times = length(p_internal))
#   temp[which(rownames(data.frame(p_nhanes)) %in%
#     rownames(data.frame(p_internal)) == T)] <- p_nhanes
#   temp[which(rownames(data.frame(p_nhanes)) %in%
#     rownames(data.frame(p_internal)) == F)] <- NA
#   p_nhanes <- temp
#   p_nhanes[which(p_nhanes == 0)] <- min(p_nhanes[which(p_nhanes > 0)], na.rm = TRUE)
#   nhanes_selection <- p_nhanes * (p_internal / (1 - p_internal))
#   ###

#   if (chop) nhanes_selection <- chopr(nhanes_selection)
#   nhanes_weight <- 1 / nhanes_selection
#   nhanes_weight <- stacked_data[dataset == dataset_name, .N] *
#     nhanes_weight /
#     sum(nhanes_weight, na.rm = TRUE)

#   ## With Cancer
#   if (cancer_factor) {
#     if (cancer_TF == TRUE) {
#       nhanes_cancer_mod <- glm(as.formula(paste0(cancer_factor_var, " ~ ", select_mod_covs)),
#         data = stacked_data[dataset == external_dataset, ],
#         weights = get(weight_outcome_var), family = quasibinomial()
#       )

#       nhanes_cancer <- predict(nhanes_cancer_mod,
#         type = "response",
#         newdata = stacked_data
#       )

#       mgi_cancer_mod <- glm(as.formula(paste0(cancer_factor_var, " ~ ", select_mod_covs)),
#         data = stacked_data[dataset == dataset_name, ],
#         family = quasibinomial()
#       )

#       mgi_cancer <- predict(mgi_cancer_mod,
#         type = "response",
#         newdata = stacked_data
#       )
#       denom <- ifelse(stacked_data[[cancer_factor_var]] == 1, mgi_cancer,
#         1 - mgi_cancer
#       )
#       num <- ifelse(stacked_data[[cancer_factor_var]] == 1, nhanes_cancer,
#         1 - nhanes_cancer
#       )
#       cancer_factor <- (num[stacked_data[, dataset] == dataset_name] /
#         denom[stacked_data[, dataset] == dataset_name])
#       if (chop) {
#         nhanes_weight <- chopr(cancer_factor * nhanes_weight)
#       } else {
#         nhanes_weight <- cancer_factor * nhanes_weight
#       }
#       nhanes_weight <- stacked_data[dataset == dataset_name, .N] *
#         nhanes_weight /
#         sum(nhanes_weight, na.rm = TRUE)
#     }
#   }

#   data.table(
#     "id"        = stacked_data[dataset == dataset_name, ][[id_var]],
#     "ip_weight" = nhanes_weight
#   )
# }

# # poststrat function -----------------------------------------------------------
# psw <- function(
#     int_data,
#     ext_data = NULL,
#     id_var = "id",
#     age_var = "AgeLastEntry",
#     age_bin = FALSE,
#     age_bin_breaks = c(0, 18, 35, 65, 80, 150),
#     covs = c("age_bin", "cad", "smoker", "diabetes", "female"),
#     not_num_vars = NULL,
#     psu_var = "psu_nhanes",
#     strata_var = "strata_nhanes",
#     weight_var = "weight_nhanes",
#     chop = TRUE) {
#   if (age_bin) {
#     cli_alert("constructing categorical age variable using breaks: {age_bin_breaks}")
#     covs <- unique(c(covs, "age_bin"))
#   }
#   cli_alert_info("estimating poststratification weights for covariates: {.field {covs}}")
#   if (chop) cli_alert_info("truncating weights to 2.5 and 97.5 percentiles")

#   # 1. load and prep nhanes data ----------------------------------------------
#   if (is.null(ext_data)) {
#     nhanes <- download_nhanes_data(datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ", "BPQ", "DPQ"))
#     phanes <- prepare_nhanes_data(nhanes)
#   } else {
#     phanes <- copy(ext_data)
#   }

#   if (age_bin == TRUE) {
#     phanes[, age_bin := cut(age, age_bin_breaks, right = FALSE)]
#   }

#   nhanes_design <- svydesign(
#     id      = ~ get(psu_var),
#     strata  = ~ get(strata_var),
#     weights = ~ get(weight_var),
#     nest    = TRUE,
#     data    = phanes
#   )

#   # 2. get nhanes proportions -------------------------------------------------
#   population_proportions <- svytable(
#     formula = as.formula(paste0("~", paste(covs, collapse = " + "))),
#     design  = nhanes_design
#   ) |>
#     prop.table() |>
#     as.data.table()
#   setnames(population_proportions, "N", "pop_prob")

#   if (!is.null(not_num_vars)) {
#     num_vars <- setdiff(covs, c(not_num_vars, "age_bin"))
#   } else {
#     num_vars <- setdiff(covs, "age_bin")
#   }

#   population_proportions[, (num_vars) := lapply(.SD, as.numeric), .SDcols = num_vars]

#   # 3. get internal proportions -----------------------------------------------
#   if (age_bin == TRUE) {
#     int_data[, age_bin := cut(get(age_var), age_bin_breaks, right = FALSE)]
#   }

#   internal_probabilities <- int_data[, ..covs] |>
#     table() |>
#     prop.table() |>
#     as.data.table()
#   setnames(internal_probabilities, "N", "int_prob")
#   internal_probabilities[, (num_vars) := lapply(.SD, as.numeric), .SDcols = num_vars]

#   # 4. merge probabilities into internal data ---------------------------------
#   sub_vars <- c(id_var, covs)

#   min_pop_prob <- min(population_proportions[, pop_prob][which(population_proportions[, pop_prob] > 0)])
#   population_proportions[pop_prob == 0, pop_prob := min_pop_prob]

#   merged <- int_data[, ..sub_vars] |>
#     merge.data.table(population_proportions, by = covs) |>
#     merge.data.table(internal_probabilities, by = covs)
#   merged[, ps_weight := pop_prob / int_prob]

#   # 5. process ----------------------------------------------------------------
#   if (chop == TRUE) merged[, ps_weight := chopr(ps_weight)]

#   merged[, ps_weight := (.N * ps_weight) / sum(ps_weight, na.rm = TRUE)]

#   return_vars <- c(id_var, "ps_weight")
#   return(
#     merged[, ..return_vars]
#   )
# }


# ###############################
# ### OLD: FOR REFERENCE ONLY ###
# ###############################
# # # calculate poststratification weights from MGI
# # # ADAPTED FROM: /net/junglebook/home/kundur/EHR/Processed Code/Weighted_using_lauren_code_bb.R
# # poststratification <- function(
# #     mgi_data,
# #     id_var             = "id",
# #     last_entry_age_var = "AgeLastEntry",
# #     cancer_var         = "cancer",
# #     chd_var            = "cad",
# #     smoke_var          = "smoking_current",
# #     diabetes_var       = "diabetes",
# #     depression_var     = "depression",
# #     female_var         = "female",
# #     hypertension_var   = "hypertension",
# #     nhw_var            = "nhw",
# #     covs               = c("cad", "smoke", "diabetes"),
# #     chop               = FALSE
# # ) {
  
# #   if (!all(covs %in% c("cad", "smoke", "diabetes", "female", "depression", "hypertension", "cancer", "nhw"))) {
# #     stop("function only supports covs 'cad', 'smoke', 'diabetes', 'female', 'depression', 'hypertension', 'cancer', and 'nhw' right now")
# #   }
  
# #   Nobs    <- nrow(mgi_data)
# #   age_vec <- as.numeric(mgi_data[[last_entry_age_var]])
# #   which_between <- function(vec, mat) {
# #     which(data.table::between(x = vec, lower = mat[["lower"]],
# #           upper = mat[["upper"]]))
# #   }
  
# #   # cancer prevalence by age (US) 
# #   # https://seer.cancer.gov/csr/1975_2016/results_merged/topic_prevalence.pdf
# #   cancer_prevalence <- age_grp_table(
# #     lower_ages   = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
# #     num_vec      = c(0.0899, 0.2023, 0.3922, 0.8989, 2.1532, 4.9326, 10.4420,
# #                       18.3168, 21.5939) / 100,
# #     num_var_name = "prevalence"
# #   )
  
# #   cancer_func_pop <- stepfun(x = cancer_prevalence[["lower"]],
# #                              y = c(0, cancer_prevalence[["prevalence"]]),
# #                              right = FALSE)
# #   b               <- aggregate(as.formula(paste0(cancer_var, " ~ ",
# #                                           last_entry_age_var)),
# #                                FUN = mean,
# #                                data = mgi_data)
# #   cancer_func_mgi <- stepfun(x = b[, 1], y = c(0, b[, 2]), right = FALSE)
  
# #   # diabetes prevalence by age (US)
# #   # https://www.cdc.gov/diabetes/pdfs/data/statistics/national-diabetes-statistics-report.pdf 
# #   # youngest group will have no people, no number provided in documentation
# #   if ("diabetes" %in% covs) {
# #     diabetes_prevalence <- age_grp_table(
# #       lower_ages   = c(0, 18, 45, 65),
# #       num_vec      = c(0.01, 0.030, 0.138, 0.214),
# #       num_var_name = "prevalence"
# #     )
    
# #     diabetes_func_pop <- stepfun(x = diabetes_prevalence[["lower"]],
# #     y = c(0, diabetes_prevalence[["prevalence"]]), right = FALSE)
# #     d_b               <- aggregate(as.formula(paste0(diabetes_var, " ~ ",
# #                                                       last_entry_age_var)),
# #                                    FUN = mean,
# #                                    data = mgi_data)
# #     diabetes_func_mgi <- stepfun(x = d_b[, 1], y = c(0, d_b[, 2]), right = FALSE)
# #   }
  
# #   # depression prevalence by age (us)
# #   # https://ourworldindata.org/grapher/prevalence-of-depression-by-age?country=~USA
# #   if ("depression" %in% covs) {
# #     depression_prevalence <- age_grp_table(
# #       lower_ages   = c(0, 10, 15, 50, 70),
# #       num_vec      = c(0.01, 0.0272, 0.0657, 0.0469, 0.0381),
# #       num_var_name = "prevalence"
# #     )
    
# #     depression_func_pop <- stepfun(x = depression_prevalence[["lower"]],
# #                                  y = c(0, depression_prevalence[["prevalence"]]), right = FALSE)
# #     dep_b                 <- aggregate(as.formula(paste0(depression_var, " ~ ",
# #                                                      last_entry_age_var)),
# #                                    FUN = mean,
# #                                    data = mgi_data)
# #     depression_func_mgi <- stepfun(x = dep_b[, 1], y = c(0, dep_b[, 2]), right = FALSE)
# #   }
  
# #   # chd prevalence by age (us)
# #   # https://www.cdc.gov/nchs/fastats/heart-disease.htm
# #   # youngest group will have no people, no number provided in documentation
# #   if ("cad" %in% covs) {
# #     chd_prevalence <- age_grp_table(
# #       lower_ages   = c(0, 18, 45, 65, 75),
# #       num_vec      = c(0.01, 0.01, 0.060, 0.155, 0.239),
# #       num_var_name = "prevalence"
# #     )
    
# #     chd_func_pop <- stepfun(x = chd_prevalence[["lower"]],
# #                             y = c(0, chd_prevalence[["prevalence"]]),
# #                             right = FALSE)
# #     chd_b        <- aggregate(as.formula(paste0(chd_var, " ~ ",
# #                                          last_entry_age_var)),
# #                               FUN = mean,
# #                               data = mgi_data)
# #     chd_func_mgi <- stepfun(x = chd_b[, 1], y = c(0, chd_b[, 2]), right = FALSE)
# #   }
  
# #   # smoking prevalence by age (us)
# #   # https://www.cdc.gov/tobacco/data_statistics/fact_sheets/adult_data/cig_smoking/index.htm
# #   # youngest group will have no people, no number provided in documentation
# #   if ("smoke" %in% covs) {
# #     smoke_prevalence <- age_grp_table(
# #       lower_ages   = c(0, 18, 25, 45, 65),
# #       num_vec      = c(0, 0.053, 0.126, 0.149, 0.083),
# #       num_var_name = "prevalence"
# #     )
    
# #     smoke_func_pop <- stepfun(x = smoke_prevalence[["lower"]],
# #                               y = c(0, smoke_prevalence[["prevalence"]]),
# #                               right = FALSE)
# #     smoke_b        <- aggregate(as.formula(paste0(smoke_var, " ~ ",
# #                                            last_entry_age_var)),
# #                                 FUN = mean, data = mgi_data)
# #     smoke_func_mgi <- stepfun(x = smoke_b[, 1], y = c(0, smoke_b[, 2]),
# #                               right = FALSE)
# #   }
  
# #   # age distribution (us)
# #   # https://www.census.gov/data/tables/2000/dec/phc-t-09.html
# #   total_population  <- 281421906
# #   male_population   <- 138053563
# #   female_population <- total_population - male_population
  
# #   low_ages <- seq(0, 85, 5)
# #   age_counts_male <- age_grp_table(
# #     lower_ages   = low_ages,
# #     num_vec      = c(9810733, 10523277, 10520197, 10391004, 9687814, 9798760,
# #                      10321769, 11318696, 11129102, 9889506, 8607724, 6508729,
# #                      5136627, 4400362, 3902912, 3044456, 1834897, 1226998),
# #     num_var_name = "counts"
# #   )
# #   age_counts_female <- age_grp_table(
# #     lower_ages   = low_ages,
# #     num_vec      = c(9365065, 10026228, 10007875, 9828886, 9276187, 9582576,
# #                      10188619, 11387968, 11312761, 10202898, 8977824, 6960508,
# #                      5668820, 5133183, 4954529, 4371357, 3110470, 3012589),
# #     num_var_name = "counts"
# #   )
# #   age_prevalence <- age_grp_table(
# #     lower_ages   = low_ages,
# #     upper_offset = 0,
# #     num_vec      = (age_counts_male[["counts"]] +
# #                       age_counts_female[["counts"]]) /
# #                       total_population,
# #     num_var_name = "prevalence"
# #   )
  
# #   age_vals      <- age_prevalence[, group]
# #   mgi_age_group <- age_vals[as.vector(unlist(apply(as.matrix(age_vec), 1,
# #                             FUN = which_between, mat = age_prevalence)))]
# #   age_func_pop  <- stepfun(x = age_prevalence[["lower"]],
# #                            y = c(0, age_prevalence[["prevalence"]]),
# #                            right = FALSE)
# #   age_func_mgi  <- stepfun(x = age_prevalence[["lower"]],
# #                            y = as.numeric(c(0, table(cut(age_vec,
# #                             breaks = c(0, age_prevalence[["upper"]]),
# #                             labels = age_prevalence[["group"]],
# #                             right = FALSE), useNA = "ifany") / Nobs)),
# #                            right = FALSE)
# #   age_func_mgi2 <- stepfun(x = age_prevalence[["lower"]],
# #                            y = c(0, as.numeric(table(factor(mgi_age_group,
# #                                  levels = age_vals))) / Nobs),
# #                            right = FALSE)
  
# #   # sex distribution by age (us)
# #   # same source as above for age
# #   if ("female" %in% covs) {
# #     female_prevalence <- age_grp_table(
# #       lower_ages   = low_ages,
# #       upper_offset = 0,
# #       num_vec      = age_counts_female[["counts"]] /
# #                        (age_counts_male[["counts"]] +
# #                        age_counts_female[["counts"]]),
# #       num_var_name = "prevalence"
# #     )
# #     female_func_pop <- stepfun(x = female_prevalence[["lower"]],
# #                                y = c(0, female_prevalence[["prevalence"]]),
# #                                right = FALSE)
# #     female_b        <- aggregate(as.formula(paste0(female_var, " ~ ",
# #                                             last_entry_age_var)),
# #                                  FUN = mean, data = mgi_data)
# #     female_func_mgi <- stepfun(x = female_b[, 1], y = c(0, female_b[, 2]),
# #                                right = FALSE)
# #   }
  
# #   # hypertension distribution by age (us)
# #   # https://www.cdc.gov/nchs/data/databriefs/db289.pdf
# #   if ("hypertension" %in% covs) {
# #     hypertension_prevalence <- age_grp_table(
# #       lower_ages   = c(0, 18, 40, 60),
# #       upper_offset = 0,
# #       num_vec      = c(0.01, 0.075, 0.332, 0.631),
# #       num_var_name = "prevalence"
# #     )
# #     hypertension_func_pop <- stepfun(x = hypertension_prevalence[["lower"]],
# #                                y = c(0, hypertension_prevalence[["prevalence"]]),
# #                                right = FALSE)
# #     hypertension_b        <- aggregate(as.formula(paste0(hypertension_var, " ~ ",
# #                                                    last_entry_age_var)),
# #                                  FUN = mean, data = mgi_data)
# #     hypertension_func_mgi <- stepfun(x = hypertension_b[, 1], y = c(0, hypertension_b[, 2]),
# #                                right = FALSE)
# #   }
  
# #   # non-hispanic white distribution by age (us)
# #   # https://data.census.gov/table?q=age,+race,+ethnicity&tid=ACSDT5Y2021.B01001H
# #   # calculated from 2021 5-year ACS estimates (tables B01001H and S0101)
# #   if ("nhw" %in% covs) {
# #     nhw_prevalence <- age_grp_table(
# #       lower_ages   = c(0, 5, 10, 15, 20, 25, 30, 35, 45, 55, 65, 75, 85),
# #       upper_offset = 0,
# #       num_vec      = c(0.48087, 0.49008, 0.49306, 0.51711, 0.53047, 0.53652, 0.55121, 
# #                        0.55932, 0.6119, 0.68936, 0.74374, 0.7696, 0.7898),
# #       num_var_name = "prevalence"
# #     )
# #     nhw_func_pop <- stepfun(x = nhw_prevalence[["lower"]],
# #                                      y = c(0, nhw_prevalence[["prevalence"]]),
# #                                      right = FALSE)
# #     nhw_b        <- aggregate(as.formula(paste0(nhw_var, " ~ ",
# #                                                          last_entry_age_var)),
# #                                        FUN = mean, data = mgi_data)
# #     nhw_func_mgi <- stepfun(x = nhw_b[, 1], y = c(0, nhw_b[, 2]),
# #                                      right = FALSE)
# #   }
  
# #   # without cancer
# #   population <- age_func_pop(age_vec)
# #   if ("diabetes" %in% covs) {
# #     population <- population * fifelse(mgi_data[[diabetes_var]] == 1,
# #                           diabetes_func_pop(age_vec),
# #                           1 - diabetes_func_pop(age_vec))
# #   }
# #   if ("cad" %in% covs) {
# #     population <- population * fifelse(mgi_data[[chd_var]] == 1,
# #                                        chd_func_pop(age_vec),
# #                                        1 - chd_func_pop(age_vec))
# #   }
# #   if ("smoke" %in% covs) {
# #     population <- population * fifelse(mgi_data[[smoke_var]] == 1,
# #                                        smoke_func_pop(age_vec),
# #                                        1 - smoke_func_pop(age_vec))
# #   }
# #   if ("female" %in% covs) {
# #     population <- population * fifelse(mgi_data[[female_var]] == 1,
# #                                        female_func_pop(age_vec),
# #                                        1 - female_func_pop(age_vec))
# #   }
# #   if ("depression" %in% covs) {
# #     population <- population * fifelse(mgi_data[[depression_var]] == 1,
# #                                        depression_func_pop(age_vec),
# #                                        1 - depression_func_pop(age_vec))
# #   }
# #   if ("hypertension" %in% covs) {
# #     population <- population * fifelse(mgi_data[[hypertension_var]] == 1,
# #                                        hypertension_func_pop(age_vec),
# #                                        1 - hypertension_func_pop(age_vec))
# #   }
# #   if ("nhw" %in% covs) {
# #     population <- population * fifelse(mgi_data[[nhw_var]] == 1,
# #                                        nhw_func_pop(age_vec),
# #                                        1 - nhw_func_pop(age_vec))
# #   }
# #   if ("cancer" %in% covs) {
# #     population <- population * fifelse(mgi_data[[cancer_var]] == 1,
# #                                        cancer_func_pop(age_vec),
# #                                        1 - cancer_func_pop(age_vec))
# #   }
  
  
# #   mgi <- age_func_mgi(age_vec)
# #   if ("diabetes" %in% covs) {
# #     mgi <- mgi* fifelse(mgi_data[[diabetes_var]] == 1,
# #                    diabetes_func_mgi(age_vec),
# #                    1 - diabetes_func_mgi(age_vec))
# #   }
# #   if ("cad" %in% covs) {
# #     mgi <- mgi * fifelse(mgi_data[[chd_var]] == 1,
# #                          chd_func_mgi(age_vec),
# #                          1 - chd_func_mgi(age_vec))
# #   }
# #   if ("smoke" %in% covs) {
# #     mgi <- mgi * fifelse(mgi_data[[smoke_var]] == 1,
# #                          smoke_func_mgi(age_vec),
# #                          1 - smoke_func_mgi(age_vec))
# #   }
# #   if ("female" %in% covs) {
# #     mgi <- mgi * fifelse(mgi_data[[female_var]] == 1,
# #                          female_func_mgi(age_vec),
# #                          1 - female_func_mgi(age_vec))
# #   }
# #   if ("depression" %in% covs) {
# #     mgi <- mgi * fifelse(mgi_data[[depression_var]] == 1,
# #                          depression_func_mgi(age_vec),
# #                          1 - depression_func_mgi(age_vec))
# #   }
# #   if ("hypertension" %in% covs) {
# #     mgi <- mgi * fifelse(mgi_data[[hypertension_var]] == 1,
# #                          hypertension_func_mgi(age_vec),
# #                          1 - hypertension_func_mgi(age_vec))
# #   }
# #   if ("nhw" %in% covs) {
# #     mgi <- mgi * fifelse(mgi_data[[nhw_var]] == 1,
# #                          nhw_func_mgi(age_vec),
# #                          1 - nhw_func_mgi(age_vec))
# #   }
# #   if ("cancer" %in% covs) {
# #     mgi <- mgi * fifelse(mgi_data[[cancer_var]] == 1,
# #                          cancer_func_mgi(age_vec),
# #                          1 - cancer_func_mgi(age_vec))
# #   }
  
# #   if (chop == TRUE) {
# #     poststrat_weights <- chopr(population / mgi)
# #   } else {
# #     poststrat_weights <- population / mgi
# #   }
  
# #   poststrat_weights <- (Nobs * poststrat_weights) / sum(poststrat_weights, na.rm = TRUE)
  
# #   return(
# #     data.table(
# #       id       = mgi_data[[id_var]],
# #       ps_weight = poststrat_weights
# #     )
# #   )
  
# # }
