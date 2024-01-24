# suppressPackageStartupMessages({
#     library(cli)
#     library(glmnet)
#     library(data.table)
#     library(doParallel)
#     library(foreach)
#     library(parallelly)
# })

# crossfit_lasso_weights <- function(
#     internal_data,
#     external_data,
#     select_vars  = c(
#         "id", "internal", "age_50", "female", "nhw",
#         "hypertension", "diabetes", "cancer", "anxiety",
#         "depression", "bmi_obese", "bmi_overweight", "bmi_underweight", "bmi_unknown"),
#     id_var       = "id",
#     internal_var = "internal",
#     ncores       = parallelly::availableCores() / 4,
#     folds         = 5
# ) {

#     stack <- rbindlist(list(
#         internal_data,
#         external_data
#     ), use.names = TRUE, fill = TRUE)

#     stack_cc <- stack[complete.cases(stack[, ..select_vars]), ]

#     # split into folds
#     split_index <- sample(folds, nrow(stack_cc), replace = TRUE)

#     # create empty list to store results
#     results <- list()

#     for (i in seq_len(folds)) {
#         cli_progress_step(paste0("Fold ", i, " of ", folds), spinner = TRUE)

#         # split into train and test
#         train <- stack_cc[split_index == i, ..select_vars]
#         test  <- stack_cc[split_index != i, ]

#         data <- copy(train)
#         data[[id_var]] <- NULL

#         # first step: using .*. for all interactions
#         f <- as.formula(paste0(internal_var, " ~ .*."))
#         y <- data[[internal_var]]

#         # create all interactions
#         x <- model.matrix(f, data)[, -1]

#         # perform cross-validation to select lambda
#         cl <- parallel::makeCluster(ncores, type = "PSOCK")
#         doParallel::registerDoParallel(cl)
#         cv.lambda.lasso <- cv.glmnet(
#             x = x, y = y,
#             alpha = 1, family = "binomial",
#             parallel = TRUE
#         )
#         stopCluster(cl)

#         # use lambda.1se to select model with parsimony
#         l.lasso.1se <- cv.lambda.lasso$lambda.1se
#         lasso.model <- glmnet(
#             x = x, y = y,
#             alpha = 1, family = "binomial",
#             lambda = l.lasso.1se
#         )

#         reg <- coef(lasso.model)
#         use_these_terms <- rownames(reg)[reg[, 1] != 0]
#         if ("(Intercept)" %in% use_these_terms) {
#             use_these_terms <- use_these_terms[use_these_terms != "(Intercept)"]
#         }
#         pretty_terms <- gsub(":", "*", use_these_terms)

#         # get weights
#         results[[i]] <- ipw(
#             stacked_data       = stack_cc[split_index != i, ],
#             weight_outcome_var = "WTFA_A",
#             samp_var           = "samp_WTFA_A",
#             external_dataset   = "NHIS",
#             dataset_name       = "MGI",
#             cancer_factor      = FALSE,
#             covs               = pretty_terms
#         )

#     }
#     cli_progress_done()

#     # combine results
#     results <- rbindlist(results, use.names = TRUE, fill = TRUE)

#     # average over predictions
#     results[, .(ip_weight = mean(ip_weight)), by = "id"]

# }
