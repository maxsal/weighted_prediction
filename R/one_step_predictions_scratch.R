
time_threshold <- 1

mgi_data <- qread(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "{opt$outcome}/time_restricted_phenomes/mgi_",
  "{opt$mgi_version}_{opt$outcome}_",
  "t{time_threshold}_pim_",
  "r{opt$matching_ratio}.qs"
))
mgi_ids <- mgi_data[, .(id, case)]

ukb_data <- qread(glue(
  "data/private/ukb/{opt$ukb_version}/",
  "{opt$outcome}/time_restricted_phenomes/ukb_",
  "{opt$ukb_version}_{opt$outcome}_",
  "t{time_threshold}_pim_",
  "r{opt$matching_ratio}.qs"
))
ukb_ids <- ukb_data[, .(id, case)]

p <- ms::pheinfox


mgi_train_obs <- sample(seq_len(nrow(mgi_data)), size = round(nrow(mgi_data) * opt$split_prop))

mgi_train <- mgi_data[mgi_train_obs, ]
mgi_test  <- mgi_data[!mgi_train_obs, ]

mgi_train_covs <- mgi_train[, !c("id", "case", "length_followup", "unique_encounters", "unique_phecodes", "age", "sex")]
mgi_train_y    <- mgi_train[, case]

phecode_vars <- names(mgi_data)[names(mgi_data) %in% ms::pheinfox[, phecode]]
dials::grid_regular(
  dials::finalize(dials::mtry(),mgi_data[, ..phecode_vars]),
  levels = 5
)

mtry_seq <- round(floor(sqrt(ncol(mgi_data[, ..phecode_vars]))) * c(0.25, 0.5, 1, 2, 4))

learners <- create.Learner("SL.ranger", tune = list(mtry = mtry_seq))

options(mc.cores = 12)
set.seed(1, "L'Ecuyer-CMRG")
SL_library <- list(
  "SL.ranger", "SL.glm", c("SL.glmnet", "screen.corP")
)
cv_sl <- CV.SuperLearner(
  Y          = y,
  X          = x,
  SL.library = SL_library,
  family     = binomial(),
  parallel   = "multicore",
  cvControl  = list(V = 3)
)




out <- rbindlist(list(
  out_ridge,
  out_enet,
  out_rf,
  out_lasso
))

out[, value := unlist(value)]

# save results -----------------------------------------------------------------
fwrite(
  x = out,
  file = paste0("results/mgi/", opt$mgi_version, "/", opt$outcome, "/mgi_", opt$mgi_version, "_", opt$outcome, "_t0_hyperparameters.csv")
)

libri(caret)

preProcess(
  mgi_tr_merged[[1]][, ..in_x],
  method = c("center", "scale", "nzv")
)

nearZeroVar(mgi_tr_merged[[1]][, ..in_x])

findCorrelation(mgi_tr_merged[[1]][, ..in_x], cutoff = 0.9)


########
ukb_version <- "20221117"
aou_version <- "20230309"
outcome     <- "CA_101.8"
time_threshold <- 1

aou_models <- qread(paste0("results/aou/", aou_version, "/", outcome, "/aou_20230309_CA_101.8_t1_models.qs"))

ukb_data <- qread(glue(
    "data/private/ukb/{opt$ukb_version}/",
    "{opt$outcome}/time_restricted_phenomes/ukb_",
    "{opt$ukb_version}_{opt$outcome}_",
    "t1_pim_",
    "r2.qs"
  ))

drop_and_add <- function(
  newx, model_names) {
  newx_names <- names(newx)

  not_in_model <- newx_names[!newx_names %in% model_names]
  not_in_newx  <- model_names[!model_names %in% newx_names]

  out <- data.table::copy(newx)

  # drop columns not in model and add columns not in newx
  if (length(not_in_model) > 0) {
    set(out, , not_in_model, NULL)
  }
  if (length(not_in_newx) > 0) {
    set(out, , not_in_newx, 0)
  }

  return(as.matrix(out[]))
}

ukb_data

ukb_fitted <- data.table(
  id   = ukb_data$id,
  case = ukb_data$case
)

for (i in 1:6) {
  ukb_prepped <- drop_and_add(
    newx        = ukb_data,
    model_names = rownames(aou_models[[i]]$beta)
  )
  var_name <- names(aou_models)[i]
  ukb_fitted[[var_name]] <- predict(
      aou_models[[i]],
      newx = ukb_prepped,
      type = "response"
    )[, 1]
}


ukb_data |>
  drop_and_add(
    model_names = rownames(aou_models[[1]]$beta)
  ) |>
  (\(x) {
    x |>
    dplyr::mutate(
      predicted = predict(
        aou_models[[1]],
        newx = as.matrix(x),
        type = "response"
      )
    )
  })()

  predict(
    aou_models[[1]],
    newx = _,
    type = "response"
  )

tmp <- ukb_data |>
mutate(predicted = predict(
  aou_models[[1]],
  newx = as.matrix(ukb_sub),
  type = "response"
)[, 1])

pROC::auc(tmp$case, tmp$predicted)

map(
  names(ukb_fitted)[-c(1, 2)],
  ~ pROC::auc(ukb_fitted$case, ukb_fitted[[.x]])
)

drop_and_add(
  newx        = ukb_data,
  model_names = aou_models[["unweighted_rf"]][["independent.variable.names"]]
)

rf_pred <- predict(
  aou_models[["unweighted_rf"]],
  data = drop_and_add(
  newx        = ukb_data,
  model_names = aou_models[["unweighted_rf"]][["independent.variable.names"]]
)
)

ukb_fitted[["unweighted_rf"]] <- rf_pred$predictions

ukb_fitted[["superlearner"]] <- predict(
  aou_models[["superlearner"]][["AllSL"]][["3"]],
  data = drop_and_add(
    newx        = ukb_data,
    model_names = aou_names
  )
)$pred[, 1]

aou_names <- c()
for (i in 1:6) {
  aou_names <- unique(c(aou_names, rownames(aou_models[[i]]$beta)))
}
drop_and_add(
  newx        = ukb_data,
  model_names = aou_names
)


aou_models[["superlearner"]][["AllSL"]][["1"]] |> names()

