
names(aou_models)
aou_models[[1]]
names(aou_models[[1]])

tmp_data <- ukb_tr_merged[["t2_threshold"]]

tmp_data[, rfpred_bin := cut_interval(rfpred, n = 10)]

tmp_data[, .N, c("rfpred_bin", "case")][case == 1, ] |>
    ggplot(aes(x = rfpred_bin, y =N)) +
    geom_bar(stat = "identity")

tmp_data |>
    ggplot(aes(x = rfpred_bin, y = case)) +
    geom_bar(stat = "identity")

tmp_data[, .N, c("rfpred_bin", "case")] |>
    dcast(rfpred_bin~case) |>
    mutate(prop = `1`/(`0` + `1`)) |>
    ggplot(aes(x = rfpred_bin, y = prop)) +
    geom_bar(stat = "identity")

# ridge
## unweighted
rownames(aou_models[["unweighted_ridge"]][["beta"]])
urdata <- data0(tmp_data, rownames(aou_models[["unweighted_ridge"]][["beta"]]))

tmp2 <- predict(aou_models[["unweighted_ridge"]], newx = as.matrix(urdata)) |> scale()
lapply(aou_models, class)

wrdata <- data0(tmp_data, rownames(aou_models[["weighted_ridge"]][["beta"]]))
tmp3 <- predict(aou_models[["weighted_ridge"]], newx = as.matrix(wrdata)) |> scale()
tmp_data[, wrpred := tmp3]

wedata <- data0(tmp_data, rownames(aou_models[["weighted_enet"]][["beta"]]))
tmp3 <- predict(aou_models[["weighted_enet"]], newx = as.matrix(wedata)) |> scale()
tmp_data[, wepred := tmp3]

tmp_rf_pred <- predict(aou_models[["unweighted_rf"]], data0(tmp_data, aou_models[["unweighted_rf"]]$independent.variable.names))$predictions

tmp_data[, rfpred := tmp_rf_pred]
glm(case ~ rfpred, data = tmp_data, family = "binomial") |> summary()
pROC::roc(case ~ rfpred, data = tmp_data)

tmp_data[, urpred:= tmp2]
glm(case ~ wepred, data = tmp_data, family = "binomial") |> summary()
pROC::roc(case ~ wepred, data = tmp_data)

aimTwo::get_bin_diagnostics(tmp_data, outcome = "case", exposure = "rfpred")

# prep
test_data     <- ukb_tr_merged[["t2_threshold"]]
fitted_models <- aou_models

# glmnet
glmnet_models <- names(fitted_models)[sapply(fitted_models, function(x) "glmnet" %in% class(x))]
for (i in seq_along(glmnet_models)) {
    vars <- rownames(fitted_models[[i]]$beta)
    tmp2 <- predict(fitted_models[[i]], newx = as.matrix(data0(test_data, vars)))
    test_data[, paste0("phers_", glmnet_models[i]) := tmp2]
}

# rf
rf_vars <- fitted_models[["unweighted_rf"]]$independent.variable.names
test_data[, phers_unweighted_rf := predict(fitted_models[["unweighted_rf"]], data = data0(test_data, rf_vars))$predictions]


test_data_phers_diagnostics <- map(
    phers_vars,
    \(x) {
        aimTwo::get_bin_diagnostics(test_data, outcome = "case", exposure = x)
    }
) |> set_names(phers_vars)


one_step_phers <- function(
    test_data,
    fitted_models,
    glmnet_models,
    rf_models,
    id_var      = "id",
    outcome_var = "case"
) {
    data <- data.table::copy(test_data)
    # glmnet
    # glmnet_models <- names(fitted_models)[sapply(fitted_models, function(x) "glmnet" %in% class(x))]
    for (i in seq_along(glmnet_models)) {
        vars <- rownames(fitted_models[[i]]$beta)
        tmp2 <- predict(fitted_models[[i]], newx = as.matrix(data0(data, vars)))
        data[, paste0("phers_", glmnet_models[i]) := tmp2]
    }

    # rf
    for (i in seq_along(rf_models)) {
        vars <- fitted_models[[rf_models[i]]]$independent.variable.names
        data[, paste0("phers_", rf_models[i]) := predict(fitted_models[[rf_models[i]]], data = data0(data, vars))$predictions]
    }

    phers_vars <- names(test_data)[startsWith(names(test_data), "phers_")]
    phers_id_vars <- c(id_var, outcome_var, phers_vars)

    data_phers_diagnostics <- map(
        phers_vars,
        \(x) {
            aimTwo::get_bin_diagnostics(data, outcome = "case", exposure = x)
        }
    ) |> set_names(phers_vars)

    return(
        list(
            data = data[, ..phers_id_vars],
            diagnostics = data_phers_diagnostics
        )
    )

}

one_step_res <- one_step_phers(
    test_data     = test_data,
    fitted_models = fitted_models,
    glmnet_models = names(fitted_models)[sapply(fitted_models, function(x) "glmnet" %in% class(x))],
    rf_models     = "unweighted_rf"
)

