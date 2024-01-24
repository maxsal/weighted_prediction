suppressPackageStartupMessages({
    library(data.table)
    library(glmnet)
    library(ranger)
    library(wlasso)
    library(doParallel)
    library(parallelly)
    library(foreach)
})

# enet function - for use in `get_hyperparameters()`
cva.glmnet <- function(
    x,
    y,
    type.measure = "deviance",
    family       = "binomial",
    try_alpha    = seq(0, 1, by = 0.05),
    nfolds        = 10,
    ...
) {
    cli_progress_step("Fitting models...")
    list_of_fits <- list()
    cli_progress_bar(
        name   = "alpha tuning",
        format = "alpha: {try_alpha[i]} {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}",
        total  = length(try_alpha)
    )
    out <- data.table()

    for (i in seq_along(try_alpha)) {
        fit_name <- paste0("alpha", try_alpha[i])

        tmp.enet.fit <- cv.glmnet(x, y, type.measure = type.measure, alpha = try_alpha[i], family = family, nfolds = nfolds, ...)

        out <- rbindlist(list(
            out,
            data.table(
                alpha      = try_alpha[i],
                lambda.min = tmp.enet.fit$lambda.min,
                lambda.1se = tmp.enet.fit$lambda.1se,
                measure    = tmp.enet.fit$cvm[which(tmp.enet.fit$lambda == tmp.enet.fit$lambda.min)]
            )
        ))
        cli_progress_update()
    }

    cli_progress_done()

    cli_alert_info("optimal alpha: {out[which.min(out$measure), 'alpha']}; lambda.1se: {out[which.min(out$measure), 'lambda.1se']}")

    return(out)
}

get_hyperparameters <- function(
    x,        # exposures as a data.table
    y,        # outcome as a vector
    w             = NULL, # weight as a vector
    methods       = c("ridge", "lasso", "enet", "rf"),
    parallel      = TRUE,
    n_cores       = 10,
    core_prop     = 1/8,
    n_node_size   = 5,
    n_trees       = 500
) {
    if (parallel) {
        cl <- parallel::makeCluster(n_cores, type = "PSOCK")
        doParallel::registerDoParallel(cl)
    }

    if (is.null(w)) {
        comb <- cbind(
            data.table(case = y),
            as.data.table(x)
        )
        comb <- comb[complete.cases(comb), ]
    } else {
        comb <- cbind(
            data.table(case = y),
            as.data.table(x),
            data.table(weight = w)
        )
        comb <- comb[complete.cases(comb[, !c("weight")]), ]
        wcomb <- comb[complete.cases(comb), ]
    }

    out <- data.table()

    if ("ridge" %in% methods) {
        # ridge
        cli_progress_step("Fitting ridge...")
        ## unweighted
        ridge.fit <- cv.glmnet(as.matrix(x), y, type.measure = "deviance", alpha = 0, family = "binomial", parallel = parallel) # set type = 'deviance' for logistic regression
        out <- rbindlist(list(
            out,
            data.table(
                "parameter" = c("ridge.lambda.min", "ridge.lambda.1se"),
                "value"     = c(ridge.fit$lambda.min, ridge.fit$lambda.1se)
            )
        ), use.names = TRUE, fill = TRUE)
        ## weighted
        if (!is.null(w)) {
            wridge.fit <- cv.glmnet(x = as.matrix(wcomb[, !c("case")]), y = wcomb[, case],
                type.measure   = "deviance", alpha = 0, family = "binomial",
                penalty.factor = as.numeric(names(wcomb[, !c("case")]) != "weight"),
                parallel = parallel
            )
            out <- rbindlist(list(
                out,
                data.table(
                    "parameter" = c("wridge.lambda.min", "wridge.lambda.1se"),
                    "value"     = c(wridge.fit$lambda.min, wridge.fit$lambda.1se)
                )
            ), use.names = TRUE, fill = TRUE)
        }

    }

    if ("lasso" %in% methods) {
        # lasso
        cli_progress_step("Fitting lasso...")
        ## unweighted
        lasso.fit <- cv.glmnet(as.matrix(x), y, type.measure = "deviance", alpha = 1, family = "binomial", parallel = parallel)
        out <- rbindlist(list(
            out,
            data.table(
                "parameter" = c("lasso.lambda.min", "lasso.lambda.1se"),
                "value"     = c(lasso.fit$lambda.min, lasso.fit$lambda.1se)
            )
        ), use.names = TRUE, fill = TRUE)
        ## weighted - using wlasso
        if (!is.null(w)) {
            cli_progress_step("Fitting weighted lasso...")
            wlasso.fit <- wlasso(
                data    = wcomb,
                col.y   = "case",
                col.x   = names(comb)[!(names(comb) %in% c("case", "weight"))],
                family  = "binomial",
                weights = "weight",
                method  = "dCV", k = 10
            )
            out <- rbindlist(list(
                out,
                data.table(
                    "parameter" = "wlasso.lambda.min",
                    "value"     = wlasso.fit$lambda.min
                )
            ), use.names = TRUE, fill = TRUE)
        }
    }

    if ("enet" %in% methods) {
        # elastic net
        cli_progress_step("Fitting elastic net...")
        # unweighted
        enet.fit <- cva.glmnet(x = as.matrix(x), y = y, type.measure = "deviance", fam = "binomial", nfolds = 10, parallel = parallel)
        out <- rbindlist(list(
            out,
            data.table(
                "parameter" = c("enet.alpha", "enet.lambda.min", "enet.lambda.1se"),
                "value"     = c(
                    enet.fit[which.min(enet.fit$measure), "alpha"],
                    enet.fit[which.min(enet.fit$measure), "lambda.min"],
                    enet.fit[which.min(enet.fit$measure), "lambda.1se"]
                )
            )
        ), use.names = TRUE, fill = TRUE)
        if (!is.null(w)) {
            # weighted
            wenet.fit <- cva.glmnet(
                x = as.matrix(wcomb[, !c("case")]), y = wcomb[, case],
                type.measure = "deviance", fam = "binomial", nfolds = 10,
                penalty.factor = as.numeric(names(wcomb[, !c("case")]) != "weight"), parallel = parallel
            )
            out <- rbindlist(list(
                out,
                data.table(
                    "parameter" = c("wenet.alpha", "wenet.lambda.min", "wenet.lambda.1se"),
                    "value" = c(
                        wenet.fit[which.min(wenet.fit$measure), "alpha"],
                        wenet.fit[which.min(wenet.fit$measure), "lambda.min"],
                        wenet.fit[which.min(wenet.fit$measure), "lambda.1se"]
                    )
                )
            ), use.names = TRUE, fill = TRUE)
        }
    }

    if ("rf" %in% methods) {
        # random forest
        cli_progress_step("Fitting random forest...")
        rf <- expand.grid(
            mtry        = round(floor(sqrt(ncol(comb[, !c("case", "weight")]))) * c(0.5, 1, 2)),
            node_size   = round(seq(3, round(nrow(comb) / 100), length.out = n_node_size)),
            oob_rmse    = 0
        ) |> as.data.table()

        cli_alert("Trying 3 values of mtry, {n_node_size} of node size ({nrow(rf)} combinations)")

        cli_progress_bar(name = "grid search...", total = nrow(rf))
        for (i in 1:nrow(rf)) {
            model <- ranger(
                formula         = case ~ .,
                data            = comb[, !c("weight")],
                num.trees       = n_trees,
                mtry            = rf$mtry[i],
                min.node.size   = rf$node_size[i],
                num.threads     = parallelly::availableCores() * core_prop
            )
            rf$oob_rmse[i] <- sqrt(model$prediction.error)
            cli_progress_update()
        }
        cli_progress_done()
        out <- rbindlist(list(
            out,
            data.table(
                "parameter" = c("rf.mtry", "rf.node_size"),
                "value"     = c(
                    rf[which.min(rf$oob_rmse), mtry],
                    rf[which.min(rf$oob_rmse), node_size]
                )
            )
        ), use.names = TRUE, fill = TRUE)
    }
    if (parallel) parallel::stopCluster(cl)
    return(out)

}
