ms::libri(ms, qs, glmnet, data.table, cli, ggplot2, ranger, wlasso)

set.seed(42)

# pick some settings
n_mtry        <- 6
n_node_size   <- 5
n_sample_frac <- 4
n_trees       <- 500
n_cores       <- parallelly::availableCores() - 2
weight_var    <- "weight"


# make data
n      <- 1000
p      <- 5000
real_p <- 15

x <- matrix(rnorm(n*p), nrow = n, ncol = p)
y <- as.numeric(apply(x[, 1:real_p], 1, sum) + rnorm(n) > 0)
w <- runif(n)

# enet function
cva.glmnet <- function(x, y, type.measure = "deviance", family = "binomial", try_alpha = seq(0, 1, by = 0.05), nfolds = 10, ...) {
    cli_progress_step("Fitting models...")
    list_of_fits <- list()
    cli_progress_bar(
        name   = "alpha tuning",
        format = "alpha: {try_alpha[i]} {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}",
        total  = length(try_alpha)
    )
    for (i in seq_along(try_alpha)) {
        fit_name <- paste0("alpha", try_alpha[i])

        list_of_fits[[fit_name]] <- cv.glmnet(x, y, type.measure = type.measure, alpha = try_alpha[i], family = family, nfolds = nfolds, ...)

        cli_progress_update()
    }
    cli_progress_done()

    cli_progress_step("Estimating MSEs...")
    results <- data.table()
    cli_progress_bar(
        name   = "getting results",
        format = "alpha: {try_alpha[i]} {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}",
        total  = length(try_alpha)
    )
    for (i in seq_along(try_alpha)) {
        fit_name <- paste0("alpha", try_alpha[i])

        predicted <- predict(list_of_fits[[fit_name]], newx = x.test, s = "lambda.1se")

        temp <- data.table(
            alpha      = try_alpha[i],
            lambda.min = list_of_fits[[fit_name]]$lambda.min,
            lambda.1se = list_of_fits[[fit_name]]$lambda.1se,
            measure    = list_of_fits[[fit_name]]$cvm[which(list_of_fits[[fit_name]]$lambda == list_of_fits[[fit_name]]$lambda.min)]
        )
        results <- rbind(results, temp)
        cli_progress_update()
    }
    cli_progress_done()

    cli_alert_info("optimal alpha: {results[which.min(results$measure), 'alpha']}; lambda.1se: {results[which.min(results$measure), 'lambda.1se']}")

    return(results)
}


get_hyperparameters <- function(
    x,
    y,
    w = NULL,
    methods = c("ridge", "lasso", "enet", "rf")
) {

    if (is.null(w)) {
        comb <- cbind(
            data.table(case = y),
            as.data.table(x)
        )
    } else {
        comb <- cbind(
            data.table(case = y),
            as.data.table(x),
            data.table(weight = w)
        )
        wcomb <- comb[complete.cases(comb), ]
    }

    if ("ridge" %in% methods) {
        # ridge
        cli_progress_step("Fitting ridge...")
        ## unweighted
        ridge.fit <- cv.glmnet(x, y, type.measure = "deviance", alpha = 0, family = "binomial") # set type = 'deviance' for logistic regression
        ## weighted
        if (!is.null(w)) {
            wridge.fit <- cv.glmnet(wcomb[, !c("y")], y,
                type.measure = "deviance", alpha = 0, family = "binomial",
                penalty.factor = c(rep(1, ncol(x) - 1), 0)
            )
        }
    }

    if ("lasso" %in% methods) {
        # lasso
        cli_progress_step("Fitting lasso...")
        ## unweighted
        lasso.fit <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "binomial")
        ## weighted - using wlasso
        if (!is.null(w)) {
            cli_progress_step("Fitting weighted lasso...")
            wlasso.fit <- wlasso(
                data = wcomb,
                col.y = "case", col.x = names(comb)[!(names(comb) %in% c("case", weight_var))],
                family = "binomial",
                weights = weight_var,
                method = "dCV", k = 10
            )
        }
    }

    if ("enet" %in% methods) {
        # elastic net
        cli_progress_step("Fitting elastic net...")
        # unweighted
        enet.fit <- cva.glmnet(x = x, y = y, type.measure = "deviance", fam = "binomial", nfolds = 10)
        if (!is.null(w)) {
            # weighted
            cli_alert_danger("WEIGHTED ENET IS NOT READY - DO NOT USE - SKIPPING")
            # wenet.fit <- cva.glmnet(x = x, y = y, type.measure = "deviance", fam = "binomial", nfolds = 10)
        }
    }

    if ("rf" %in% methods) {
        # random forest
        cli_progress_step("Fitting random forest...")
        rf <- expand.grid(
            mtry        = round(seq(5, round((ncol(comb) - 1) / 2), length.out = n_mtry)),
            node_size   = round(seq(3, round(nrow(comb) / 100), length.out = n_node_size)),
            sample_frac = seq(0.5, 0.8, length.out = n_sample_frac),
            oob_rmse    = 0
        ) |> as.data.table()

        cli_alert("Trying {n_mtry} values of mtry, {n_node_size} of node size, {n_sample_frac} values of sample fraction ({nrow(rf)} combinations)")

        cli_progress_bar(name = "grid search...", total = nrow(rf))
        for (i in 1:nrow(rf)) {
            model <- ranger(
                formula         = case ~ .,
                data            = comb,
                num.trees       = n_trees,
                mtry            = rf$mtry[i],
                min.node.size   = rf$node_size[i],
                sample.fraction = rf$sample_frac[i],
                num.threads     = n_cores
            )
            rf$oob_rmse[i] <- sqrt(model$prediction.error)
            cli_progress_update()
        }
        cli_progress_done()
    }

    # OUTPUT HYPERPARAMETERS

}

comb <- cbind(
    data.table(case = y),
    as.data.table(x),
    data.table(weight = w)
)
wcomb <- comb[complete.cases(comb), ]


# ridge
cli_progress_step("Fitting ridge...")
## unweighted
ridge.fit <- cv.glmnet(x, y, type.measure = "deviance", alpha = 0, family = "binomial") # set type = 'deviance' for logistic regression
## weighted
wridge.fit <- cv.glmnet(wcomb[, !c("y")], y, type.measure = "deviance", alpha = 0, family = "binomial",
                        penalty.factor = c(rep(1, ncol(x) - 1), 0))

# lasso
cli_progress_step("Fitting lasso...")
## unweighted
lasso.fit <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "binomial")
## weighted - using wlasso
cli_progress_step("Fitting weighted lasso...")
wlasso.fit <- wlasso(
    data  = wcomb,
    col.y = "case", col.x = names(comb)[!(names(comb) %in% c("case", weight_var))],
    family = "binomial",
    weights = weight_var,
    method = "dCV", k = 10
)

# elastic net
cli_progress_step("Fitting elastic net...")
# unweighted
enet.fit <- cva.glmnet(x = x, y = y, type.measure = "deviance", fam = "binomial", nfolds = 10)
# weighted
enet.fit <- cva.glmnet(x = x, y = y, type.measure = "deviance", fam = "binomial", nfolds = 10)

# random forest
cli_progress_step("Fitting random forest...")
rf <- expand.grid(
    mtry        = round(seq(5, round((ncol(comb) - 1) / 2), length.out = n_mtry)),
    node_size   = round(seq(3, round(nrow(comb) / 100), length.out = n_node_size)),
    sample_frac = seq(0.5, 0.8, length.out = n_sample_frac),
    oob_rmse    = 0
) |> as.data.table()

cli_alert("Trying {n_mtry} values of mtry, {n_node_size} of node size, {n_sample_frac} values of sample fraction ({nrow(rf)} combinations)")

cli_progress_bar(name = "grid search...", total = nrow(rf))
for (i in 1:nrow(rf)) {
    model <- ranger(
        formula         = case ~ .,
        data            = comb,
        num.trees       = n_trees,
        mtry            = rf$mtry[i],
        min.node.size   = rf$node_size[i],
        sample.fraction = rf$sample_frac[i],
        num.threads     = n_cores
    )
    rf$oob_rmse[i] <- sqrt(model$prediction.error)
    cli_progress_update()
}
cli_progress_done()

data.table(
    "parameter" = c(
        "ridge.lambda.min",
        "ridge.lambda.1se",
        "lasso.lambda.min",
        "lasso.lambda.1se",
        "wlasso.lambda.min",
        "enet.alpha",
        "enet.lambda.min",
        "enet.lambda.1se",
        "rf.mtry",
        "rf.node_size",
        "rf.sample_frac"
    ),
    "value" = c(
        ridge.fit$lambda.min,
        ridge.fit$lambda.1se,
        lasso.fit$lambda.min,
        lasso.fit$lambda.1se,
        wlasso.fit$lambda.min,
        enet.fit[which.min(enet.fit$measure), "alpha"],
        enet.fit[which.min(enet.fit$measure), "lambda.min"],
        enet.fit[which.min(enet.fit$measure), "lambda.1se"],
        rf[which.min(rf$oob_rmse), mtry],
        rf[which.min(rf$oob_rmse), node_size],
        rf[which.min(rf$oob_rmse), sample_frac]
    )
)
