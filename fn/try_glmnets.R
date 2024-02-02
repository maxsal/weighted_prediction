try_glmnets <- function(
    x,
    y,
    weights = NULL,
    alpha = 1,
    lambda,
    alt_lambda,
    family = "binomial",
    maxit = 1e5,
    cor_cutoff = 0.25,
    ...
) {

    out <- list(
        model = NA,
        diag = data.table::data.table()
    )

    # Class 1: primary lambda, model weights, no screen
    out <- tryCatch(
        {
            model <- glmnet(
                x = x,
                y = y,
                weights = weights,
                alpha = alpha,
                lambda = lambda,
                family = family,
                maxit = maxit,
                ...
            )
            if (model$dev.ratio < 0 | model$df == 0 | is.null(model$lambda)) {
                warning(
                    paste0("dev.ratio: ", model$dev.ratio, ", df: ", model$df, "; lambda: ", model$lambda, ".")
                )
            } else {
                list(
                    model = model,
                    diag = rbindlist(list(
                        out[["diag"]],
                        data.table::data.table(
                            class = 1,
                            lambda = "primary",
                            weight = "model",
                            screen = FALSE,
                            dev.ratio = model$dev.ratio,
                            df = model$df,
                            n0_beta = sum(as.matrix(model$beta)[, 1] != 0),
                            converge = ifelse(is.null(model$lambda) | is.infinite(model$lambda), FALSE, TRUE)
                        )
                    ), use.names = TRUE, fill = TRUE)
                )
            }
        },
        warning = function(w) {
            warning(w)
            model <- glmnet(
                x = x,
                y = y,
                weights = weights,
                alpha = alpha,
                lambda = lambda,
                family = family,
                maxit = maxit,
                ...
            )
            list(
                model = NA,
                diag = rbindlist(list(
                    out[["diag"]],
                    data.table::data.table(
                        class = 1,
                        lambda = "primary",
                        weight = "model",
                        screen = FALSE,
                        dev.ratio = model$dev.ratio,
                        df = model$df,
                        converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                    )
                ), use.names = TRUE, fill = TRUE)
            )
        }
    )

    # Class 2: alternate lambda, model weights, no screen
    if (!("glmnet" %in% class(out$model))) { 
        out <- tryCatch(
            {
                model <- glmnet(
                    x = x,
                    y = y,
                    weights = weights,
                    alpha = alpha,
                    lambda = lambda,
                    family = family,
                    maxit = maxit,
                    ...
                )
                if (model$dev.ratio < 0 | model$df == 0 | is.null(model$lambda) | is.infinite(model$lambda)) {
                    warning(
                        paste0("dev.ratio: ", model$dev.ratio, ", df: ", model$df, "; lambda: ", model$lambda, ".")
                    )
                } else {
                    list(
                        model = model,
                        diag = rbindlist(list(
                            out[["diag"]], data.table(
                                class = 2,
                                lambda = "alternate",
                                weight = "model",
                                screen = FALSE,
                                dev.ratio = model$dev.ratio,
                                df = model$df,
                                n0_beta = sum(as.matrix(model$beta)[, 1] != 0),
                                converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                            )
                        ), use.names = TRUE, fill = TRUE)
                    )
                }
            },
            warning = function(w) {
                warning(w)
                model <- glmnet(
                    x = x,
                    y = y,
                    weights = weights,
                    alpha = alpha,
                    lambda = lambda,
                    family = family,
                    maxit = maxit,
                    ...
                )
                list(
                    model = NA,
                    diag = rbindlist(list(
                        out[["diag"]], data.table(
                            class = 2,
                            lambda = "alternate",
                            weight = "model",
                            screen = FALSE,
                            dev.ratio = model$dev.ratio,
                            df = model$df,
                            converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                        )
                    ), use.names = TRUE, fill = TRUE)
                )
            }
        )
    }

    # Class 3: primary lambda, model weights, screen
    if (!("glmnet" %in% class(out$model))) {
        tmp_x <- as.data.table(x)
        multi_values <- more_than_one_unique(tmp_x, names(tmp_x))
        remove_these <- caret::findCorrelation(cor(tmp_x[, ..multi_values]), cutoff = cor_cutoff, names = TRUE)
        use_these <- setdiff(multi_values, remove_these)
        new_x <- as.matrix(tmp_x[, ..use_these])
        out <- tryCatch(
            {
                model <- glmnet(
                    x = new_x,
                    y = y,
                    weights = weights,
                    alpha = alpha,
                    lambda = lambda,
                    family = family,
                    maxit = maxit,
                    ...
                )
                if (model$dev.ratio < 0 | model$df == 0 | is.infinite(model$lambda) | is.null(model$lambda)) {
                    warning(
                        paste0("dev.ratio: ", model$dev.ratio, ", df: ", model$df, "; lambda: ", model$lambda, ".")
                    )
                } else {
                    list(
                        model = model,
                        diag = rbindlist(list(
                            out[["diag"]], data.table(
                                class = 3,
                                lambda = "primary",
                                weight = "model",
                                screen = TRUE,
                                cor_cutoff = cor_cutoff,
                                dev.ratio = model$dev.ratio,
                                df = model$df,
                                n0_beta = sum(as.matrix(model$beta)[, 1] != 0),
                                converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                            )
                        ), use.names = TRUE, fill = TRUE)
                    )
                }
            },
            warning = function(w) {
                warning(w)
                model <- glmnet(
                    x = new_x,
                    y = y,
                    weights = weights,
                    alpha = alpha,
                    lambda = lambda,
                    family = family,
                    maxit = maxit,
                    ...
                )
                list(
                    model = NA,
                    diag = rbindlist(list(
                        out[["diag"]], data.table(
                            class = 3,
                            lambda = "primary",
                            weight = "model",
                            screen = TRUE,
                            cor_cutoff = cor_cutoff,
                            dev.ratio = model$dev.ratio,
                            df = model$df,
                            converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                        )
                    ), use.names = TRUE, fill = TRUE)
                )
            }
        )
    }

    # Class 4: alternate lambda, model weights, screen
    if (!("glmnet" %in% class(out$model))) {
        out <- tryCatch(
            {
                model <- glmnet(
                    x = new_x,
                    y = y,
                    weights = weights,
                    alpha = alpha,
                    lambda = alt_lambda,
                    family = family,
                    maxit = maxit,
                    ...
                )
                if (model$dev.ratio < 0 | model$df == 0 | is.infinite(model$lambda) | is.null(model$lambda)) {
                    warning(
                        paste0("dev.ratio: ", model$dev.ratio, ", df: ", model$df, "; lambda: ", model$lambda, ".")
                    )
                } else {
                    list(
                        model = model,
                        diag = rbindlist(list(
                            out[["diag"]], data.table(
                                class = 4,
                                lambda = "alternate",
                                weight = "model",
                                screen = TRUE,
                                cor_cutoff = cor_cutoff,
                                dev.ratio = model$dev.ratio,
                                df = model$df,
                                n0_beta = sum(as.matrix(model$beta)[, 1] != 0),
                                converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                            )
                        ), use.names = TRUE, fill = TRUE)
                    )
                }
            },
            warning = function(w) {
                warning(w)
                model <- glmnet(
                    x = new_x,
                    y = y,
                    weights = weights,
                    alpha = alpha,
                    lambda = alt_lambda,
                    family = family,
                    maxit = maxit,
                    ...
                )
                list(
                    model = NA,
                    diag = rbindlist(list(
                        out[["diag"]], data.table(
                            class = 4,
                            lambda = "alternate",
                            weight = "model",
                            screen = TRUE,
                            cor_cutoff = cor_cutoff,
                            dev.ratio = model$dev.ratio,
                            df = model$df,
                            converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                        )
                    ), use.names = TRUE, fill = TRUE)
                )
            }
        )
    }

    # Class 5: primary lambda, weights as predictors, screen
    if (!("glmnet" %in% class(out$model))) {
        out <- tryCatch(
            {
                new_x_w <- cbind(as.matrix(x), weight = weights)
                model <- glmnet(
                    x = new_x_w,
                    y = y,
                    alpha = alpha,
                    lambda = lambda,
                    family = family,
                    maxit = maxit,
                    penalty.factor = as.numeric(unlist(dimnames(new_x_w)) != "weight"),
                    ...
                )
                if (model$dev.ratio < 0 | model$df == 0 | is.infinite(model$lambda) | is.null(model$lambda) |
                    (length(names(which(as.matrix(model$beta)[, 1] != 0, useNames = TRUE))) == 1 &
                        "weight" %in% names(which(as.matrix(model$beta)[, 1] != 0, useNames = TRUE)))) {
                    warning(
                        paste0("dev.ratio: ", model$dev.ratio, ", df: ", model$df, "; lambda: ", model$lambda, ".")
                    )
                } else {
                    list(
                        model = model,
                        diag = rbindlist(list(
                            out[["diag"]], data.table(
                                class = 5,
                                lambda = "primary",
                                weight = "as_predictor",
                                screen = TRUE,
                                cor_cutoff = cor_cutoff,
                                dev.ratio = model$dev.ratio,
                                df = model$df,
                                n0_beta = sum(as.matrix(model$beta)[, 1] != 0),
                                converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                            )
                        ), use.names = TRUE, fill = TRUE)
                    )
                }
            },
            warning = function(w) {
                warning(w)
                new_x_w <- cbind(as.matrix(x), weight = weights)
                model <- glmnet(
                    x = new_x_w,
                    y = y,
                    alpha = alpha,
                    lambda = lambda,
                    family = family,
                    maxit = maxit,
                    penalty.factor = as.numeric(unlist(dimnames(new_x_w)) != "weight"),
                    ...
                )
                list(
                    model = NA,
                    diag = rbindlist(list(
                        out[["diag"]], data.table(
                            class = 5,
                            lambda = "primary",
                            weight = "as_predictor",
                            screen = TRUE,
                            cor_cutoff = cor_cutoff,
                            dev.ratio = model$dev.ratio,
                            df = model$df,
                            converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE),
                            weight_only_predictor = TRUE
                        )
                    ), use.names = TRUE, fill = TRUE)
                )
            }
        )
    }

    # Class 6: alternate lambda, weights as predictors, screen
    if (!("glmnet" %in% class(out$model))) {
        out <- tryCatch(
            {
                new_x_w <- cbind(as.matrix(x), weight = weights)
                model <- glmnet(
                    x = new_x_w,
                    y = y,
                    alpha = alpha,
                    lambda = alt_lambda,
                    family = family,
                    maxit = maxit,
                    penalty.factor = as.numeric(unlist(dimnames(new_x_w)) != "weight"),
                    ...
                )
                if (model$dev.ratio < 0 | model$df == 0 | is.infinite(model$lambda) | is.null(model$lambda) |
                    (length(names(which(as.matrix(model$beta)[, 1] != 0, useNames = TRUE))) == 1 &
                        "weight" %in% names(which(as.matrix(model$beta)[, 1] != 0, useNames = TRUE)))) {
                    warning(
                        paste0("dev.ratio: ", model$dev.ratio, ", df: ", model$df, "; lambda: ", model$lambda, ".")
                    )
                } else {
                    list(
                        model = model,
                        diag = rbindlist(list(
                            out[["diag"]], data.table(
                                class = 6,
                                lambda = "alternate",
                                weight = "as_predictor",
                                screen = TRUE,
                                cor_cutoff = cor_cutoff,
                                dev.ratio = model$dev.ratio,
                                df = model$df,
                                n0_beta = sum(as.matrix(model$beta)[, 1] != 0),
                                converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE)
                            )
                        ), use.names = TRUE, fill = TRUE)
                    )
                }
            },
            warning = function(w) {
                warning(w)
                new_x_w <- cbind(as.matrix(x), weight = weights)
                model <- glmnet(
                    x = new_x_w,
                    y = y,
                    alpha = alpha,
                    lambda = alt_lambda,
                    family = family,
                    maxit = maxit,
                    penalty.factor = as.numeric(unlist(dimnames(new_x_w)) != "weight"),
                    ...
                )
                list(
                    model = NA,
                    diag = rbindlist(list(
                        out[["diag"]], data.table(
                            class = 6,
                            lambda = "alternate",
                            weight = "as_predictor",
                            screen = TRUE,
                            cor_cutoff = cor_cutoff,
                            dev.ratio = model$dev.ratio,
                            df = model$df,
                            converge = ifelse(is.infinite(model$lambda) | is.null(model$lambda), FALSE, TRUE),
                            weight_only_predictor = TRUE
                        )
                    ), use.names = TRUE, fill = TRUE)
                )
            }
        )
    }

    return(out)
}
