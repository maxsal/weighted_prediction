library(glmnet)

set.seed(42)

n      <- 1000
p      <- 5000
real_p <- 15

x <- matrix(rnorm(n*p), nrow = n, ncol = p)
y <- as.numeric(apply(x[, 1:real_p], 1, sum) + rnorm(n) > 0)

train_rows <- sample(1:n, 0.66 * n)
x.train <- x[train_rows, ]
x.test  <- x[-train_rows, ]

y.train <- y[train_rows]
y.test  <- y[-train_rows]

# ridge
alpha0.fit <- cv.glmnet(x.train, y.train, type.measure = "deviance", alpha = 0, family = "binomial") # set type = 'deviance' for logistic regression

alpha0.predicted <- predict(alpha0.fit, newx = x.test, s = "lambda.1se")

###
# from cv.lognet.R: https://github.com/cran/glmnet/blob/master/R/cv.lognet.R
# how (capped) binomial deviance is calculated
# deviance <- {
#     predmat <- pmin(pmax(predmat, prob_min), prob_max)
#     lp <- y[, 1] * log(1 - predmat) + y[, 2] * log(predmat)
#     ly <- log(y)
#     ly[y == 0] <- 0
#     ly <- drop((y * ly) %*% c(1, 1))
#     2 * (ly - lp)
# }
###
mean((y.test - alpha0.predicted)^2)

# lasso
alpha1.fit <- cv.glmnet(x.train, y.train, type.measure = "deviance", alpha = 1, family = "binomial")
plot(alpha1.fit)

alpha1.predicted <- predict(alpha1.fit, newx = x.test, s = "lambda.1se")

mean((y.test - alpha1.predicted)^2)

# elastic net - set alpha
alpha0.5.fit <- cv.glmnet(x.train, y.train, type.measure = "deviance", alpha = 0.5, family = "binomial")

alpha0.5.predicted <- predict(alpha0.5.fit, newx = x.test, s = "lambda.1se")

deviance(alpha0.5.fit$glmnet.fit)
mean((y.test - alpha0.5.predicted)^2)

fitted <- glmnet(x.train, y.train, type.measure = "deviance", alpha = 0.5, family = "binomial", lambda = alpha0.5.fit$lambda.1se)
deviance(fitted)

# elastic net - try alpha
try_alpha <- seq(0, 1, by = 0.05)
list.of.fits <- list()
cli_progress_bar(
    name = "alpha tuning",
    format = "alpha: {try_alpha[i]} {cli::pb_bar} {cli::pb_percent} {cli::pb_eta}",
    total = length(try_alpha)
)
for (i in seq_along(try_alpha)) {
    fit.name <- paste0("alpha", try_alpha[i])

    list.of.fits[[fit.name]] <- cv.glmnet(x.train, y.train, type.measure = "mse", alpha = try_alpha[i], family = "gaussian")

    cli_progress_update()
}
cli_progress_done()

results <- data.table()
for (i in seq_along(try_alpha)) {
    cli_alert("alpha = {try_alpha[i]}")
    fit.name <- paste0("alpha", try_alpha[i])

    predicted <- predict(list.of.fits[[fit.name]], newx = x.test, s = "lambda.1se")
    mse       <- mean((y.test - predicted)^2)

    temp      <- data.table(alpha = try_alpha[i], mse = mse, fit.name = fit.name)
    results   <- rbind(results, temp)
}

results[which.min(results$mse), ]

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

        temp    <- data.table(
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

results <- cva.glmnet(x = x.train, y = y.train, type.measure = "deviance", fam = "binomial", nfolds = 10)

results |>
    ggplot(aes(x = alpha, y = measure)) +
    geom_vline(xintercept = results[which.min(results$mse), alpha], color = "red") +
    geom_line() +
    geom_point() +
    theme_ms()
