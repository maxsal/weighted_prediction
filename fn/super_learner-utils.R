# source: https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html#fit-ensemble-with-external-cross-validation
# Review meta-weights (coefficients) from a CV.SuperLearner object
review_weights <- function(cv_sl) {
  meta_weights <- coef(cv_sl)
  means        <- colMeans(meta_weights)
  sds          <- apply(meta_weights, MARGIN = 2,  FUN = sd)
  mins         <- apply(meta_weights, MARGIN = 2, FUN = min)
  maxs         <- apply(meta_weights, MARGIN = 2, FUN = max)
  # Combine the stats into a single matrix.
  sl_stats <- cbind("mean(weight)" = means, "sd" = sds, "min" = mins, "max" = maxs)
  # Sort by decreasing mean weight.
  sl_stats[order(sl_stats[, 1], decreasing = TRUE), ]
}