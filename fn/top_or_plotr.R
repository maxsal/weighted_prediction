suppressPackageStartupMessages({
  library(ggpubr)
  library(ggnewscale)
  library(scales)
  library(logistf)
})

## HELPER FUNCTIONS
prep_top_or_data <- function(.phers_data, .phers_var = "phers", .probs = c(0.25, 0.1, 0.01), .cut_points) {
  
  cases    <- .phers_data[case == 1, ]
  controls <- .phers_data[case == 0, ]
  
  case_density    <- density(cases[, phers])
  control_density <- density(controls[, phers])
  
  case_out <- data.table(
    x    = case_density$x,
    y    = case_density$y,
    case = "Cases"
  )
  control_out <- data.table(
    x    = control_density$x,
    y    = control_density$y,
    case = "Controls"
  )
  
  out <- rbindlist(list(
    case_out,
    control_out
  ))
  
  out[, `:=` (
    t = cut(x, breaks = c(.cut_points, max(x, na.rm = TRUE)), labels = glue("Top {label_percent()(.probs)}"))
  )][]
  
}

approximate_density <- function(.input_dat, .probs = c(0.25, 0.1, 0.01), .cut_points) {
  # helper
  approximatr <- function(xs, ys, b) {
    out <- lapply(b,
                  \(x) approx(x = xs, y = ys, xout = x)$y)
    return(unlist(out))
  }
  
  out <- data.table(
    new_x = seq(min(.input_dat[, x], na.rm = TRUE), max(.input_dat[, x], na.rm = TRUE), length.out = 500)
  )[, `:=` (
    case_y    = approximatr(xs = .input_dat[case == "Cases", x], ys = .input_dat[case == "Cases", y], b = new_x),
    control_y = approximatr(xs = .input_dat[case == "Controls", x], ys = .input_dat[case == "Controls", y], b = new_x)
  )][, `:=` (
    max_y = fcase(
      is.na(case_y), control_y,
      is.na(control_y), case_y,
      case_y > control_y, case_y,
      control_y > case_y, control_y
    ),
    t = cut(new_x, breaks = c(.cut_points, max(new_x, na.rm = TRUE)), labels = glue("Top {label_percent()(.probs)}")))
  ][]
  
  return(out)
}

make_top_pct_indicators <- function(.phers_data, .phers_var = "phers", .probs, .cut_points) {
  out <- .phers_data
  for (i in seq_along(.cut_points)) {
    out[[paste0("top_", .probs[i])]] <- fifelse(.phers_data[[.phers_var]] >= .cut_points[i], 1, 0)
  }
  return(out[])
}

top_or_extractr <- function(x, p, r = 1) {
  mod <- logistf(glue("case ~ top_{p}"), data = x, control = logistf.control(maxit = 1000))
  suppressMessages(y <- confint(mod))
  data.table(
    or_est = exp(coef(mod)[[glue("top_{p}")]]),
    or_lo = exp(y[glue("top_{p}"), 1]),
    or_hi = exp(y[glue("top_{p}"), 2])
  )[, print := paste0(format(round(or_est, r), big.mark = ",", nsmall = r), " (", format(round(or_lo, r), big.mark = ",", nsmall = r), ", ", format(round(or_hi, r), big.mark = ",", nsmall = r), ")")][]
}

top_or_plotr <- function(phers_data, phers_var = "phers", probs = c(0.25, 0.1, 0.01),
                         .title = "PheRS distribution by case status",
                         .subtitle = "In DATA SAMPLE",
                         .caption = "") {
  
  cut_points <- quantile(phers_data[[phers_var]], 1 - probs)
  
  if (length(unique(cut_points)) != length(cut_points)) {
    cli_alert_warning("quantile cut points not unique - keeping lowest quantile of each unique cut point value")
    ind <- match(unique(cut_points), cut_points)
    cut_points <- cut_points[ind]
    probs      <- probs[ind]
  }
  
  prep_data <- prep_top_or_data(.phers_data = phers_data, .phers_var = phers_var, .probs = probs, .cut_points = cut_points)
  prep_data2 <- make_top_pct_indicators(.phers_data = phers_data, .probs = probs, .cut_points = cut_points)
  
  density_data <- approximate_density(.input_dat = prep_data, .probs = probs, .cut_points = cut_points)
  
  ors <- lapply(probs,
                \(x) {top_or_extractr(x = prep_data2, p = x)[, t := glue("Top {label_percent()(x)}")]}) |>
    rbindlist()
  
  or_table <- ggtexttable(ors[, .("Percentile" = t, "OR (95% CI)" = print)],
                          rows = NULL,
                          theme = ttheme(
                            colnames.style = colnames_style(fill = "white"),
                            tbody.style = tbody_style(color = "black",
                                                      fill = paste0(palette_OkabeIto[seq_along(cut_points)], "98"))))
  
  top_or_plot <- ggplot() +
    geom_ribbon(
      data = density_data[!is.na(t), ],
      aes(x = new_x, ymin = 0, ymax = max_y, fill = t), alpha = 0.4
    ) +
    geom_vline(
      data = density_data[!is.na(t), .SD[which.min(new_x)], by = t],
      aes(xintercept = new_x, color = t), linewidth = 1
    ) +
    annotation_custom(
      ggplotGrob(or_table),
      xmax = max(prep_data[, x], na.rm = TRUE),
      xmin = max(prep_data[, x], na.rm = TRUE) - 2, ymin = max(prep_data[, y], na.rm = TRUE) - 0.2,
      ymax = max(prep_data[, y], na.rm = TRUE)
    ) +
    scale_fill_OkabeIto() +
    scale_color_OkabeIto() +
    new_scale_color() +
    geom_line(data = prep_data, aes(x = x, y = y, group = factor(case), color = factor(case)), linewidth = 1) +
    scale_color_manual(values = c("red", "black")) +
    labs(
      x = "PheRS (mean standardized)",
      y = "Density",
      title = .title,
      subtitle = .subtitle,
      caption = .caption
    ) +
    cowplot::theme_half_open() +
    guides(color = guide_legend(order = 1)) +
    theme(
      legend.position = "top",
      legend.title = element_blank()
    )
  
  return(top_or_plot)
  
}
