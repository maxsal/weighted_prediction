# plot variable importance for top n
suppressPackageStartupMessages({
  require(vip)
  require(cowplot)
  require(stringr)
})

plot_vip_with_description <- function(tidy_fit, pheinfo, n = 20, relative = TRUE) {
  
  # get top n variables
  top_n <- tidy_fit |>
    extract_fit_parsnip() |>
    vip::vi() |>
    mutate(`Relative importance` = Importance / max(Importance, na.rm = TRUE)) |>
    head(n)
  
  # merge in descriptions
  plt_data <- top_n |>
    left_join(pheinfo[, .(phecode, description, color, group)],
              by = c("Variable" = "phecode")) |>
    mutate(Variable = paste0(Variable, ": ", description))
  
  cols <- unique(pheinfo[group %in% plt_data[["group"]], .(group, color)])
  colors <- cols[, color]
  names(colors) <- cols[, group]

  if (relative == TRUE) {  
    plt_data |>
      ggplot(aes(x = reorder(Variable, `Relative importance`), y = `Relative importance`)) +
      geom_bar(aes(fill = group), stat = "identity") +
      scale_fill_manual(values = colors) +
      scale_x_discrete(labels = \(x) stringr::str_wrap(x, width = 50)) +
      labs(
        title = "Variable importance plot",
        x = ""
      ) +
      guides(fill=guide_legend(ncol = 3, byrow = TRUE)) +
      coord_flip() +
      cowplot::theme_minimal_grid() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )
  } else {
    plt_data |>
      ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_bar(aes(fill = group), stat = "identity") +
      scale_fill_manual(values = colors) +
      scale_x_discrete(labels = \(x) stringr::str_wrap(x, width = 50)) +
      labs(
        title = "Variable importance plot",
        x = ""
      ) +
      guides(fill=guide_legend(ncol = 3, byrow = TRUE)) +
      coord_flip() +
      cowplot::theme_minimal_grid() +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )
  }
  
}