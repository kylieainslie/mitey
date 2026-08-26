
plot_convergence_diagnostics <- function(history_data, 
                                         target_scenario, 
                                         tol_limit = 1e-6, 
                                         save_plot = FALSE, 
                                         output_dir = "notes/Simulation_results",
                                         file_name = NULL) {
  

  plot_data_long <- history_data %>%
    filter(scenario_id == target_scenario) %>%
    select(rep, iteration, mu_change, sigma_change, w_change) %>%
    pivot_longer(
      cols = ends_with("change"),
      names_to = "parameter",
      values_to = "change_magnitude"
    ) %>%
    filter(!is.na(change_magnitude))

  if (nrow(plot_data_long) == 0) {
    warning(sprintf("No data found for scenario %d. Returning NULL.", target_scenario))
    return(NULL)
  }


  summary_data <- plot_data_long %>%
    group_by(parameter, iteration) %>%
    summarise(
      median_change = median(change_magnitude, na.rm = TRUE),
      lower_bound = quantile(change_magnitude, 0.05, na.rm = TRUE),
      upper_bound = quantile(change_magnitude, 0.95, na.rm = TRUE),
      .groups = "drop"
    )


  convergence_rates <- plot_data_long %>%
    group_by(parameter, rep) %>%
    summarise(

      indiv_converged = last(change_magnitude, order_by = iteration) < tol_limit,
      .groups = "drop"
    ) %>%
    group_by(parameter) %>%
    summarise(
      success_count = sum(indiv_converged), 
      total_reps = n(),                     
      .groups = "drop"
    ) %>%
   
    mutate(label_text = paste0("Converged: ", success_count, " / ", total_reps))
  

  p <- ggplot() +
    geom_ribbon(
      data = summary_data,
      aes(x = iteration, ymin = lower_bound, ymax = upper_bound),
      fill = "steelblue", alpha = 0.3
    ) +
    geom_line(
      data = summary_data,
      aes(x = iteration, y = median_change),
      color = "#2c3e50", linewidth = 1.2
    ) +
    geom_hline(yintercept = tol_limit, linetype = "dashed", color = "#e74c3c", linewidth = 0.8) +
    geom_text(
      data = convergence_rates,
      aes(x = max(summary_data$iteration) * 0.9, y = 1e-1, label = label_text),
      hjust = 1, vjust = 1, size = 4.5, fontface = "bold", color = "#2c3e50"
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    facet_wrap(~ parameter, scales = "fixed", ncol = 1) +
    labs(
      title = paste("EM Algorithm Convergence Diagnostics (Scenario", target_scenario, ")"),
      subtitle = "Solid line: Median trajectory | Shaded area: 5th - 95th percentile range",
      x = "Iteration Number",
      y = "Parameter Change Magnitude (Log10 Scale)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
    

  if (save_plot) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    if (is.null(file_name)) {
      file_name <- sprintf("Convergence_Diagnostics_Scenario_%d.png", target_scenario)
    }
    
    full_path <- file.path(output_dir, file_name)
    ggsave(filename = full_path, plot = p, width = 8, height = 10, dpi = 300)
    message("✔ Plot successfully saved to: ", full_path)
  }
  
  return(p)
}


