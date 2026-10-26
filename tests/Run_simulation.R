# Run simulation

source("tests/simu_full_function.R")

mitey_pkg_root <- r"(C:\Users\Administrator.DESKTOP-JQG2K7N\桌面\2025 melb\master_project\mitey\mitey)"

my_output_dir <- file.path(
  mitey_pkg_root,
  "notes",
  "Simulation_results"
)

# N  = 200
results <- run_full_simulation(
  pkg_root = mitey_pkg_root,
  results_path = my_output_dir, # Custom output path
  N_sample = 200,
  reps = 100,
  n_iter = 300,
  seed = 2026
)

# N = 100
results <- run_full_simulation(
  pkg_root = mitey_pkg_root,
  results_path = my_output_dir, # Custom output path
  N_sample = 100,
  reps = 100,
  n_iter = 300,
  seed = 2026
)

# N= 50
results <- run_full_simulation(
  pkg_root = mitey_pkg_root,
  results_path = my_output_dir, # Custom output path
  N_sample = 50,
  reps = 100,
  n_iter = 300,
  seed = 2026
)

# N = 20
results <- run_full_simulation(
  pkg_root = mitey_pkg_root,
  results_path = my_output_dir, # Custom output path
  N_sample = 20,
  reps = 100,
  n_iter = 300,
  seed = 2026
)

# N  = 1000
results <- run_full_simulation(
  pkg_root = mitey_pkg_root,
  results_path = my_output_dir, # Custom output path
  N_sample = 1000,
  reps = 5,
  n_iter = 200,
  seed = 2026
)





source(file.path(mitey_pkg_root, "tests", "plot&format_results.R"))

combined_plot_dir <- file.path(my_output_dir, "Combined_summary_panels")

simulation_workbooks <- list.files(
  my_output_dir,
  pattern = "^Simu_results_s[0-9]+_r[0-9]+_i[0-9]+\\.xlsx$",
  full.names = TRUE
)

for (input_file in simulation_workbooks) {
  generate_simulation_report(
    input_path = input_file,
    output_dir = my_output_dir,
    generate_tables = FALSE,
    combined_output_dir = combined_plot_dir,
    save_individual_plots = FALSE
  )
}
