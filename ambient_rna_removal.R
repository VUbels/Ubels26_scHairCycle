###############################################################
# AMBIENT RNA REMOVAL USING CELLBENDER
###############################################################

remove_ambient_rna <- function(input_folder, output_folder = "./", epochs = 150, cellbender_learning_rate = 0.0001) {
  
  datasets <- list.files(input_folder, pattern = "\\.h5$", recursive = FALSE)
  
  for (dataset in datasets) {
    
    # Construct paths
    cellbender_input_path <- file.path(input_folder, dataset)
    dataset_base <- sub("\\.h5$", "", dataset)
    
    # Create output directory
    cellbender_output_dir <- file.path(output_folder, "preprocessing/", paste0(dataset_base, "_cellbender_results"))
    dir.create(cellbender_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Output file inside that directory
    output_name <- paste0(dataset_base, "_postcellbender.h5")
    cellbender_output_path <- file.path(cellbender_output_dir, output_name)
    
    # Build CellBender command
    cmd <- sprintf(
      "%s remove-background --input %s --output %s --cuda --epochs %d --learning-rate %f",
      conda_info_env$cellbender_bin,
      shQuote(cellbender_input_path),
      shQuote(cellbender_output_path),
      epochs,
      cellbender_learning_rate)
    
    # Run CellBender
    message("Processing: ", dataset)
    exit_code <- system(cmd)
    
    if (exit_code != 0) {
      warning("CellBender failed for: ", dataset)
    } else {
      message("Completed: ", output_name)
    }
  }
}