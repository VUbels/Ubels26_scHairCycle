###############################################################
# AMBIENT RNA REMOVAL USING CELLBENDER
###############################################################

remove_ambient_rna <- function(input_folder, output_folder = "./", epochs = 150, cellbender_learning_rate = 0.0001) {
  
  cellbender <- reticulate::import("cellbender")
  
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

###############################################################
# COMPRESS CELLBENDER OUTPUT FOR SEURAT
###############################################################
make_seurat_compatible <- function(output_folder) {
  
  # Find all CellBender filtered output files
  cellbender_output_dirs <- list.dirs(
    file.path(output_folder, "preprocessing"),
    recursive = FALSE
  )
  
  for (cb_dir in cellbender_output_dirs) {
    
    # Look for the _filtered.h5 file
    filtered_files <- list.files(
      cb_dir,
      pattern = "_filtered\\.h5$",
      full.names = TRUE
    )
    
    if (length(filtered_files) == 0) {
      warning("No filtered.h5 file found in: ", cb_dir)
      next
    }
    
    for (input_file in filtered_files) {
      
      # Create output filename
      output_file <- sub("_filtered\\.h5$", "_filtered_seurat.h5", input_file)
      
      # Build ptrepack command
      cmd <- sprintf(
        "ptrepack --complevel 5 %s:/matrix %s:/matrix",
        shQuote(input_file),
        shQuote(output_file)
      )
      
      # Run compression
      message("Compressing: ", basename(input_file))
      exit_code <- system(cmd)
      
      if (exit_code != 0) {
        warning("ptrepack failed for: ", basename(input_file))
      } else {
        message("Created: ", basename(output_file))
      }
    }
  }
}